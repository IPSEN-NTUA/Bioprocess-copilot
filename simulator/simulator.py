"""
Fed-batch xylitol bioreactor simulator
======================================

Process simulator for the BIO40 RAG / KG / Hybrid benchmark. It generates the
*true* state of a yeast xylitol cultivation (batch, fed-batch, pulsed batch) from
which all synthetic data sources (HPLC, online logs, off-gas, lab notes, ...)
are later rendered.

Core kinetics
-------------
VTT kinetic **model 31** with the **BC18 refined parameters** (realdata_summary.json),
kept exactly as fitted:

    mu1  = mumax1 * S1 / (KS1 + S1)                                   (growth on glucose)
    qP31 = mumax2 * S2 / ((KP + S2 + S2^2/KSI2) * (1 + P/KPI2))        (xylitol formation)
    dX/dt  = mu1 X - kd X - (F/V) X     (+ inhibitor death to the non-viable pool Xd)
    dS1/dt = -(1/YXS1) mu1 X + (F/V)(S1f - S1)
    dS2/dt = -(1/YPS2) qP31 X + (F/V)(S2f - S2)
    dP/dt  = qP31 X - (F/V) P

YPS2 = 1.24 g/g is kept as fitted even though it exceeds the stoichiometric
maximum (1.013 g/g); see `carbon_balance()` for how the excess is reported.

Extensions (all reduce to model 31 at the BC18 reference conditions:
30 degC, pH 5.5, aerobic, glucose-limited, no inhibitors)
----------------------------------------------------------------------------
1. Crabtree overflow: an extra glucose uptake channel that only opens at high
   glucose (Hill, n=2) and goes to ethanol, glycerol, CO2, acetate, succinate,
   pyruvate (fermentative C-mol split) with a small biomass yield.
2. Pasteur effect: the model-31 (oxidative) glucose channel needs O2; under O2
   limitation the non-respired part is diverted to the fermentative channel.
3. Diauxic re-consumption of ethanol, glycerol, acetate and pyruvate, repressed
   by glucose through a slow induction state (enzyme adaptation, tau ~ 1.5 h).
4. Cofactor coupling of XR: xylitol formation is multiplied by
   phi = A_C / (A_C + K_AC), A_C = carbon uptake of external co-substrates
   (glucose, ethanol, glycerol, acetate, pyruvate, formate via FDH), capped
   by the total electron supply (soft-min, so O2 uptake can never go negative).
5. Oxygen: dissolved-O2 balance with kLa(stirrer, airflow), temperature-dependent
   saturation, optional DO-stirrer cascade. O2 uptake and CO2 evolution follow
   from degree-of-reduction and carbon balances of every flux, so off-gas,
   HPLC and biomass are mutually consistent.
6. Lignocellulosic hydrolysate (CIMV C5 fraction): formate & acetate toxicity via
   the undissociated acid (pH dependent), furfural / HMF inhibition and their
   NAD(P)H-consuming reduction to furan alcohols, FDH formate oxidation,
   arabinose (not consumed), xylo-oligomers released by added enzyme, and a
   tracer for "% CIMV in the reactor".
7. Temperature (CTMI) and pH (CPM) cardinal factors, normalised to 1 at 30 degC
   and pH 5.5.
8. Reactor physics: feed, base and acid addition, sampling, pulses, water
   evaporation through the condenser, ethanol stripping, pH control by base/acid
   computed from proton production (secreted acids, NH4+ assimilation).

Units: concentrations g/L (DO in mmol/L), volume L, time h, specific rates g/g/h
unless noted (mmol/g/h for O2/CO2/electrons, C-mmol/g/h for carbon fluxes).

Strain design (NRA on a large-scale kinetic model) will later perturb the Strain
parameters; this module only needs a Strain object.
"""
from __future__ import annotations

import dataclasses as dc
import math
from typing import Callable

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Physical constants and compound data
# ---------------------------------------------------------------------------
MW = dict(Glc=180.156, Xyl=150.130, Xol=152.146, Ara=150.130, Olig=150.130,
          EtOH=46.068, Gly=92.094, Ac=60.052, Suc=118.088, Pyr=88.062,
          For=46.025, Fur=96.085, HMF=126.110, FurOH=98.101, HMFOH=128.127,
          CO2=44.010, O2=31.998)
NC = dict(Glc=6, Xyl=5, Xol=5, Ara=5, Olig=5, EtOH=2, Gly=3, Ac=2, Suc=4, Pyr=3,
          For=1, Fur=5, HMF=6, FurOH=5, HMFOH=6)
# degree of reduction per C-mol (NH3 as N reference)
GAMMA = dict(Glc=4.0, Xyl=4.0, Xol=4.4, Ara=4.0, Olig=4.0, EtOH=6.0, Gly=14 / 3,
             Ac=4.0, Suc=3.5, Pyr=10 / 3, For=2.0, Fur=4.0, HMF=4.0, FurOH=4.4,
             HMFOH=13 / 3, X=4.2)
G_PER_CMOL_X = 24.6          # g DW per C-mol biomass (CH1.8O0.5N0.2 + ash)
N_PER_CMOL_X = 0.2           # mol N per C-mol biomass
PKA = dict(For=3.75, Ac=4.76, Pyr=2.50, Suc1=4.21, Suc2=5.64)


def g_per_cmol(sp: str) -> float:
    return MW[sp] / NC[sp]


def undissociated(pH: float, pKa: float) -> float:
    return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def succinate_charge(pH: float) -> float:
    """Mean negative charge of succinate at pH (protons released per mol secreted)."""
    h = 10.0 ** -pH
    k1, k2 = 10.0 ** -PKA["Suc1"], 10.0 ** -PKA["Suc2"]
    d = h * h + k1 * h + k1 * k2
    return (k1 * h + 2 * k1 * k2) / d


def o2_saturation_mmolL(T: float) -> float:
    """Air-saturated dissolved O2 in water (mmol/L), Benson-Krause polynomial fit."""
    mgL = 14.62 - 0.3898 * T + 0.006969 * T ** 2 - 5.897e-5 * T ** 3
    return mgL / MW["O2"]


def water_vapour_fraction(T: float) -> float:
    """Saturated water vapour mole fraction at 1 atm (Antoine, T in degC)."""
    p_mmHg = 10 ** (8.07131 - 1730.63 / (233.426 + T))
    return p_mmHg / 760.0


def cardinal(x: float, xmin: float, xopt: float, xmax: float) -> float:
    """Rosso cardinal model (CTMI / CPM), 0 outside (xmin, xmax), 1 at xopt."""
    if x <= xmin or x >= xmax:
        return 0.0
    num = (x - xmax) * (x - xmin) ** 2
    den = (xopt - xmin) * ((xopt - xmin) * (x - xopt) - (xopt - xmax) * (xopt + xmin - 2 * x))
    return max(num / den, 0.0)


# ---------------------------------------------------------------------------
# Strain
# ---------------------------------------------------------------------------
@dc.dataclass
class Strain:
    """Kinetic description of a production strain.

    The model-31 block carries the VTT fit for BC18 (H1274, refined parameters).
    The other blocks are physiology extensions with literature-typical values for
    S. cerevisiae; they are the parameters the NRA strain-design step will perturb.
    """
    name: str = "H1274"
    species: str = "Saccharomyces cerevisiae"
    genotype: str = "VTT B-67002 + plasmid XR (Scheffersomyces stipitis XYL1)"

    # --- model 31: BC18 refined fit (realdata_summary.json), re-tuned within +/-20 %
    #     against BC18 with the extended simulator (tune_bc18.py). Originals in MODEL31_REFINED.
    mumax1: float = 0.395485  # 1/h (+19.5 % vs refined 0.3308)
    mumax2: float = 4.79508  # g xylitol / g DW / h (+20.0 % vs refined 3.997)
    YXS1: float = 0.768728  # g DW / g glucose (+11.5 % vs refined 0.6896)
    YPS2: float = 1.08357  # g xylitol / g xylose (-12.6 % vs refined 1.24)
    KS1: float = 0.0813307  # g/L (-20.0 % vs refined 0.1017)
    KP: float = 32.001  # g/L (-20.0 % vs refined 40)
    kd: float = 0.0252444  # 1/h (+20.0 % vs refined 0.02104)
    KSI2: float = 29.0019  # g/L (-14.3 % vs refined 33.84)
    KPI2: float = 13.7017  # g/L (+19.9 % vs refined 11.43)

    # --- 1. Crabtree overflow --------------------------------------------------
    qGof_max: float = 1.6      # g glucose / g / h
    KS_of: float = 1.5         # g/L (half-saturation of overflow channel)
    n_of: float = 2.0          # Hill coefficient (overflow opens only at high glucose)
    # fermentative channel, C-mol fractions of glucose carbon (sum = 1)
    ferm_split: dict = dc.field(default_factory=lambda: dict(
        X=0.12, EtOH=0.53, CO2=0.287, Gly=0.045, Ac=0.010, Suc=0.004, Pyr=0.004))
    # oxidative (model 31) channel: succinate leak taken from the CO2 share
    ox_suc_frac: float = 0.003

    # --- 2. Oxygen --------------------------------------------------------------
    K_O2: float = 0.0015       # mmol/L (~0.6 % air saturation), glucose respiration
    K_O2_E: float = 0.04       # mmol/L, O2 affinity of ethanol/glycerol/acetate oxidation
                               # (lower affinity -> glucose gets O2 priority under limitation)

    # --- 3. Re-consumption of C2/C3 by-products (glucose repressed) -------------
    qE_max: float = 0.15       # g EtOH / g / h
    K_E: float = 0.10          # g/L
    qGly_max: float = 0.012    # g/g/h
    K_Gly: float = 0.5
    qAc_max: float = 0.02
    K_Ac: float = 0.2
    qPyr_max: float = 0.005
    K_Pyr: float = 0.5
    Ki_rep: float = 0.10       # g/L glucose giving half repression
    tau_ind: float = 1.5       # h, induction time constant (diauxic lag)
    # C-mol fates of re-consumed substrates
    YX_E: float = 0.60         # C-mol X / C-mol EtOH (~0.64 g/g)
    Ac_from_E: float = 0.05    # C-mol acetate / C-mol EtOH (Ald leak)
    YX_Gly: float = 0.40
    YX_Ac: float = 0.30
    YX_Pyr: float = 0.30
    E_max: float = 90.0        # g/L ethanol stopping all activity

    # --- 4. Cofactor coupling of XR --------------------------------------------
    K_AC: float = 0.49         # C-mmol / g / h co-substrate uptake giving half XR activation
                               # (phi ~ 0.71 in the BC18 fed phase, as with the earlier e- form)
    theta_e: float = 0.90      # max share of total electron supply usable by reductases
    f_endo: float = 0.60       # share of decayed biomass that is respired (rest lysed)
    k_lys_dead: float = 0.005  # 1/h lysis of non-viable (inhibitor-killed) cells, half-life ~6 d
    kd_tox_max: float = 0.30   # 1/h ceiling of inhibitor-induced death (half-life >= 2.3 h)

    # --- 6. Hydrolysate inhibitors & detox --------------------------------------
    Ki_ForU: float = 0.030     # g/L undissociated formic acid (Hill n=2)
    Ki_AcU: float = 1.0        # g/L undissociated acetic acid
    Ki_Fur: float = 1.0        # g/L furfural
    Ki_HMF: float = 2.0        # g/L HMF
    Kd_ForU: float = 0.05      # g/L undissociated formic acid doubling death (n=2)
    Kd_Fur: float = 2.0        # g/L furfural doubling death
    fdh: float = 0.0           # 0 = native only, 1 = heterologous FDH expressed
    qFor_native: float = 0.002 # g/g/h
    qFor_FDH: float = 0.25     # g/g/h at fdh = 1
    K_For: float = 1.0         # g/L
    qFur_max: float = 0.10     # g/g/h furfural -> furfuryl alcohol
    K_Fur: float = 0.2
    qHMF_max: float = 0.04     # g/g/h HMF -> HMF alcohol
    K_HMF: float = 0.3

    # --- 8. Production-phase modifiers (strain design, see strain_design.py) -----
    # Applied only in the glucose-derepressed production state: each multiplier m acts as
    # m**Ind, so at Ind = 0 (batch, glucose excess) the parent kinetics are untouched.
    prod_mumax1_mult: float = 1.0      # growth capacity on glucose
    prod_YXS1_mult: float = 1.0        # biomass yield on glucose
    prod_suc_mult: float = 1.0         # succinate leak of the oxidative channel
    prod_split: dict | None = None     # C-mol split of O2-limited (Pasteur) glucose; None = ferm_split
    prod_div_mult: float = 1.0         # scales the share of glucose diverted to fermentation under O2 limit
    design_note: str = ""

    # --- 7. Cardinal temperature / pH ------------------------------------------
    T_card: tuple = (5.0, 32.0, 40.0)
    pH_card: tuple = (2.5, 5.0, 8.5)
    T_ref: float = 30.0
    pH_ref: float = 5.5

    def copy(self, **changes) -> "Strain":
        return dc.replace(self, **changes)


MODEL31_REFINED = dict(mumax1=0.3308406077952456, mumax2=3.9973119444823175, YXS1=0.6895977214979163, YPS2=1.240136074448464, KS1=0.10166337963741057, KP=40.0, kd=0.021038494365469833, KSI2=33.84249597051714, KPI2=11.429154090271716)

H1274 = Strain()


# ---------------------------------------------------------------------------
# Reactor
# ---------------------------------------------------------------------------
@dc.dataclass
class Reactor:
    name: str = "Labfors 1"
    vendor: str = "Infors HT"
    V_total: float = 1.0          # L vessel
    V_work_max: float = 0.80      # L, feed and pulses are blocked above this volume
    kLa_ref: float = 88.0         # 1/h at N_ref, vvm_ref, V_ref, no antifoam (BC18 pO2 crash 8-17 h)
    N_ref: float = 900.0          # rpm
    vvm_ref: float = 2.33         # 0.7 L/min in 0.3 L
    a_N: float = 1.2              # kLa ~ N^a (P/V ~ N^3, van 't Riet exponent 0.4)
    c_V: float = 0.4              # kLa ~ V^-c at fixed N (P/V falls as broth volume grows)
    V_ref: float = 0.30           # L, broth volume at which kLa_ref holds
    b_Q: float = 0.5              # kLa ~ superficial gas velocity^b
    T_condenser: float = 4.0      # degC off-gas condenser outlet
    K_AF: float = 100.0           # ppm antifoam halving kLa (silicone/PPG antifoams cut kLa 30-60 %)
    fixed_kLa: float | None = None  # 1/h, e.g. shake flasks (no stirrer / sparger correlation)
    evap_L_h: float | None = None   # L/h, overrides the condenser evaporation model (flasks)
    k_strip_ref: float = 0.002    # 1/h ethanol stripping at vvm_ref

    def kLa(self, N: float, Q_Lmin: float, V: float) -> float:
        vvm = Q_Lmin / max(V, 1e-6)
        return (self.kLa_ref * (N / self.N_ref) ** self.a_N * (V / self.V_ref) ** -self.c_V
                * (vvm / self.vvm_ref) ** self.b_Q)


LABFORS = Reactor()
MULTIFORS = Reactor(name="Multifors 2", V_total=0.75, kLa_ref=80.0)


# ---------------------------------------------------------------------------
# Experiment recipe
# ---------------------------------------------------------------------------
SPECIES = ["Glc", "Xyl", "Xol", "Ara", "Olig", "EtOH", "Gly", "Ac", "Suc", "Pyr",
           "For", "Fur", "HMF", "FurOH", "HMFOH"]
STATES = (["V", "X", "Xd"] + SPECIES +
          ["Trc", "Enz", "AF", "DO", "Ind",
           "cumFeed", "cumBase", "cumAcid", "cumEvap", "cumO2", "cumCO2",
           "cumStripE", "cumLysis", "cumXol"])
IDX = {s: i for i, s in enumerate(STATES)}
# states that are diluted / concentrated by volume changes
LIQ = ["X", "Xd"] + SPECIES + ["Trc", "Enz", "AF"]


@dc.dataclass
class FeedSegment:
    t_start: float
    t_end: float
    rate: float                      # L/h at t_start
    mu_exp: float = 0.0              # exponential feeding: F = rate*exp(mu_exp*(t-t_start))

    def F(self, t: float) -> float:
        if self.t_start <= t < self.t_end:
            return self.rate * math.exp(self.mu_exp * (t - self.t_start))
        return 0.0


@dc.dataclass
class Pulse:
    t: float
    volume: float                    # L
    composition: dict                # g/L of SPECIES (+ "Trc" in %)
    label: str = "pulse"


@dc.dataclass
class Experiment:
    run_id: str
    strain: Strain = dc.field(default_factory=Strain)
    reactor: Reactor = dc.field(default_factory=Reactor)
    V0: float = 0.300                # L
    X0: float = 0.9                  # g/L DW after inoculation
    medium: dict = dc.field(default_factory=lambda: {"Glc": 20.0})
    N_source: str = "YP"             # "YP" or "NH4"
    pH_control: bool = True          # False: no base/acid dosing (buffered shake flasks)
    T_set: float = 30.0
    T_schedule: list = dc.field(default_factory=list)   # [(t, T)]
    pH_set: float = 5.5
    base_M: float = 2.0              # mol/L OH- (KOH) or NH4OH
    acid_M: float = 2.0              # mol/L H+ (1 M H2SO4)
    stirrer: float = 900.0           # rpm (fixed, or minimum when cascade on)
    do_cascade: dict | None = None   # {"setpoint": 30, "N_min": 300, "N_max": 1200, "gain": 40}
    airflow: float = 0.7             # L/min
    feed_segments: list = dc.field(default_factory=list)
    feed_composition: dict = dc.field(default_factory=dict)  # g/L (+ "Trc" %)
    feed_density: float = 1.12       # kg/L, for "feed vessel (g)" bookkeeping
    feed_vessel0: float = 800.0      # g on balance at start
    pulses: list = dc.field(default_factory=list)
    enzyme_additions: list = dc.field(default_factory=list)  # [(t, U/L)]
    antifoam_additions: list = dc.field(default_factory=list)  # [(t, ppm in broth)]
    k_enz: float = 0.02              # L/(U h) hydrolysis rate constant
    K_olig: float = 5.0              # g/L
    k_enz_deact: float = 0.01        # 1/h
    sample_times: list = dc.field(default_factory=list)
    sample_volume: float = 0.007     # L
    t_end: float = 150.0
    notes: str = ""

    def F_feed(self, t: float) -> float:
        return sum(seg.F(t) for seg in self.feed_segments)

    def T_at(self, t: float) -> float:
        T = self.T_set
        for ts, Ts in sorted(self.T_schedule):
            if t >= ts:
                T = Ts
        return T


# ---------------------------------------------------------------------------
# Kinetics
# ---------------------------------------------------------------------------
def _hill(S: float, K: float, n: float) -> float:
    S = max(S, 0.0)
    return S ** n / (K ** n + S ** n) if S > 0 else 0.0


def _mm(S: float, K: float) -> float:
    S = max(S, 0.0)
    return S / (K + S)


def rates(t: float, y: np.ndarray, ex: Experiment, full: bool = False):
    """Right-hand side. With full=True also returns a dict of diagnostics."""
    s = ex.strain
    r = ex.reactor
    c = {k: max(y[IDX[k]], 0.0) for k in STATES}
    V = max(y[IDX["V"]], 1e-6)
    X = c["X"]
    T = ex.T_at(t)
    pH = ex.pH_set

    # ---- environment factors ---------------------------------------------
    fT = cardinal(T, *s.T_card) / cardinal(s.T_ref, *s.T_card)
    fpH = cardinal(pH, *s.pH_card) / cardinal(s.pH_ref, *s.pH_card)
    ForU = c["For"] * undissociated(pH, PKA["For"])
    AcU = c["Ac"] * undissociated(pH, PKA["Ac"])
    f_inh = (1.0 / (1.0 + (ForU / s.Ki_ForU) ** 2)
             * 1.0 / (1.0 + AcU / s.Ki_AcU)
             * s.Ki_Fur / (s.Ki_Fur + c["Fur"])
             * s.Ki_HMF / (s.Ki_HMF + c["HMF"])
             * max(0.0, 1.0 - c["EtOH"] / s.E_max))
    f_env = fT * fpH * f_inh
    DOsat = o2_saturation_mmolL(T)
    fO = _mm(c["DO"], s.K_O2)

    # ---- production-phase (derepressed) state and strain-design modifiers --
    ind = min(max(c["Ind"], 0.0), 1.0)
    mumax1 = s.mumax1 * s.prod_mumax1_mult ** ind
    YXS1 = s.YXS1 * s.prod_YXS1_mult ** ind
    YXS1 = min(YXS1, 0.95 * G_PER_CMOL_X / g_per_cmol("Glc"))       # keep e- supply >= 0
    suc_ox = s.ox_suc_frac * s.prod_suc_mult ** ind
    fs = s.ferm_split
    ps = s.prod_split or fs
    split_div = {k: (1.0 - ind) * fs[k] + ind * ps[k] for k in fs}   # Pasteur (O2-limited) split

    # ---- glucose: model-31 oxidative channel + overflow + Pasteur diversion -
    qG31 = (mumax1 / YXS1) * _mm(c["Glc"], s.KS1) * f_env            # g/g/h
    qG_div = qG31 * min(1.0, (1.0 - fO) * s.prod_div_mult ** ind)
    qG_ox = qG31 - qG_div
    qG_of = s.qGof_max * _hill(c["Glc"], s.KS_of, s.n_of) * f_env
    qG_f = qG_of + qG_div
    mu_ox = YXS1 * qG_ox                                              # model 31 growth

    # ---- re-consumption (glucose repressed, induction state Ind) ----------
    resp = ind * _mm(c["DO"], s.K_O2_E) * f_env
    qE = s.qE_max * _mm(c["EtOH"], s.K_E) * resp
    qGly = s.qGly_max * _mm(c["Gly"], s.K_Gly) * resp
    qAc = s.qAc_max * _mm(c["Ac"], s.K_Ac) * resp
    qPyr = s.qPyr_max * _mm(c["Pyr"], s.K_Pyr) * resp

    # C-mmol/g/h fluxes
    C = lambda q, sp: q / g_per_cmol(sp) * 1000.0
    cG_ox, cG_f = C(qG_ox, "Glc"), C(qG_f, "Glc")
    cG_of, cG_div = C(qG_of, "Glc"), C(qG_div, "Glc")
    cE, cGly, cAc, cPyr = C(qE, "EtOH"), C(qGly, "Gly"), C(qAc, "Ac"), C(qPyr, "Pyr")
    YX_ox_C = YXS1 * g_per_cmol("Glc") / G_PER_CMOL_X
    alloc = {k: fs[k] * cG_of + split_div[k] * cG_div for k in fs}   # fermentative C-mmol/g/h

    # biomass growth (C-mmol/g/h -> 1/h)
    cX = (YX_ox_C * cG_ox + alloc["X"] + s.YX_E * cE + s.YX_Gly * cGly
          + s.YX_Ac * cAc + s.YX_Pyr * cPyr)
    mu_extra = (cX - YX_ox_C * cG_ox) * G_PER_CMOL_X / 1000.0
    mu = mu_ox + mu_extra

    # ---- death -------------------------------------------------------------
    # basal decay (kd, model 31) is endogenous respiration + lysis; inhibitor-induced
    # death (formic acid, furfural) turns viable cells into non-viable but intact cells
    # (Xd), which still count in dry weight / OD and lyse slowly.
    z = s.kd * ((ForU / s.Kd_ForU) ** 2 + c["Fur"] / s.Kd_Fur)
    kd_tox = z / (1.0 + z / s.kd_tox_max)                              # ~linear when mild, saturates
    kd_eff = s.kd + kd_tox
    c_decay = s.kd / G_PER_CMOL_X * 1000.0                            # C-mmol/g/h (basal only)

    # ---- detox and FDH -----------------------------------------------------
    qFor = (s.qFor_native + s.fdh * s.qFor_FDH) * _mm(c["For"], s.K_For) * fT
    qFur = s.qFur_max * _mm(c["Fur"], s.K_Fur) * fT
    qHMF = s.qHMF_max * _mm(c["HMF"], s.K_HMF) * fT

    # ---- electron (NAD(P)H) supply -----------------------------------------
    gX = GAMMA["X"]
    e_ox = cG_ox * (GAMMA["Glc"] - YX_ox_C * gX - suc_ox * GAMMA["Suc"])
    e_f = cG_f * GAMMA["Glc"] - sum(alloc[k] * (gX if k == "X" else GAMMA[k])
                                    for k in alloc if k != "CO2")
    e_E = cE * (GAMMA["EtOH"] - s.YX_E * gX - s.Ac_from_E * GAMMA["Ac"])
    e_Gly = cGly * (GAMMA["Gly"] - s.YX_Gly * gX)
    e_Ac = cAc * (GAMMA["Ac"] - s.YX_Ac * gX)
    e_Pyr = cPyr * (GAMMA["Pyr"] - s.YX_Pyr * gX)
    e_For = qFor / MW["For"] * 1000.0 * 2.0
    e_endo = s.f_endo * c_decay * gX
    A_ext = e_ox + e_E + e_Gly + e_Ac + e_Pyr + e_For                 # co-substrate e- supply
    A_C = cG_ox + cE + cGly + cAc + cPyr + qFor / MW["For"] * 1000.0  # co-substrate C uptake
    E_tot = A_ext + e_f + e_endo

    # ---- xylitol (model 31 x environment x cofactor) -----------------------
    qP31 = s.mumax2 * c["Xyl"] / ((s.KP + c["Xyl"] + c["Xyl"] ** 2 / s.KSI2) * (1.0 + c["Xol"] / s.KPI2))
    phi = A_C / (A_C + s.K_AC) if A_C > 0 else 0.0      # XR needs active co-substrate catabolism
    qP_dem = qP31 * f_env * phi
    D = (2.0 * qP_dem / MW["Xol"] + 2.0 * qFur / MW["Fur"] + 2.0 * qHMF / MW["HMF"]) * 1000.0
    cap = s.theta_e * E_tot
    if D > 0:
        psi = (1.0 + (D / cap) ** 4) ** -0.25 if cap > 0 else 0.0
    else:
        psi = 1.0
    qP = qP_dem * psi
    qFur_a, qHMF_a = qFur * psi, qHMF * psi
    D_used = D * psi

    # ---- O2 and CO2 --------------------------------------------------------
    qO2 = max(E_tot - D_used, 0.0) / 4.0                               # mmol/g/h
    cCO2 = (cG_ox * (1.0 - YX_ox_C - suc_ox) + alloc["CO2"]
            + cE * (1.0 - s.YX_E - s.Ac_from_E) + cGly * (1.0 - s.YX_Gly)
            + cAc * (1.0 - s.YX_Ac) + cPyr * (1.0 - s.YX_Pyr)
            + qFor / MW["For"] * 1000.0 + s.f_endo * c_decay)          # mmol/g/h
    OUR = qO2 * X                                                      # mmol/L/h
    CER = cCO2 * X

    # ---- species production rates (g/L/h, before dilution) ----------------
    gc = g_per_cmol
    rX = (mu - kd_eff) * X
    rXd = kd_tox * X - s.k_lys_dead * c["Xd"]
    rGlc = -(qG_ox + qG_f) * X
    rXyl = -(qP / s.YPS2) * X
    rXol = qP * X
    rEtOH = (alloc["EtOH"] * gc("EtOH") / 1000.0 - qE) * X
    rGly = (alloc["Gly"] * gc("Gly") / 1000.0 - qGly) * X
    rAc = ((alloc["Ac"] + s.Ac_from_E * cE) * gc("Ac") / 1000.0 - qAc) * X
    rSuc = (alloc["Suc"] + suc_ox * cG_ox) * gc("Suc") / 1000.0 * X
    rPyr = (alloc["Pyr"] * gc("Pyr") / 1000.0 - qPyr) * X
    rFor = -qFor * X
    rFur = -qFur_a * X
    rHMF = -qHMF_a * X
    rFurOH = qFur_a * X * MW["FurOH"] / MW["Fur"]
    rHMFOH = qHMF_a * X * MW["HMFOH"] / MW["HMF"]

    # enzymatic release of xylose from oligomers (xylose equivalents)
    r_enz = ex.k_enz * c["Enz"] * c["Olig"] / (ex.K_olig + c["Olig"])

    # ---- protons -> base / acid (mmol/L/h) --------------------------------
    H = ((rAc + qAc * X) / MW["Ac"] * (1 - undissociated(pH, PKA["Ac"]))        # secreted
         - qAc * X / MW["Ac"] * (1 - undissociated(pH, PKA["Ac"]))             # consumed
         + rSuc / MW["Suc"] * succinate_charge(pH)
         + (rPyr + qPyr * X) / MW["Pyr"] - qPyr * X / MW["Pyr"]
         - qFor * X / MW["For"] * (1 - undissociated(pH, PKA["For"])))
    H *= 1000.0
    n_frac = 1.0 if ex.N_source == "NH4" else 0.15                        # YP: mostly amino-N
    H += n_frac * N_PER_CMOL_X * max(cX, 0.0) * X                          # NH4+ assimilation
    F_base = max(H, 0.0) * V / (ex.base_M * 1000.0) if ex.pH_control else 0.0   # L/h
    F_acid = max(-H, 0.0) * V / (ex.acid_M * 1000.0) if ex.pH_control else 0.0

    # ---- reactor flows -----------------------------------------------------
    F_feed = ex.F_feed(t) / (1.0 + math.exp(min((V - r.V_work_max) / 0.002, 50.0)))  # overflow guard
    Q = ex.airflow
    N = ex.stirrer
    if ex.do_cascade:
        dc_ = ex.do_cascade
        pO2 = c["DO"] / DOsat * 100.0
        N = min(max(dc_.get("N_min", ex.stirrer) + dc_.get("gain", 40.0) * (dc_["setpoint"] - pO2),
                    dc_.get("N_min", ex.stirrer)), dc_.get("N_max", 1200.0))
    kLa = (r.fixed_kLa if r.fixed_kLa is not None else r.kLa(N, Q, V)) / (1.0 + c["AF"] / r.K_AF)
    mol_air = Q * 60.0 / 24.06                                            # mol/h at 20 degC
    F_evap = (r.evap_L_h if r.evap_L_h is not None
              else mol_air * water_vapour_fraction(r.T_condenser) * 18.015 / 1000.0)  # L/h
    k_strip = r.k_strip_ref * (Q / V) / r.vvm_ref

    F_in = F_feed + F_base + F_acid
    dV = F_in - F_evap

    dy = np.zeros(len(STATES))
    dy[IDX["V"]] = dV
    prod = dict(X=rX, Xd=rXd, Glc=rGlc, Xyl=rXyl + r_enz, Xol=rXol, Ara=0.0, Olig=-r_enz,
                EtOH=rEtOH, Gly=rGly, Ac=rAc, Suc=rSuc, Pyr=rPyr, For=rFor, Fur=rFur,
                HMF=rHMF, FurOH=rFurOH, HMFOH=rHMFOH, Trc=0.0, Enz=-ex.k_enz_deact * c["Enz"], AF=0.0)
    for k in LIQ:
        conc = y[IDX[k]]
        feed_in = F_feed * ex.feed_composition.get(k, 0.0)
        dy[IDX[k]] = prod[k] + (feed_in - F_in * conc + F_evap * conc) / V
    dy[IDX["EtOH"]] -= k_strip * c["EtOH"]
    dy[IDX["DO"]] = kLa * (DOsat - y[IDX["DO"]]) - OUR
    target = s.Ki_rep ** 2 / (s.Ki_rep ** 2 + c["Glc"] ** 2)
    dy[IDX["Ind"]] = (target - y[IDX["Ind"]]) / s.tau_ind
    dy[IDX["cumFeed"]] = F_feed
    dy[IDX["cumBase"]] = F_base
    dy[IDX["cumAcid"]] = F_acid
    dy[IDX["cumEvap"]] = F_evap
    dy[IDX["cumO2"]] = OUR * V
    dy[IDX["cumCO2"]] = CER * V
    dy[IDX["cumStripE"]] = k_strip * c["EtOH"] * V                       # g
    dy[IDX["cumLysis"]] = ((1 - s.f_endo) * s.kd * X + s.k_lys_dead * c["Xd"]) * V   # g DW lysed
    dy[IDX["cumXol"]] = rXol * V                                          # g formed

    if not full:
        return dy
    diag = dict(T=T, pH=pH, fT=fT, fpH=fpH, f_inh=f_inh, fO2=fO, ForU=ForU, AcU=AcU,
                mu=mu, mu_model31=mu_ox, kd_eff=kd_eff, qG_ox=qG_ox, qG_overflow=qG_of,
                qG_pasteur=qG_div, qE=qE, qGly=qGly, qAc=qAc, qXol=qP, qXol_model31=qP31,
                phi_cofactor=phi, psi_electron=psi, e_supply=E_tot, e_cosubstrate=A_ext, C_cosubstrate=A_C,
                qO2=qO2, qCO2=cCO2, OUR=OUR, CER=CER,
                RQ=(CER / OUR) if OUR > 1e-9 else np.nan,
                kLa=kLa, stirrer=N, DOsat=DOsat, pO2=c["DO"] / DOsat * 100.0,
                F_feed=F_feed, F_base=F_base, F_acid=F_acid, F_evap=F_evap, H_prod=H)
    return dy, diag


# ---------------------------------------------------------------------------
# Simulation driver
# ---------------------------------------------------------------------------
@dc.dataclass
class Result:
    experiment: Experiment
    truth: pd.DataFrame          # states + diagnostics on a fine grid
    samples: pd.DataFrame        # true state at each sampling event (pre-sample)
    events: pd.DataFrame         # discrete event log

    def offgas(self, y_O2_in: float = 0.2095, y_CO2_in: float = 0.0004) -> pd.DataFrame:
        """Dry off-gas composition from OUR/CER (inert-gas balance)."""
        tr = self.truth
        mol_air = self.experiment.airflow * 60.0 / 24.06                  # mol/h
        inert = mol_air * (1 - y_O2_in - y_CO2_in)
        nO2 = mol_air * y_O2_in - tr["OUR"] * tr["V"] / 1000.0
        nCO2 = mol_air * y_CO2_in + tr["CER"] * tr["V"] / 1000.0
        tot = inert + nO2 + nCO2
        return pd.DataFrame(dict(t=tr["t"], O2_pct=100 * nO2 / tot, CO2_pct=100 * nCO2 / tot,
                                 OUR=tr["OUR"], CER=tr["CER"], RQ=tr["RQ"]))

    def carbon_balance(self) -> dict:
        """Carbon closure (C-mol). The xylitol term is split into the
        stoichiometric part and the excess implied by YPS2 > 1.013 g/g."""
        ex, tr, sm = self.experiment, self.truth, self.samples
        last = tr.iloc[-1]
        cm = lambda g, sp: g / g_per_cmol(sp)
        sp_c = [s for s in SPECIES if s != "Olig"] + ["Olig"]
        def liquid_C(row):
            tot = (row["X"] + row["Xd"]) * row["V"] / G_PER_CMOL_X
            for s_ in sp_c:
                tot += cm(row[s_] * row["V"], s_)
            return tot
        C0 = liquid_C(tr.iloc[0])
        C_feed = sum(cm(ex.feed_composition.get(s_, 0.0), s_) for s_ in sp_c) * last["cumFeed"]
        C_pulse = sum(sum(cm(p.composition.get(s_, 0.0), s_) for s_ in sp_c) * p.volume for p in ex.pulses)
        C_samples = sum(liquid_C(row) * ex.sample_volume / row["V"] for _, row in sm.iterrows())
        C_end = liquid_C(last)
        C_CO2 = last["cumCO2"] / 1000.0
        C_strip = cm(last["cumStripE"], "EtOH")
        C_lysis = last["cumLysis"] / G_PER_CMOL_X
        xol_formed = last["cumXol"]
        xyl_stoich = xol_formed / self.experiment.strain.YPS2 * MW["Xol"] / MW["Xyl"]
        C_excess = cm(xol_formed, "Xol") - cm(xol_formed / self.experiment.strain.YPS2, "Xyl")
        inputs = C0 + C_feed + C_pulse
        outputs = C_end + C_samples + C_CO2 + C_strip + C_lysis
        return dict(C_in=inputs, C_out=outputs, closure=outputs / inputs,
                    closure_without_YPS2_excess=(outputs - C_excess) / inputs,
                    xylitol_excess_Cmol=C_excess)


def simulate(ex: Experiment, dt_out: float = 1 / 60, rtol: float = 1e-6,
             atol: float = 1e-8) -> Result:
    y = np.zeros(len(STATES))
    y[IDX["V"]] = ex.V0
    y[IDX["X"]] = ex.X0
    for k, v in ex.medium.items():
        y[IDX[k]] = v
    y[IDX["DO"]] = o2_saturation_mmolL(ex.T_at(0.0))
    y[IDX["Ind"]] = ex.strain.Ki_rep ** 2 / (ex.strain.Ki_rep ** 2 + ex.medium.get("Glc", 0.0) ** 2)

    # breakpoints
    events = []
    for ts in ex.sample_times:
        events.append((ts, "sample", None))
    for p in ex.pulses:
        events.append((p.t, "pulse", p))
    for ta, ppm in ex.antifoam_additions:
        events.append((ta, "antifoam", ppm))
    for te, U in ex.enzyme_additions:
        events.append((te, "enzyme", U))
    bps = {0.0, ex.t_end}
    bps.update(t for t, _, _ in events)
    for seg in ex.feed_segments:
        bps.update([seg.t_start, seg.t_end])
    bps.update(t for t, _ in ex.T_schedule)
    bps = sorted(b for b in bps if 0.0 <= b <= ex.t_end)

    rows, sample_rows, log = [], [], []
    log.append(dict(t=0.0, event="inoculation", detail=f"X0={ex.X0} g/L, V0={ex.V0*1000:.0f} mL"))
    for seg in ex.feed_segments:
        log.append(dict(t=seg.t_start, event="feed_start", detail=f"{seg.rate*1000:.2f} mL/h, mu_exp={seg.mu_exp}"))
        log.append(dict(t=seg.t_end, event="feed_stop", detail=""))
    for ts, Ts in ex.T_schedule:
        log.append(dict(t=ts, event="temperature_setpoint", detail=f"{Ts} degC"))

    fun = lambda t, yy: rates(t, yy, ex)
    for a, b in zip(bps[:-1], bps[1:]):
        # events at a (apply before integrating the interval)
        for te, kind, payload in sorted([e for e in events if abs(e[0] - a) < 1e-9], key=lambda e: e[1] != "sample"):
            y = _apply_event(y, ex, te, kind, payload, sample_rows, log)
        if b - a <= 1e-12:
            continue
        t_eval = np.arange(math.ceil(a / dt_out - 1e-9) * dt_out, b, dt_out)
        t_eval = np.unique(np.concatenate([t_eval[(t_eval >= a) & (t_eval < b)], [b]]))
        sol = solve_ivp(fun, (a, b), y, method="LSODA", t_eval=t_eval, rtol=rtol, atol=atol,
                        first_step=1e-4, max_step=0.5)
        if not sol.success:
            raise RuntimeError(f"{ex.run_id}: integration failed in [{a}, {b}]: {sol.message}")
        for i, tt in enumerate(sol.t[:-1]):
            rows.append((tt, sol.y[:, i].copy()))
        y = sol.y[:, -1].copy()
    for te, kind, payload in [e for e in events if abs(e[0] - ex.t_end) < 1e-9]:
        y = _apply_event(y, ex, te, kind, payload, sample_rows, log)
    rows.append((ex.t_end, y.copy()))

    recs = []
    for tt, yy in rows:
        _, diag = rates(tt, yy, ex, full=True)
        rec = {"t": tt}
        rec.update({k: yy[IDX[k]] for k in STATES})
        rec.update(diag)
        recs.append(rec)
    truth = pd.DataFrame(recs).drop_duplicates("t", keep="last").reset_index(drop=True)
    truth["feed_vessel_g"] = ex.feed_vessel0 - truth["cumFeed"] * ex.feed_density * 1000.0
    for k in ["X", "Xd"] + SPECIES:
        truth[k] = truth[k].clip(lower=0.0)
    truth["X_total"] = truth["X"] + truth["Xd"]                         # what dry weight / OD see
    if (truth["V"] > ex.reactor.V_work_max - 0.005).any():
        t_full = float(truth.loc[truth["V"] > ex.reactor.V_work_max - 0.005, "t"].iloc[0])
        log.append(dict(t=t_full, event="max_volume_reached",
                        detail=f"feed blocked at V={ex.reactor.V_work_max*1000:.0f} mL"))
    return Result(ex, truth, pd.DataFrame(sample_rows), pd.DataFrame(log).sort_values("t").reset_index(drop=True))


def _apply_event(y, ex: Experiment, t, kind, payload, sample_rows, log):
    y = y.copy()
    V = y[IDX["V"]]
    if kind == "sample":
        rec = {"t": t, "sample_volume_L": ex.sample_volume}
        rec.update({k: y[IDX[k]] for k in STATES})
        sample_rows.append(rec)
        y[IDX["V"]] = V - ex.sample_volume
        log.append(dict(t=t, event="sample", detail=f"{ex.sample_volume*1000:.1f} mL removed"))
    elif kind == "pulse":
        p: Pulse = payload
        Vn = V + p.volume
        for k in LIQ:
            y[IDX[k]] = (y[IDX[k]] * V + p.composition.get(k, 0.0) * p.volume) / Vn
        y[IDX["V"]] = Vn
        log.append(dict(t=t, event=p.label, detail=f"{p.volume*1000:.1f} mL of {p.composition}"))
    elif kind == "antifoam":
        y[IDX["AF"]] += payload
        log.append(dict(t=t, event="antifoam_addition", detail=f"{payload} ppm"))
    elif kind == "enzyme":
        y[IDX["Enz"]] += payload
        log.append(dict(t=t, event="enzyme_addition", detail=f"{payload} U/L"))
    return y


# ---------------------------------------------------------------------------
# Reference recipe: BC18 (VTT, Jan 2011)
# ---------------------------------------------------------------------------
def bc18_experiment(strain: Strain | None = None) -> Experiment:
    """H1274 in YP + 20 g/L glucose, 30 degC, pH 5.5, 900 rpm, 0.7 L/min,
    fed from 7.66 h at ~2 mL/h with 95 g/L glucose + 317 g/L xylose."""
    return Experiment(
        run_id="BC18-sim", strain=strain or Strain(), reactor=Reactor(),
        V0=0.300, X0=0.888, medium={"Glc": 19.8, "Gly": 0.577}, N_source="YP",
        T_set=30.0, pH_set=5.5, stirrer=900.0, airflow=0.7,
        feed_segments=[FeedSegment(7.66, 143.0, 0.00204)],
        feed_composition={"Glc": 95.1, "Xyl": 317.0}, feed_vessel0=780.0,
        sample_times=[0.05, 4.3, 7.65, 23.4, 29.1, 45.6, 73.2, 97.7, 120.0],
        sample_volume=0.007, t_end=143.0,
        # INFERRED, not in the VTT records: an antifoam dose when the feed took over is the most
        # plausible cause of the ~60 % kLa drop between the ethanol phase and the fed phase.
        antifoam_additions=[(16.7, 200.0)],
        notes="Reference recipe reconstructed from data_for_fluxes.xlsx and online log BC18")


if __name__ == "__main__":
    res = simulate(bc18_experiment())
    tr = res.truth
    cols = ["t", "V", "X", "Glc", "Xyl", "Xol", "EtOH", "Gly", "Ac", "pO2", "OUR", "CER", "RQ", "phi_cofactor"]
    print(tr.loc[tr["t"].round(2).isin([0.0, 4.3, 7.65, 12.0, 23.4, 45.6, 73.2, 97.7, 120.0, 143.0]), cols].round(3).to_string())
    print(res.carbon_balance())
