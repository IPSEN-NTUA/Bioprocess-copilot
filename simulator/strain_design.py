"""
Strain-design layer: genetic modifications -> enzyme activities -> NRA (log-linear) flux
response of the ground-truth kinetic model -> production-phase parameters of the macro model.

Ground truth
------------
Kinetic model 978,2 (index 922) of the A_oxidative ORACLE ensemble (939 models), stored in
strain_design/ground_truth_978_2.npz: flux / concentration control coefficients (283 x 283,
190 x 283) and the reference fluxes and concentrations of reference state 978. Chosen because
it is among the most representative models of the ensemble (weighted distance of its control
coefficients on the nine boundary fluxes to the ensemble median: 0.35 vs ensemble median 0.76),
all top-10 effects per boundary flux have the ensemble sign, and its xylitol control pattern
correlates 0.95 with the ensemble median (XR control 0.37, all NADPH handles present).

Pipeline
--------
1. Each modification produces a realised activity fold on one enzyme:
   * overexpression / downregulation / targeted point mutation:
         fold = F_nominal * u,  u ~ U(0.5, 1.5)
     with F_nominal in (1, 2] (up) or [0.5, 1) (down), and u truncated so that the realised
     fold stays in [1, 2] for an overexpression and [0.5, 1] for a downregulation: the
     direction is never reversed and the expression change never exceeds 2-fold.
     The draw is per transformant, so sister clones differ.
   * deletion: fold = DELETION_FOLD = 0.10 (deterministic, 10 % residual activity).
   * random mutation (mutagenesis / evolution, or an untargeted point mutation in a gene):
     fold drawn from a distribution of fitness effects skewed to loss of function; the
     enzyme is drawn at random from the design handles if no target is given.
2. Activities enter NRA in log space: dln e = ln(fold). Modifications along a lineage add.
3. Flux response (linear, the NRA formulation): dln J = C^J . dln e, using the FCC of the
   ground-truth model. Large responses are compressed smoothly,
         dln J_eff = L tanh(dln J / L),  L = ln 3,
   which is ~1:1 below |dln J| ~ 0.3 (the linearity limit used in the NRA notebook) and
   saturates at a 3-fold change of any flux.
4. Mapping to the macro model (production phase only; each multiplier m acts as m**Ind,
   Ind = glucose-derepression state, so the batch phase is untouched):
       XYLTt (xylitol export)     -> mumax2            x exp(dlnJ)
       biomass lump               -> prod_mumax1_mult  = exp(dlnJ_growth)
       biomass / GLCt1            -> prod_YXS1_mult    = exp(dlnJ_growth - dlnJ_GLC)
       ETOHt, GLYCt, SUCCt, ACtr  -> prod_div_mult + prod_split: the share of O2-limited
                                     glucose diverted to fermentation and its product split,
                                     solved so that every by-product flux changes by exactly
                                     exp(dlnJ) while carbon and electrons still balance
       SUCCt                      -> prod_suc_mult (succinate leak of the oxidative channel)
       O2t, CO2t                  -> not mapped (the macro model computes them from balances);
                                     reported for consistency checks.
   Traits outside the kinetic network (heterologous FDH, inhibitor tolerance, temperature
   optimum) stay macro-level and are set directly on the Strain.
"""
from __future__ import annotations

import dataclasses as dc
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from simulator import GAMMA, G_PER_CMOL_X, Strain

HERE = Path(__file__).resolve().parent
GROUND_TRUTH = HERE / "strain_design" / "ground_truth_978_2.npz"

L_COMPRESS = math.log(3.0)
DELETION_FOLD = 0.10            # deletion = 10 % residual activity (log-linear NRA needs a finite fold)
U_RANGE = (0.5, 1.5)
MAX_EXPRESSION_FOLD = 2.0     # an expression change never exceeds 2x up or 0.5x down

# reactions that are not genetic design handles (transport, biomass lump, maintenance)
TRANSPORT = set("""2DDA7Ptm 2OBUTtm 2OXOADPTm 2PHETOHt 2PHETOHtm 3C3HMPtm 3MOPtm ACALDtm ACtr AKGCITtm AKGMAL
AKGMALtm AKGt ALAtmi AMPtrm ASPGLU2m ATPtm ATPtrm CA2t3ec CITtam CITtbm CITtcm CMPt6 CO2t CO2tgm CO2tm
COAtrm CTPtrm CU2TRA CYSt2r Clt DHPTtm D_LACt2m D_LACtm E4Ptm ETOHt ETOHtm FE2t FRDPtm GCALDt GLCt1 GLUt7m
GLYCt GLYtm H2Ot H2Otce H2Oter H2Otm H2Otrm Ht Htce Htg Htgm Htm Htr Htrm Htv ID3ACALDt ID3ACALDtm ILEtmi
IND3ETHt IND3ETHtm Kt2r MALt MALtm MGt2ec MNt2 NADHtrm NADPHtru NADPtru NADtrm NAt3_1 NH4t O2t O2ter O2tm
O2trm OAAt2m ORNt3m PACALDt PACALDtm PIt2m PIt2r PItrm PPItm PPItrm PYRt2m SER__Ltrm SO4ti SUCCt SUCCtm
SUCFUMtm VALt5m XYLTt XYLt r2073_1""".split())
NON_HANDLES = TRANSPORT | {"LMPD_s_0450_c_1_256", "ATPM", "ATPM_g", "ATPM_v"}

# boundary fluxes that drive the macro model
TARGETS = dict(xylitol="XYLTt", xylose="XYLt", growth="LMPD_s_0450_c_1_256", glucose="GLCt1",
               O2="O2t", CO2="CO2t", ethanol="ETOHt", glycerol="GLYCt", succinate="SUCCt",
               acetate="ACtr")

# S. cerevisiae gene names for the main handles (approximate one-gene-one-reaction mapping,
# used for realistic strain names and lab notes; the reaction ID is what the model uses)
GENE_TO_REACTION = {
    "XYL1": "XYLR", "ZWF1": "G6PDH2r", "SOL3": "PGL", "GND1": "GND", "PGI1": "PGI",
    "GDH1": "GLUDy", "PDC1": "PYRDC", "ADH1": "ALCD2ir", "GPD1": "G3PD1", "GPP1": "G3PT",
    "ALD6": "ALDD2y", "HXK2": "HEX1", "TPI1": "TPI", "SDH1": "SUCD2_u6m", "LSC1": "SUCOASm",
    "PYC1": "PC", "AMD1": "AMPDA", "IDH1": "ICDHxm", "COX4": "FERCOXOXI", "TKL1": "TKT1",
    "TAL1": "TALA", "RPE1": "RPE",
}

# distribution of fitness (activity) effects for random / untargeted mutations
#   (probability, low fold, high fold) - drawn log-uniformly inside each class
DFE = [(0.45, 0.97, 1.03),    # neutral / silent
       (0.30, 0.50, 0.97),    # mild loss of function
       (0.15, 0.05, 0.50),    # strong loss of function
       (0.05, 0.005, 0.05),   # near null
       (0.05, 1.03, 2.00)]    # gain of function (rare)


# ---------------------------------------------------------------------------
class GroundTruth:
    def __init__(self, path: Path = GROUND_TRUTH):
        z = np.load(path, allow_pickle=True)
        self.fcc = z["fcc"]
        self.ccc = z["ccc"]
        self.reactions = [str(r) for r in z["reactions"]]
        self.mets = [str(m) for m in z["mets"]]
        self.ln_ref_flux = z["ln_ref_flux"]
        self.model_key = str(z["model_key"])
        self.ri = {r: i for i, r in enumerate(self.reactions)}
        self.mi = {m: i for i, m in enumerate(self.mets)}
        self.handles = [r for r in self.reactions if r not in NON_HANDLES]

    def resolve(self, target: str) -> str:
        r = GENE_TO_REACTION.get(target.upper(), target)
        if r not in self.ri:
            raise KeyError(f"unknown enzyme / gene '{target}'")
        if r in NON_HANDLES:
            raise ValueError(f"'{r}' is a transport/lumped reaction, not a design handle")
        return r

    def response(self, dlne: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Linear NRA response (dln J, dln x) to a vector of enzyme log-fold changes."""
        return self.fcc @ dlne, self.ccc @ dlne

    def ref_ratio(self, a: str, b: str) -> float:
        return math.exp(self.ln_ref_flux[self.ri[a]] - self.ln_ref_flux[self.ri[b]])


GT = None


def ground_truth() -> GroundTruth:
    global GT
    if GT is None:
        GT = GroundTruth()
    return GT


# ---------------------------------------------------------------------------
@dc.dataclass
class Modification:
    target: str | None                 # gene or reaction ID; None for random mutation (random enzyme)
    kind: str                          # overexpression | downregulation | deletion |
                                       # point_mutation (targeted, needs fold) | random_mutation
    fold: float | None = None          # nominal fold (overexpression > 1, downregulation < 1)
    label: str = ""                    # e.g. "pTEF1-XYL1", "3x XYL1 integrated", "GPD1 K->A"

    # filled when realised
    reaction: str | None = None
    u: float | None = None
    realized_fold: float | None = None
    dfe_class: str | None = None


def _truncated_u(F: float, rng: np.random.Generator) -> float:
    lo, hi = U_RANGE
    cap = MAX_EXPRESSION_FOLD
    if F > 1:          # overexpression: 1 <= F*u <= cap
        lo, hi = max(lo, 1.0 / F), min(hi, cap / F)
    elif F < 1:        # downregulation: 1/cap <= F*u <= 1
        lo, hi = max(lo, 1.0 / (cap * F)), min(hi, 1.0 / F)
    return float(rng.uniform(lo, hi))


def realise(mod: Modification, rng: np.random.Generator) -> Modification:
    gt = ground_truth()
    m = dc.replace(mod)
    if m.kind in ("overexpression", "downregulation", "point_mutation"):
        if m.fold is None:
            raise ValueError(f"{m.kind} needs a nominal fold")
        if m.kind == "overexpression" and not 1 < m.fold <= MAX_EXPRESSION_FOLD:
            raise ValueError(f"overexpression needs 1 < fold <= {MAX_EXPRESSION_FOLD}")
        if m.kind == "downregulation" and not 1 / MAX_EXPRESSION_FOLD <= m.fold < 1:
            raise ValueError(f"downregulation needs {1 / MAX_EXPRESSION_FOLD} <= fold < 1")
        if m.kind == "point_mutation" and not 1 / MAX_EXPRESSION_FOLD <= m.fold <= MAX_EXPRESSION_FOLD:
            raise ValueError(f"targeted point mutation needs fold in [{1 / MAX_EXPRESSION_FOLD}, {MAX_EXPRESSION_FOLD}]")
        m.reaction = gt.resolve(m.target)
        m.u = _truncated_u(m.fold, rng)
        m.realized_fold = m.fold * m.u
    elif m.kind == "deletion":
        m.reaction = gt.resolve(m.target)
        m.realized_fold = DELETION_FOLD
    elif m.kind == "random_mutation":
        m.reaction = gt.resolve(m.target) if m.target else str(rng.choice(gt.handles))
        p = np.array([d[0] for d in DFE])
        k = int(rng.choice(len(DFE), p=p / p.sum()))
        lo, hi = DFE[k][1], DFE[k][2]
        m.realized_fold = float(math.exp(rng.uniform(math.log(lo), math.log(hi))))
        m.dfe_class = ["neutral", "mild_loss", "strong_loss", "near_null", "gain"][k]
    else:
        raise ValueError(f"unknown modification kind '{m.kind}'")
    return m


# ---------------------------------------------------------------------------
def compress(x: np.ndarray | float, L: float = L_COMPRESS):
    return L * np.tanh(np.asarray(x) / L)


def _production_split(base: dict, dE: float, dG: float, dS: float, dA: float) -> tuple[dict, float]:
    """Diversion multiplier m and product split of fermentative glucose such that ethanol,
    glycerol, succinate and acetate fluxes change by exactly exp(d) while carbon and electrons
    balance (biomass yield of the channel and its electron slack kept)."""
    gX = GAMMA["X"]
    slack = 4.0 - sum(base[k] * (gX if k == "X" else GAMMA[k]) for k in base if k != "CO2")
    f = dict(EtOH=base["EtOH"] * math.exp(dE), Gly=base["Gly"] * math.exp(dG),
             Suc=base["Suc"] * math.exp(dS), Ac=base["Ac"] * math.exp(dA), Pyr=base["Pyr"])
    e_products = sum(f[k] * GAMMA[k] for k in f)
    m = e_products / (4.0 - slack - base["X"] * gX)
    split = {k: v / m for k, v in f.items()}
    split["X"] = base["X"]
    split["CO2"] = 1.0 - sum(split.values())
    if split["CO2"] < 0:
        raise ValueError("by-product changes exceed the carbon available in the fermentative channel")
    return split, m


@dc.dataclass
class StrainRecord:
    name: str
    parent: str | None
    modifications: list                 # realised Modification objects (whole lineage)
    dlne: np.ndarray                    # cumulative enzyme log-fold changes (283)
    strain: Strain                      # macro-model strain
    dlnJ: dict                          # compressed boundary-flux log changes
    dlnJ_linear: dict                   # uncompressed (pure NRA) values

    def summary(self) -> pd.DataFrame:
        rows = [dict(flux=k, reaction=TARGETS[k], linear_pct=100 * (math.exp(self.dlnJ_linear[k]) - 1),
                     applied_pct=100 * (math.exp(self.dlnJ[k]) - 1)) for k in TARGETS]
        return pd.DataFrame(rows).set_index("flux").round(1)

    def genotype(self) -> str:
        parts = []
        for m in self.modifications:
            tag = m.label or f"{m.target or m.reaction}"
            if m.kind == "deletion":
                parts.append(f"{tag}Δ")
            elif m.kind == "random_mutation":
                parts.append(f"{tag}* ({m.dfe_class})")
            else:
                parts.append(f"{tag} ({m.realized_fold:.2f}x)")
        return ", ".join(parts) if parts else "parent"


def parent_record(base: Strain | None = None, name: str = "H1274") -> StrainRecord:
    gt = ground_truth()
    base = base or Strain()
    z = {k: 0.0 for k in TARGETS}
    return StrainRecord(name, None, [], np.zeros(len(gt.reactions)), base.copy(name=name), z, dict(z))


def engineer(parent: StrainRecord, mods: list[Modification], name: str,
             rng: np.random.Generator, root: Strain | None = None) -> StrainRecord:
    """Apply modifications to a parent strain. The macro strain is always rebuilt from the
    root (unengineered) parameters with the cumulative dln e of the whole lineage, so the
    response stays linear in dln e and compression is applied once."""
    gt = ground_truth()
    root = root or Strain()
    realised = [realise(m, rng) for m in mods]
    dlne = parent.dlne.copy()
    for m in realised:
        dlne[gt.ri[m.reaction]] += math.log(m.realized_fold)
    dJ, _ = gt.response(dlne)
    lin = {k: float(dJ[gt.ri[r]]) for k, r in TARGETS.items()}
    eff = {k: float(compress(v)) for k, v in lin.items()}

    split, m_div = _production_split(root.prod_split or root.ferm_split, eff["ethanol"],
                                     eff["glycerol"], eff["succinate"], eff["acetate"])
    strain = root.copy(
        name=name,
        mumax2=root.mumax2 * math.exp(eff["xylitol"]),
        prod_mumax1_mult=math.exp(eff["growth"]),
        prod_YXS1_mult=math.exp(eff["growth"] - eff["glucose"]),
        prod_suc_mult=math.exp(eff["succinate"]),
        prod_div_mult=m_div,
        prod_split=split,
        fdh=parent.strain.fdh,
        design_note=", ".join(m.label or m.reaction for m in parent.modifications + realised))
    return StrainRecord(name, parent.name, parent.modifications + realised, dlne, strain, eff, lin)


def design_effect(targets: dict[str, float]) -> pd.DataFrame:
    """Deterministic (u = 1) NRA effect of a set of enzyme folds, e.g. {'XYLR': 2, 'PGI': 0.5}."""
    gt = ground_truth()
    dlne = np.zeros(len(gt.reactions))
    for t, f in targets.items():
        dlne[gt.ri[gt.resolve(t)]] += math.log(f)
    dJ, _ = gt.response(dlne)
    return pd.DataFrame({k: [100 * (math.exp(dJ[gt.ri[r]]) - 1),
                             100 * (math.exp(float(compress(dJ[gt.ri[r]]))) - 1)]
                         for k, r in TARGETS.items()}, index=["linear_pct", "applied_pct"]).T.round(1)


if __name__ == "__main__":
    gt = ground_truth()
    print(f"ground truth: model {gt.model_key}, {len(gt.reactions)} reactions, {len(gt.handles)} design handles")
    print("\nnotebook design GND+ PGI- XYLR+ (2x / 0.5x / 2x, u = 1):")
    print(design_effect({"GND": 2, "PGI": 0.5, "XYLR": 2}).to_string())
