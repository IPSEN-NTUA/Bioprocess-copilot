"""
Measurement layer: true simulator states -> what the instruments and people record.

All noise levels are module constants so they can be reported and changed in one place.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

def resample(df: pd.DataFrame, dt_h: float, cols: list) -> pd.DataFrame:
    """Values on an exact time grid (events add off-grid rows to the simulator output)."""
    t = np.arange(0.0, df["t"].iloc[-1] + 1e-9, dt_h)
    return pd.DataFrame({"t": t, **{c: np.interp(t, df["t"].values, df[c].values) for c in cols}})


# ---------------------------------------------------------------------------
# HPLC
# ---------------------------------------------------------------------------
HPLC_ANALYTES_V1 = ["Glc", "Xyl", "Ara", "Xol", "Gly", "Ac", "EtOH", "For", "Pyr"]
HPLC_ANALYTES_V2 = HPLC_ANALYTES_V1 + ["Suc", "Fur", "HMF"]          # SOP v2 from 2025-06-01
HPLC_NAMES = dict(Glc="Glucose", Xyl="Xylose", Ara="Arabinose", Xol="Xylitol", Gly="Glycerol",
                  Ac="Acetate", EtOH="Ethanol", For="Formate", Pyr="Pyruvate", Suc="Succinate",
                  Fur="Furfural", HMF="HMF")
HPLC_CV = 0.025                 # relative repeatability (1 sd)
HPLC_ABS_SD = 0.01              # g/L absolute noise in the vial
HPLC_SEQ_BIAS_SD = 0.015        # per-sequence, per-analyte calibration offset (relative, 1 sd)
HPLC_LOD = dict(Glc=0.05, Xyl=0.05, Ara=0.05, Xol=0.05, Gly=0.05, EtOH=0.05,
                Ac=0.03, For=0.03, Pyr=0.02, Suc=0.03, Fur=0.01, HMF=0.01)   # g/L in the vial
HPLC_MAX_VIAL = 10.0            # g/L, top of calibration -> sample is diluted
DILUTIONS = [1, 2, 5, 10, 20, 50]
# response factors (area units per g/L in the vial) for the raw exports
RF = dict(Glc=101800, Xyl=98400, Ara=97100, Xol=104900, Gly=88300, Ac=41200, EtOH=46100,
          For=35600, Pyr=127000, Suc=52800, Fur=611000, HMF=702000)
RT = dict(Pyr=8.1, Glc=8.9, Xyl=9.5, Ara=10.0, Xol=10.6, Suc=11.3, Gly=13.0, For=13.6, Ac=14.7,
          EtOH=21.4, HMF=29.0, Fur=43.6)                                  # retention time (min)


def choose_dilution(true: dict) -> int:
    top = max(true.values())
    for d in DILUTIONS:
        if top / d <= HPLC_MAX_VIAL:
            return d
    return DILUTIONS[-1]


def hplc_inject(true: dict, analytes: list, seq_bias: dict, rng, extra_bias: dict | None = None) -> dict:
    """One injection. Returns per analyte dict(area, vial, conc, flag) and the dilution."""
    d = choose_dilution({a: true.get(a, 0.0) for a in analytes})
    out = {}
    for a in analytes:
        vial_true = max(true.get(a, 0.0), 0.0) / d
        bias = 1.0 + seq_bias.get(a, 0.0) + (extra_bias or {}).get(a, 0.0)
        vial = vial_true * bias * (1.0 + rng.normal(0, HPLC_CV)) + rng.normal(0, HPLC_ABS_SD)
        flag = ""
        if vial_true <= 0.0:          # compound absent: no peak (random draws kept so the stream is unchanged)
            vial = 0.0
        if vial < HPLC_LOD[a]:
            flag = "<LOD"
        area = max(vial, 0.0) * RF[a] * (1.0 + rng.normal(0, 0.002))
        out[a] = dict(area=round(area), vial=vial, conc=vial * d, flag=flag)
    return dict(dilution=d, results=out)


def sequence_bias(analytes: list, rng) -> dict:
    return {a: float(rng.normal(0, HPLC_SEQ_BIAS_SD)) for a in analytes}


# ---------------------------------------------------------------------------
# Biomass
# ---------------------------------------------------------------------------
DW_TUBE_CV = 0.04          # dry-weight tube-to-tube (1 sd, relative)
DW_TUBES = 3
DW_VOLUME_ML = 2.0
OD_PER_DW = 2.5            # OD600 units per g DW/L for this host
OD_CV = 0.05


def dry_weight(X: float, rng) -> dict:
    tubes = []
    for _ in range(DW_TUBES):
        tare = rng.uniform(1.0050, 1.0250)
        cells_g = max(X, 0) * (1 + rng.normal(0, DW_TUBE_CV)) * DW_VOLUME_ML / 1000.0
        tubes.append(dict(tare=round(tare, 4), gross=round(tare + cells_g, 4)))
    vals = [(t["gross"] - t["tare"]) * 1000.0 / DW_VOLUME_ML for t in tubes]
    return dict(tubes=tubes, values=vals, mean=float(np.mean(vals)),
                sem=float(np.std(vals, ddof=1) / math.sqrt(len(vals))))


def od600(X: float, rng) -> dict:
    od = max(X, 0) * OD_PER_DW * (1 + rng.normal(0, OD_CV))
    dil = 1 if od < 0.8 else (10 if od < 8 else (20 if od < 16 else 50))
    reading = od / dil
    return dict(dilution=dil, reading=round(reading, 3), od=round(reading * dil, 2))


# ---------------------------------------------------------------------------
# Online sensors (reactors)
# ---------------------------------------------------------------------------
ONLINE_DT_MIN = 10
T_SD = 0.02                # degC
PH_SD = 0.008
PH_OFFSET_SD = 0.03        # per-run probe calibration offset
PO2_SD = 0.4               # % air saturation
PO2_DRIFT_SD = 0.02        # %/h per-run calibration drift (normal runs)
STIRRER_SD = 0.6           # rpm
PUMP_ML_PER_COUNT = dict(acid=0.05, base=0.05, feed=0.10)
BALANCE_SD = 0.05          # g, feed balance


def online_signals(truth: pd.DataFrame, rng, pH_set: float, T_set: float, drift_pct_per_h: float | None = None,
                   gaps: list | None = None, feed_log_factor: float = 1.0) -> pd.DataFrame:
    """10-min records from the 1-min truth grid. feed_log_factor: logged/actual feed (pump
    calibration faults). gaps: [(t0, t1)] hours with no records (logger off)."""
    tr = resample(truth, ONLINE_DT_MIN / 60.0, ["stirrer", "T", "pO2", "cumAcid", "cumBase", "cumFeed",
                                                 "feed_vessel_g"])
    off = rng.normal(0, PH_OFFSET_SD)
    drift = drift_pct_per_h if drift_pct_per_h is not None else rng.normal(0, PO2_DRIFT_SD)
    n = len(tr)
    rec = pd.DataFrame(dict(t=tr["t"].values))
    rec["stirrer"] = np.round(tr["stirrer"].values + rng.normal(0, STIRRER_SD, n))
    rec["temp"] = np.round(tr["T"].values + rng.normal(0, T_SD, n), 2)
    rec["pH"] = np.round(pH_set + off + rng.normal(0, PH_SD, n), 3)
    po2 = tr["pO2"].values * (1 + drift * tr["t"].values / 100.0) + rng.normal(0, PO2_SD, n)
    rec["pO2"] = np.round(np.clip(po2, 0, 110), 1)
    rec["acid_counts"] = np.floor(tr["cumAcid"].values * 1000 / PUMP_ML_PER_COUNT["acid"]).astype(int)
    rec["base_counts"] = np.floor(tr["cumBase"].values * 1000 / PUMP_ML_PER_COUNT["base"]).astype(int)
    rec["feed_counts"] = np.floor(tr["cumFeed"].values * 1000 * feed_log_factor / PUMP_ML_PER_COUNT["feed"]).astype(int)
    rec["feed_balance_g"] = np.round(tr["feed_vessel_g"].values + rng.normal(0, BALANCE_SD, n), 2)
    rec["airflow"] = np.nan
    for g0, g1 in (gaps or []):
        rec = rec[(rec["t"] < g0) | (rec["t"] >= g1)]
    return rec.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Off-gas analyser
# ---------------------------------------------------------------------------
OFFGAS_CYCLE_MIN = 30
OFFGAS_O2_SD = 0.010       # % abs
OFFGAS_CO2_SD = 0.006      # % abs
OFFGAS_O2_DRIFT = -0.004   # % abs per day (span drift)


def offgas_cycles(offgas: pd.DataFrame, rng, t_shift: float) -> pd.DataFrame:
    """Analyser cycles every 30 min; the analyser clock starts t_shift h before inoculation."""
    og = resample(offgas, OFFGAS_CYCLE_MIN / 60.0, ["O2_pct", "CO2_pct"])
    n = len(og)
    days = og["t"].values / 24
    out = pd.DataFrame(dict(t=og["t"].values))
    out["O2"] = np.round(og["O2_pct"].values + OFFGAS_O2_DRIFT * days + rng.normal(0, OFFGAS_O2_SD, n), 4)
    out["CO2"] = np.round(np.clip(og["CO2_pct"].values + rng.normal(0, OFFGAS_CO2_SD, n), 0, None), 4)
    out["N2"] = np.round(100 - out["O2"] - out["CO2"] - 0.93 + rng.normal(0, 0.02, n), 3)
    out["Ar"] = np.round(0.93 + rng.normal(0, 0.004, n), 4)
    out["rel_time_s"] = np.round((out["t"] + t_shift) * 3600, 1)
    return out
