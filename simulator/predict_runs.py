"""Blind prediction of VTT runs BC12, BC15, BC19 with the current simulator (no parameter changes).

Recipes are rebuilt from data_for_fluxes.xlsx (initial medium, feed composition, feed start and
piecewise feed rates, sampling times) and the online logs (stirrer, pH, inoculation offset).
Usage (from simulator/):  python predict_runs.py   ->  prediction_BC12_BC15_BC19.png
"""
from __future__ import annotations

import dataclasses as dc
from pathlib import Path

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd

from simulator import Experiment, FeedSegment, Reactor, Strain, simulate

HERE = Path(__file__).resolve().parent
# partner data are not in the repository; see data/README.md
VTT = Path(os.environ.get("PARTNER_DATA_DIR", HERE.parent / "data" / "partner"))
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
SIM, SIM_AF, DATA = "#2a78d6", "#1baf7a", "#eb6834"

RUNS = {
    # inoculation offset of the online log (h) and pH setpoint read from the logs;
    # BC12's log started ~18.6 h before inoculation and has two empty tab columns.
    "BC12": dict(log_offset=18.6, pH=5.4, log_cols=["Date", "Time", "Stirrer", "Temp", "pH", "Acid_Pump",
                                                    "_1", "Base_Pump", "_2", "pO2", "Feed_Pump"],
                 note="H1274 multi-copy, YP, 63 g/L xylose in batch medium"),
    "BC15": dict(log_offset=0.0, pH=5.45, log_cols=None, note="H1274 multi-copy, YP"),
    "BC19": dict(log_offset=0.0, pH=5.45, log_cols=None,
                 note="H1274 multi-copy, YP + corn steep solids (CSS effect not modelled)"),
}


def read_flux_sheet() -> list[tuple]:
    ws = openpyxl.load_workbook(VTT / "data_for_fluxes.xlsx", data_only=True).active
    rows = list(ws.iter_rows(values_only=True))
    col = {}
    for i, name in enumerate(rows[1]):
        if name and name not in col:
            col[name] = i
    return rows, col


def offline(run: str):
    rows, col = read_flux_sheet()
    num = lambda v: float(v) if isinstance(v, (int, float)) else np.nan
    samples, feed = [], None
    for r in rows[2:]:
        nm = r[col["Sample name"]]
        if not isinstance(nm, str) or not nm.startswith(run + "_"):
            continue
        tag = nm[len(run) + 1:]
        if tag == "Feed":
            feed = dict(Glc=num(r[col["Glucose"]]), Xyl=num(r[col["Xylose"]]))
        elif tag.isdigit():
            samples.append(dict(t=num(r[col["Time (h)"]]), V=num(r[col["Reactor volume (mL)"]]),
                                X=num(r[col["DW (g/L)"]]), Glc=num(r[col["Glucose"]]),
                                Xyl=num(r[col["Xylose"]]), Xol=num(r[col["Xylitol"]]),
                                Gly=num(r[col["Glycerol"]]), Ac=num(r[col["Acetate"]]),
                                EtOH=num(r[col["Ethanol"]]), feed_added=num(r[col["Feed added (mL)"]]),
                                feed_rate=num(r[col["Feed rate (mL/h)"]])))
    return pd.DataFrame(samples), feed


def online(run: str, cfg: dict) -> pd.DataFrame:
    lines = (VTT / run).read_text(encoding="latin-1").replace("\r", "").split("\n")
    i0 = next(i for i, l in enumerate(lines) if l.startswith("Date"))
    hdr = cfg["log_cols"] or [h.strip() for h in lines[i0].split("\t")]
    recs = [dict(zip(hdr, [x.strip() for x in l.split("\t")])) for l in lines[i0 + 1:] if l.strip()]
    df = pd.DataFrame(recs)
    ts = pd.to_datetime(df["Date"] + " " + df["Time"], format="%d.%m.%Y %H:%M:%S")
    df["t"] = (ts - ts.iloc[0]).dt.total_seconds() / 3600 + 5 / 60 - cfg["log_offset"]
    df["pO2"] = pd.to_numeric(df["pO2"], errors="coerce")
    return df.loc[df["t"] >= 0, ["t", "pO2"]]


def build_experiment(run: str, cfg: dict, antifoam: bool = False) -> Experiment:
    sm, feed = offline(run)
    first = sm.iloc[0]
    fed = sm[sm["feed_added"] > 0]
    t_start = fed.iloc[0]["t"] - fed.iloc[0]["feed_added"] / fed.iloc[0]["feed_rate"]
    segs, t_prev = [], t_start
    for _, r in sm[sm["t"] > t_start].iterrows():
        rate = r["feed_rate"] if r["feed_rate"] > 0 else (segs[-1].rate * 1000 if segs else fed.iloc[0]["feed_rate"])
        segs.append(FeedSegment(t_prev, r["t"], rate / 1000.0))
        t_prev = r["t"]
    medium = {k: first[k] for k in ["Glc", "Xyl", "Gly", "EtOH", "Ac"] if first[k] > 0}
    return Experiment(
        run_id=f"{run}-pred", strain=Strain(), reactor=Reactor(),
        V0=first["V"] / 1000.0, X0=first["X"], medium=medium, N_source="YP",
        T_set=30.0, pH_set=cfg["pH"], stirrer=900.0, airflow=0.7,
        feed_segments=segs, feed_composition=feed,
        sample_times=list(sm["t"].iloc[1:-1]), sample_volume=0.007, t_end=float(sm["t"].iloc[-1]),
        antifoam_additions=[(t_start + 9.0, 200.0)] if antifoam else [],
        notes=cfg["note"])


def style(ax, title, unit):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold")
    ax.set_ylabel(unit, color=INK2, fontsize=8)
    ax.tick_params(colors=INK2, labelsize=7.5, length=0)
    ax.grid(True, color=GRID, linewidth=0.6)
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    panels = [("X", "Biomass", "g DW/L"), ("Glc", "Glucose", "g/L"), ("Xyl", "Xylose", "g/L"),
              ("Xol", "Xylitol", "g/L"), ("EtOH", "Ethanol", "g/L"), ("pO2", "Dissolved O2", "% air sat.")]
    fig, axes = plt.subplots(3, 6, figsize=(19, 10), facecolor=SURFACE)
    summary = []
    for row, (run, cfg) in enumerate(RUNS.items()):
        ex = build_experiment(run, cfg)
        tr = simulate(ex).truth
        tr_af = simulate(build_experiment(run, cfg, antifoam=True)).truth
        sm, _ = offline(run)
        on = online(run, cfg)
        rm = {}
        for col, (k, title, unit) in enumerate(panels):
            ax = axes[row, col]
            style(ax, f"{run} · {title}" if col == 0 else title, unit)
            ax.plot(tr["t"], tr[k], color=SIM, lw=2, label="predicted")
            if k == "pO2":
                ax.plot(tr_af["t"], tr_af[k], color=SIM_AF, lw=1.6, ls="--",
                        label="predicted, BC18-like antifoam")
                ax.plot(on["t"], on["pO2"], color=DATA, lw=1.0, label="measured")
                o = on[on["t"] > 20]
                rm[k] = np.sqrt(np.mean((np.interp(o.t, tr.t, tr[k]) - o.pO2) ** 2))
            else:
                d = sm[["t", k]].dropna()
                ax.plot(d["t"], d[k], "o", color=DATA, ms=5.5, mec=SURFACE, mew=1.3, label="measured")
                rm[k] = np.sqrt(np.mean((np.interp(d.t, tr.t, tr[k]) - d[k]) ** 2))
            ax.set_xlim(0, tr["t"].max() + 2)
            ax.set_ylim(bottom=0)
            if row == 2:
                ax.set_xlabel("time (h)", color=INK2, fontsize=8)
        e = tr.iloc[-1]
        d_last = sm.iloc[-1]
        summary.append(dict(run=run, t_end=round(e.t, 1),
                            Xol_pred=round(e.Xol, 1), Xol_meas=d_last.Xol,
                            Xyl_pred=round(e.Xyl, 1), Xyl_meas=d_last.Xyl,
                            X_pred=round(e.X, 1), X_meas=d_last.X,
                            **{f"RMSE_{k}": round(v, 2) for k, v in rm.items()}))
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor=INK2)
    axes[0, 5].legend(frameon=False, fontsize=7, loc="upper right", labelcolor=INK2)
    fig.suptitle("Blind prediction of BC12, BC15, BC19: simulator tuned on BC18 only, no parameter changes",
                 x=0.01, ha="left", color=INK, fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(HERE / "prediction_BC12_BC15_BC19.png", dpi=140, facecolor=SURFACE)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
