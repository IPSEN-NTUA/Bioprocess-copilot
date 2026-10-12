"""Compare the extended simulator against VTT run BC18 and plot a few what-if scenarios.

Usage (from the simulator/ folder):  python validate_bc18.py
Reads data_for_fluxes.xlsx and BC18 from the partner data folder (data/partner, or $PARTNER_DATA_DIR),
writes PNGs next to this file.
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

from simulator import (FeedSegment, Strain, bc18_experiment, simulate)

HERE = Path(__file__).resolve().parent
# partner data are not in the repository; see data/README.md
VTT = Path(os.environ.get("PARTNER_DATA_DIR", HERE.parent / "data" / "partner"))

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
SIM, DATA = "#2a78d6", "#eb6834"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def load_bc18_offline() -> pd.DataFrame:
    ws = openpyxl.load_workbook(VTT / "data_for_fluxes.xlsx", data_only=True).active
    rows = list(ws.iter_rows(values_only=True))
    hdr = rows[1]
    col = {}
    for i, name in enumerate(hdr):          # first occurrence wins (flux columns reuse names)
        if name and name not in col:
            col[name] = i
    out = []
    for r in rows[2:]:
        name = r[col["Sample name"]]
        if isinstance(name, str) and name.startswith("BC18_") and name[5:].isdigit():
            out.append(dict(t=r[col["Time (h)"]], X=r[col["DW (g/L)"]], Glc=r[col["Glucose"]],
                            Xyl=r[col["Xylose"]], Xol=r[col["Xylitol"]], Gly=r[col["Glycerol"]],
                            Ac=r[col["Acetate"]], EtOH=r[col["Ethanol"]], V=r[col["Reactor volume (mL)"]]))
    return pd.DataFrame(out).apply(pd.to_numeric, errors="coerce")


def load_bc18_online() -> pd.DataFrame:
    lines = (VTT / "BC18").read_text(encoding="latin-1").splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith("Date"))
    hdr = [h.strip() for h in lines[i0].split("\t")]
    recs = []
    for l in lines[i0 + 1:]:
        p = [x.strip() for x in l.split("\t")]
        if len(p) < 5 or not p[0]:
            continue
        recs.append(dict(zip(hdr, p)))
    df = pd.DataFrame(recs)
    ts = pd.to_datetime(df["Date"] + " " + df["Time"], format="%d.%m.%Y %H:%M:%S")
    df["t"] = (ts - ts.iloc[0]).dt.total_seconds() / 3600 + 5 / 60
    df["pO2"] = pd.to_numeric(df["pO2"], errors="coerce")
    return df[["t", "pO2"]]


def style(ax, title, unit):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold")
    ax.set_ylabel(unit, color=INK2, fontsize=8)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(True, color=GRID, linewidth=0.6)
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    res = simulate(bc18_experiment())
    tr = res.truth
    off = load_bc18_offline()
    on = load_bc18_online()

    panels = [("X", "Biomass", "g DW/L"), ("Glc", "Glucose", "g/L"), ("Xyl", "Xylose", "g/L"),
              ("Xol", "Xylitol", "g/L"), ("EtOH", "Ethanol", "g/L"), ("Gly", "Glycerol", "g/L"),
              ("pO2", "Dissolved O2 (online)", "% air sat."), ("V", "Reactor volume", "mL")]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.6), facecolor=SURFACE)
    for ax, (k, title, unit) in zip(axes.flat, panels):
        style(ax, title, unit)
        ysim = tr[k] * (1000 if k == "V" else 1)
        ax.plot(tr["t"], ysim, color=SIM, lw=2, label="simulated")
        if k == "pO2":
            ax.plot(on["t"], on["pO2"], color=DATA, lw=1.2, alpha=0.9, label="measured")
        else:
            ax.plot(off["t"], off[k], "o", color=DATA, ms=6, mec=SURFACE, mew=1.5, label="measured")
        ax.set_xlim(0, 145)
        ax.set_ylim(bottom=0)
    for ax in axes[1]:
        ax.set_xlabel("time (h)", color=INK2, fontsize=8)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK2)
    fig.suptitle("BC18 (H1274, YP, 30 °C, pH 5.5): extended model 31, parameters tuned within ±20 % of the refined fit",
                 x=0.01, ha="left", color=INK, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(HERE / "validation_bc18.png", dpi=150, facecolor=SURFACE)

    # ---- what-if scenarios ------------------------------------------------
    base = bc18_experiment()
    cimv = {"Xyl": 98.0, "Glc": 30.0, "Ara": 12.0, "Olig": 25.0, "For": 14.0, "Ac": 5.0,
            "Fur": 0.8, "HMF": 0.3, "Trc": 70.0}
    scen = [("BC18 recipe", base),
            ("xylose-only feed", dc.replace(base, feed_composition={"Xyl": 317.0})),
            ("300 rpm (O2-limited)", dc.replace(base, stirrer=300.0)),
            ("70% CIMV feed, no FDH", dc.replace(base, feed_composition=cimv)),
            ("70% CIMV feed, FDH strain", dc.replace(base, feed_composition=cimv,
                                                     strain=Strain(name="H4031", fdh=1.0)))]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), facecolor=SURFACE)
    for (name, ex), colr in zip(scen, CAT):
        r = simulate(ex).truth
        axes[0].plot(r["t"], r["Xol"], color=colr, lw=2, label=name)
        axes[1].plot(r["t"], r["X"], color=colr, lw=2, label=name)
        axes[2].plot(r["t"], r["pO2"], color=colr, lw=1.6, label=name)
    for ax, (t, u) in zip(axes, [("Xylitol", "g/L"), ("Biomass", "g DW/L"), ("Dissolved O2", "% air sat.")]):
        style(ax, t, u)
        ax.set_xlabel("time (h)", color=INK2, fontsize=8)
        ax.set_xlim(0, 145)
        ax.set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left", labelcolor=INK2)
    fig.suptitle("What-if scenarios on the BC18 recipe", x=0.01, ha="left", color=INK, fontsize=12,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(HERE / "scenarios.png", dpi=150, facecolor=SURFACE)

    cb = res.carbon_balance()
    print("carbon closure (excluding YPS2 excess): %.5f" % cb["closure_without_YPS2_excess"])
    print("xylitol C above stoichiometry: %.3f C-mol (%.1f %% of inputs)"
          % (cb["xylitol_excess_Cmol"], 100 * cb["xylitol_excess_Cmol"] / cb["C_in"]))


if __name__ == "__main__":
    main()
