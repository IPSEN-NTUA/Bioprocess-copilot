"""
Simulator checks for the experiment-design questions (answer key only, never in the corpus).

The design questions ask a system what the team should do next. The gold rubric is built from the
corpus; this script runs candidate experiments on the TRUE strains in the simulator so the paper can
also say which suggestions would actually have been informative.

All runs are noise-free and use nominal kinetics (no run-to-run variability) unless stated, so the
numbers are expected outcomes, not what a single real run would show.

    python design_checks.py [--seed 20250217] [--out ../pilot/answer_key]
"""
from __future__ import annotations

import argparse
import dataclasses as dc
import json
import math
from pathlib import Path

import numpy as np

import campaign as CP
import world as W
from simulator import Experiment, FeedSegment, simulate


def strain(C, sid, clone=None):
    s = C.strains[sid]
    c = s["chosen"] if clone is None else clone
    return s["clones"][c][0].strain


def reactor_run(st, feed, rate_mLh=2.0, stirrer=900.0, t_end=143.5, enzyme=None, cascade=None, vessel="R1",
                mult=None, ftype="A"):
    mult = mult or {}
    st = st.copy(mumax2=st.mumax2 * mult.get("mumax2", 1.0), KPI2=st.KPI2 * mult.get("KPI2", 1.0))
    r = W.VESSELS[vessel].reactor
    r = dc.replace(r, kLa_ref=r.kLa_ref * mult.get("kLa", 1.0))
    samples = CP.HYDRO_SAMPLE_H if t_end > 150 else CP.REACTOR_SAMPLE_H
    ex = Experiment(run_id="design", strain=st, reactor=r, V0=0.300, X0=0.8 * mult.get("X0", 1.0),
                    medium={"Glc": 20.0}, T_set=30.0, pH_set=5.5, stirrer=stirrer, airflow=0.7,
                    feed_segments=[FeedSegment(7.5, t_end, rate_mLh / 1000.0 * mult.get("feed", 1.0))],
                    feed_composition=dict(feed), feed_density=W.FEED_DENSITY[ftype],
                    sample_times=[t for t in samples[1:] if t < t_end], sample_volume=0.007, t_end=t_end,
                    enzyme_additions=enzyme or [], do_cascade=cascade)
    return simulate(ex)


def summary(res):
    tr = res.truth
    e = tr.iloc[-1]
    at = lambda h: float(tr.Xol[(tr.t - h).abs().idxmin()])
    fed = tr[tr.t > 30]
    return dict(Xol_end=round(float(e.Xol), 2), Xyl_end=round(float(e.Xyl), 2), Xol_48h=round(at(48), 2),
                Xol_72h=round(at(72), 2), EtOH_fed_max=round(float(fed.EtOH.max()), 2),
                pO2_fed_median=round(float(fed.pO2.median()), 1), X_end=round(float(e.X), 2),
                V_end_mL=round(float(e.V) * 1000, 0), xylitol_g=round(float(e.cumXol), 2))


def main(seed, out: Path):
    C = CP.run_pilot(seed)
    S = {k: strain(C, k) for k in ("YX-001", "YX-002", "YX-003", "YX-004", "YX-005", "YX-006")}
    R = {}

    # A/B: is the standard fed-batch a sensitive readout of strain differences? Feed rate 2 vs 3 mL/h.
    for rate in (2.0, 3.0):
        R[f"strains_feedA_{rate:.0f}mLh"] = {k: summary(reactor_run(S[k], CP.FEED_A, rate)) for k in
                                             ("YX-001", "YX-002", "YX-003", "YX-004", "YX-006")}
    # D: PN's idea (weekly report 2025-03-14): the PPP strain needs more glucose to show its advantage.
    R["strains_feedB_2mLh"] = {k: summary(reactor_run(S[k], CP.FEED_B, 2.0, ftype="B")) for k in ("YX-001", "YX-003")}
    # E: oxygen with feed B
    R["yx006_feedB_agitation"] = {f"{n:.0f} rpm": summary(reactor_run(S["YX-006"], CP.FEED_B, stirrer=n, ftype="B"))
                                  for n in (500, 700, 900)}
    R["yx006_feedB_DO_cascade_20pct"] = summary(reactor_run(S["YX-006"], CP.FEED_B, stirrer=500, ftype="B",
                                                            cascade={"setpoint": 20, "N_min": 500, "N_max": 1200, "gain": 40}))
    # F: the planned repeat - C5-2503 without enzyme, plus C5-2502 for reference
    for lot in ("C5-2502", "C5-2503"):
        f = CP.hydrolysate_feed(lot, 0.6)
        R[f"yx006_{lot}_60pct"] = {e: summary(reactor_run(S["YX-006"], f, t_end=167.5, ftype="C5",
                                                          enzyme=[(24.0, u)] if u else None))
                                   for e, u in (("no enzyme", 0), ("5 U/L", 5.0), ("20 U/L", 20.0))}
    # H: how far can the C5 fraction go with YX-006 (lot C5-2501) and with the more toxic C5-2503
    R["yx006_C5-2501_fraction"] = {f"{int(fr * 100)} %": summary(reactor_run(S["YX-006"], CP.hydrolysate_feed("C5-2501", fr),
                                                                             t_end=167.5, ftype="C5"))
                                   for fr in (0.7, 0.8, 0.9, 1.0)}
    R["yx006_C5-2503_fraction"] = {f"{int(fr * 100)} %": summary(reactor_run(S["YX-006"], CP.hydrolysate_feed("C5-2503", fr),
                                                                             t_end=167.5, ftype="C5"))
                                   for fr in (0.6, 0.7, 0.8)}
    R["yx005_vs_yx006_70pct_C5-2501"] = {k: summary(reactor_run(S[k], CP.hydrolysate_feed("C5-2501", 0.7), t_end=167.5,
                                                                ftype="C5")) for k in ("YX-005", "YX-006")}
    # I: replicates needed to resolve YX-003 vs YX-001 in the standard fed-batch (run-to-run variability as in the
    # campaign + HPLC noise ~ sqrt(2.5^2 + 1.5^2) %)
    rng = np.random.default_rng(seed + 7)
    finals = {"YX-001": [], "YX-003": []}
    for k in finals:
        for _ in range(16):
            m = {v: math.exp(rng.normal(0, sd)) for v, sd in CP.VAR_REACTOR.items()}
            x = reactor_run(S[k], CP.FEED_A, mult=m).truth.Xol.iloc[-1]
            finals[k].append(float(x * (1 + rng.normal(0, math.hypot(0.025, 0.015)))))
    d = np.mean(finals["YX-003"]) - np.mean(finals["YX-001"])
    sd = math.sqrt((np.var(finals["YX-003"], ddof=1) + np.var(finals["YX-001"], ddof=1)) / 2)
    n80 = math.ceil(2 * ((1.96 + 0.84) * sd / d) ** 2) if d > 0 else None
    R["replicates_YX003_vs_YX001"] = dict(true_difference_gL=round(d, 2), run_to_run_sd_gL=round(sd, 2),
                                          runs_per_strain_for_80pct_power=n80)

    out.mkdir(parents=True, exist_ok=True)
    (out / "design_checks.json").write_text(json.dumps(R, indent=1))
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20250217)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "pilot" / "answer_key"))
    a = ap.parse_args()
    R = main(a.seed, Path(a.out))
    print(json.dumps(R, indent=1))
