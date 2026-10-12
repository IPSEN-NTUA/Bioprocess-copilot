"""
Build the pilot corpus and its answer key.

    python build_pilot.py [--out ../pilot]

corpus/      what the benchmark systems see (instrument exports, spreadsheets, notes, records)
answer_key/  ground truth, never indexed: true states, strain truth, faults, decisions, noise settings
"""
from __future__ import annotations

import argparse
import datetime as dtm
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import campaign as CP
import data_files as DF
import measure as ms
import narrative as NR
import world as W


# ---------------------------------------------------------------------------
def first_t(samples, key, cond, after=0.0):
    for s in samples:
        if s["t"] > after and cond(s["hplc"].get(key)):
            return s["t"]
    return None


def summarise(rec) -> dict:
    spec = rec["spec"]
    bro = [s for s in rec["samples"] if s["kind"] == "broth"]
    g = lambda s, a: s["hplc"].get(a)
    et = [(g(s, "EtOH") or 0.0, s["t"]) for s in bro]
    peak, peak_t = max(et)
    sm = dict(od0=bro[0]["od"]["od"], t_end=bro[-1]["t"],
              xol_end=g(bro[-1], "Xol") or 0.0, xyl_end=g(bro[-1], "Xyl") or 0.0, for_end=g(bro[-1], "For"),
              dw_end=bro[-1]["dw"]["mean"] if "dw" in bro[-1] else None,
              dw_max=max(s["dw"]["mean"] for s in bro) if "dw" in bro[0] else None,
              etoh_peak=peak, etoh_peak_t=peak_t,
              glc_out_t=first_t(bro, "Glc", lambda v: v is None) or bro[-1]["t"],
              etoh_gone_t=first_t(bro, "EtOH", lambda v: v is None, after=peak_t) or bro[-1]["t"])
    # culture stalled (as an operator would see it): after 24 h, measured xylose climbs by > 2 g/L
    # between two samples while xylitol does not increase - substrate is fed but no longer converted
    # (must hold for two consecutive intervals, so a single odd sample does not trigger it)
    def stalled(a_, b_):
        if not all(g(x, k) is not None for x in (a_, b_) for k in ("Xyl", "Xol")):
            return False
        return g(b_, "Xyl") - g(a_, "Xyl") > 2.0 and g(b_, "Xol") - g(a_, "Xol") <= 0.5
    dead = [b_["t"] for a_, b_, c_ in zip(bro[:-2], bro[1:-1], bro[2:])
            if b_["t"] > 24 and stalled(a_, b_) and stalled(b_, c_)]
    sm["dead_t"] = dead[0] if dead else None
    # ethanol re-appearing in the fed phase (O2 limitation / overflow)
    late = [s for s in bro if s["t"] > sm["etoh_gone_t"] + 10 and (g(s, "EtOH") or 0.0) > 0.5]
    if late:
        sm.update(etoh_late_t=late[0]["t"], etoh_late=g(bro[-1], "EtOH") or 0.0)
    feed = spec.feed or {}
    sm["feed_name"] = ("A" if feed == CP.FEED_A else "B" if feed == CP.FEED_B else
                       f"C5 {int((spec.frac or 0) * 100)} %")
    on = rec.get("online")
    if on is not None:
        w = on[(on.t > 2) & (on.t < 40)]
        i = w["pO2"].idxmin()
        sm.update(po2_min=float(w.loc[i, "pO2"]), po2_min_t=float(w.loc[i, "t"]),
                  po2_fed=float(on[on.t > 30]["pO2"].median()),
                  base_mL=on["base_counts"].iloc[-1] * ms.PUMP_ML_PER_COUNT["base"],
                  acid_mL=on["acid_counts"].iloc[-1] * ms.PUMP_ML_PER_COUNT["acid"],
                  feed_cnt_mL=on["feed_counts"].iloc[-1] * ms.PUMP_ML_PER_COUNT["feed"],
                  feed_bal_mL=(on["feed_balance_g"].iloc[0] - on["feed_balance_g"].iloc[-1]) /
                  rec["experiment"].feed_density)
    return {k: (round(v, 3) if isinstance(v, float) else v) for k, v in sm.items()}


def mean_sd(vals):
    v = np.array(vals, float)
    return f"{v.mean():.1f} +/- {v.std(ddof=1):.1f}" if len(v) > 1 else f"{v.mean():.1f}"


# ---------------------------------------------------------------------------
def build(out: Path, seed: int):
    C = CP.run_pilot(seed)
    rng = np.random.default_rng(seed + 1)
    runs = {r["spec"].run_id: r for r in C.runs}
    faults = []

    # HPLC re-run bookkeeping
    for s in C.samples:
        if s.get("hplc_rerun"):
            s["hplc_rerun_date"] = CP.HPLC_BIAS_FAULT["rerun_date"]
    faults.append(dict(code="F-HPLCSTD", run="several (P3 round 1)", where="HPLC SEQ-20250321",
                       truth="xylitol standard prepared ~10 % too dilute -> xylitol results 10 % high",
                       evidence="DEV-2025-02, re-run sequence SEQ-20250409, TR results workbook comments; MD master log keeps the biased values"))
    # sample tube swap FB-016 S07 <-> FB-017 S07 (no document mentions it)
    a, b = "FB012-S07", "FB013-S07"
    sa = next(s for s in C.samples if s["sample_id"] == a)
    sb = next(s for s in C.samples if s["sample_id"] == b)
    sa["hplc"], sb["hplc"] = sb["hplc"], sa["hplc"]
    sa["dilution"], sb["dilution"] = sb["dilution"], sa["dilution"]
    for seq in C.sequences.values():
        ia = [i for i in seq["injections"] if i["sample_id"] in (a, b)]
        if len(ia) == 2:
            ia[0]["results"], ia[1]["results"] = ia[1]["results"], ia[0]["results"]
            ia[0]["dilution"], ia[1]["dilution"] = ia[1]["dilution"], ia[0]["dilution"]
    faults.append(dict(code="F-SWAP", run="FB-012 / FB-013", where="sample S07 (72 h)",
                       truth="tubes swapped before HPLC; the HPLC values of the two samples are exchanged",
                       evidence="only in the data: S07 breaks the trend of both runs"))
    for rec in C.runs:
        for f in rec["spec"].faults:
            faults.append(dict(code=f["code"], run=rec["spec"].run_id, where=rec["spec"].vessel, truth=f["what"],
                               evidence=dict(**{"F-INOC": "MD lab note", "F-FEEDCAL": "EV lab note, DEV-2025-01, maintenance log, feed balance vs pump counter",
                                                "F-DODRIFT": "online pO2 only, maintenance log (membrane replaced)",
                                                "F-POWER": "gap in online log, MD note, DEV-2025-03, maintenance log",
                                                "F-ANTIFOAM": "MD note only; pO2 step in the online log",
                                                "F-WRONGLOT": "MD media sheet (lot used) vs run plan/notes, feed-sample formate vs CoA, DEV-2025-04"}).get(f["code"], "")))

    summ = {rid: summarise(r) for rid, r in runs.items() if r["spec"].kind == "reactor"}
    for rid, sm in summ.items():
        codes = {f["code"] for f in runs[rid]["spec"].faults}
        sm["md_note"] = ("power cut ~49h" if "F-POWER" in codes else "AF added 30h" if "F-ANTIFOAM" in codes
                         else "old preculture" if "F-INOC" in codes else "")

    ratio5 = summ["FB-005"]["feed_bal_mL"] / summ["FB-005"]["feed_cnt_mL"]
    rr = [s["hplc"]["Xol"] / s["hplc_rerun"]["Xol"] for s in C.samples if s.get("hplc_rerun") and s["hplc_rerun"].get("Xol")]
    rr_pct = [100 * (1 - 1 / x) for x in rr]
    rr_mean, rr_min, rr_max, rr_n = float(np.mean(rr_pct)), float(np.min(rr_pct)), float(np.max(rr_pct)), len(rr)
    for rid, sm in summ.items():
        sm["rerun_samples"] = [s["sample_id"].split("-")[-1] for s in runs[rid]["samples"] if s.get("hplc_rerun")]

    # ---------------- measured comparisons used in decisions ----------------
    def xol(rid, rerun=True):
        r = runs[rid]
        s = [x for x in r["samples"] if x["kind"] == "broth"][-1]
        v = s.get("hplc_rerun", {}).get("Xol") if rerun and s.get("hplc_rerun") else s["hplc"]["Xol"]
        return v or 0.0
    p1 = [xol(r) for r in ("FB-001", "FB-002", "FB-003")]
    p3 = {"YX-001": [xol("FB-004"), xol("FB-010")], "YX-002": [xol("FB-008")],
          "YX-003": [xol("FB-006"), xol("FB-009")], "YX-004": [xol("FB-007"), xol("FB-011")]}
    ctrl = np.mean(p3["YX-001"])
    gain = {k: 100 * (np.mean(v) / ctrl - 1) for k, v in p3.items()}
    sd_p1 = float(np.std(p1, ddof=1))
    lsd = 2 * np.sqrt(2) * sd_p1                                   # ~95 % band for a difference of two single runs
    status = {"YX-004": "dropped Apr-2025" if gain["YX-004"] < 0 else "on hold",
              "YX-002": "parent of YX-003/004 (kept)"}

    # measured flask values (what PN sees)
    def fl(rid, hour, key):
        r = runs[rid]
        s = min((x for x in r["samples"] if x["kind"] == "broth"), key=lambda x: abs(x["t"] - hour))
        return s["od"]["od"] if key == "OD" else (s["hplc"].get(key) or 0.0)
    def flasks(phase, **kw):
        return [rid for rid, r in runs.items() if r["spec"].phase == phase and all(
            (getattr(r["spec"], k) == v) for k, v in kw.items())]
    P2 = flasks("P2")
    od24 = [fl(r, 24, "OD") for r in P2]
    od0 = [fl(r, 0.1, "OD") for r in P2]
    txt_0307 = (f"Started the clone screen on Tuesday ({len(P2)} flasks). "
                + ("All flasks grew" if min(od24) > 2 * max(od0) else "Not all flasks grew")
                + (f"; OD at 24 h looked similar between clones ({min(od24):.1f}-{max(od24):.1f})."
                   if max(od24) / min(od24) < 1.25 else f"; OD at 24 h ranged from {min(od24):.1f} to {max(od24):.1f}."))
    y1_70 = flasks("P4", strain_id="YX-001", frac=0.7)[0]
    fdh_50 = [r for r in flasks("P4", frac=0.5, buffered=True) if runs[r]["spec"].strain_id != "YX-001"]
    sick = fl(y1_70, 24, "OD") < 0.5 * np.median([fl(r, 24, "OD") for r in fdh_50])
    txt_0425 = ("C5 tolerance flasks started on Tuesday with lot C5-2501. "
                + (f"The 70 % flask without FDH looks sick already at 24 h (OD {fl(y1_70, 24, 'OD'):.1f})."
                   if sick else "Hard to tell differences at 24 h from OD alone."))
    y1_50b = flasks("P4", strain_id="YX-001", frac=0.5, buffered=True)[0]
    unb = flasks("P4", buffered=False)
    worse = all(fl(u, 72, "Xol") < 0.5 * fl(y1_50b, 72, "Xol") for u in unb)
    unb_x = ", ".join((f"{fl(u, 72, 'Xol'):.1f}" if fl(u, 72, 'Xol') > 0 else "n.d.") for u in unb)
    txt_0502 = (f"Unbuffered flasks ({', '.join(runs[u]['spec'].strain_id for u in unb)}) were clearly worse than buffered ones at the same C5 level"
                f" (xylitol at 72 h {unb_x} g/L without CaCO3 vs "
                f"{fl(y1_50b, 72, 'Xol'):.1f} g/L for buffered YX-001) - probably because more formate is in the acid form at low pH."
                if worse else "Buffered and unbuffered flasks were not very different.")
    fdh_70 = [r for r in flasks("P4", frac=0.7) if runs[r]["spec"].strain_id != "YX-001"]
    x170, od170, od170_0 = fl(y1_70, 72, 'Xol'), fl(y1_70, 72, 'OD'), fl(y1_70, 0.1, 'OD')
    if od170 < 1.5 * od170_0:
        head = ("At 70 % C5 the parent YX-001 made " + (f"{x170:.1f} g/L xylitol" if x170 > 0 else "no detectable xylitol")
                + f" and did not grow (OD {od170:.1f} at 72 h, inoculated at OD {od170_0:.1f}); ")
    else:
        head = f"At 70 % C5 the parent YX-001 made {x170:.1f} g/L xylitol (OD {od170:.1f} at 72 h); "
    p4_min = head + "; ".join(f"{runs[r]['spec'].strain_id} ({r}) reached OD {fl(r, 72, 'OD'):.1f} and "
                              f"{fl(r, 72, 'Xol'):.1f} g/L xylitol" for r in fdh_70) + "."

    out_c, out_k = out / "corpus", out / "answer_key"
    if out.exists():
        shutil.rmtree(out)
    for p in ["online", "offgas", "hplc/raw", "hplc/processed", "lims", "spreadsheets", "strains", "materials",
              "equipment", "deviations", "sops", "notes", "meetings", "reports"]:
        (out_c / p).mkdir(parents=True, exist_ok=True)
    (out_k / "timeseries").mkdir(parents=True, exist_ok=True)

    # ---------------- instrument data ----------------
    for rid, rec in runs.items():
        if rec["spec"].kind != "reactor":
            continue
        v = W.VESSELS[rec["spec"].vessel]
        if v.log_format == "bioctrl_tsv":
            DF.online_bioctrl(rec, out_c / "online" / f"{rid}_{v.vessel_id}.txt", rng)
        else:
            DF.online_vesselview(rec, out_c / "online" / f"{rid.replace('-', '')}_{v.vessel_id}.csv", rng)
        if "offgas" in rec:
            DF.offgas_xlsx(rec, out_c / "offgas" / f"GA4_{rid.replace('-', '').lower()}_{v.vessel_id}.xlsx")
    for d, seq in sorted(C.sequences.items()):
        DF.hplc_sequence_csv(d, seq, out_c / "hplc/raw" / f"SEQ-{d:%Y%m%d}.csv")
    bymonth = {}
    for s in C.samples:
        bymonth.setdefault(s["hplc_date"].strftime("%Y-%m"), []).append(s)
    for m, ss in bymonth.items():
        DF.hplc_processed_xlsx(sorted(ss, key=lambda x: x["sample_id"]), m, out_c / "hplc/processed" / f"HPLC_results_{m}.xlsx", C.sequences)
    DF.sample_register_csv(C, out_c / "lims" / "sample_register.csv")

    # ---------------- people's spreadsheets ----------------
    for rid, rec in runs.items():
        if rec["spec"].kind == "reactor" and rec["spec"].operator == "EV":
            sm = summ[rid]
            calc = {"Final xylitol (g/L)": round(sm["xol_end"], 1), "Residual xylose (g/L)": round(sm["xyl_end"], 1),
                    "Max DW (g/L)": round(sm["dw_max"], 2), "Feed added by pump (mL)": (round(sm["feed_cnt_mL"], 1) if rec["spec"].vessel in ("R1", "R2")
                                                                       else "not logged (VesselView)"),
                    "Feed bottle start (g)": round(float(rec["online"]["feed_balance_g"].iloc[0]), 1),
                    "Feed bottle end (g)": round(float(rec["online"]["feed_balance_g"].iloc[-1]), 1),
                    "Feed added by balance (mL)": round(sm["feed_bal_mL"], 1),
                    "Xylitol productivity (g/L/h)": round(sm["xol_end"] / sm["t_end"], 3)}
            DF.ev_run_workbook(rec, out_c / "spreadsheets" / f"EV_{rid}.xlsx", calc)
    DF.md_master_log(C, out_c / "spreadsheets" / "MD_fermentation_master_2025.xlsx", summ)
    p2 = [r for r in C.runs if r["spec"].phase == "P2"]
    p4 = [r for r in C.runs if r["spec"].phase == "P4"]
    DF.pn_wide_workbook(p2, "Clone screen", [
        "Shake-flask clone screen of YX-002, YX-003 and YX-004 (3 clones each, duplicate flasks) against the parent YX-001.",
        "Medium: YP + 20 g/L glucose + 50 g/L xylose, CaCO3 buffered. Sheets named by strain; flask IDs SF-xxx.",
        "OD600 measured after dissolving CaCO3 in HCl. HPLC by TR (SOP-ANA-003 v1)."],
        [["Strains", "YX-001 parent; YX-002 c1-c3; YX-003 c1-c3 (built on YX-002 c1); YX-004 c1-c3 (built on YX-002 c1)"],
         ["Culture", "50 mL in 250 mL baffled flasks, 30 C, 200 rpm, inoculated to OD 2 from YPD preculture"],
         ["Replicates", "two flasks per clone"]], out_c / "spreadsheets" / "PN_clone_screen_Mar2025.xlsx", flask=True)
    DF.pn_wide_workbook(p4, "C5 tolerance", [
        "Tolerance of YX-001, YX-005 (FDH) and YX-006 (PPP + FDH) to the C5 fraction lot C5-2501 at 30, 50 and 70 %.",
        "FDH clone screen at 50 % (3 clones each); chosen clones then tested at 30 and 70 % and without CaCO3.",
        "Medium: YP + 20 g/L glucose + C5 fraction; CaCO3 buffered unless noted."],
        [["Culture", "50 mL in 250 mL baffled flasks, 30 C, 200 rpm, OD0 2"],
         ["Note", "C5 fraction pH adjusted to ~5.5 with NaOH before mixing; not sterilised"]],
        out_c / "spreadsheets" / "PN_C5_tolerance_Apr2025.xlsx", flask=True)
    DF.pn_wide_workbook([runs["FB-022"], runs["FB-023"]], "Lot comparison",
                        ["FB-022 and FB-023: YX-006 on C5 lot C5-2502 (60 %), FB-023 with enzyme at 24 h.",
                         "FB-024/FB-025 (lot C5-2503) were run by MD - see his master log."], [],
                        out_c / "spreadsheets" / "PN_FB022-023_lot_comparison.xlsx", flask=False)
    DF.strain_registry_xlsx(C, out_c / "strains" / "strain_registry.xlsx", status)

    # ---------------- notes ----------------
    for rid, rec in runs.items():
        if rec["spec"].kind == "reactor":
            (out_c / "notes" / f"{rid}_{rec['spec'].operator}.md").write_text(NR.note_reactor(rec, summ[rid], rng))
    dec = {d["strain"]: d for d in C.decisions}

    def flask_sums(rs, hour, label_fn):
        out = {}
        for r in rs:
            s = min((x for x in r["samples"] if x["kind"] == "broth"), key=lambda x: abs(x["t"] - hour))
            out.setdefault(label_fn(r), []).append(s["hplc"]["Xol"])
        return {k: (None if all(x is None for x in v) else float(np.mean([x or 0.0 for x in v])), len(v))
                for k, v in out.items()}
    lab = lambda r: f"{r['spec'].strain_id}" + (f" c{r['clone'] + 1}" if len(C.strains[r['spec'].strain_id]['clones']) > 1 else "") + \
        (f", {int(r['spec'].frac * 100)} % C5" if r["spec"].frac else "") + ("" if r["spec"].buffered else ", no CaCO3")
    def scatter(rs, hour):
        by = {}
        for r in rs:
            by.setdefault(lab(r), []).append(fl(r["spec"].run_id, hour, "Xol"))
        dup = [np.std(v, ddof=1) for v in by.values() if len(v) > 1]
        return float(np.std([np.mean(v) for v in by.values()], ddof=1)), float(np.sqrt(np.mean(np.square(dup))))
    (out_c / "notes" / "SF_screen_P2_PN.md").write_text(NR.note_flask_screen(
        "P2", p2, flask_sums(p2, 48, lab), [dec["YX-002"], dec["YX-003"], dec["YX-004"]], rng, scatter(p2, 48), hour=48))
    (out_c / "notes" / "SF_screen_P4_PN.md").write_text(NR.note_flask_screen(
        "P4", p4, flask_sums(p4, 72, lab), [dec["YX-005"], dec["YX-006"]], rng, hour=72))

    # ---------------- strains, materials, SOPs, maintenance ----------------
    for plan in W.STRAIN_PLAN:
        (out_c / "strains" / f"construction_{plan.strain_id}.md").write_text(
            NR.construction_record(plan, C.strains[plan.strain_id], dec.get(plan.strain_id)))
    coa_vals = {}
    for lot in W.LOTS:
        (out_c / "materials" / f"CoA_{lot}.md").write_text(NR.coa(lot, rng, coa_vals))
    (out_c / "sops" / "SOP-ANA-003_v1.md").write_text(NR.sop_hplc(1))
    (out_c / "sops" / "SOP-ANA-003_v2.md").write_text(NR.sop_hplc(2))
    (out_c / "sops" / "SOP-FB-001_v3.md").write_text(NR.sop_fedbatch())
    fb10 = runs["FB-010"]["spec"].start + dtm.timedelta(hours=49.5)
    maint = [f"# {W.SYNTHETIC_NOTICE}", "date,equipment,event,by",
             "2025-01-20,HPLC-2,new guard column,TR",
             "2025-02-28,R1-R4,pH probes two-point calibration (pH 4/7),MD",
             f"2025-03-26,R2 feed pump,recalibrated: delivered {2.0 * ratio5:.2f} mL/h at 2.00 mL/h setpoint - tubing worn; replaced,EV",
             f"{fb10:%Y-%m-%d},Lab UPS,UPS tripped at {fb10:%H:%M}; R4 controller and feed pump off ~3.5 h,MD",
             "2025-04-07,R3 pO2 probe,slow response and low reading; membrane and electrolyte replaced,MD",
             "2025-05-27,GA-4,span calibration with 20.95 % O2 / 2.00 % CO2 test gas,TR",
             "2025-06-01,HPLC-2,UV detector installed in series (SOP-ANA-003 v2),TR",
             "2025-06-16,R2,exhaust filter replaced (wet),MD",
             "2025-07-01,R1-R4,annual service by vendor,MD"]
    (out_c / "equipment" / "maintenance_log.csv").write_text("\n".join(maint) + "\n")
    (out_c / "equipment" / "reactor_setup.md").write_text(NR.reactor_setup())

    # ---------------- deviations ----------------
    f24 = next(s for s in runs["FB-024"]["samples"] if s["kind"] == "feed")["hplc"]["For"]
    exp03 = coa_vals["C5-2503"]["For"] * 0.6                     # what the author can compute: CoA x 60 %
    exp02 = coa_vals["C5-2502"]["For"] * 0.6
    devs = {
        "DEV-2025-01": NR.deviation("DEV-2025-01", "2025-03-26", "EV", "FB-005 (R2)", "Feed pump over-delivery",
            f"Feed balance at the end of FB-005 showed {summ['FB-005']['feed_bal_mL']:.0f} mL consumed while the pump counter gave "
            f"{summ['FB-005']['feed_cnt_mL']:.0f} mL.", f"Worn pump tubing on the R2 feed pump; the pump delivered about {100 * (ratio5 - 1):.0f} % more than its setpoint.",
            "FB-005 received more glucose and xylose than planned. Exclude FB-005 from the strain comparison.",
            "Tubing replaced and pump recalibrated (maintenance log 2025-03-26). Repeat run FB-008."),
        "DEV-2025-02": NR.deviation("DEV-2025-02", "2025-04-04", "TR", "HPLC sequence SEQ-20250321", "Xylitol standard error",
            "The xylitol check standard in SEQ-20250321 read low against the new standard lot.",
            "Xylitol stock standard prepared from a partially hydrated batch, so the xylitol calibration read low and sample results high.",
            f"Xylitol values in SEQ-20250321 (P3 round 1 samples) are biased high: on re-analysis they were on average "
            f"{rr_mean:.0f} % lower (range {rr_min:.0f}-{rr_max:.0f} %, n = {rr_n}). Other analytes unaffected.",
            "Retained samples re-analysed in SEQ-20250409 with fresh standards; corrected values in HPLC_results_2025-03.xlsx (comment column). "
            "Run owners to update their spreadsheets."),
        "DEV-2025-03": NR.deviation("DEV-2025-03", "2025-04-03", "MD", "FB-010 (R4)", "Power interruption",
            f"The lab UPS tripped on {fb10:%d %b} at about {fb10:%H:%M} (about 49.5 h into FB-010). R4 controller and feed pump were off "
            "for about 3.5 h; temperature and pH control resumed automatically, feed was restarted by hand.",
            "UPS overload during building maintenance.", "No feed for ~3.5 h; online data missing for the same period. Run continued.",
            "Building services informed; R4 moved to a different UPS circuit."),
        "DEV-2025-04": NR.deviation("DEV-2025-04", "2025-07-09", "MD", "FB-024 (R3)", "Wrong C5 lot in feed",
            f"FB-024 was planned with lot C5-2503. The formate in the FB-024 feed sample was {f24:.1f} g/L, close to the "
            f"{exp02:.1f} g/L expected for C5-2502 at 60 % and well below the {exp03:.1f} g/L expected for C5-2503.",
            "Feed for R3 prepared from the C5-2502 container (both containers in the same cold-room shelf).",
            "FB-024 is a C5-2502 run, not C5-2503. The lot comparison has only one C5-2503 run (FB-025, with enzyme).",
            "Containers relabelled and stored separately. Media sheet already shows the lot actually used."),
    }
    for k, v in devs.items():
        (out_c / "deviations" / f"{k}.md").write_text(v)

    # ---------------- minutes and weekly reports ----------------
    s = summ
    x900 = np.mean([s["FB-018"]["xol_end"], s["FB-020"]["xol_end"]]); x500 = np.mean([s["FB-019"]["xol_end"], s["FB-021"]["xol_end"]])
    dB = s["FB-021"]["xol_end"] - s["FB-020"]["xol_end"]
    p6_action = ("Keep 900 rpm as default." if x900 >= x500 else
                 "500 rpm was not worse on average; keep 900 rpm anyway for feed B.") + (
                 f" At 500 rpm feed B lost {-dB:.0f} g/L xylitol vs 900 rpm." if dB < -3 else "")
    de = s["FB-023"]["xol_end"] - s["FB-022"]["xol_end"]
    ff = {r: next(x for x in runs[r]["samples"] if x["kind"] == "feed")["hplc"]["For"] for r in ("FB-022", "FB-024", "FB-025")}
    txt_0704 = ("All four P7 runs finished. Enzyme effect on lot C5-2502 " +
                (f"looks small ({de:+.1f} g/L, FB-023 vs FB-022). " if abs(de) < 0.05 * s["FB-022"]["xol_end"] else f"is {de:+.1f} g/L (FB-023 vs FB-022). ") +
                (f"The FB-024 feed sample has {ff['FB-024']:.1f} g/L formate, much closer to FB-022 ({ff['FB-022']:.1f}) than to FB-025 "
                 f"({ff['FB-025']:.1f}), although FB-024 and FB-025 should both be C5-2503 - asked Tomas to double-check."
                 if abs(ff['FB-024'] - ff['FB-022']) < abs(ff['FB-024'] - ff['FB-025']) else ""))
    mins = {
        "2025-02-14": NR.minutes("2025-02-14", ["EV", "MD", "PN", "TR"], [
            f"Baseline fed-batches FB-001 to FB-003 (YX-001): final xylitol {mean_sd(p1)} g/L at ~143 h.",
            f"FB-002 started slowly (old preculture, OD0 {s['FB-002']['od0']:.1f}) but "
            + ("ended close to the others." if abs(p1[1] / np.mean([p1[0], p1[2]]) - 1) < 0.05 else
               f"ended {100 * (p1[1] / np.mean([p1[0], p1[2]]) - 1):+.0f} % from the others."),
            "HPLC turnaround fine with the Tuesday/Friday batches."],
            ["PN to finish YX-002; YX-003 (ZWF1/GND1) and YX-004 (gpd1) to be built on the first confirmed YX-002 clone.",
             "Clone screen in flasks once all three strains are ready."]),
        "2025-03-14": NR.minutes("2025-03-14", ["EV", "MD", "PN", "TR"], [
            "Flask screen done (see PN note). Clones picked: " + ", ".join(f"{k} {d['chosen']}" for k, d in dec.items() if k in ("YX-002", "YX-003", "YX-004")) + ".",
            "YX-003 and YX-004 were built on YX-002 c1 before the screen; " + (
                "c1 also came out best, so nothing to redo." if dec["YX-002"]["chosen"] == "c1" else
                f"the screen preferred YX-002 {dec['YX-002']['chosen']}, but we keep c1 as the background for now.")],
            ["Bioreactor comparison (P3) in two rounds from 17 Mar: YX-001 control, YX-002, YX-003, YX-004."]),
        "2025-04-11": NR.minutes("2025-04-11", ["EV", "MD", "PN", "TR"], [
            "FB-005 excluded (feed pump over-delivery, DEV-2025-01). Round-1 xylitol values corrected after DEV-2025-02.",
            f"Final xylitol (mean of runs): control {np.mean(p3['YX-001']):.1f}, YX-002 {np.mean(p3['YX-002']):.1f}, "
            f"YX-003 {np.mean(p3['YX-003']):.1f}, YX-004 {np.mean(p3['YX-004']):.1f} g/L.",
            f"YX-004 change vs control {gain['YX-004']:+.2f} %; YX-003 {gain['YX-003']:+.2f} %.",
            "FB-010 lost ~3.5 h of feed in the power cut (DEV-2025-03).",
            f"Run-to-run scatter is about {sd_p1:.1f} g/L (sd of FB-001 to FB-003), so differences below ~{lsd:.0f} g/L "
            "between single runs are not meaningful. Final titre is mostly set by the xylose fed, so it is a weak readout "
            "of strain differences.",
            "FB-009 pO2 drifted low in the second half; R3 membrane and electrolyte were replaced on 7 Apr."],
            [("Drop YX-004 (no gain over control)." if status["YX-004"].startswith("dropped") else
              "Keep YX-004 on hold: the gain rests on two runs per strain and the control includes FB-010 (power cut); "
              "without GPD1 the strain may also be less robust in the inhibitor-rich C5 media."),
             "C5 work continues with the FDH strains already built: YX-005 (FDH in YX-001) and YX-006 (FDH in YX-003, "
             "the cofactor-supply PPP background). Flask screen from 22 Apr.",
             "Treat FB-009 pO2 after ~70 h with caution."]),
        "2025-05-02": NR.minutes("2025-05-02", ["EV", "MD", "PN"], [
            f"FDH clone screen at 50 % C5: picked YX-005 {dec['YX-005']['chosen']} and YX-006 {dec['YX-006']['chosen']}.",
            p4_min],
            ["P5 fed-batches on lot C5-2501: YX-003 (no FDH) at 70 % as negative control, YX-005 and YX-006 at 50 and 70 %."]),
        "2025-05-30": NR.minutes("2025-05-30", ["EV", "MD", "PN", "TR"], [
            f"FB-012 (YX-003, no FDH, 70 % C5) failed: DW {s['FB-012']['dw_end']:.1f} g/L, formate {s['FB-012']['for_end'] or 0:.1f} g/L at the end.",
            "FDH strains reached " + ", ".join(f"{r} {s[r]['xol_end']:.0f} g/L" for r in ("FB-013", "FB-014", "FB-015", "FB-016", "FB-017")) + ".",
            "Foaming in FB-014 at ~30 h, antifoam added.",
            "HPLC v2 method from 1 June adds succinate, furfural and HMF."],
            ["P6: YX-006 at 900 vs 500 rpm with feed A and feed B (more glucose).", "P7: compare C5 lots 2502 and 2503, with and without enzyme."]),
        "2025-06-20": NR.minutes("2025-06-20", ["EV", "MD", "PN", "TR"], [
            f"P6 final xylitol: 900 rpm A {s['FB-018']['xol_end']:.0f}, 500 rpm A {s['FB-019']['xol_end']:.0f}, "
            f"900 rpm B {s['FB-020']['xol_end']:.0f}, 500 rpm B {s['FB-021']['xol_end']:.0f} g/L.",
            f"pO2 in the fed phase: {s['FB-018']['po2_fed']:.0f} / {s['FB-019']['po2_fed']:.0f} / {s['FB-020']['po2_fed']:.0f} / {s['FB-021']['po2_fed']:.0f} %."],
            [p6_action]),
        "2025-07-11": NR.minutes("2025-07-11", ["EV", "MD", "PN", "TR"], [
            "DEV-2025-04: FB-024 was fed with C5-2502, not C5-2503.",
            f"Final xylitol: FB-022 {s['FB-022']['xol_end']:.0f}, FB-023 (+enzyme) {s['FB-023']['xol_end']:.0f}, "
            f"FB-024 {s['FB-024']['xol_end']:.0f}, FB-025 (C5-2503 +enzyme) {s['FB-025']['xol_end']:.0f} g/L."],
            ["Repeat C5-2503 without enzyme in August."]),
    }
    for d, t in mins.items():
        (out_c / "meetings" / f"minutes_{d}.md").write_text(t)
    rep = {
        "2025-03-07": [txt_0307,
                       "Next week: HPLC results and clone picks."],
        "2025-03-14": ["Clone screen results are in (see note). Picks: " + ", ".join(f"{k} {dec[k]['chosen']}" for k in ("YX-002", "YX-003", "YX-004")) + ".",
                       "I think the PPP strain might need more glucose to show its NADPH advantage - worth testing later."],
        "2025-04-25": [txt_0425],
        "2025-05-02": ["Finished the C5 flasks. FDH clones picked: " + ", ".join(f"{k} {dec[k]['chosen']}" for k in ("YX-005", "YX-006")) + ".",
                       txt_0502],
        "2025-06-27": ["Started FB-022 and FB-023 (lot C5-2502); Mark has FB-024/025 on C5-2503. Enzyme goes into FB-023 and FB-025 tomorrow (24 h)."],
        "2025-07-04": [txt_0704],
    }
    for d, ps in rep.items():
        (out_c / "reports" / f"weekly_PN_{d}.md").write_text(NR.weekly_report(d, ps))

    # ---------------- answer key ----------------
    exp_rows, ts_cols = [], ["t", "V", "X", "Xd", "Glc", "Xyl", "Xol", "Ara", "Olig", "EtOH", "Gly", "Ac", "Suc", "Pyr", "For",
                              "Fur", "HMF", "DO", "pO2", "OUR", "CER", "RQ", "stirrer", "cumFeed", "cumBase", "cumAcid", "AF"]
    for rec in C.runs:
        sp = rec["spec"]
        st = rec["experiment"].strain
        exp_rows.append(dict(run_id=sp.run_id, phase=sp.phase, kind=sp.kind, vessel=sp.vessel, operator=sp.operator,
                             start=sp.start.isoformat(), strain=sp.strain_id, clone=f"c{rec['clone'] + 1}",
                             purpose=sp.purpose, lot_recorded=sp.lot, lot_actual=sp.lot_actual, frac=sp.frac,
                             stirrer=sp.stirrer, pH=sp.pH, buffered=sp.buffered, enzyme_h=sp.enzyme_h,
                             faults=";".join(f["code"] for f in sp.faults),
                             **{f"draw_{k}": round(v, 4) for k, v in rec["draws"].items()},
                             mumax2_used=st.mumax2, KPI2_used=st.KPI2, fdh_used=st.fdh,
                             final_Xol_true=rec["result"].truth["Xol"].iloc[-1], final_Xyl_true=rec["result"].truth["Xyl"].iloc[-1],
                             final_X_true=rec["result"].truth["X"].iloc[-1],
                             final_Xd_true=rec["result"].truth["Xd"].iloc[-1]))
        tr = ms.resample(rec["result"].truth, 1 / 6, [c for c in ts_cols if c != "t"])
        tr.round(4).to_csv(out_k / "timeseries" / f"{sp.run_id}.csv", index=False)
    pd.DataFrame(exp_rows).to_csv(out_k / "experiments.csv", index=False)
    srows = []
    for sid, sd in C.strains.items():
        for i, (r, fdh) in enumerate(sd["clones"]):
            srows.append(dict(strain=sid, clone=f"c{i + 1}", chosen=(sd["chosen"] == i), built_on=sd["built_on"],
                              genotype=r.genotype(), fdh=fdh, mumax2=r.strain.mumax2,
                              **{f"NRA_{k}_pct": round(100 * (np.exp(v) - 1), 2) for k, v in r.dlnJ.items()}))
    pd.DataFrame(srows).to_csv(out_k / "strains_truth.csv", index=False)
    pd.DataFrame([dict(sample_id=s["sample_id"], run_id=s["run_id"], kind=s["kind"], t=s["t"],
                       **{f"true_{k}": round(v, 4) for k, v in s["truth"].items()},
                       **{f"meas_{k}": (None if v is None else round(v, 4)) for k, v in s["hplc"].items()})
                  for s in C.samples]).to_csv(out_k / "samples_truth.csv", index=False)
    pd.DataFrame(faults).to_csv(out_k / "faults.csv", index=False)
    pd.DataFrame(C.decisions).to_csv(out_k / "clone_decisions.csv", index=False)
    noise = {k: getattr(ms, k) for k in dir(ms) if k.isupper() and isinstance(getattr(ms, k), (int, float, dict, list))}
    noise.update(VAR_REACTOR=CP.VAR_REACTOR, VAR_FLASK=CP.VAR_FLASK)
    (out_k / "noise_and_variability.json").write_text(json.dumps(noise, indent=1, default=str))
    (out_k / "README.md").write_text(f"# Answer key (do not index)\n\n{W.SYNTHETIC_NOTICE}\n\nGround truth for the pilot corpus: "
                                     "experiments.csv (recipes, variability draws, faults), timeseries/ (true states every 10 min), "
                                     "samples_truth.csv (true vs measured), strains_truth.csv (clone truth, NRA flux changes), "
                                     "faults.csv (fault, truth, where the evidence is), clone_decisions.csv, noise_and_variability.json.\n")
    return C, summ, faults, gain


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "pilot"))
    ap.add_argument("--seed", type=int, default=20250217)
    a = ap.parse_args()
    C, summ, faults, gain = build(Path(a.out), a.seed)
    n = sum(1 for _ in Path(a.out, "corpus").rglob("*") if _.is_file())
    print(f"corpus files: {n}; strain gains vs control (measured, P3): " + ", ".join(f"{k} {v:+.1f}%" for k, v in gain.items()))
