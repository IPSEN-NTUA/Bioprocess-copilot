"""
Independent check of the benchmark gold answers.

questions.py derives the answers mostly from processed spreadsheets and notes. This script re-derives the
key facts by a different route (raw HPLC sequence files, raw online logs, LIMS register) and checks:

  A. every number quoted in a gold answer is either found in one of the cited evidence files or is a
     value this script re-derives independently (means, differences, gaps, ...);
  B. re-derived values agree with the gold answers;
  C. the deliberate faults referenced by questions are the ones in answer_key/faults.csv;
  D. how far the measured gold values are from the simulator truth (information only).

    python verify_questions.py [--pilot ../pilot]
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

OUT = []


def ok(cond, msg):
    OUT.append(("PASS" if cond else "FAIL", msg))


def raw_hplc(corpus: Path) -> pd.DataFrame:
    """Sample concentrations from the raw sequence files: vial conc x dilution; re-runs flagged by sequence."""
    rows = []
    for f in sorted(glob.glob(str(corpus / "hplc/raw/SEQ-*.csv"))):
        head = open(f).read().split("\n")[1]
        d = pd.read_csv(f, comment="#")
        d["seq"] = Path(f).stem
        d["rerun"] = "re-analysis" in head
        rows.append(d)
    d = pd.concat(rows)
    d["conc"] = np.where(d.Flag.fillna("") == "<LOD", np.nan, d.Conc_vial_gL * d.Dilution)
    return d


def file_text(corpus_root: Path, rel: str) -> str:
    """Plain text of an evidence path (file or directory); xlsx flattened cell by cell."""
    p = corpus_root.parent / rel
    files = [p] if p.is_file() else [Path(x) for x in glob.glob(str(p / "**/*"), recursive=True) if Path(x).is_file()]
    out = []
    for f in files:
        if f.suffix == ".xlsx":
            wb = openpyxl.load_workbook(f, data_only=True)
            for ws in wb.worksheets:
                for r in ws.iter_rows(values_only=True):
                    out.append(" ".join("" if v is None else str(v) for v in r))
        else:
            out.append(f.read_text(errors="ignore"))
    return "\n".join(out)


def numbers(s: str):
    return [x for x in re.findall(r"(?<![\w-])\d+(?:\.\d+)?(?![\w-])", s)]


def main(pilot: Path):
    corpus, key = pilot / "corpus", pilot / "answer_key"
    Q = {q["id"]: q for q in json.loads((pilot / "questions/questions.json").read_text())}
    raw = raw_hplc(corpus)
    lims = pd.read_csv(corpus / "lims/sample_register.csv", comment="#")
    broth = lims[lims.type == "broth"]
    last = broth.sort_values("process_time_h").groupby("run_id").sample_id.last()

    def raw_val(sample, analyte, rerun=True):
        d = raw[(raw.SampleID == sample) & (raw.Analyte == analyte)]
        if rerun and d.rerun.any():
            d = d[d.rerun]
        else:
            d = d[~d.rerun]
        return float(d.conc.iloc[0])

    finals = {r: raw_val(s, "Xylitol") for r, s in last.items()}
    derived = set()                                   # numbers re-derived here (count as supported)

    def add(*vals, nd=(0, 1, 2)):
        for v in vals:
            for k in nd:
                derived.add(f"{v:.{k}f}")

    # ---- B: re-derive key answers ------------------------------------------------
    # baseline mean/sd from raw HPLC
    b = [finals[r] for r in ("FB-001", "FB-002", "FB-003")]
    ok(abs(np.mean(b) - Q["L3-03"]["value"]) < 0.1, f"L3-03 baseline mean raw {np.mean(b):.2f} vs gold {Q['L3-03']['value']:.2f}")
    add(np.mean(b), np.std(b, ddof=1), *b)
    # ranking of P5 by raw HPLC
    p5 = ["FB-012", "FB-013", "FB-014", "FB-015", "FB-016", "FB-017"]
    rank = sorted(p5, key=finals.get, reverse=True)
    ok(rank == Q["L3-02"]["value"], f"L3-02 ranking from raw HPLC {rank}")
    add(*[finals[r] for r in p5])
    # campaign maximum
    top = sorted(finals, key=finals.get, reverse=True)
    ok(top[0] in Q["L5-06"]["key_points"] and top[1] in Q["L5-06"]["gold_answer"], f"L5-06 top runs {top[:3]}")
    add(*[finals[r] for r in top[:4]])
    # P3 means
    grp = {"YX-001": ["FB-004", "FB-010"], "YX-002": ["FB-008"], "YX-003": ["FB-006", "FB-009"], "YX-004": ["FB-007", "FB-011"]}
    mean = {k: np.mean([finals[r] for r in v]) for k, v in grp.items()}
    best = max([k for k in grp if k != "YX-001"], key=mean.get)
    gain = 100 * (mean[best] / mean["YX-001"] - 1)
    ok(best in Q["L5-01"]["key_points"] and abs(gain - Q["L5-01"]["value"]) < 0.2, f"L5-01 best {best} {gain:+.2f} %")
    add(*mean.values(), gain, 100 * (mean[best] / finals["FB-004"] - 1), mean[best] - mean["YX-001"], finals["FB-010"],
        finals["FB-004"])
    # FB-022 vs FB-023
    d23 = finals["FB-023"] - finals["FB-022"]
    ok(abs(d23 - Q["L3-07"]["value"]) <= 0.15, f"L3-07 enzyme difference raw {d23:+.2f} vs gold {Q['L3-07']['value']:+.1f}")
    add(finals["FB-022"], finals["FB-023"], d23)
    # FB005-S07 original vs re-run
    o, c = raw_val("FB005-S07", "Xylitol", rerun=False), raw_val("FB005-S07", "Xylitol")
    ok(f"{o:.2f}" in Q["L4-05"]["gold_answer"] and f"{c:.2f}" in Q["L4-05"]["gold_answer"], f"L4-05 raw {o:.2f} -> {c:.2f}")
    # feed-sample formate (FB-024 lot)
    ff = {r: raw_val(f"{r.replace('-', '')}-F", "Formate") for r in ("FB-022", "FB-024", "FB-025")}
    ok(abs(ff["FB-024"] - ff["FB-022"]) < abs(ff["FB-024"] - ff["FB-025"]) and f"{ff['FB-024']:.1f}" in Q["L4-04"]["gold_answer"],
       f"L4-04 feed formate raw {ff}")
    add(*ff.values())
    # FB-010 gap from raw timestamps
    on = pd.read_csv(corpus / "online/FB010_R4.csv", comment="#")
    gap = on.Elapsed_h.diff().max()
    ok(abs(gap - Q["L4-02"]["value"]) < 0.05, f"L4-02 gap raw {gap:.2f} h")
    # S07 swap visible in raw
    x = {r: [raw_val(f"{r}-S{i:02d}", "Xylitol") for i in (6, 7, 8)] for r in ("FB012", "FB013")}
    ok(not (x["FB012"][0] < x["FB012"][1] < x["FB012"][2]) or x["FB012"][1] > x["FB012"][2],
       f"L5-05 FB-012 S06-S08 xylitol {np.round(x['FB012'], 2)} not monotone")
    # off-gas runs = R1/R2 per LIMS
    og = sorted("FB-" + re.search(r"fb(\d+)", f).group(1) for f in glob.glob(str(corpus / "offgas/*.xlsx")))
    v12 = sorted(lims[lims.vessel.isin(["R1", "R2"])].run_id.unique())
    ok(og == v12 and og == Q["L3-05"]["key_points"][:-2], "L3-05 off-gas runs = LIMS R1/R2 runs")
    # P6 pO2 medians from raw online files (fed phase > 30 h)
    # (BioCtrl: hours from the Inoculation header; VesselView: agitation step)
    # handled by questions.py; here only the ranking: FB-021 lowest
    # FB-021 late ethanol from raw HPLC
    e = [raw_val(f"FB021-S{i:02d}", "Ethanol") for i in (8, 9, 10)]
    ok(all(v > 0.5 for v in e), f"L5-03 FB-021 late ethanol raw {np.round(e, 2)}")
    # ---- raw online logs (independent parser) ----
    import datetime as dtm

    def online(run):
        f = glob.glob(str(corpus / f"online/{run}_*.txt"))
        if f:
            L = open(f[0]).read().splitlines()
            inoc = dtm.datetime.strptime(next(l for l in L if l.startswith("Inoculation")).split(":", 1)[1].strip(),
                                         "%d %b %Y %H:%M:%S")
            i = next(k for k, l in enumerate(L) if l.startswith("Date"))
            d = pd.DataFrame([l.split("\t") for l in L[i + 1:] if l.strip()], columns=L[i].split("\t"))
            t = (pd.to_datetime(d.Date + " " + d.Time, format="%d.%m.%Y %H:%M:%S") - inoc).dt.total_seconds() / 3600
            return pd.DataFrame(dict(t=t, pO2=d["pO2[%]"].astype(float)))
        d = pd.read_csv(glob.glob(str(corpus / f"online/{run.replace('-', '')}_*.csv"))[0], comment="#")
        # inoculation = LIMS S01 time minus 0.1 h (independent of the agitation step used in questions.py)
        s01 = pd.Timestamp(lims[lims.sample_id == f"{run.replace('-', '')}-S01"].sampled_at.iloc[0]) - pd.Timedelta(hours=0.1)
        t = (pd.to_datetime(d.Timestamp) - s01).dt.total_seconds() / 3600
        return pd.DataFrame(dict(t=t, pO2=d.DO_pct))

    def near(a, b, tol):
        return abs(float(a) - float(b)) <= tol

    o5 = online("FB-005"); w = o5[(o5.t > 2) & (o5.t < 40)]
    ok(near(w.pO2.min(), Q["L2-06"]["value"], 0.05), f"L2-06 FB-005 min pO2 raw {w.pO2.min():.1f}")
    o14 = online("FB-014")
    pre, post = o14[(o14.t > 27) & (o14.t < 30)].pO2.mean(), o14[(o14.t > 31) & (o14.t < 34)].pO2.mean()
    ok(f"about {pre:.0f} %" in Q["L4-07"]["gold_answer"] and f"about {post:.0f} %" in Q["L4-07"]["gold_answer"],
       f"L4-07 FB-014 pO2 raw {pre:.1f} -> {post:.1f}")
    o9, o6 = online("FB-009"), online("FB-006")
    at = lambda o, h: o[(o.t > h - 1) & (o.t < h + 1)].pO2.median()
    ok(all(f"{v:.0f} %" in Q["L4-03"]["gold_answer"] or f"{v:.0f}-" in Q["L4-03"]["gold_answer"] or f"-{v:.0f} %" in Q["L4-03"]["gold_answer"]
           for v in (at(o9, 40), at(o9, 140), at(o6, 40), at(o6, 140))),
       f"L4-03 pO2 raw FB-009 {at(o9, 40):.0f}->{at(o9, 140):.0f}, FB-006 {at(o6, 40):.0f}->{at(o6, 140):.0f}")
    med = {r: online(r).query("t > 30").pO2.median() for r in ("FB-018", "FB-019", "FB-020", "FB-021")}
    ok(all(f"{v:.0f}" in Q["L5-03"]["gold_answer"] for v in med.values()) and min(med, key=med.get) == "FB-021",
       f"L5-03 fed-phase pO2 medians raw {', '.join(f'{k} {v:.0f}' for k, v in med.items())}")
    for r in ("FB-018", "FB-019", "FB-020", "FB-021"):
        ok(f"{finals[r]:.1f}" in Q["L5-03"]["gold_answer"], f"L5-03 final xylitol {r} raw {finals[r]:.2f}")
    # FB-012 (L5-02) from raw HPLC + LIMS
    s12 = last["FB-012"]
    ok(f"{raw_val(s12, 'Xylitol'):.1f} g/L xylitol" in Q["L5-02"]["gold_answer"] and
       f"formate {raw_val(s12, 'Formate'):.1f}" in Q["L5-02"]["gold_answer"], "L5-02 FB-012 final xylitol / formate raw")
    for r in ("FB-015", "FB-016", "FB-017"):
        ok(f"{r} {finals[r]:.1f}" in Q["L5-02"]["gold_answer"], f"L5-02 {r} raw {finals[r]:.2f}")
    # flask: parent at 70 % C5 (PN tolerance workbook, independent of the minutes)
    ws = openpyxl.load_workbook(corpus / "spreadsheets/PN_C5_tolerance_Apr2025.xlsx", data_only=True)["YX-001"]
    fl = [r for r in ws.iter_rows(min_row=3, values_only=True) if r[1]]
    first = {}
    for r in fl:
        first.setdefault(r[1], r)
    rows70 = [r for r in fl if r[1] == next(k for k in first if any("70 %" in str(x[0]) for x in fl if x[1] == k))]
    od0, od72 = rows70[0][3], rows70[-1][3]
    ok(f"OD {od72:.1f} at 72 h" in Q["L5-02"]["gold_answer"] and od72 < 1.5 * od0,
       f"L5-02 YX-001 70 % flask OD {od0} -> {od72} (no growth)")
    # FB-002 slow start (L4-08) from raw sources: OD/DW in workbooks, all three runs
    # FB-005 feed: EV note and DEV agree
    n5 = (corpus / "notes/FB-005_EV.md").read_text()
    m = re.search(r"about (\d+) mL by balance vs (\d+) mL by pump", n5)
    ok(m and m.group(1) in Q["L4-01"]["gold_answer"] and m.group(2) in Q["L4-01"]["gold_answer"], "L4-01 EV note = DEV numbers")
    # glucose LOD time FB-003 from raw
    g3 = [(float(broth.set_index("sample_id").loc[sid, "process_time_h"]), raw_val(sid, "Glucose"))
          for sid in broth[broth.run_id == "FB-003"].sample_id]
    t_lod = min(t for t, v in g3 if np.isnan(v) and t > 0.2)
    ok(t_lod == Q["L2-01"]["value"], f"L2-01 FB-003 glucose <LOD first at {t_lod} h (raw)")
    # MD transcription slip
    ws = openpyxl.load_workbook(corpus / "spreadsheets/MD_fermentation_master_2025.xlsx", data_only=True)["Samples"]
    rows = list(ws.iter_rows(values_only=True))
    md = pd.DataFrame(rows[1:], columns=rows[0])
    slip, stale = [], []
    for _, r in md.iterrows():
        sid = f"{r.Run.replace('FB', 'FB0')}-{r.Smpl}"
        if not isinstance(r.Xol, (int, float)):
            continue
        d = raw[(raw.SampleID == sid) & (raw.Analyte == "Xylitol")]
        first, rer = d[~d.rerun].conc.iloc[0], (d[d.rerun].conc.iloc[0] if d.rerun.any() else None)
        if abs(r.Xol / first - 10) < 0.05:
            slip.append(sid)
        elif rer is not None and abs(r.Xol - first) < 0.011:
            stale.append(sid)
    g = Q["L4-06"]["gold_answer"]
    ok(len(slip) == 1 and f"{slip[0][:5].replace('FB0', 'FB-0')} {slip[0][-3:]}" in g, f"L4-06 slip {slip} (raw)")
    ok(f"{len(stale)} values" in g and all(x[-3:] in g for x in stale), f"L4-06 stale re-run values {stale} (raw)")

    # ---- U: the information really is absent ---------------------------------------
    alltext = file_text(corpus, "corpus").lower()
    ok("dissolved co2" not in alltext and "pco2" not in alltext, "U-01 no dissolved-CO2 measurement anywhere")
    ok(raw[(raw.SampleID.str.startswith("FB006")) & (raw.Analyte == "Succinate")].empty, "U-02 no succinate for FB-006")
    ok(not glob.glob(str(corpus / "offgas/*fb006*")), "U-03 no off-gas for FB-006")
    ok("fb-026" not in alltext and "fb026" not in alltext, "U-04 FB-026 never mentioned")
    ok(not re.search(r"\b(price|cost|eur|usd|invoice)\b|[$\u20ac\u00a3]", alltext), "U-06 no prices/costs anywhere")
    reg = openpyxl.load_workbook(corpus / "strains/strain_registry.xlsx", data_only=True)["Strains"]
    st = {r[0]: r[9] for r in reg.iter_rows(min_row=2, values_only=True)}
    ok(not str(st["YX-004"]).startswith("dropped"), f"U-05 YX-004 status '{st['YX-004']}'")
    lots = pd.read_csv(key / "experiments.csv")
    ok(set(lots[(lots.lot_actual == "C5-2503")].enzyme_h.notna()) == {True}, "U-07 every true C5-2503 run had enzyme")

    # ---- C: faults referenced -----------------------------------------------------
    faults = pd.read_csv(key / "faults.csv")
    codes = dict(zip(faults.code, faults.run))
    ok(codes.get("F-FEEDCAL") == "FB-005", "fault F-FEEDCAL is FB-005 (L4-01, L5-06)")
    ok(codes.get("F-POWER") == "FB-010", "fault F-POWER is FB-010 (L4-02, L5-01)")
    ok(codes.get("F-DODRIFT") == "FB-009", "fault F-DODRIFT is FB-009 (L4-03)")
    ok(codes.get("F-WRONGLOT") == "FB-024", "fault F-WRONGLOT is FB-024 (L4-04, L5-04, U-07)")
    ok(codes.get("F-ANTIFOAM") == "FB-014", "fault F-ANTIFOAM is FB-014 (L4-07)")
    ok(codes.get("F-INOC") == "FB-002", "fault F-INOC is FB-002 (L4-08)")
    ok("F-SWAP" in codes, "fault F-SWAP exists (L5-05)")

    # ---- IO / CA / XD: independent re-derivations -----------------------------------
    import openpyxl as ox

    def nearest(o, h):
        return float(o.pO2[(o.t - h).abs().idxmin()])
    # IO-04: VesselView FB-021 at 47.5 h, inoculation from the agitation step (route differs from questions_extra)
    vv = pd.read_csv(glob.glob(str(corpus / "online/FB021_*.csv"))[0], comment="#")
    k0 = vv.index[vv.AGIT_rpm > 400][0]
    vv["t"] = vv.Elapsed_h - vv.Elapsed_h[k0]
    v_ag = float(vv.DO_pct[(vv.t - 47.5).abs().idxmin()])
    ok(abs(v_ag - Q["IO-04"]["value"]) <= Q["IO-04"]["tolerance"], f"IO-04 FB-021 pO2 at 47.5 h via agitation step {v_ag:.1f}")
    # IO-05: CO2 peak, inoculation from LIMS (not the BioCtrl header)
    rows = list(ox.load_workbook(glob.glob(str(corpus / "offgas/GA4_fb001_*.xlsx"))[0], data_only=True)["GA4"].iter_rows(values_only=True))
    g = pd.DataFrame(rows[5:], columns=rows[4])
    s01 = pd.Timestamp(lims[lims.sample_id == "FB001-S01"].sampled_at.iloc[0]) - pd.Timedelta(hours=0.1)
    tpk = (pd.Timestamp(g.Date[g["CO2 [%]"].idxmax()] + " " + g.Time[g["CO2 [%]"].idxmax()]) - s01).total_seconds() / 3600
    ok(abs(tpk - Q["IO-05"]["value"]) <= 0.5, f"IO-05 CO2 peak {tpk:.1f} h after inoculation (LIMS route)")
    add(g["RelTime[s]"][g["CO2 [%]"].idxmax()] / 3600)
    derived.add("3600")
    # IO-01 / IO-07 / IO-06 from raw HPLC
    ok(abs(raw_val("FB006-S08", "Xylitol") - Q["IO-01"]["value"]) < 0.02, "IO-01 FB006-S08 xylitol raw")
    ok(abs(raw_val("FB006-S05", "Xylitol") - Q["IO-07"]["value"]) < 0.02 and
       abs(raw_val("FB006-S05", "Xylitol", rerun=False) - Q["IO-07"]["value"]) > 1, "IO-07 re-run value raw, differs from original")
    ok(abs(raw_val(last["FB-016"], "Arabinose") - Q["IO-06"]["value"]) < 0.02, "IO-06 FB-016 final arabinose raw")
    ok(np.isnan(raw_val("FB002-S04", "Glycerol")), "IO-02 FB002-S04 glycerol <LOD in raw")
    ok(str(lims[lims.sample_id == "FB002-S01"].sampled_at.iloc[0]).startswith(Q["IO-03"]["value"]), "IO-03 date = LIMS")
    # IO-10 base, CA-04 feed pump: own parser for BioCtrl
    def bio_last(run):
        L_ = [l for l in open(glob.glob(str(corpus / f"online/{run}_*.txt"))[0]).read().splitlines() if l and not l.startswith("#")]
        h = L_[[i for i, l in enumerate(L_) if l.startswith("Date")][0]].split("\t")
        return dict(zip(h, L_[-1].split("\t")))
    setup = (corpus / "equipment/reactor_setup.md").read_text()
    mlc = {k: float(re.search(rf"{k} pump: ([\d.]+)", setup).group(1)) for k in ("Base", "Feed")}
    b2 = int(bio_last("FB-002")["Base_Pump[cnt]"]) * mlc["Base"]
    b6 = float(pd.read_csv(glob.glob(str(corpus / "online/FB006_*.csv"))[0], comment="#").BASE_mL.iloc[-1])
    ok(abs(b2 - Q["IO-10"]["value"][0]) < 0.05 and abs(b6 - Q["IO-10"]["value"][1]) < 0.01, f"IO-10 base {b2:.1f} / {b6:.2f} mL")
    sop = (corpus / "sops/SOP-FB-001_v3.md").read_text()
    dA = float(re.search(r"Feed A ([\d.]+) g/mL", sop).group(1))
    v6 = pd.read_csv(glob.glob(str(corpus / "online/FB006_*.csv"))[0], comment="#")
    f6 = (v6.FEED_BAL_g.max() - v6.FEED_BAL_g.iloc[-1]) / dA
    ok(abs(f6 - Q["CA-03"]["value"]) <= Q["CA-03"]["tolerance"], f"CA-03 FB-006 feed {f6:.1f} mL")
    cw = {r[0]: r[1] for r in ox.load_workbook(corpus / "spreadsheets/EV_FB-005.xlsx", data_only=True)["Calculated"].iter_rows(values_only=True)}
    over = 100 * ((cw["Feed bottle start (g)"] - cw["Feed bottle end (g)"]) / dA / (int(bio_last("FB-005")["Feed_Pump[cnt]"]) * mlc["Feed"]) - 1)
    ok(abs(over - Q["CA-04"]["value"]) < 0.5, f"CA-04 over-delivery {over:.1f} %")
    # CA-01 growth rate from OD instead of DW
    w1 = ox.load_workbook(corpus / "spreadsheets/EV_FB-001.xlsx", data_only=True)["Offline data"]
    r1 = [r for r in w1.iter_rows(values_only=True)][1:3]
    mu_od = np.log(r1[1][4] / r1[0][4]) / (r1[1][1] - r1[0][1])
    ok(abs(mu_od - Q["CA-01"]["value"]) <= Q["CA-01"]["tolerance"], f"CA-01 mu from OD {mu_od:.3f} vs DW-based gold")
    ok(abs(raw_val("FB018-S06", "Xylitol") - Q["CA-02"]["value"]) < 0.02, "CA-02 raw recompute")
    # CA-06 RQ with the measured N2+Ar as inert (alternative to 100 - O2 - CO2)
    rows = list(ox.load_workbook(glob.glob(str(corpus / "offgas/GA4_fb018_*.xlsx"))[0], data_only=True)["GA4"].iter_rows(values_only=True))
    g = pd.DataFrame(rows[5:], columns=rows[4])
    s01 = pd.Timestamp(lims[lims.sample_id == "FB018-S01"].sampled_at.iloc[0]) - pd.Timedelta(hours=0.1)
    g["t"] = (pd.to_datetime(g.Date + " " + g.Time) - s01).dt.total_seconds() / 3600
    inert_out = (g["N2 [%]"] + g["Ar [%]"]) / 100
    f_ = (1 - 0.2095 - 0.0004) / inert_out
    rq = ((g["CO2 [%]"] / 100 * f_ - 0.0004) / (0.2095 - g["O2 [%]"] / 100 * f_))[g.t > 30].median()
    ok(abs(rq - Q["CA-06"]["value"]) <= Q["CA-06"]["tolerance"], f"CA-06 RQ with measured inert {rq:.2f}")
    # CA-07 undissociated formic acid, formate from raw HPLC
    und = raw_val(last["FB-012"], "Formate") / (1 + 10 ** 1.75)
    ok(abs(und - Q["CA-07"]["value"]) < 0.005, f"CA-07 undissociated formic acid {und:.3f} g/L")
    # CA-08 yield vs simulator truth
    ts = pd.read_csv(key / "timeseries/FB-018.csv")
    sm = pd.read_csv(key / "samples_truth.csv")
    s18 = sm[(sm.run_id == "FB-018") & (sm.kind == "broth")]
    fed_x = ts.cumFeed.iloc[-1] * 300.0
    left = ts.Xyl.iloc[-1] * ts.V.iloc[-1] + (s18.true_Xyl * 0.007).iloc[:-1].sum()
    prod_ = ts.Xol.iloc[-1] * ts.V.iloc[-1] + (s18.true_Xol * 0.007).iloc[:-1].sum()
    y_true = prod_ / (fed_x - left)
    ok(abs(y_true - Q["CA-08"]["value"]) <= Q["CA-08"]["tolerance"], f"CA-08 measured yield {Q['CA-08']['value']} vs truth {y_true:.3f}")
    # CA-09 residual xylose from raw
    rx = {k: np.mean([raw_val(last[r], "Xylose") for r in v]) for k, v in grp.items()}
    ok(min(rx, key=rx.get) == Q["CA-09"]["key_points"][0], f"CA-09 lowest residual xylose {min(rx, key=rx.get)} (raw)")
    sd_b = np.std(b, ddof=1)
    ok(abs(np.ceil(2 * (2.8 * sd_b / 5) ** 2) - Q["CA-10"]["value"]) <= 1, "CA-10 power calculation (raw baseline sd)")
    # XD: simulator checks present and rubric evidence runs exist
    dcj = json.loads((key / "design_checks.json").read_text())
    for qid in ("XD-01", "XD-02", "XD-03"):
        rb = Q[qid]["rubric"]
        ok(all("not available" not in r["simulator_check"] for r in rb), f"{qid}: simulator checks filled")
        runs_cited = set(re.findall(r"FB-0\d\d", " ".join(r["evidence"] for r in rb)))
        ok(runs_cited <= set(finals), f"{qid}: cited runs exist")
    ok(dcj["strains_feedA_3mLh"]["YX-003"]["Xol_end"] - dcj["strains_feedA_3mLh"]["YX-001"]["Xol_end"] >
       dcj["strains_feedA_2mLh"]["YX-003"]["Xol_end"] - dcj["strains_feedA_2mLh"]["YX-001"]["Xol_end"],
       "XD R3/D2: higher feed rate enlarges the strain difference (simulator)")

    # ---- A: every quoted number is supported (core and IO sets) ---------------------
    for qid, q in Q.items():
        if q.get("set", "core") not in ("core", "IO"):
            continue
        ev = "\n".join(file_text(corpus, e) for e in q["evidence"])
        evn = set(numbers(ev))
        # numbers in evidence may be written with comma decimals (MD notes) or more decimals
        evn |= {n.replace(",", ".") for n in re.findall(r"\d+,\d+", ev)}
        evf = {round(float(n), k) for n in evn for k in (0, 1, 2)}
        miss = []
        for n in numbers(q["gold_answer"]):
            if n in evn or n in derived or round(float(n), 2) in evf:
                continue
            if re.fullmatch(r"20\d\d", n) or float(n) in (0, 1, 2, 3, 4, 5, 10, 60, 70, 50, 30, 900, 500, 15, 3.5):
                continue                               # years, small counts, design constants
            miss.append(n)
        ok(not miss, f"{qid}: numbers supported by evidence" + (f" - unsupported: {miss}" if miss else ""))

    # ---- D: measured gold vs truth (information) -------------------------------------
    exp = pd.read_csv(key / "experiments.csv").set_index("run_id")
    dev = {r: 100 * (finals[r] / exp.loc[r, "final_Xol_true"] - 1) for r in finals}
    info = (f"measured vs true final xylitol: median |error| {np.median([abs(v) for r, v in dev.items() if r.startswith('FB')]):.1f} % (reactors), "
            f"max {max(dev, key=lambda r: abs(dev[r]))} {dev[max(dev, key=lambda r: abs(dev[r]))]:+.1f} %")
    tm = {k: np.mean([exp.loc[r, 'final_Xol_true'] for r in v]) for k, v in grp.items()}
    info2 = "true P3 means: " + ", ".join(f"{k} {v:.1f}" for k, v in tm.items())

    n = {k: sum(1 for r in OUT if r[0] == k) for k in ("PASS", "FAIL")}
    for s, m in OUT:
        if s != "PASS":
            print(s, m)
    print(f"\n{n['PASS']} PASS, {n['FAIL']} FAIL")
    print("INFO", info)
    print("INFO", info2)
    return n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default=str(Path(__file__).resolve().parents[1] / "pilot"))
    main(Path(ap.parse_args().pilot))
