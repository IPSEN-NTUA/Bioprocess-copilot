"""
Extra question sets for the pilot corpus:

  IO  semantic interoperability: the answer needs heterogeneous sources to be reconciled
      (identifiers, below-LOD encodings, number/date formats, time bases, vocabularies, units,
      superseded values, implicit context)
  CA  calculations on the recorded data (unit conversion, gas balance, mass balance, statistics)
  XD  experiment design: what to do next, graded with a rubric of corpus-grounded suggestions

Like questions.py, every gold answer is computed from the corpus files only. The XD rubric items
also carry a `simulator_check` (from answer_key/design_checks.json, produced by design_checks.py)
that says whether the suggested experiment would really have been informative - for analysis only,
not for grading.
"""
from __future__ import annotations

import datetime as dtm
import glob
import json
import math
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

MD = "corpus/spreadsheets/MD_fermentation_master_2025.xlsx"
REG = "corpus/strains/strain_registry.xlsx"
SETUP = "corpus/equipment/reactor_setup.md"
SOPFB = "corpus/sops/SOP-FB-001_v3.md"

# interoperability tags used across all question sets
INTEROP_TAGS = {
    "id_resolution": "the same entity is named differently across sources (FB02 / FB-002 / FB002-S01, YX-002-c1 / 'YX-002 (clone c3)')",
    "lod_encoding": "below-detection results are written as '<0.05', '<LOD', 'n.d.' or a blank cell",
    "number_date_format": "comma decimals and dd/mm/yy or dd.mm.yyyy dates in some sources",
    "time_alignment": "process time must be reconstructed from wall-clock stamps, logger start or analyser start",
    "vocabulary": "different names for the same quantity (Xol / Xylitol / 'Xylitol (g/L)', formic acid / formate)",
    "units_conversion": "pump counts, balance grams and mL must be converted with documented factors",
    "versioning": "superseded results, method versions and corrected values",
    "source_authority": "sources disagree; the answer must pick the authoritative one and say why",
    "implicit_context": "a fact is not written in the record and must be inferred from another document",
}

CORE_TAGS = {
    "L1-01": ["id_resolution"], "L2-01": ["lod_encoding"], "L2-06": ["time_alignment"],
    "L3-01": ["id_resolution", "implicit_context"], "L3-02": ["id_resolution"], "L3-03": ["id_resolution"],
    "L3-06": ["id_resolution"], "L3-08": ["number_date_format"], "L4-02": ["time_alignment"],
    "L4-03": ["time_alignment"], "L4-04": ["vocabulary", "source_authority"], "L4-05": ["versioning"],
    "L4-06": ["versioning", "source_authority", "number_date_format"], "L4-07": ["time_alignment"],
    "L4-08": ["number_date_format"], "L5-01": ["versioning"], "L5-05": ["id_resolution"],
    "L5-06": ["versioning", "source_authority"], "U-02": ["versioning"],
}


def offgas(C, run):
    f = glob.glob(str(C.p(f"offgas/GA4_fb{run[3:]}_*.xlsx")))[0]
    rows = list(openpyxl.load_workbook(f, data_only=True)["GA4"].iter_rows(values_only=True))
    g = pd.DataFrame(rows[5:], columns=rows[4])
    g["ts"] = pd.to_datetime(g.Date + " " + g.Time)
    return g, rows, Path(f).relative_to(C.root.parent)


def bioctrl_inoculation(C, run):
    t = open(glob.glob(str(C.p(f"online/{run}_*.txt")))[0]).read()
    return pd.Timestamp(dtm.datetime.strptime(re.search(r"Inoculation\t: (.+)", t).group(1).strip(), "%d %b %Y %H:%M:%S"))


def bioctrl_last(C, run):
    lines = [l for l in open(glob.glob(str(C.p(f"online/{run}_*.txt")))[0]) if not l.startswith("#")]
    i = next(k for k, l in enumerate(lines) if l.startswith("Date"))
    hdr = lines[i].strip().split("\t")
    return dict(zip(hdr, lines[-1].strip().split("\t")))


def vesselview(C, run):
    f = glob.glob(str(C.p(f"online/{run.replace('-', '')}_*.csv")))[0]
    return pd.read_csv(f, comment="#"), Path(f).relative_to(C.root.parent)


def build_extra(C, Q) -> None:
    lims = C.lims()
    L = lims.set_index("sample_id")
    tr = C.tr()
    mds, mdr = C.md("Samples"), C.md("Runs")
    reg = C.registry()
    notes = C.notes()
    setup = C.text("equipment/reactor_setup.md")
    sop = C.text("sops/SOP-FB-001_v3.md")
    ml_per = {k: float(re.search(rf"{k} pump: ([\d.]+) mL per count", setup).group(1)) for k in ("Acid", "Base", "Feed")}
    dens = {k: float(re.search(rf"{k} ([\d.]+) g/mL", sop).group(1)) for k in ("Feed A", "Feed B", "C5 feeds")}
    process_t = lambda sid: float(L.loc[sid, "process_time_h"])
    inoc_lims = lambda run: pd.Timestamp(L.loc[f"{run.replace('-', '')}-S01", "sampled_at"]) - pd.Timedelta(hours=0.1)

    def q(set_, level, category, question, answer, evidence, grading, value=None, tol=None, key_points=None,
          hops=1, tags=None, trap=None, rubric=None, penalties=None):
        n = sum(1 for x in Q if x["id"].startswith(set_ + "-")) + 1
        Q.append(dict(id=f"{set_}-{n:02d}", set=set_, tier=level, category=category, question=question,
                      gold_answer=answer, value=value, tolerance=tol, grading=grading, key_points=key_points or [],
                      evidence=evidence, hops=hops, interop=tags or [], trap=trap, rubric=rubric, penalties=penalties))

    # =================== IO: semantic interoperability ===================
    v6 = C.tr_value("FB006-S08", "Xylitol")
    m6 = mds[(mds.Run == "FB06") & (mds.Smpl == "S08")].Xol.iloc[0]
    assert process_t("FB006-S08") == 96.0 and abs(v6 - m6) < 0.011
    q("IO", "L2", "lookup", "What xylitol concentration was measured in FB-006 at 96 h?",
      f"{v6:.2f} g/L (sample FB006-S08; written as run 'FB06', sample 'S08', column 'Xol' in the MD master log and as "
      "'FB006-S08' / 'Xylitol' in the HPLC results).", [MD, "corpus/hplc/processed/HPLC_results_2025-03.xlsx",
                                                          "corpus/lims/sample_register.csv"],
      "numeric", value=round(v6, 2), tol=0.02, hops=2, tags=["id_resolution", "vocabulary"])

    g2 = mds[(mds.Run == "FB02") & (mds.Smpl == "S04")].Gly.iloc[0]
    t2 = tr[(tr.SampleID == "FB002-S04") & (tr.Analyte == "Glycerol")].iloc[0]
    assert g2 == "<0.05" and t2.Flag == "<LOD" and process_t("FB002-S04") == 23.5
    q("IO", "L2", "lookup", "What glycerol concentration was measured in FB-002 at 23.5 h?",
      "Below the detection limit: the MD master log shows '<0.05' and the HPLC results flag it '<LOD' (LOD 0.02-0.05 g/L in "
      "the vial). It was measured, so it is not missing, and it should not be reported as 0 or as a number.",
      [MD, "corpus/hplc/processed/HPLC_results_2025-02.xlsx", "corpus/sops/SOP-ANA-003_v1.md"], "key_points",
      key_points=["below detection", "<0.05"], hops=2, tags=["lod_encoding", "id_resolution"],
      trap="0 g/L, or 'not measured'")

    dmd = re.search(r"- (\d\d/\d\d/\d\d) inoc (\d\d:\d\d)", notes["FB-002_MD"]).groups()
    dl = L.loc["FB002-S01", "sampled_at"]
    d_iso = dtm.datetime.strptime(dmd[0], "%d/%m/%y").date().isoformat()
    assert dl.startswith(d_iso)
    q("IO", "L2", "lookup", "On which date and at what time was FB-002 inoculated?",
      f"{dtm.date.fromisoformat(d_iso):%d %B %Y}, {dmd[1]} ({d_iso}). MD writes dates as dd/mm/yy ('{dmd[0]}' in his note "
      "and master log), which must not be read as month/day; the LIMS register confirms the date.",
      ["corpus/notes/FB-002_MD.md", MD, "corpus/lims/sample_register.csv"], "key_points",
      key_points=[d_iso, dmd[1]], value=d_iso, hops=2, tags=["number_date_format"],
      trap=f"{dtm.datetime.strptime(dmd[0], '%m/%d/%y'):%d %B %Y} (month/day read the wrong way round)")

    vv, vp = vesselview(C, "FB-021")
    vv["t"] = (pd.to_datetime(vv.Timestamp) - inoc_lims("FB-021")).dt.total_seconds() / 3600
    i, j = (vv.t - 47.5).abs().idxmin(), (vv.Elapsed_h - 47.5).abs().idxmin()
    near = vv[(vv.t > 47.5 - 0.2) & (vv.t < 47.5 + 0.2)].DO_pct
    lead = float(vv.t.iloc[0])
    q("IO", "L3", "time_series", "What was the dissolved-oxygen (pO2) reading in FB-021 at 47.5 h after inoculation?",
      f"About {vv.DO_pct[i]:.1f} % (record {vv.Timestamp[i]}). In the VesselView export Elapsed_h counts from the logger "
      f"start, which was {-lead:.1f} h before inoculation, so the row with Elapsed_h = 47.5 is actually "
      f"{vv.t[j]:.1f} h after inoculation (pO2 {vv.DO_pct[j]:.1f} %). Inoculation time from the LIMS S01 sample / run sheet.",
      [str(vp), "corpus/lims/sample_register.csv", ev_path_("FB-021"), SETUP], "numeric", value=round(float(vv.DO_pct[i]), 1),
      tol=round(max(1.0, float(near.max() - near.min())), 1), hops=3, tags=["time_alignment"],
      trap=f"{vv.DO_pct[j]:.1f} % (row with Elapsed_h = 47.5)")

    g, rows, gp = offgas(C, "FB-001")
    ino = bioctrl_inoculation(C, "FB-001")
    g["t"] = (g.ts - ino).dt.total_seconds() / 3600
    k = g["CO2 [%]"].idxmax()
    q("IO", "L3", "time_series", "How many hours after inoculation did the exhaust CO2 of FB-001 peak, and what was the peak value?",
      f"At about {g.t[k]:.1f} h after inoculation ({g.ts[k]:%d.%m.%Y %H:%M}), CO2 {g['CO2 [%]'][k]:.2f} %, during the "
      f"glucose batch phase. Inoculation is in the BioCtrl header ({ino:%d %b %Y %H:%M}); the analyser RelTime counts "
      f"from the analyser start ({rows[2][1]}), so RelTime/3600 = {g['RelTime[s]'][k] / 3600:.1f} h is not process time.",
      [str(gp), "corpus/online/FB-001_R1.txt", SETUP], "key_points", value=round(float(g.t[k]), 1), tol=0.5,
      key_points=[f"{g.t[k]:.1f} h", f"{g['CO2 [%]'][k]:.2f} %"], hops=3, tags=["time_alignment"],
      trap=f"{g['RelTime[s]'][k] / 3600:.1f} h (analyser RelTime)")

    a16 = C.final("FB-016", "Arabinose")
    assert "Ara" not in mds.columns
    q("IO", "L3", "lookup", "What was the arabinose concentration at the end of FB-016?",
      f"{a16:.2f} g/L (FB016-{C.broth_samples('FB-016')[-1].split('-')[1]}, HPLC results). The MD master log has no "
      "arabinose column, so the value is only in the HPLC results workbook; arabinose comes in with the C5 fraction and is "
      "not consumed.", ["corpus/hplc/processed/HPLC_results_2025-05.xlsx", MD], "numeric", value=round(a16, 2), tol=0.02,
      hops=2, tags=["source_authority", "id_resolution"], trap="'not measured' (because the MD log lacks the column)")

    c5 = C.tr_value("FB006-S05", "Xylitol")
    o5 = C.tr_value("FB006-S05", "Xylitol", rerun=False)
    m5 = mds[(mds.Run == "FB06") & (mds.Smpl == "S05")].Xol.iloc[0]
    assert abs(m5 - o5) < 0.011 and abs(c5 - o5) > 1
    q("IO", "L3", "provenance", "What is the valid xylitol concentration of sample FB006-S05?",
      f"{c5:.2f} g/L from the re-run sequence SEQ-20250409. The first result ({o5:.2f} g/L, SEQ-20250321) was superseded "
      f"after DEV-2025-02 (biased xylitol standard); the MD master log still shows the old value ({m5}).",
      ["corpus/hplc/processed/HPLC_results_2025-03.xlsx", "corpus/deviations/DEV-2025-02.md", MD], "numeric",
      value=round(c5, 2), tol=0.02, hops=3, tags=["versioning", "source_authority"], trap=f"{o5:.2f} g/L (superseded)")

    p3 = re.search(r"\*\*Parent:\*\* (\S+)", C.text("strains/construction_YX-003.md")).group(1)
    ev5 = dict(r[:2] for r in C.xlsx_sheet("spreadsheets/EV_FB-005.xlsx", "Run info"))["Strain"]
    ev8 = dict(r[:2] for r in C.xlsx_sheet("spreadsheets/EV_FB-008.xlsx", "Run info"))["Strain"]
    c58 = re.search(r"clone (c\d)", ev5).group(1)
    assert c58 == re.search(r"clone (c\d)", ev8).group(1) and p3.endswith("c1") and c58 != "c1"
    q("IO", "L3", "provenance", "Was YX-003 built on the same YX-002 clone that was tested in the bioreactor runs FB-005 and FB-008?",
      f"No. YX-003 was built on {p3} (construction record, registry parent), while FB-005 and FB-008 used YX-002 clone "
      f"{c58} (EV workbooks: '{ev5}'). The same clone is written 'YX-002-c1', 'YX-002 c1' and 'clone c1' in different documents.",
      ["corpus/strains/construction_YX-003.md", REG, "corpus/spreadsheets/EV_FB-005.xlsx", "corpus/spreadsheets/EV_FB-008.xlsx"],
      "key_points", key_points=["no", "c1", c58], hops=3, tags=["id_resolution"])

    dw16 = re.search(r"DW (\d+,\d)", notes["FB-016_MD"].split("end")[1]).group(1)
    dwl = float(mds[mds.Run == "FB16"].DW.iloc[-1])
    assert abs(float(dw16.replace(",", ".")) - dwl) < 0.06
    q("IO", "L2", "lookup", "What final dry weight did Mark Doherty write in his FB-016 note, and does it agree with his master log?",
      f"{dw16.replace(',', '.')} g/L - written '{dw16}' with a decimal comma in the note; the master log gives {dwl:.2f} g/L, "
      "so they agree.", ["corpus/notes/FB-016_MD.md", MD], "key_points",
      key_points=[dw16.replace(",", "."), "agree"], value=float(dw16.replace(",", ".")), tol=0.05, hops=2,
      tags=["number_date_format"], trap="33 g/L (comma read as a thousands separator)")

    last2 = bioctrl_last(C, "FB-002")
    b2 = int(last2["Base_Pump[cnt]"]) * ml_per["Base"]
    vv6, vp6 = vesselview(C, "FB-006")
    b6 = float(vv6.BASE_mL.iloc[-1])
    q("IO", "L3", "calculation", "How much base (2 M KOH) was added in FB-002 and in FB-006?",
      f"FB-002 (R2, BioCtrl): {int(last2['Base_Pump[cnt]'])} base-pump counts x {ml_per['Base']:.2f} mL/count = {b2:.1f} mL. "
      f"FB-006 (R3, VesselView): {b6:.2f} mL, logged directly in mL (BASE_mL). The two controllers log base differently; "
      "the pump calibration is in the reactor set-up document.",
      ["corpus/online/FB-002_R2.txt", str(vp6), SETUP], "key_points", key_points=[f"{b2:.1f}", f"{b6:.2f}"],
      value=[round(b2, 1), round(b6, 2)], tol=0.2, hops=3, tags=["units_conversion", "vocabulary"],
      trap=f"{int(last2['Base_Pump[cnt]'])} mL (counts taken as mL)")

    an = tr[(tr.SampleID == "FB015-S08") & (tr.Analyte == "Xylitol")].iloc[0]
    q("IO", "L4", "provenance", "The MD master log has a row 'FB15 / S08'. Which LIMS sample is this, when was it taken, and in "
      "which HPLC sequence was it analysed?",
      f"LIMS sample FB015-S08 (run FB-015, {process_t('FB015-S08'):g} h), sampled {L.loc['FB015-S08', 'sampled_at']}, "
      f"analysed on {an.Analysed} in sequence SEQ-{str(an.Analysed).replace('-', '')} (method {an.SOP}). This is also the "
      "row with the decimal-point slip in the master log.",
      [MD, "corpus/lims/sample_register.csv", "corpus/hplc/processed/HPLC_results_2025-05.xlsx", "corpus/hplc/raw/"],
      "key_points", key_points=["FB015-S08", str(L.loc["FB015-S08", "sampled_at"])[:10],
                                f"SEQ-{str(an.Analysed).replace('-', '')}"], hops=3, tags=["id_resolution"])

    t4 = C.text("notes/SF_screen_P4_PN.md")
    c6 = reg.loc["YX-006", "Selected clone"]
    lab = f"YX-006 {c6}, 70 % C5"
    x = float(re.search(rf"\| {re.escape(lab)} \| ([\d.]+) \|", t4).group(1))
    rows6 = [r for r in C.xlsx_sheet("spreadsheets/PN_C5_tolerance_Apr2025.xlsx", "YX-006")[2:] if r[1]]
    fl = {}
    for r in rows6:
        fl.setdefault(r[1], []).append(r)
    sf = [k for k, v in fl.items() if any("70 %" in str(rr[0]) for rr in v)]
    assert len(sf) == 1
    sf = sf[0]
    od72 = [r for r in fl[sf] if r[2] == 72][0][3]
    x72 = [r for r in fl[sf] if r[2] == 72][0][6]
    assert abs(x72 - x) < 0.06
    q("IO", "L3", "provenance", f"PN's April flask note lists '{lab}' with {x:.1f} g/L xylitol. Which flask is that, and what OD600 "
      "did it reach at 72 h?",
      f"Flask {sf} (sheet YX-006 of PN_C5_tolerance_Apr2025.xlsx, 'Tolerance, 70 % C5'); OD600 {od72} at 72 h. The minutes "
      f"of 2025-05-02 refer to the same flask by its ID.",
      ["corpus/notes/SF_screen_P4_PN.md", "corpus/spreadsheets/PN_C5_tolerance_Apr2025.xlsx", "corpus/meetings/minutes_2025-05-02.md"],
      "key_points", key_points=[sf, f"{od72}"], value=od72, tol=0.1, hops=2, tags=["id_resolution"])

    s13 = mdr[mdr.Run == "FB13"].Strain.iloc[0]
    sel5 = reg.loc["YX-005", "Selected clone"]
    pick = re.search(r"Clone (c\d) selected after the flask screen on (\S+)\.", C.text("strains/construction_YX-005.md"))
    assert s13 == "YX-005" and pick.group(1) == sel5
    start13 = dtm.datetime.strptime(mdr[mdr.Run == "FB13"].Start.iloc[0], "%d/%m/%y").date()
    assert dtm.date.fromisoformat(pick.group(2)) < start13
    q("IO", "L4", "inference", "Which clone of YX-005 was used in FB-013?",
      f"Clone {sel5}. MD's records only say 'YX-005'; the clone follows from the strain registry / construction record "
      f"(clone {sel5} selected after the flask screen on {pick.group(2)}), which was before FB-013 started on {start13}.",
      [MD, "corpus/notes/FB-013_MD.md", REG, "corpus/strains/construction_YX-005.md"], "key_points",
      key_points=[sel5], value=sel5, hops=3, tags=["implicit_context", "number_date_format"])

    # =================== CA: calculations ===================
    ev1 = C.ev("FB-001")[0]
    b1 = ev1[ev1.Sample.str.contains("-S")].iloc[:2]
    mu = math.log(b1["DW (g/L)"].iloc[1] / b1["DW (g/L)"].iloc[0]) / (float(b1["Time (h)"].iloc[1]) - float(b1["Time (h)"].iloc[0]))
    q("CA", "L2", "calculation", "Estimate the specific growth rate of FB-001 in the batch phase from the dry weight between the first "
      "two samples.",
      f"mu = ln({b1['DW (g/L)'].iloc[1]:.2f}/{b1['DW (g/L)'].iloc[0]:.2f}) / ({float(b1['Time (h)'].iloc[1]):g} - "
      f"{float(b1['Time (h)'].iloc[0]):g} h) = {mu:.2f} 1/h (OD600 gives a similar value).",
      [ev_path_("FB-001")], "numeric", value=round(mu, 3), tol=0.03, hops=1, tags=[])

    rw = None
    for f in glob.glob(str(C.p("hplc/raw/SEQ-*.csv"))):
        d = pd.read_csv(f, comment="#")
        d = d[(d.SampleID == "FB018-S06") & (d.Analyte == "Xylitol")]
        if len(d):
            rw, rf = d.iloc[0], Path(f).name
    xr = rw.Conc_vial_gL * rw.Dilution
    assert abs(xr - C.tr_value("FB018-S06", "Xylitol")) < 0.01
    q("CA", "L2", "calculation", "Recompute the xylitol concentration of sample FB018-S06 from the raw HPLC sequence file.",
      f"{rw.Conc_vial_gL:.4f} g/L in the vial x dilution {int(rw.Dilution)} = {xr:.2f} g/L ({rf}); this matches the processed "
      "result.", ["corpus/hplc/raw/" + rf, "corpus/hplc/processed/HPLC_results_2025-06.xlsx"], "numeric",
      value=round(xr, 2), tol=0.02, hops=2, tags=["units_conversion"])

    fv = (vv6.FEED_BAL_g.iloc[0] - vv6.FEED_BAL_g.iloc[-1]) / dens["Feed A"]
    q("CA", "L3", "calculation", "How much feed (mL) was added to FB-006?",
      f"Balance loss {vv6.FEED_BAL_g.iloc[0]:.2f} - {vv6.FEED_BAL_g.iloc[-1]:.2f} = "
      f"{vv6.FEED_BAL_g.iloc[0] - vv6.FEED_BAL_g.iloc[-1]:.1f} g; divided by the Feed A density {dens['Feed A']:.2f} g/mL "
      f"(SOP-FB-001) = {fv:.0f} mL. VesselView has no feed pump counter, so the balance is the only record.",
      [str(vp6), SOPFB, SETUP, MD], "numeric", value=round(fv, 1), tol=2.0, hops=3, tags=["units_conversion"])

    c5w = dict(C.xlsx_sheet("spreadsheets/EV_FB-005.xlsx", "Calculated"))
    l5 = bioctrl_last(C, "FB-005")
    pump5 = int(l5["Feed_Pump[cnt]"]) * ml_per["Feed"]
    bal5 = (c5w["Feed bottle start (g)"] - c5w["Feed bottle end (g)"]) / dens["Feed A"]
    q("CA", "L3", "calculation", "Check the claim in DEV-2025-01 that the R2 feed pump over-delivered by about 15 % in FB-005, "
      "using the raw records.",
      f"Pump: {int(l5['Feed_Pump[cnt]'])} counts x {ml_per['Feed']:.2f} mL = {pump5:.1f} mL. Balance: "
      f"({c5w['Feed bottle start (g)']} - {c5w['Feed bottle end (g)']}) g / {dens['Feed A']:.2f} g/mL = {bal5:.1f} mL. "
      f"Ratio {bal5 / pump5:.3f}, i.e. {100 * (bal5 / pump5 - 1):.0f} % more delivered than counted - the claim holds.",
      ["corpus/online/FB-005_R2.txt", "corpus/spreadsheets/EV_FB-005.xlsx", SOPFB, SETUP, "corpus/deviations/DEV-2025-01.md"],
      "numeric", value=round(100 * (bal5 / pump5 - 1), 1), tol=1.0, hops=4, tags=["units_conversion"])

    olig = C.coa("C5-2502")["Xylo-oligosaccharides (xylose eq.)"]
    q("CA", "L3", "calculation", "At most, how much extra xylose (g per litre of feed) could the enzyme release from the FB-023 feed, "
      "if it hydrolysed all xylo-oligosaccharides? Compare it with the xylose in that feed.",
      f"The feed is 60 % of lot C5-2502, whose CoA lists {olig:.2f} g/L xylo-oligosaccharides (xylose equivalents): "
      f"0.6 x {olig:.2f} = {0.6 * olig:.1f} g/L of feed, about {100 * 0.6 * olig / 140:.0f} % on top of the 140 g/L xylose the "
      "feed is topped up to (SOP). Oligomers are not measured by the HPLC method, so the actual release was not verified.",
      ["corpus/materials/CoA_C5-2502.md", SOPFB, "corpus/notes/FB-023_PN.md", "corpus/sops/SOP-ANA-003_v2.md"],
      "numeric", value=round(0.6 * olig, 1), tol=0.5, hops=3, tags=["vocabulary"])

    g18, _, gp18 = offgas(C, "FB-018")
    i18 = bioctrl_inoculation(C, "FB-018")
    g18["t"] = (g18.ts - i18).dt.total_seconds() / 3600
    yin = re.search(r"O2 ([\d.]+) %, CO2 ([\d.]+) %", str(_[2][4])).groups()
    o2in, co2in = float(yin[0]) / 100, float(yin[1]) / 100
    o2, co2 = g18["O2 [%]"] / 100, g18["CO2 [%]"] / 100
    r_ = (1 - o2in - co2in) / (1 - o2 - co2)
    g18["RQ"] = (co2 * r_ - co2in) / (o2in - o2 * r_)
    fed = g18[g18.t > 30]
    rq = float(fed.RQ.median())
    q("CA", "L4", "calculation", "Calculate the respiratory quotient (RQ) of FB-018 during the fed phase (after 30 h) from the off-gas data.",
      f"With an inert-gas balance (inlet {o2in * 100:.2f} % O2 / {co2in * 100:.2f} % CO2; outlet dry gas): "
      f"RQ = (yCO2,out*f - yCO2,in) / (yO2,in - yO2,out*f), f = (1 - yO2,in - yCO2,in)/(1 - yO2,out - yCO2,out). Median over "
      f"the fed phase = {rq:.2f} (interquartile range {fed.RQ.quantile(.25):.2f}-{fed.RQ.quantile(.75):.2f}). The O2 "
      "difference over the reactor is only ~0.1 %, so single cycles are noisy. RQ > 1 is expected: xylose reductase takes "
      "NAD(P)H that would otherwise be reoxidised with O2.",
      [str(gp18), "corpus/online/FB-018_R1.txt", SETUP], "numeric", value=round(rq, 2), tol=0.15, hops=3,
      tags=["time_alignment"])

    f12 = float(mds[mds.Run == "FB12"].For.iloc[-1])
    und = f12 / (1 + 10 ** (5.5 - 3.75))
    q("CA", "L4", "calculation", "What concentration of undissociated formic acid was present at the end of FB-012?",
      f"Formate at the end was {f12} g/L (MD master log); at the controlled pH 5.5 (SOP) and pKa 3.75 the undissociated "
      f"fraction is 1/(1 + 10^(5.5 - 3.75)) = {1 / (1 + 10 ** 1.75):.4f}, so about {und:.2f} g/L. The pKa of formic acid is "
      "general chemistry, not in the corpus.",
      [MD, SOPFB], "numeric", value=round(und, 3), tol=0.01, hops=2, tags=["vocabulary"])

    ev18, calc18, _i = C.ev("FB-018")
    b18 = ev18[ev18.Sample.str.contains("-S")].copy()
    for c_ in ("Time (h)", "Volume (mL)", "Xylose (g/L)", "Xylitol (g/L)"):
        b18[c_] = b18[c_].astype(float).fillna(0.0) if c_ != "Volume (mL)" else b18[c_].astype(float)
    feed_ml = (calc18["Feed bottle start (g)"] - calc18["Feed bottle end (g)"]) / dens["Feed A"]
    xf = C.tr_value("FB018-F", "Xylose")
    smp = b18.iloc[:-1]
    end = b18.iloc[-1]
    prod = end["Xylitol (g/L)"] * end["Volume (mL)"] / 1000 + (smp["Xylitol (g/L)"] * 7 / 1000).sum()
    fedx = feed_ml * xf / 1000
    left = end["Xylose (g/L)"] * end["Volume (mL)"] / 1000 + (smp["Xylose (g/L)"] * 7 / 1000).sum()
    yld = prod / (fedx - left)
    assert yld < 1.06
    q("CA", "L5", "calculation", "Estimate the xylitol yield on consumed xylose (g/g) for FB-018 with a mass balance, and say whether "
      "it is physically plausible.",
      f"Xylose fed: {feed_ml:.0f} mL feed (balance loss / {dens['Feed A']:.2f} g/mL) x {xf:.1f} g/L (feed sample FB018-F) = "
      f"{fedx:.1f} g. Xylose left: {end['Xylose (g/L)']:.2f} g/L x {end['Volume (mL)']:.0f} mL in the reactor plus "
      f"{(smp['Xylose (g/L)'] * 7 / 1000).sum():.2f} g removed in 7 mL samples = {left:.1f} g. Xylitol formed: "
      f"{end['Xylitol (g/L)']:.1f} g/L x {end['Volume (mL)']:.0f} mL + {(smp['Xylitol (g/L)'] * 7 / 1000).sum():.2f} g in "
      f"samples = {prod:.1f} g. Yield = {prod:.1f}/{fedx - left:.1f} = {yld:.2f} g/g. Plausible: the stoichiometric maximum "
      "is 1.01 g/g (152.15/150.13), and the estimate is within the measurement error of it (volumes are recorded to 10 mL, "
      "HPLC ~3 %).",
      [ev_path_("FB-018"), "corpus/hplc/processed/HPLC_results_2025-06.xlsx", SOPFB], "numeric", value=round(yld, 3),
      tol=0.05, hops=4, tags=["units_conversion"])

    grp = {"YX-001": ["FB-004", "FB-010"], "YX-002": ["FB-008"], "YX-003": ["FB-006", "FB-009"], "YX-004": ["FB-007", "FB-011"]}
    rx = {k: float(np.mean([C.final(r, "Xylose") for r in v])) for k, v in grp.items()}
    lo = min(rx, key=rx.get)
    q("CA", "L3", "calculation", "In the P3 comparison (excluding FB-005), which strain left the least residual xylose at the end "
      "(mean of its runs)?",
      f"{lo} ({rx[lo]:.1f} g/L); " + ", ".join(f"{k} {v:.1f}" for k, v in rx.items() if k != lo) + " g/L. Residual xylose "
      "reflects how much of the fed xylose a strain converted, so it is a more direct readout of conversion capacity than "
      "the final titre.", ["corpus/hplc/processed/", "corpus/lims/sample_register.csv", "corpus/deviations/DEV-2025-01.md"],
      "key_points", key_points=[lo, f"{rx[lo]:.1f}"], value=round(rx[lo], 1), tol=0.2, hops=3, tags=["id_resolution"])

    base = [C.ev("FB-001")[1]["Final xylitol (g/L)"], float(mds[mds.Run == "FB02"].Xol.iloc[-1]),
            C.ev("FB-003")[1]["Final xylitol (g/L)"]]
    sd = float(np.std(base, ddof=1))
    n = 2 * ((1.96 + 0.84) * sd / 5.0) ** 2
    q("CA", "L4", "statistics", "Using the run-to-run scatter of the three baseline runs, how many runs per strain would be needed to "
      "detect a 5 g/L difference in final xylitol between two strains (two-sided alpha 0.05, 80 % power)?",
      f"sd = {sd:.1f} g/L (FB-001 to FB-003, n = 3). n = 2((1.96 + 0.84) sd / 5)^2 = {n:.1f}, i.e. about {math.ceil(n)} runs per "
      "strain. With only 3 baseline runs the sd itself is very uncertain, so this is a lower bound; the P3 comparison used "
      "1-2 runs per strain.", [ev_path_("FB-001"), ev_path_("FB-003"), MD], "numeric", value=math.ceil(n), tol=1, hops=3,
      tags=["id_resolution"])

    # =================== XD: experiment design ===================
    dc_ = json.loads((C.root.parent / "answer_key" / "design_checks.json").read_text()) if \
        (C.root.parent / "answer_key" / "design_checks.json").exists() else {}

    def sim(fn):
        try:
            return fn(dc_)
        except Exception:
            return "not available (run design_checks.py)"
    fin = {r: C.final(r) for r in ["FB-004", "FB-010", "FB-006", "FB-009", "FB-007", "FB-011", "FB-008",
                                    "FB-012", "FB-013", "FB-014", "FB-015", "FB-016", "FB-017", "FB-020", "FB-021",
                                    "FB-022", "FB-023"]}
    p6o = {r: float(C.online(r)[0].query("t > 30").pO2.median()) for r in ("FB-020", "FB-021")}
    xyl_feed = re.search(r"Feed A: \d+ g/L glucose \+ (\d+) g/L xylose", sop).group(1)
    c5_top = re.search(r"xylose topped up to (\d+) g/L", sop).group(1)
    fdh_for = {r: mds[mds.Run == r.replace("FB-0", "FB")].For.iloc[-1] for r in ("FB-014", "FB-015", "FB-017")}
    fdh_txt = ("below detection (<0.05 g/L) in FB-014, FB-015 and FB-017" if all(str(v).startswith("<") for v in fdh_for.values())
               else ", ".join(f"{k} {v}" for k, v in fdh_for.items()) + " g/L")
    rubric = [
        dict(id="R1", suggestion="Repeat the C5-2503 run without enzyme (ideally with replicates).",
             evidence="FB-024 was fed C5-2502 by mistake (DEV-2025-04), so lot C5-2503 has no no-enzyme run; the "
                      "2025-07-11 minutes plan a repeat.",
             simulator_check=sim(lambda d: "Expected enzyme effect at 5 U/L on C5-2503 only "
                                 f"{d['yx006_C5-2503_60pct']['5 U/L']['Xol_end'] - d['yx006_C5-2503_60pct']['no enzyme']['Xol_end']:+.1f} g/L, "
                                 "below single-run scatter: one repeat run will not settle it.")),
        dict(id="R2", suggestion="Re-run the strain comparison with more replicates per strain (balanced over vessels/days) "
                                 "and a more sensitive readout (residual xylose, xylitol at 48-72 h, or a higher feed rate).",
             evidence=f"Differences between strains (control {np.mean([fin['FB-004'], fin['FB-010']]):.1f}, best "
                      f"{max(np.mean([fin['FB-007'], fin['FB-011']]), np.mean([fin['FB-006'], fin['FB-009']])):.1f} g/L) are of "
                      f"the order of the run-to-run sd ({sd:.1f} g/L); only 1-2 runs per strain; YX-002 has a single valid run "
                      "(FB-008); the control includes the power-cut run FB-010.",
             simulator_check=sim(lambda d: f"True YX-003 vs YX-001 difference {d['replicates_YX003_vs_YX001']['true_difference_gL']} g/L "
                                 f"with run-to-run sd {d['replicates_YX003_vs_YX001']['run_to_run_sd_gL']} g/L: about "
                                 f"{d['replicates_YX003_vs_YX001']['runs_per_strain_for_80pct_power']} runs per strain needed.")),
        dict(id="R3", suggestion="Increase the xylose supply (higher feed rate or more concentrated xylose) - titre is "
                                 "limited by the xylose fed.",
             evidence=f"Feed A carries {xyl_feed} g/L xylose and gives ~140-150 g/L xylitol; the C5 feeds are topped up to only "
                      f"{c5_top} g/L xylose and give ~62-75 g/L, even with FDH strains whose formate is fully removed; residual "
                      "xylose stays low (~7-18 g/L), so most fed xylose is converted.",
             simulator_check=sim(lambda d: "Feed A at 3 instead of 2 mL/h: YX-006 "
                                 f"{d['strains_feedA_2mLh']['YX-006']['Xol_end']:.0f} -> {d['strains_feedA_3mLh']['YX-006']['Xol_end']:.0f} g/L, "
                                 "and the YX-003 vs YX-001 difference grows from "
                                 f"{d['strains_feedA_2mLh']['YX-003']['Xol_end'] - d['strains_feedA_2mLh']['YX-001']['Xol_end']:.1f} to "
                                 f"{d['strains_feedA_3mLh']['YX-003']['Xol_end'] - d['strains_feedA_3mLh']['YX-001']['Xol_end']:.1f} g/L.")),
        dict(id="R4", suggestion="With the glucose-rich feed B keep oxygen non-limiting: >= 700-900 rpm or a pO2-controlled "
                                 "stirrer cascade.",
             evidence=f"FB-021 (500 rpm, feed B) had a fed-phase pO2 median of {p6o['FB-021']:.0f} % and ethanol re-appeared "
                      f"late; FB-020 (900 rpm) stayed at {p6o['FB-020']:.0f} %.",
             simulator_check=sim(lambda d: f"Feed B: 500 rpm {d['yx006_feedB_agitation']['500 rpm']['Xol_end']:.0f} g/L with "
                                 f"{d['yx006_feedB_agitation']['500 rpm']['EtOH_fed_max']:.1f} g/L ethanol; 700 rpm "
                                 f"{d['yx006_feedB_agitation']['700 rpm']['Xol_end']:.0f}, 900 rpm "
                                 f"{d['yx006_feedB_agitation']['900 rpm']['Xol_end']:.0f}, pO2 cascade at 20 % "
                                 f"{d['yx006_feedB_DO_cascade_20pct']['Xol_end']:.0f} g/L.")),
        dict(id="R5", suggestion="Push the C5 fraction above 70 % (up to undiluted) with YX-006, and test lot C5-2503 "
                                 "(more formic acid) at 70 %.",
             evidence=f"With FDH, 70 % C5 (FB-015 {fin['FB-015']:.1f}, FB-017 {fin['FB-017']:.1f} g/L) performs like 50 % "
                      f"(FB-014 {fin['FB-014']:.1f}) and formate is removed (end formate {fdh_txt}); "
                      "C5-2503 has the highest formic acid in its CoA.",
             simulator_check=sim(lambda d: "YX-006 on C5-2501: " + ", ".join(f"{k} {v['Xol_end']:.0f}" for k, v in
                                                                          d['yx006_C5-2501_fraction'].items())
                                 + " g/L - no penalty up to 100 %.")),
        dict(id="R6", suggestion="Do an enzyme dose-response and measure the oligomers directly (e.g. acid hydrolysis), "
                                 "since the HPLC method does not see them.",
             evidence=f"The +{fin['FB-023'] - fin['FB-022']:.1f} g/L of FB-023 vs FB-022 is a single-run difference within "
                      f"scatter; the CoAs list 29-43 g/L oligomers, i.e. up to ~{0.6 * olig:.0f} g/L extra xylose in a 60 % feed.",
             simulator_check=sim(lambda d: "C5-2502: no enzyme " f"{d['yx006_C5-2502_60pct']['no enzyme']['Xol_end']:.0f}, "
                                 f"5 U/L {d['yx006_C5-2502_60pct']['5 U/L']['Xol_end']:.0f}, 20 U/L "
                                 f"{d['yx006_C5-2502_60pct']['20 U/L']['Xol_end']:.0f} g/L.")),
        dict(id="R7", suggestion="Test PN's idea that the PPP strain (YX-003) needs more glucose to show its NADPH "
                                 "advantage: YX-001 vs YX-003 on feed B.",
             evidence="PN weekly report 2025-03-14; the P6 feed-B runs used only YX-006.",
             simulator_check=sim(lambda d: "Informative but negative: YX-003 advantage on feed B "
                                 f"{d['strains_feedB_2mLh']['YX-003']['Xol_end'] - d['strains_feedB_2mLh']['YX-001']['Xol_end']:+.1f} g/L "
                                 "vs on feed A "
                                 f"{d['strains_feedA_2mLh']['YX-003']['Xol_end'] - d['strains_feedA_2mLh']['YX-001']['Xol_end']:+.1f} g/L.")),
        dict(id="R8", suggestion="Close analytical gaps: analyse C5 runs with HPLC v2 (succinate, furfural, HMF), add a "
                                 "viability measurement, and repeat the single YX-002 run.",
             evidence="P5 C5 runs were analysed with v1 (no furans); FB-012's death is only inferred from xylose "
                      "accumulation; YX-002 rests on FB-008 alone.",
             simulator_check="Not simulated (measurement gap, not a process change)."),
    ]
    penalties = [
        "Uses FB-005 as valid data (excluded, DEV-2025-01).",
        "Uses superseded SEQ-20250321 xylitol values or the 490.7 g/L slip in the MD log.",
        "States that a strain is clearly better than the parent from the P3 data without addressing the scatter.",
        "States that the enzyme clearly increases xylitol, or treats FB-024 as a C5-2503 run.",
        "Recommends unbuffered flasks or a strain without FDH for >= 50 % C5.",
    ]
    q("XD", "L5", "experiment_design", "Based on all the data so far, what should the team try next? Propose up to five experiments, "
      "and for each give the evidence that motivates it and what result would change a decision.",
      "Any well-argued subset of the rubric items; at least three grounded in the corpus evidence, with no penalised claim.",
      ["corpus/"], "rubric", key_points=[r["id"] for r in rubric], hops=6, tags=["source_authority", "versioning"],
      rubric=rubric, penalties=penalties)

    d01 = np.mean([fin["FB-006"], fin["FB-009"]]) - np.mean([fin["FB-004"], fin["FB-010"]])
    q("XD", "L5", "experiment_design", "The team wants to show in the bioreactor that the PPP strain YX-003 makes more xylitol than "
      "the parent YX-001. Design that experiment (conditions, readout, number of runs, controls) and justify each choice from "
      "the data.",
      "Key elements: (1) replicates - the measured difference so far is "
      f"{d01:+.1f} g/L against a run-to-run sd of {sd:.1f} g/L, so ~{math.ceil(2 * ((1.96 + 0.84) * sd / max(d01, 1)) ** 2)} "
      "runs per strain (power calculation), not 1-2; (2) make xylose conversion, not xylose supply, limiting - higher feed "
      "rate or xylose concentration - and use residual xylose and xylitol at 48-72 h as readouts besides final titre; "
      "(3) balance strains over vessels R1-R4 and days, same feed batch, and keep pO2 non-limiting (900 rpm); "
      "(4) exclude or repeat runs with deviations (feed pump, power cut, probe drift) and use re-run-corrected HPLC values.",
      ["corpus/"], "rubric", key_points=["replicates", "feed rate", "residual xylose", "balanced vessels"], hops=5,
      tags=["source_authority"],
      rubric=[dict(id="D1", suggestion="replicates from a power calculation", evidence=f"sd {sd:.1f} g/L, difference {d01:+.1f} g/L",
                   simulator_check=sim(lambda d: f"{d['replicates_YX003_vs_YX001']['runs_per_strain_for_80pct_power']} runs per strain "
                                       "at standard conditions")),
              dict(id="D2", suggestion="higher xylose supply so that conversion capacity limits",
                   evidence="titre follows the xylose fed; residual xylose low",
                   simulator_check=sim(lambda d: "difference grows from "
                                       f"{d['strains_feedA_2mLh']['YX-003']['Xol_end'] - d['strains_feedA_2mLh']['YX-001']['Xol_end']:.1f} "
                                       f"to {d['strains_feedA_3mLh']['YX-003']['Xol_end'] - d['strains_feedA_3mLh']['YX-001']['Xol_end']:.1f} g/L at 3 mL/h")),
              dict(id="D3", suggestion="residual xylose / intermediate xylitol as readouts",
                   evidence=f"P3 residual xylose: lowest for {lo}", simulator_check=sim(lambda d: "residual xylose at 3 mL/h: YX-001 "
                                                                                  f"{d['strains_feedA_3mLh']['YX-001']['Xyl_end']:.1f} vs YX-003 "
                                                                                  f"{d['strains_feedA_3mLh']['YX-003']['Xyl_end']:.1f} g/L")),
              dict(id="D4", suggestion="balanced vessel/day assignment, non-limiting O2, exclude deviation runs",
                   evidence="vessel kLa differ, FB-005/FB-010/FB-009 deviations", simulator_check="Not simulated.")],
      penalties=["Relies on final titre from 1-2 runs per strain.", "Uses FB-005 or superseded HPLC values."])

    q("XD", "L4", "experiment_design", "Is it worth running YX-006 on undiluted (100 %) C5 fraction? What does the existing data say, "
      "and what risk remains untested?",
      f"Probably yes: with FDH, 70 % C5 gave the same titre as 50 % (FB-015 {fin['FB-015']:.1f}, FB-017 {fin['FB-017']:.1f} vs "
      f"FB-014 {fin['FB-014']:.1f} g/L) and formate was removed (end formate {fdh_txt} in the "
      f"master log), while the strain without FDH (FB-012) failed at 70 %. Untested: fractions above 70 %, furfural/HMF and "
      "acetate at higher fractions (P5 runs were analysed with HPLC v1, which has no furans), and lot C5-2503, which has more "
      "formic acid. Titre will stay limited by the xylose in the feed unless the xylose top-up is raised.",
      [MD, "corpus/meetings/minutes_2025-05-30.md", "corpus/sops/SOP-ANA-003_v1.md", "corpus/materials/CoA_C5-2503.md", SOPFB],
      "rubric", key_points=["70 % = 50 %", "formate removed", "furans not measured", "C5-2503"], hops=5,
      tags=["versioning", "vocabulary"],
      rubric=[dict(id="U1", suggestion="worth it: no penalty 50 -> 70 % with FDH", evidence="FB-013..017",
                   simulator_check=sim(lambda d: "YX-006: " + ", ".join(f"{k} {v['Xol_end']:.0f}" for k, v in
                                                                       d['yx006_C5-2501_fraction'].items()) + " g/L")),
              dict(id="U2", suggestion="risk: furans/acetate not measured in P5", evidence="v1 method", simulator_check="Not simulated."),
              dict(id="U3", suggestion="risk: C5-2503 formic acid", evidence="CoA",
                   simulator_check=sim(lambda d: "YX-006 on C5-2503: " + ", ".join(f"{k} {v['Xol_end']:.0f}" for k, v in
                                                                                d['yx006_C5-2503_fraction'].items()) + " g/L")),
              dict(id="U4", suggestion="titre stays xylose-limited", evidence=f"C5 feeds topped up to {c5_top} g/L",
                   simulator_check="See design_checks: titre flat with C5 fraction at constant xylose top-up.")],
      penalties=["Claims the C5 fraction is the main limit on titre.", "Ignores that FB-012 lacked FDH."])


def ev_path_(r):
    return f"corpus/spreadsheets/EV_{r}.xlsx"
