"""
Audit of a generated pilot: physics, measurement sanity and cross-document consistency.
Reads ONLY the corpus files and the answer key (not generator internals).

    python audit_pilot.py [--pilot ../pilot]
Prints PASS / FAIL / NOTE lines; deliberate faults are reported as NOTE with their fault code.
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import openpyxl
import pandas as pd

RESULTS = []


def check(ok, msg, note=False):
    RESULTS.append(("NOTE" if note else ("PASS" if ok else "FAIL"), msg))


def tr_results(corpus: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(str(corpus / "hplc/processed/*.xlsx"))):
        ws = openpyxl.load_workbook(f, data_only=True)["Results"]
        it = ws.iter_rows(values_only=True)
        hdr = next(it)
        rows += [dict(zip(hdr, r)) for r in it]
    return pd.DataFrame(rows)


def main(pilot: Path):
    corpus, key = pilot / "corpus", pilot / "answer_key"
    exp = pd.read_csv(key / "experiments.csv")
    smp = pd.read_csv(key / "samples_truth.csv")

    # ---------------- physics ----------------
    for _, e in exp[exp.kind == "flask"].iterrows():
        d = smp[(smp.run_id == e.run_id) & (smp.kind == "broth")]
        first, last = d.iloc[0], d.iloc[-1]
        cons = first.true_Xyl * first.true_V_mL - last.true_Xyl * last.true_V_mL      # mg, mass basis
        if cons > 250:
            y = last.true_Xol * last.true_V_mL / cons
            check(y <= 1.013, f"{e.run_id}: xylitol formed / xylose consumed {y:.3f} g/g (mass basis, max 1.013)")
    check(exp.final_X_true.max() < 25, f"max true biomass {exp.final_X_true.max():.1f} g/L plausible (<25)")
    check(exp[exp.kind == 'reactor'].final_Xol_true.max() < 200, "reactor xylitol titres < 200 g/L")
    for f in glob.glob(str(key / "timeseries/FB-*.csv")):
        ts = pd.read_csv(f)
        check(ts.V.max() < 0.8, f"{Path(f).stem}: max volume {ts.V.max()*1000:.0f} mL < 800 mL working limit")
        check((ts[["Glc", "Xyl", "Xol", "EtOH", "X"]] >= -1e-6).all().all(), f"{Path(f).stem}: no negative concentrations")

    # ---------------- measurement error ----------------
    for a in ["Xol", "Xyl", "Glc", "EtOH", "For"]:
        d = smp[(smp[f"true_{a}"] > 2) & smp[f"meas_{a}"].notna()]
        err = d[f"meas_{a}"] / d[f"true_{a}"] - 1
        big = d[err.abs() > 0.15]
        check(err.abs().median() < 0.04, f"HPLC {a}: median |error| {100*err.abs().median():.1f} %")
        for _, r in big.iterrows():
            check(False, f"HPLC {a} {r.sample_id}: {100*(r[f'meas_{a}']/r[f'true_{a}']-1):+.0f} % vs truth",
                  note=r.sample_id in ("FB012-S07", "FB013-S07") or a == "Xol")

    # ---------------- cross-document consistency ----------------
    T = tr_results(corpus)
    reg = pd.read_csv(corpus / "lims/sample_register.csv", comment="#")
    check(set(reg.sample_id) == set(T.SampleID), f"LIMS register ({len(reg)}) and HPLC results cover the same samples")
    raw_ids = set()
    for f in glob.glob(str(corpus / "hplc/raw/*.csv")):
        raw_ids |= set(pd.read_csv(f, comment="#").SampleID)
    check(set(reg.sample_id) <= raw_ids, "every registered sample appears in a raw HPLC sequence")
    T["Analysed"] = pd.to_datetime(T.Analysed)
    T["Sampled"] = pd.to_datetime(T.Sampled)
    check((T.Analysed > T.Sampled).all(), "every HPLC analysis date is after its sampling date")
    v2 = T[T.SOP == "v2"]
    check((v2.Analysed >= "2025-06-01").all() and (T[T.SOP == "v1"].Analysed < "2025-06-01").all(),
          "SOP v1/v2 used according to the effective date")
    xol = T[(T.Analyte == "Xylitol")]
    first = xol[~xol.Comment.fillna("").str.startswith("re-run")].set_index("SampleID").Result_gL
    rerun = xol[xol.Comment.fillna("").str.startswith("re-run")].set_index("SampleID").Result_gL

    # EV workbooks == TR (re-run values where they exist)
    for f in glob.glob(str(corpus / "spreadsheets/EV_FB-*.xlsx")):
        ws = openpyxl.load_workbook(f, data_only=True)["Offline data"]
        it = ws.iter_rows(values_only=True)
        hdr = next(it)
        bad = 0
        for r in it:
            d = dict(zip(hdr, r))
            sid = d["Sample"]
            if sid not in first.index or d.get("Xylitol (g/L)") is None:
                continue
            ref = rerun.get(sid, first[sid])
            if ref is None or abs(d["Xylitol (g/L)"] - ref) > 0.011:
                bad += 1
        check(bad == 0, f"{Path(f).name}: xylitol matches HPLC results (re-run where available)")

    # MD master log == TR first values, except the known slip
    ws = openpyxl.load_workbook(corpus / "spreadsheets/MD_fermentation_master_2025.xlsx", data_only=True)["Samples"]
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    for r in it:
        d = dict(zip(hdr, r))
        sid = f"{d['Run'].replace('FB', 'FB0')}-{d['Smpl']}"
        v = d["Xol"]
        if isinstance(v, (int, float)) and sid in first.index and first[sid] is not None:
            if abs(v - first[sid]) > 0.011:
                check(False, f"MD master log {sid}: xylitol {v} vs HPLC {first[sid]:.2f}",
                      note=abs(v / first[sid] - 10) < 0.05)

    # vessels never double-booked
    re_ = exp[exp.kind == "reactor"].copy()
    re_["start"] = pd.to_datetime(re_.start)
    tsend = {Path(f).stem: pd.read_csv(f).t.iloc[-1] for f in glob.glob(str(key / "timeseries/FB-*.csv"))}
    re_["end"] = re_.start + pd.to_timedelta(re_.run_id.map(tsend), unit="h")
    for v, g in re_.groupby("vessel"):
        g = g.sort_values("start")
        ok = (g.start.values[1:] > g.end.values[:-1]).all()
        check(ok, f"vessel {v}: {len(g)} runs, no overlap")

    # clone choices: registry vs answer key vs minutes
    dec = pd.read_csv(key / "clone_decisions.csv")
    ws = openpyxl.load_workbook(corpus / "strains/strain_registry.xlsx", data_only=True)["Strains"]
    rows = list(ws.iter_rows(values_only=True))
    regmap = {r[0]: r[7] for r in rows[1:]}
    for _, d in dec.iterrows():
        check(regmap.get(d.strain) == d.chosen, f"registry clone for {d.strain} = {d.chosen}")
    mins = " ".join(Path(f).read_text() for f in glob.glob(str(corpus / "meetings/*.md")))
    for _, d in dec.iterrows():
        check(f"{d.strain} {d.chosen}" in mins, f"minutes mention {d.strain} {d.chosen}")

    # notes name the strain/clone actually used
    for _, e in re_.iterrows():
        f = glob.glob(str(corpus / f"notes/{e.run_id}_*.md"))
        txt = Path(f[0]).read_text() if f else ""
        check(e.strain in txt, f"note {e.run_id} names strain {e.strain}")

    # lot recorded in MD media sheet = lot actually used
    ws = openpyxl.load_workbook(corpus / "spreadsheets/MD_fermentation_master_2025.xlsx", data_only=True)["Media & feeds"]
    media = {r[5]: r[3] for r in list(ws.iter_rows(values_only=True))[1:]}
    for _, e in re_[re_.operator == "MD"].iterrows():
        if isinstance(e.lot_actual, str):
            check(media.get(e.run_id) == e.lot_actual, f"MD media sheet lot for {e.run_id} = {e.lot_actual}")

    # ---------------- report ----------------
    n = {k: sum(1 for r in RESULTS if r[0] == k) for k in ("PASS", "FAIL", "NOTE")}
    for s, m in RESULTS:
        if s != "PASS":
            print(f"{s:4s} {m}")
    print(f"\n{n['PASS']} PASS, {n['FAIL']} FAIL, {n['NOTE']} NOTE (deliberate faults)")
    return n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default=str(Path(__file__).resolve().parents[1] / "pilot"))
    main(Path(ap.parse_args().pilot))
