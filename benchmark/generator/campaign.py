"""
Pilot campaign: storyline, run recipes, run-to-run variability, faults, clone selection.

Every number that ends up in a document is produced here: the simulator gives the true process,
measure.py turns it into recorded values. Decisions (which clone, which strain to drop) are taken
on the *measured* values, as a lab would, so a noisy measurement can lead to a wrong pick.
"""
from __future__ import annotations

import dataclasses as dc
import datetime as dtm
import math

import numpy as np
import pandas as pd

import measure as ms
import world as W
from simulator import Experiment, FeedSegment, Pulse, simulate

# ---------------------------------------------------------------------------
# Run-to-run variability (biological + process), lognormal sd on the multiplier
# ---------------------------------------------------------------------------
VAR_REACTOR = dict(mumax2=0.07, KPI2=0.07, kLa=0.08, X0=0.08, feed=0.02)
VAR_FLASK = dict(mumax2=0.05, KPI2=0.05, kLa=0.12, X0=0.08)

REACTOR_SAMPLE_H = [0.1, 4.0, 8.0, 23.5, 30.0, 47.5, 72.0, 96.0, 120.0, 143.5]
HYDRO_SAMPLE_H = REACTOR_SAMPLE_H + [167.5]
FLASK_SAMPLE_H = [0.1, 24.0, 48.0, 72.0]

FEED_A = {"Glc": 100.0, "Xyl": 300.0}         # standard synthetic feed (g/L)
FEED_B = {"Glc": 150.0, "Xyl": 300.0}         # higher glucose co-feed


def hydrolysate_feed(lot: str, frac: float, glc_target: float = 30.0, xyl_target: float = 140.0) -> dict:
    """Feed = C5 fraction diluted to `frac`, glucose supplemented to glc_target and xylose topped up
    to xyl_target (g/L). Trc carries % C5 fraction."""
    L = W.LOTS[lot]
    f = {k: L[k] * frac for k in ["Xyl", "Glc", "Ara", "Olig", "For", "Ac", "Fur", "HMF"]}
    f["Glc"] = max(f["Glc"], glc_target)
    f["Xyl"] = max(f["Xyl"], xyl_target)
    f["Trc"] = 100.0 * frac
    return f


@dc.dataclass
class RunSpec:
    run_id: str
    phase: str
    kind: str                    # reactor | flask
    vessel: str                  # R1..R4 or flask position
    operator: str
    start: dtm.datetime          # inoculation
    strain_id: str
    clone: int | None = None     # None -> the currently chosen clone
    purpose: str = ""
    feed: dict | None = None
    feed_rate_mLh: float = 2.0
    feed_start_h: float = 7.5
    medium: dict = dc.field(default_factory=lambda: {"Glc": 20.0})
    lot: str | None = None       # lot recorded as used (planned)
    lot_actual: str | None = None
    frac: float | None = None
    stirrer: float = 900.0
    pH: float = 5.5
    buffered: bool = True
    enzyme_h: float | None = None
    t_end: float = 143.5
    faults: list = dc.field(default_factory=list)
    replicate_of: str | None = None


def _d(s, hm="09:30"):
    return dtm.datetime.fromisoformat(f"{s} {hm}")


# ---------------------------------------------------------------------------
# Storyline
# ---------------------------------------------------------------------------
def reactor_plan() -> list[RunSpec]:
    R = []
    # P1 baseline
    R += [RunSpec("FB-001", "P1", "reactor", "R1", "EV", _d("2025-02-03"), "YX-001", purpose="Baseline fed-batch, parent strain"),
          RunSpec("FB-002", "P1", "reactor", "R2", "MD", _d("2025-02-03", "10:15"), "YX-001", purpose="Baseline fed-batch, parent strain",
                  faults=[dict(code="F-INOC", what="old preculture, low inoculum", X0_factor=0.55)]),
          RunSpec("FB-003", "P1", "reactor", "R3", "EV", _d("2025-02-03", "11:00"), "YX-001", purpose="Baseline fed-batch, parent strain")]
    # P3 engineered strains, two rounds
    R += [RunSpec("FB-004", "P3", "reactor", "R1", "EV", _d("2025-03-17"), "YX-001", purpose="Control"),
          RunSpec("FB-005", "P3", "reactor", "R2", "EV", _d("2025-03-17", "10:00"), "YX-002", purpose="Second XYL1 copy",
                  faults=[dict(code="F-FEEDCAL", what="feed pump delivered 15 % more than logged", feed_factor=1.15)]),
          RunSpec("FB-006", "P3", "reactor", "R3", "MD", _d("2025-03-17", "10:40"), "YX-003", purpose="ZWF1/GND1 overexpression"),
          RunSpec("FB-007", "P3", "reactor", "R4", "MD", _d("2025-03-17", "11:10"), "YX-004", purpose="gpd1 deletion"),
          RunSpec("FB-008", "P3", "reactor", "R1", "EV", _d("2025-03-31"), "YX-002", purpose="Second XYL1 copy, repeat", replicate_of="FB-005"),
          RunSpec("FB-009", "P3", "reactor", "R3", "MD", _d("2025-03-31", "10:20"), "YX-003", purpose="ZWF1/GND1, repeat", replicate_of="FB-006",
                  faults=[dict(code="F-DODRIFT", what="pO2 membrane fouling, signal drifts down", drift=-0.35)]),
          RunSpec("FB-010", "P3", "reactor", "R4", "MD", _d("2025-03-31", "10:50"), "YX-001", purpose="Control, repeat", replicate_of="FB-004",
                  faults=[dict(code="F-POWER", what="UPS failure, logger and feed pump off", t0=49.5, t1=53.0)]),
          RunSpec("FB-011", "P3", "reactor", "R2", "EV", _d("2025-03-31", "11:30"), "YX-004", purpose="gpd1 deletion, repeat", replicate_of="FB-007")]
    # P5 hydrolysate fed-batch (lot C5-2501)
    for rid, v, s, frac, day, hm, flt in [
            ("FB-012", "R1", "YX-003", 0.70, "2025-05-05", "09:30", []),
            ("FB-013", "R2", "YX-005", 0.50, "2025-05-05", "10:05", []),
            ("FB-014", "R3", "YX-006", 0.50, "2025-05-05", "10:40", [dict(code="F-ANTIFOAM", what="foaming, antifoam added", t=30.0, ppm=150.0)]),
            ("FB-015", "R4", "YX-006", 0.70, "2025-05-05", "11:15", []),
            ("FB-016", "R1", "YX-005", 0.70, "2025-05-19", "09:30", []),
            ("FB-017", "R2", "YX-006", 0.70, "2025-05-19", "10:05", [])]:
        R.append(RunSpec(rid, "P5", "reactor", v, "MD", _d(day, hm), s, purpose=f"{int(frac*100)} % C5 fraction feed",
                         feed=hydrolysate_feed("C5-2501", frac), lot="C5-2501", lot_actual="C5-2501", frac=frac,
                         t_end=167.5, faults=flt))
    R[-1].replicate_of = "FB-015"
    # P6 O2 and feed-ratio study (synthetic feeds)
    R += [RunSpec("FB-018", "P6", "reactor", "R1", "EV", _d("2025-06-09"), "YX-006", purpose="900 rpm, feed A"),
          RunSpec("FB-019", "P6", "reactor", "R2", "EV", _d("2025-06-09", "10:00"), "YX-006", purpose="500 rpm, feed A", stirrer=500.0),
          RunSpec("FB-020", "P6", "reactor", "R3", "EV", _d("2025-06-09", "10:30"), "YX-006", purpose="900 rpm, feed B (more glucose)", feed=FEED_B),
          RunSpec("FB-021", "P6", "reactor", "R4", "EV", _d("2025-06-09", "11:00"), "YX-006", purpose="500 rpm, feed B", stirrer=500.0, feed=FEED_B)]
    # P7 lot comparison, +/- enzyme
    R += [RunSpec("FB-022", "P7", "reactor", "R1", "PN", _d("2025-06-23"), "YX-006", purpose="Lot C5-2502, 60 %",
                  feed=hydrolysate_feed("C5-2502", 0.6), lot="C5-2502", lot_actual="C5-2502", frac=0.6, t_end=167.5),
          RunSpec("FB-023", "P7", "reactor", "R2", "PN", _d("2025-06-23", "10:00"), "YX-006", purpose="Lot C5-2502, 60 %, + enzyme",
                  feed=hydrolysate_feed("C5-2502", 0.6), lot="C5-2502", lot_actual="C5-2502", frac=0.6, enzyme_h=24.0, t_end=167.5),
          RunSpec("FB-024", "P7", "reactor", "R3", "MD", _d("2025-06-23", "10:30"), "YX-006", purpose="Lot C5-2503, 60 %",
                  feed=hydrolysate_feed("C5-2502", 0.6), lot="C5-2503", lot_actual="C5-2502", frac=0.6, t_end=167.5,
                  faults=[dict(code="F-WRONGLOT", what="feed made from C5-2502 instead of C5-2503")]),
          RunSpec("FB-025", "P7", "reactor", "R4", "MD", _d("2025-06-23", "11:00"), "YX-006", purpose="Lot C5-2503, 60 %, + enzyme",
                  feed=hydrolysate_feed("C5-2503", 0.6), lot="C5-2503", lot_actual="C5-2503", frac=0.6, enzyme_h=24.0, t_end=167.5)]
    for r in R:
        if r.feed is None:
            r.feed = dict(FEED_A)
        if r.t_end > 150:
            pass
    return R


def flask_plan_p2() -> list[RunSpec]:
    F, i = [], 1
    for s, clones in [("YX-001", [None]), ("YX-002", [0, 1, 2]), ("YX-003", [0, 1, 2]), ("YX-004", [0, 1, 2])]:
        for c in clones:
            for rep in range(2):
                F.append(RunSpec(f"SF-{i:03d}", "P2", "flask", f"shelf A{i}", "PN", _d("2025-03-04", "10:00"), s, clone=c,
                                 purpose="Clone screen, YP glucose/xylose", medium={"Glc": 20.0, "Xyl": 50.0},
                                 pH=6.0, t_end=72.0))
                i += 1
    return F


def flask_plan_p4(start_index: int) -> list[RunSpec]:
    F, i = [], start_index
    def add(s, c, frac, buffered=True, purpose=""):
        nonlocal i
        L = W.LOTS["C5-2501"]
        med = {k: L[k] * frac for k in ["Xyl", "Glc", "Ara", "Olig", "For", "Ac", "Fur", "HMF"]}
        med["Glc"] = med["Glc"] + 20.0
        med["Trc"] = frac * 100
        F.append(RunSpec(f"SF-{i:03d}", "P4", "flask", f"shelf B{i - start_index + 1}", "PN", _d("2025-04-22", "10:00"), s,
                         clone=c, purpose=purpose, medium=med, lot="C5-2501", lot_actual="C5-2501", frac=frac,
                         pH=6.0 if buffered else 4.6, buffered=buffered, t_end=72.0))
        i += 1
    for frac in (0.3, 0.5, 0.7):
        add("YX-001", None, frac, purpose=f"Tolerance, {int(frac*100)} % C5")
    for s in ("YX-005", "YX-006"):
        for c in (0, 1, 2):
            add(s, c, 0.5, purpose="FDH clone screen, 50 % C5")
    for s in ("YX-005", "YX-006"):
        for frac in (0.3, 0.7):
            add(s, "chosen", frac, purpose=f"Tolerance, {int(frac*100)} % C5")
    add("YX-001", None, 0.5, buffered=False, purpose="Unbuffered, 50 % C5")
    add("YX-006", "chosen", 0.5, buffered=False, purpose="Unbuffered, 50 % C5")
    return F


# ---------------------------------------------------------------------------
# Simulation of one run
# ---------------------------------------------------------------------------
def lognorm(rng, sd):
    return float(math.exp(rng.normal(0, sd)))


def build_experiment(spec: RunSpec, strain, rng) -> tuple[Experiment, dict]:
    var = VAR_FLASK if spec.kind == "flask" else VAR_REACTOR
    draws = {k: lognorm(rng, v) for k, v in var.items()}
    st = strain.copy(mumax2=strain.mumax2 * draws["mumax2"], KPI2=strain.KPI2 * draws["KPI2"])
    X0 = 0.8 * draws["X0"]
    for f in spec.faults:
        X0 *= f.get("X0_factor", 1.0)
    if spec.kind == "flask":
        reactor = dc.replace(W.FLASK, fixed_kLa=W.FLASK.fixed_kLa * draws["kLa"])
        ex = Experiment(run_id=spec.run_id, strain=st, reactor=reactor, V0=0.050, X0=X0, medium=dict(spec.medium),
                        T_set=30.0, pH_set=spec.pH, pH_control=False, stirrer=200.0, airflow=0.0,
                        sample_times=FLASK_SAMPLE_H[1:], sample_volume=0.001, t_end=spec.t_end)
        return ex, draws
    v = W.VESSELS[spec.vessel].reactor
    reactor = dc.replace(v, kLa_ref=v.kLa_ref * draws["kLa"])
    feed_factor = draws["feed"]
    for f in spec.faults:
        feed_factor *= f.get("feed_factor", 1.0)
    rate = spec.feed_rate_mLh / 1000.0 * feed_factor
    segs = [FeedSegment(spec.feed_start_h, spec.t_end, rate)]
    antifoam = []
    for f in spec.faults:
        if f["code"] == "F-POWER":
            segs = [FeedSegment(spec.feed_start_h, f["t0"], rate), FeedSegment(f["t1"], spec.t_end, rate)]
        if f["code"] == "F-ANTIFOAM":
            antifoam.append((f["t"], f["ppm"]))
    samples = HYDRO_SAMPLE_H if spec.t_end > 150 else REACTOR_SAMPLE_H
    ftype = "A" if spec.feed == FEED_A else "B" if spec.feed == FEED_B else "C5"
    ex = Experiment(run_id=spec.run_id, strain=st, reactor=reactor, V0=0.300, X0=X0, medium=dict(spec.medium),
                    feed_density=W.FEED_DENSITY[ftype],
                    T_set=30.0, pH_set=spec.pH, stirrer=spec.stirrer, airflow=0.7,
                    feed_segments=segs, feed_composition=dict(spec.feed), feed_vessel0=float(rng.uniform(760, 820)),
                    sample_times=[t for t in samples[1:] if t < spec.t_end], sample_volume=0.007,
                    t_end=spec.t_end, antifoam_additions=antifoam,
                    enzyme_additions=[(spec.enzyme_h, 5.0)] if spec.enzyme_h else [])
    draws["feed_factor_total"] = feed_factor
    return ex, draws


# ---------------------------------------------------------------------------
# HPLC sequences: samples are frozen and run in batches (Tue / Fri)
# ---------------------------------------------------------------------------
SOP_V2_DATE = dtm.date(2025, 6, 1)
HPLC_BIAS_FAULT = dict(code="F-HPLCSTD", seq_date=dtm.date(2025, 3, 21), analyte="Xol", bias=0.10,
                       rerun_date=dtm.date(2025, 4, 9))         # extra Wednesday sequence


def hplc_batch_date(sample_dt: dtm.datetime) -> dtm.date:
    d = sample_dt.date()
    while True:
        d += dtm.timedelta(days=1)
        if d.weekday() in (1, 4):          # Tuesday, Friday
            return d


# ---------------------------------------------------------------------------
class Campaign:
    def __init__(self, seed: int = 20250203):
        self.rng = np.random.default_rng(seed)
        self.strains = {}            # strain_id -> dict(clones=[(rec, fdh)], chosen=int, built_on=str, plan)
        self.runs = []               # dicts with spec, experiment, draws, result, measurements
        self.samples = []            # measured sample records (all runs)
        self.sequences = {}          # date -> dict(bias, analytes, version, injections=[])
        self.decisions = []
        self.events = []             # lab-level events (maintenance, deviations ...)

    # -- strains -----------------------------------------------------------
    def construct(self, plan: W.StrainPlan, parent_id: str | None, parent_clone: int | None):
        if parent_id is None:
            rec = W.root_record()
            self.strains[plan.strain_id] = dict(plan=plan, clones=[(rec, 0.0)], chosen=0, built_on=None)
            return
        p = self.strains[parent_id]
        pc = p["chosen"] if parent_clone is None else parent_clone
        prec, pfdh = p["clones"][pc]
        clones = W.make_clones(plan, prec, pfdh, self.rng)
        self.strains[plan.strain_id] = dict(plan=plan, clones=clones, chosen=None,
                                            built_on=f"{parent_id}-c{pc + 1}" if len(p["clones"]) > 1 else parent_id)

    def strain_for(self, spec: RunSpec):
        s = self.strains[spec.strain_id]
        c = spec.clone
        if c is None or c == "chosen":
            c = s["chosen"] if s["chosen"] is not None else 0
        return s["clones"][c][0].strain, c

    # -- one run -----------------------------------------------------------
    def run(self, spec: RunSpec):
        strain, clone = self.strain_for(spec)
        ex, draws = build_experiment(spec, strain, self.rng)
        res = simulate(ex)
        rec = dict(spec=spec, clone=clone, experiment=ex, draws=draws, result=res)
        self.runs.append(rec)
        self._measure(rec)
        return rec

    def _measure(self, rec):
        spec, res = rec["spec"], rec["result"]
        rng = self.rng
        smp = res.samples.copy()
        # t0 sample (after inoculation) from truth
        t0row = res.truth.iloc[(res.truth["t"] - 0.1).abs().argmin()]
        rows = [("S01", float(t0row["t"]), t0row)]
        for i, (_, r) in enumerate(smp.iterrows()):
            rows.append((f"S{i + 2:02d}", float(r["t"]), r))
        last = res.truth.iloc[-1]
        if abs(rows[-1][1] - last["t"]) > 0.05:
            rows.append((f"S{len(rows) + 1:02d}", float(last["t"]), last))
        medium = {k: v for k, v in rec["experiment"].medium.items()}
        out = []
        prefix = spec.run_id.replace("-", "")
        # medium (and feed) samples
        out.append(self._sample(spec, f"{prefix}-M", -0.5, {**medium, "X": 0.0}, kind="medium"))
        if spec.kind == "reactor":
            out.append(self._sample(spec, f"{prefix}-F", -0.5, {**spec.feed, "X": 0.0}, kind="feed"))
        for sid, t, r in rows:
            true = {k: float(r[k]) for k in ms.HPLC_ANALYTES_V2 + ["X", "Xd"]}
            true["V_mL"] = float(r["V"]) * 1000
            out.append(self._sample(spec, f"{prefix}-{sid}", t, true, kind="broth"))
        rec["samples"] = out
        # online + off-gas for reactors
        if spec.kind == "reactor":
            drift = None
            gaps = []
            feed_log = 1.0 / rec["draws"]["feed_factor_total"] * rec["draws"]["feed"]   # logged = setpoint-calibrated
            for f in spec.faults:
                if f["code"] == "F-DODRIFT":
                    drift = f["drift"]
                if f["code"] == "F-POWER":
                    gaps.append((f["t0"], f["t1"]))
            rec["online"] = ms.online_signals(res.truth, rng, spec.pH, 30.0, drift, gaps, feed_log_factor=feed_log)
            if W.VESSELS[spec.vessel].vessel_id in ("R1", "R2"):
                rec["offgas"] = ms.offgas_cycles(res.offgas(), rng, t_shift=float(rng.uniform(0.5, 18)))

    def _sample(self, spec, sample_id, t, true, kind):
        rng = self.rng
        when = spec.start + dtm.timedelta(hours=max(t, -0.5))
        seq_date = hplc_batch_date(when)
        version = 2 if seq_date >= SOP_V2_DATE else 1
        analytes = ms.HPLC_ANALYTES_V2 if version == 2 else ms.HPLC_ANALYTES_V1
        if seq_date not in self.sequences:
            self.sequences[seq_date] = dict(bias=ms.sequence_bias(ms.HPLC_ANALYTES_V2, rng), version=version,
                                            injections=[])
        seq = self.sequences[seq_date]
        extra = {}
        if seq_date == HPLC_BIAS_FAULT["seq_date"]:
            extra = {HPLC_BIAS_FAULT["analyte"]: HPLC_BIAS_FAULT["bias"]}
        inj = ms.hplc_inject(true, analytes, seq["bias"], rng, extra)
        seq["injections"].append(dict(sample_id=sample_id, run_id=spec.run_id, **inj))
        s = dict(sample_id=sample_id, run_id=spec.run_id, kind=kind, t=round(t, 2), datetime=when,
                 hplc_date=seq_date, sop_version=version, dilution=inj["dilution"],
                 hplc={a: (None if v["flag"] else v["conc"]) for a, v in inj["results"].items()},
                 truth=true)
        if kind == "broth":
            s["dw"] = ms.dry_weight(true["X"] + true["Xd"], rng)      # viable + non-viable intact cells
            s["od"] = ms.od600(true["X"] + true["Xd"], rng)
            s["V_mL"] = true["V_mL"]
        self.samples.append(s)
        return s

    # -- HPLC re-run of the biased sequence -----------------------------------
    def hplc_rerun(self):
        bad = self.sequences.get(HPLC_BIAS_FAULT["seq_date"])
        if not bad:
            return
        rd = HPLC_BIAS_FAULT["rerun_date"]
        assert rd not in self.sequences, "re-run must not overwrite a routine sequence"
        self.sequences[rd] = dict(bias=ms.sequence_bias(ms.HPLC_ANALYTES_V2, self.rng), version=1, injections=[],
                                  rerun_of=HPLC_BIAS_FAULT["seq_date"])
        for inj in bad["injections"]:
            s = next(x for x in self.samples if x["sample_id"] == inj["sample_id"])
            if s["kind"] != "broth":
                continue
            new = ms.hplc_inject(s["truth"], ms.HPLC_ANALYTES_V1, self.sequences[rd]["bias"], self.rng)
            self.sequences[rd]["injections"].append(dict(sample_id=s["sample_id"], run_id=s["run_id"], **new))
            s["hplc_rerun"] = {a: (None if v["flag"] else v["conc"]) for a, v in new["results"].items()}

    # -- clone choice from measured data ----------------------------------------
    def choose_clone(self, strain_id: str, metric_hour: float, phase_runs: list, decided_on: str, by: str):
        vals = {}
        for rec in phase_runs:
            if rec["spec"].strain_id != strain_id:
                continue
            s = min((x for x in rec["samples"] if x["kind"] == "broth"), key=lambda x: abs(x["t"] - metric_hour))
            v = s["hplc"]["Xol"] or 0.0
            vals.setdefault(rec["clone"], []).append(v)
        means = {c: float(np.mean(v)) for c, v in vals.items()}
        chosen = max(means, key=means.get)
        true_best = max(range(len(self.strains[strain_id]["clones"])),
                        key=lambda c: self.strains[strain_id]["clones"][c][0].strain.mumax2
                        * (1 + self.strains[strain_id]["clones"][c][1]))
        self.strains[strain_id]["chosen"] = chosen
        self.decisions.append(dict(date=decided_on, by=by, kind="clone_choice", strain=strain_id,
                                   chosen=f"c{chosen + 1}", measured_means=means,
                                   true_best=f"c{true_best + 1}", correct=chosen == true_best))
        return chosen


# ---------------------------------------------------------------------------
def run_pilot(seed: int = 20250203) -> Campaign:
    C = Campaign(seed)
    plans = {p.strain_id: p for p in W.STRAIN_PLAN}
    # strains available before P2: YX-002 built on YX-001; YX-003/004 built on YX-002 c1 (before screening)
    C.construct(plans["YX-001"], None, None)
    C.construct(plans["YX-002"], "YX-001", 0)
    C.construct(plans["YX-003"], "YX-002", 0)
    C.construct(plans["YX-004"], "YX-002", 0)
    rplan = {r.run_id: r for r in reactor_plan()}
    # P1
    for rid in ("FB-001", "FB-002", "FB-003"):
        C.run(rplan[rid])
    # P2 clone screen
    p2 = [C.run(s) for s in flask_plan_p2()]
    for sid in ("YX-002", "YX-003", "YX-004"):
        C.choose_clone(sid, 48.0, p2, "2025-03-14", "PN")
    # P3
    for rid in ("FB-004", "FB-005", "FB-006", "FB-007", "FB-008", "FB-009", "FB-010", "FB-011"):
        C.run(rplan[rid])
    C.hplc_rerun()
    # FDH strains: YX-005 on YX-001, YX-006 on the chosen YX-003 clone
    C.construct(plans["YX-005"], "YX-001", 0)
    C.construct(plans["YX-006"], "YX-003", None)
    # P4: FDH clone screen at 50 % first, then the chosen clones at 30/70 % and unbuffered
    p4_specs = flask_plan_p4(len(p2) + 1)
    first = [s for s in p4_specs if s.clone != "chosen"]
    p4 = [C.run(s) for s in first]
    for sid in ("YX-005", "YX-006"):
        C.choose_clone(sid, 72.0, [r for r in p4 if r["spec"].purpose.startswith("FDH")], "2025-04-29", "PN")
    p4 += [C.run(s) for s in p4_specs if s.clone == "chosen"]
    # P5, P6, P7
    for rid in [r for r in rplan if rplan[r].phase in ("P5", "P6", "P7")]:
        C.run(rplan[rid])
    return C


if __name__ == "__main__":
    import time
    t0 = time.time()
    C = run_pilot()
    print(f"{len(C.runs)} runs, {len(C.samples)} samples, {len(C.sequences)} HPLC sequences in {time.time()-t0:.0f} s")
    for d in C.decisions:
        print(d)
    for r in C.runs:
        if r["spec"].kind == "reactor":
            e = r["result"].truth.iloc[-1]
            print(f"{r['spec'].run_id} {r['spec'].strain_id}-c{r['clone']+1} {r['spec'].purpose:32s} Xol {e.Xol:6.1f} Xyl {e.Xyl:5.1f} X {e.X:5.1f} EtOHmax {r['result'].truth.EtOH.max():4.1f} Gly {e.Gly:4.2f} For {e.For:5.2f}")
