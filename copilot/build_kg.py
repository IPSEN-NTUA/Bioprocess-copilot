"""
Build the knowledge graph from benchmark/pilot/corpus.

Deterministic: no LLM is involved, and nothing outside corpus/ is read (never answer_key/).
Every source keeps its own statements; conflicting values from different documents are all kept, each
with its provenance, and are never silently merged. What the builder does reconcile is the *form* of
the records, which is where the interoperability work lies:

  identifiers   FB02 / FB-002 / fb002 / FB002-S08 / ('FB15', 'S08')  -> one run and one sample node
  vocabulary    Xol / Xylitol, For / Formic acid / Formate           -> one analyte node
  numbers       decimal commas, '<0.05', '<LOD', 'n.d.'              -> value or belowLOD flag
  dates         dd/mm/yy, dd.mm.yyyy, '23 Jun 2025', ISO             -> xsd:date / xsd:dateTime
  time bases    logger Elapsed_h, analyser RelTime, wall clock       -> processTimeH (h after inoculation)
  units         pump counts                                          -> mL (calibration read from the corpus)
  versions      'superseded by re-run' comments                      -> supersededBy links, isCurrent flag

Usage:  python -m copilot.build_kg [--corpus DIR] [--store DIR]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dtm
import gc
import re
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import openpyxl
from pyoxigraph import DefaultGraph, Literal, NamedNode, Quad, RdfFormat, Store

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "benchmark" / "pilot" / "corpus"
STORE = ROOT / "copilot" / ".cache" / "kg_store"
ONTOLOGY = Path(__file__).resolve().parent / "ontology" / "bioprocess.ttl"

BP = "http://example.org/bioprocess#"
RES = "http://example.org/bioprocess/id/"
GRAPH = "http://example.org/bioprocess/graph/"
XSD = "http://www.w3.org/2001/XMLSchema#"
RDF_TYPE = NamedNode("http://www.w3.org/1999/02/22-rdf-syntax-ns#type")
RDFS_LABEL = NamedNode("http://www.w3.org/2000/01/rdf-schema#label")
SKOS_ALT = NamedNode("http://www.w3.org/2004/02/skos/core#altLabel")

PREFIXES = f"""PREFIX bp: <{BP}>
PREFIX res: <{RES}>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX xsd: <{XSD}>
"""

# Analyte vocabulary: preferred label -> the other spellings found in the corpus.
ANALYTES = {
    "Glucose": ["Glc"], "Xylose": ["Xyl", "Xylose (monomeric)"], "Arabinose": ["Ara"], "Xylitol": ["Xol"],
    "Glycerol": ["Gly"], "Acetate": ["Ac", "Acetic acid"], "Ethanol": ["EtOH"],
    "Formate": ["For", "Formic acid"], "Pyruvate": ["Pyr"], "Succinate": ["Suc"], "Furfural": ["Fur"],
    "HMF": ["5-HMF"], "Xylo-oligosaccharides": ["Xylo-oligosaccharides (xylose eq.)"],
    "OD600": ["OD"], "Dry weight": ["DW"],
}
_ANALYTE_LOOKUP = {name.lower(): pref for pref, alts in ANALYTES.items() for name in [pref, *alts]}


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
def norm_run(text) -> str | None:
    """FB02, FB-002, fb002, FB002, SF-7 -> FB-002 / SF-007."""
    m = re.search(r"\b(FB|SF)[-_ ]?0*(\d{1,3})(?!\d)", str(text), re.I)
    return f"{m.group(1).upper()}-{int(m.group(2)):03d}" if m else None


def norm_sample(text, run=None) -> str | None:
    """'FB005-S07', ('S08' with run 'FB15') -> 'FB015-S08'. Suffix is Snn, M (medium) or F (feed)."""
    m = re.search(r"\b(FB|SF)[-_ ]?0*(\d{1,3})[-_ /]+(S\d{1,2}|M|F)\b", str(text), re.I)
    if m:
        run_id, suffix = f"{m.group(1).upper()}-{int(m.group(2)):03d}", m.group(3).upper()
    elif run and re.fullmatch(r"S\d{1,2}|M|F", str(text).strip(), re.I):
        run_id, suffix = norm_run(run), str(text).strip().upper()
    else:
        return None
    if suffix.startswith("S"):
        suffix = f"S{int(suffix[1:]):02d}"
    return f"{run_id.replace('-', '')}-{suffix}"


def norm_analyte(text) -> str | None:
    key = re.sub(r"\s*\(g/L\)\s*$", "", str(text)).strip().lower()
    return _ANALYTE_LOOKUP.get(key)


def to_number(v):
    """Returns (value, below_lod, lod_limit). Handles decimal commas, '<0.05', '<LOD', 'n.d.', blanks."""
    if v is None:
        return None, False, None
    if isinstance(v, (int, float)):
        return float(v), False, None
    s = str(v).strip()
    if s == "":
        return None, False, None
    if s.startswith("<") or s.lower() in ("n.d.", "nd", "<lod"):
        m = re.search(r"<\s*([\d.,]+)", s)
        return None, True, float(m.group(1).replace(",", ".")) if m else None
    try:
        return float(s.replace(",", ".")), False, None
    except ValueError:
        return None, False, None


_DATE_FORMATS = ["%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
                 "%d %b %Y %H:%M:%S", "%d %b %Y", "%d.%m.%Y %H:%M:%S", "%d.%m.%Y"]


def to_datetime(text, day_first_short=False) -> dtm.datetime | None:
    """Parse the date styles of the corpus. dd/mm/yy (day first) only when the source is known to use it."""
    s = str(text).strip()
    formats = (["%d/%m/%y %H:%M", "%d/%m/%y"] if day_first_short else []) + _DATE_FORMATS
    for f in formats:
        try:
            return dtm.datetime.strptime(s, f)
        except ValueError:
            pass
    return None


def initials(name: str) -> str:
    return "".join(w[0] for w in name.split()).upper()


# ---------------------------------------------------------------------------
# Quad collector
# ---------------------------------------------------------------------------
class KG:
    def __init__(self):
        self.quads: list[Quad] = []
        self.graph = DefaultGraph()         # all data in the default graph; provenance is explicit (reportedIn)
        self.notes: list[str] = []          # data-quality observations made while building

    def source(self, rel_path: str) -> NamedNode:
        """Document node of a source file."""
        return self.node("doc", rel_path)

    @staticmethod
    def node(kind: str, ident: str) -> NamedNode:
        return NamedNode(f"{RES}{kind}/{quote(str(ident), safe='/')}")

    @staticmethod
    def lit(v):
        if isinstance(v, bool):
            return Literal("true" if v else "false", datatype=NamedNode(XSD + "boolean"))
        if isinstance(v, int):
            return Literal(str(v), datatype=NamedNode(XSD + "integer"))
        if isinstance(v, float):
            return Literal(repr(round(v, 6)), datatype=NamedNode(XSD + "double"))
        if isinstance(v, dtm.datetime):
            return Literal(v.strftime("%Y-%m-%dT%H:%M:%S"), datatype=NamedNode(XSD + "dateTime"))
        if isinstance(v, dtm.date):
            return Literal(v.isoformat(), datatype=NamedNode(XSD + "date"))
        return Literal(str(v))

    def add(self, s: NamedNode, p: str, o):
        if o is None or o == "":
            return
        pred = p if isinstance(p, NamedNode) else NamedNode(BP + p)
        obj = o if isinstance(o, NamedNode) else self.lit(o)
        self.quads.append(Quad(s, pred, obj, self.graph))

    def typed(self, kind: str, ident: str, cls: str, label: str | None = None) -> NamedNode:
        n = self.node(kind, ident)
        self.quads.append(Quad(n, RDF_TYPE, NamedNode(BP + cls), self.graph))
        if label:
            self.quads.append(Quad(n, RDFS_LABEL, Literal(label), self.graph))
        return n

    # -- entity helpers -----------------------------------------------------
    def run(self, run_id: str) -> NamedNode:
        n = self.typed("run", run_id, "BioreactorRun" if run_id.startswith("FB") else "FlaskRun", run_id)
        self.quads.append(Quad(n, RDF_TYPE, NamedNode(BP + "Run"), self.graph))
        self.add(n, "runId", run_id)
        return n

    def sample(self, sample_id: str) -> NamedNode:
        n = self.typed("sample", sample_id, "Sample", sample_id)
        self.add(n, "sampleId", sample_id)
        run_id = norm_run(sample_id)
        self.add(n, "ofRun", self.node("run", run_id))
        return n

    def person(self, name_or_initials: str) -> NamedNode | None:
        s = str(name_or_initials).strip()
        if not s:
            return None
        ini = s if re.fullmatch(r"[A-Z]{2}", s) else initials(s)
        n = self.typed("person", ini, "Person")
        self.add(n, "initials", ini)
        if ini != s:
            self.add(n, "name", s)
            self.add(n, RDFS_LABEL, s)
        return n

    def strain(self, strain_id: str) -> NamedNode:
        n = self.typed("strain", strain_id, "Strain", strain_id)
        self.add(n, "strainId", strain_id)
        return n

    def clone(self, strain_id: str, clone: str) -> NamedNode:
        cid = f"{strain_id}-{clone.lower()}"
        n = self.typed("clone", cid, "Clone", cid)
        self.add(n, "cloneId", cid)
        self.add(n, "cloneLabel", clone.lower())
        self.add(n, "cloneOf", self.strain(strain_id))
        return n

    def strain_ref(self, subject: NamedNode, text: str, strain_p="usesStrain", clone_p="usesClone"):
        """'YX-002 (clone c3)', 'YX-002-c1', 'YX-001 c1', 'YX-005' -> strain link and, if stated, clone link."""
        m = re.search(r"(YX-\d{3})(?:\s*\(clone\s*|[-\s]+)?(c\d)?", str(text))
        if not m:
            return
        self.add(subject, strain_p, self.strain(m.group(1)))
        if m.group(2):
            self.add(subject, clone_p, self.clone(m.group(1), m.group(2)))

    def vessel(self, vessel_id: str) -> NamedNode:
        n = self.typed("vessel", vessel_id, "Vessel", vessel_id)
        self.add(n, "vesselId", vessel_id)
        return n

    def lot(self, lot_id: str) -> NamedNode:
        n = self.typed("lot", lot_id, "Lot", lot_id)
        self.add(n, "lotId", lot_id)
        return n

    def analyte(self, pref: str) -> NamedNode:
        return self.node("analyte", pref.replace(" ", "_"))

    def measurement(self, key: str, sample_id: str, analyte: str, raw, *, source_type: str, doc: NamedNode,
                    unit="g/L", **props) -> NamedNode | None:
        value, below, lod = to_number(raw)
        if value is None and not below:
            return None
        m = self.typed("meas", key, "Measurement")
        self.add(m, "ofSample", self.sample(sample_id))
        self.add(m, "analyte", self.analyte(analyte))
        self.add(m, "value", value)
        self.add(m, "belowLOD", below)
        self.add(m, "lodLimit", lod)
        self.add(m, "unit", unit)
        self.add(m, "sourceType", source_type)
        self.add(m, "reportedIn", doc)
        for p, v in props.items():
            self.add(m, p, v)
        return m


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8").replace("\r\n", "\n")


def md_fields(text: str) -> dict:
    """'**Author:** Elena Varga' header fields of the markdown documents."""
    return {m.group(1).strip(): m.group(2).strip() for m in re.finditer(r"^\*\*(.+?):\*\*\s*(.*)$", text, re.M)}


def md_sections(text: str) -> dict:
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith("## "):
            name = line[3:].strip()
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "\n".join(v).strip() for k, v in out.items()}


def sheet_rows(path: Path) -> dict:
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    return {ws.title: [list(r) for r in ws.iter_rows(values_only=True)] for ws in wb.worksheets}


def csv_rows(path: Path) -> tuple[list[str], list[dict]]:
    lines = read_text(path).splitlines()
    comments = [ln for ln in lines if ln.startswith("#")]
    return comments, list(csv.DictReader(ln for ln in lines if not ln.startswith("#")))


_MENTION_PATTERNS = [
    ("run", r"\b(?:FB|SF)[-_ ]?\d{1,3}\b", norm_run),
    ("strain", r"\bYX-\d{3}\b", lambda s: s),
    ("lot", r"\bC5-\d{4}\b", lambda s: s),
    ("deviation", r"\bDEV-\d{4}-\d{2}\b", lambda s: s),
    ("sequence", r"\bSEQ-\d{8}\b", lambda s: s),
    ("vessel", r"\bR[1-4]\b", lambda s: s),
]


def document(kg: KG, rel: str, doc_type: str, text: str = "") -> NamedNode:
    doc = kg.source(rel)
    kg.quads.append(Quad(doc, RDF_TYPE, NamedNode(BP + "Document"), kg.graph))
    kg.add(doc, "path", rel)
    kg.add(doc, "docType", doc_type)
    if text:
        title = re.search(r"^#\s+(.+)$", text, re.M)
        kg.add(doc, "title", title.group(1).strip() if title else None)
        f = md_fields(text)
        if f.get("Author"):
            kg.add(doc, "author", kg.person(f["Author"]))
        d = to_datetime(f.get("Date", ""))
        kg.add(doc, "date", d.date() if d else None)
        seen = set()
        for kind, pat, norm in _MENTION_PATTERNS:
            for m in re.finditer(pat, text):
                ident = norm(m.group(0))
                if ident and (kind, ident) not in seen:
                    seen.add((kind, ident))
                    kg.add(doc, "mentions", kg.node(kind, ident))
    return doc


# ---------------------------------------------------------------------------
# Source parsers
# ---------------------------------------------------------------------------
def parse_lims(kg: KG, corpus: Path, ctx: dict):
    rel = "lims/sample_register.csv"
    doc = document(kg, rel, "lims_register")
    _, rows = csv_rows(corpus / rel)
    first_sample = {}
    for r in rows:
        run_id, sid = norm_run(r["run_id"]), norm_sample(r["sample_id"])
        run = kg.run(run_id)
        kg.add(run, "describedIn", doc)
        if re.fullmatch(r"R\d", r["vessel"]):
            kg.add(run, "inVessel", kg.vessel(r["vessel"]))
        else:
            kg.add(run, "location", r["vessel"])
        s = kg.sample(sid)
        when = to_datetime(r["sampled_at"])
        kg.add(s, "sampledAt", when)
        kg.add(s, "sampleType", r["type"])
        kg.add(s, "volumeML", to_number(r["volume_mL"])[0])
        kg.add(s, "requestedAnalyses", r["requested"])
        kg.add(s, "storage", r["stored"])
        t = to_number(r["process_time_h"])[0]
        if t is not None:
            kg.add(s, "processTimeH", t)
            ctx["sample_at"][(run_id, round(t, 1))] = sid
            if run_id not in first_sample or t < first_sample[run_id][0]:
                first_sample[run_id] = (t, when)
    # Inoculation is not a LIMS field; it follows from any broth sample: sampled_at - process_time_h.
    for run_id, (t, when) in first_sample.items():
        ctx["inoc_lims"][run_id] = when - dtm.timedelta(hours=t)


def parse_hplc_raw(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "hplc" / "raw").glob("SEQ-*.csv")):
        rel = f"hplc/raw/{path.name}"
        doc = document(kg, rel, "hplc_raw_sequence")
        comments, rows = csv_rows(path)
        header = " ".join(comments)
        seq_id = re.search(r"SEQ-\d{8}", header).group(0)
        seq = kg.typed("sequence", seq_id, "HPLCSequence", seq_id)
        seq_date = dtm.datetime.strptime(seq_id[4:], "%Y%m%d").date()
        kg.add(seq, "sequenceId", seq_id)
        kg.add(seq, "runDate", seq_date)
        kg.add(seq, "reportedIn", doc)
        ver = re.search(r"SOP-ANA-003 v(\d)", header)
        kg.add(seq, "methodVersion", f"v{ver.group(1)}" if ver else None)
        rerun = re.search(r"re-analysis of (SEQ-\d{8})", header)
        if rerun:
            kg.add(seq, "rerunOf", kg.node("sequence", rerun.group(1)))
        ctx["sequences"][seq_date] = seq_id
        for r in rows:
            sid, an = norm_sample(r["SampleID"]), norm_analyte(r["Analyte"])
            vial, dil = float(r["Conc_vial_gL"]), float(r["Dilution"])
            below = r["Flag"].strip() == "<LOD"
            m = kg.measurement(f"raw/{seq_id}/{sid}/{an}", sid, an, "<LOD" if below else vial * dil,
                               source_type="hplc_raw", doc=doc, kind="HPLC", fromSequence=seq,
                               dilution=dil, analysedOn=seq_date, sopVersion=f"v{ver.group(1)}" if ver else None,
                               retentionTimeMin=float(r["RT_min"]), peakArea=float(r["Area"]))
            kg.add(m, "vialConcentration", vial)


def parse_hplc_processed(kg: KG, corpus: Path, ctx: dict):
    superseded = []                       # (measurement node, sample, analyte, re-run date)
    by_key = {}
    for path in sorted((corpus / "hplc" / "processed").glob("*.xlsx")):
        rel = f"hplc/processed/{path.name}"
        doc = document(kg, rel, "hplc_results_workbook")
        rows = sheet_rows(path)["Results"]
        col = {name: i for i, name in enumerate(rows[0])}
        for r in rows[1:]:
            if not r or r[col["SampleID"]] is None:
                continue
            sid, an = norm_sample(r[col["SampleID"]]), norm_analyte(r[col["Analyte"]])
            analysed = to_datetime(r[col["Analysed"]]).date()
            comment = r[col["Comment"]] or ""
            raw = "<LOD" if (r[col["Flag"]] or "").strip() == "<LOD" else r[col["Result_gL"]]
            seq_id = ctx["sequences"].get(analysed)
            m = kg.measurement(f"hplc/{sid}/{an}/{analysed:%Y%m%d}", sid, an, raw, source_type="hplc_processed",
                               doc=doc, kind="HPLC", analysedOn=analysed, sopVersion=r[col["SOP"]],
                               dilution=float(r[col["Dilution"]]), comment=comment or None,
                               fromSequence=kg.node("sequence", seq_id) if seq_id else None)
            if m is None:
                continue
            by_key[(sid, an, analysed)] = m
            sup = re.match(r"superseded by re-run (\d{4}-\d{2}-\d{2})", comment)
            if sup:
                superseded.append((m, sid, an, to_datetime(sup.group(1)).date()))
            else:
                ctx["current"].append(m)
            if comment.startswith("re-run"):
                kg.add(m, "isRerun", True)
    for m, sid, an, date in superseded:
        kg.add(m, "isCurrent", False)
        if (sid, an, date) in by_key:
            kg.add(m, "supersededBy", by_key[(sid, an, date)])
        else:
            kg.notes.append(f"HPLC {sid} {an}: marked superseded by a re-run on {date}, which is not in the results")
    for m in ctx["current"]:
        kg.add(m, "isCurrent", True)


def parse_ev_workbooks(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "spreadsheets").glob("EV_*.xlsx")):
        rel = f"spreadsheets/{path.name}"
        doc = document(kg, rel, "run_workbook")
        sheets = sheet_rows(path)
        info = {r[0]: r[1] for r in sheets["Run info"] if r and r[0]}
        run_id = norm_run(info["Run"])
        run = kg.run(run_id)
        kg.add(run, "describedIn", doc)
        kg.add(run, "inVessel", kg.vessel(re.match(r"R\d", info["Vessel"]).group(0)))
        kg.add(run, "operatedBy", kg.person(info["Operator"]))
        kg.strain_ref(run, info["Strain"])
        kg.add(run, "purpose", info.get("Purpose"))
        inoc = to_datetime(info["Inoculation"])
        kg.add(run, "inoculatedAt", inoc)
        ctx["inoc"][run_id] = inoc
        kg.add(run, "batchMedium", info.get("Batch medium"))
        kg.add(run, "feedDescription", info.get("Feed"))
        feed_setup(kg, run, info.get("Feed", ""))
        m = re.match(r"([\d.]+) h / ([\d.]+) mL/h", info.get("Feed start / rate", ""))
        if m:
            kg.add(run, "feedStartH", float(m.group(1)))
            kg.add(run, "feedRateSetpointMLPerH", float(m.group(2)))
        m = re.search(r"/ ([\d.]+) / (\d+) rpm", info.get("T / pH / agitation / air", ""))
        if m:
            kg.add(run, "pHSetpoint", float(m.group(1)))
            kg.add(run, "stirrerRpm", float(m.group(2)))
        kg.add(run, "note", info.get("Note"))
        rows = sheets["Offline data"]
        head = rows[0]
        for r in rows[1:]:
            if not r or r[0] is None:
                continue
            sid = norm_sample(r[0])
            s = kg.sample(sid)
            if isinstance(r[1], (int, float)):
                kg.add(s, "reactorVolumeML", to_number(r[3])[0])
            for i, name in enumerate(head):
                if name == "DW sem (g/L)" or i < 4:
                    continue
                an = norm_analyte(name)
                if an is None:
                    continue
                unit = "AU" if an == "OD600" else "g/L"
                m_ = kg.measurement(f"ev/{sid}/{an}", sid, an, r[i], source_type="run_workbook", doc=doc,
                                    unit=unit, kind="HPLC" if an not in ("OD600", "Dry weight") else an)
                if m_ is not None and an == "Dry weight":
                    kg.add(m_, "standardError", to_number(r[head.index("DW sem (g/L)")])[0])
        for r in sheets.get("Calculated", []):
            if r and r[0] is not None:
                reported_value(kg, run_id, doc, str(r[0]), r[1])


def reported_value(kg: KG, run_id: str, doc: NamedNode, label: str, raw, at_h=None, analyte=None):
    """A summary number that someone wrote down for a run (not a per-sample measurement)."""
    value, below, _ = to_number(raw)
    if value is None and not below:
        return
    key = f"{run_id}/{quote(str(doc.value).rsplit('/', 1)[-1])}/{re.sub(r'[^A-Za-z0-9]+', '_', label)}"
    n = kg.typed("reported", key, "ReportedValue", label)
    kg.add(n, "ofRun", kg.node("run", run_id))
    kg.add(n, "valueLabel", label)
    kg.add(n, "value", value)
    kg.add(n, "belowLOD", below)
    kg.add(n, "atProcessTimeH", at_h)
    kg.add(n, "reportedIn", doc)
    if analyte:
        kg.add(n, "analyte", kg.analyte(analyte))


def feed_setup(kg: KG, run: NamedNode, text: str):
    """Feed type, C5 fraction and stated lot from the feed descriptions used by the different people."""
    t = str(text)
    m = re.search(r"(\d+)\s*%(?:\s*of the)?\s*C5", t) or re.search(r"C5\s*(\d+)\s*%", t)
    if m:
        kg.add(run, "feedType", "C5")
        kg.add(run, "c5FractionPct", float(m.group(1)))
    elif re.fullmatch(r"\s*(feed\s+)?[AB]\s*", t, re.I):
        kg.add(run, "feedType", t.strip()[-1].upper())
    elif "Glc 100 g/L" in t:
        kg.add(run, "feedType", "A")
    elif "Glc 150 g/L" in t:
        kg.add(run, "feedType", "B")
    lot = re.search(r"C5-\d{4}", t)
    if lot:
        kg.add(run, "statedLot", kg.lot(lot.group(0)))


def parse_md_master(kg: KG, corpus: Path, ctx: dict):
    rel = "spreadsheets/MD_fermentation_master_2025.xlsx"
    doc = document(kg, rel, "master_log")
    sheets = sheet_rows(corpus / rel)
    owner = kg.person(Path(rel).name.split("_")[0])            # file naming: <initials>_...
    for r in sheets["Runs"][1:]:
        if not r or r[0] is None:
            continue
        run_id = norm_run(r[0])
        run = kg.run(run_id)
        kg.add(run, "describedIn", doc)
        kg.add(run, "operatedBy", owner)
        kg.add(run, "inVessel", kg.vessel(r[1]))
        kg.strain_ref(run, r[2])
        d = to_datetime(r[3], day_first_short=True)             # this log writes dd/mm/yy
        kg.add(run, "startDate", d.date() if d else None)
        kg.add(run, "batchMedium", r[4])
        kg.add(run, "feedDescription", r[5])
        feed_setup(kg, run, r[5])
        if r[6]:
            kg.add(run, "statedLot", kg.lot(r[6]))
        kg.add(run, "stirrerRpm", to_number(r[7])[0])
        kg.add(run, "note", r[8])
    head = sheets["Samples"][0]
    for r in sheets["Samples"][1:]:
        if not r or r[0] is None:
            continue
        sid = norm_sample(r[1], run=r[0])
        for i, name in enumerate(head):
            an = norm_analyte(name)
            if an is None:
                continue
            unit = "AU" if an == "OD600" else "g/L"
            kg.measurement(f"md/{sid}/{an}", sid, an, r[i], source_type="operator_log", doc=doc, unit=unit,
                           kind="HPLC" if an not in ("OD600", "Dry weight") else an, comment=r[head.index("Comment")])
    for r in sheets["Media & feeds"][1:]:
        if not r or r[0] is None:
            continue
        fb = kg.typed("feedbatch", r[1], "FeedBatch", r[1])
        kg.add(fb, "batchId", r[1])
        d = to_datetime(r[0], day_first_short=True)
        kg.add(fb, "preparedOn", d.date() if d else None)
        kg.add(fb, "contents", r[2])
        if r[3] and re.fullmatch(r"C5-\d{4}", str(r[3])):
            kg.add(fb, "lotUsed", kg.lot(r[3]))
        kg.add(fb, "volumeML", to_number(r[4])[0])
        kg.add(fb, "forRun", kg.node("run", norm_run(r[5])))
        kg.add(fb, "preparedBy", kg.person(r[6]))
        kg.add(fb, "reportedIn", doc)


def parse_pn_workbooks(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "spreadsheets").glob("PN_*.xlsx")):
        rel = f"spreadsheets/{path.name}"
        sheets = sheet_rows(path)
        desc = " ".join(str(r[0]) for r in sheets.get("Description", []) if r and r[0])
        doc = document(kg, rel, "study_workbook")
        kg.add(doc, "title", desc.split(".")[0] if desc else None)
        kg.add(doc, "description", desc)
        owner = kg.person(path.name.split("_")[0])
        for name, rows in sheets.items():
            if name in ("About", "Description", "Methods") or len(rows) < 3:
                continue
            head = rows[1]
            flask = "Flask" in head
            id_col = head.index("Flask" if flask else "Sample name")
            t_col = head.index("Time (h)")
            notes = defaultdict(list)
            for r in rows[2:]:
                if not r or r[id_col] is None:
                    continue
                run_id = norm_run(r[id_col])
                t = to_number(r[t_col])[0]
                sid = ctx["sample_at"].get((run_id, round(t, 1))) if flask else norm_sample(r[id_col])
                if r[0]:
                    notes[run_id].append(str(r[0]))
                if sid is None:
                    kg.notes.append(f"{rel} sheet {name}: no LIMS sample for {run_id} at {t} h")
                    continue
                s = kg.sample(sid)
                for i, col_name in enumerate(head):
                    if col_name == "Feed added (mL)":
                        kg.add(s, "feedAddedML", to_number(r[i])[0])
                    elif col_name == "Reactor volume (mL)":
                        kg.add(s, "reactorVolumeML", to_number(r[i])[0])
                    an = norm_analyte(col_name)
                    if an is None or i <= t_col:
                        continue
                    unit = "AU" if an == "OD600" else "g/L"
                    m_ = kg.measurement(f"pn/{sid}/{an}", sid, an, r[i], source_type="run_workbook", doc=doc,
                                        unit=unit, kind="HPLC" if an not in ("OD600", "Dry weight") else an)
                    if m_ is not None and an == "Dry weight" and "sem" in head:
                        kg.add(m_, "standardError", to_number(r[head.index("sem")])[0])
            for run_id, lines in notes.items():           # column A: strain/clone first, then the conditions
                run = kg.run(run_id)
                kg.add(run, "describedIn", doc)
                kg.add(run, "operatedBy", owner)
                kg.strain_ref(run, lines[0])
                cond = "; ".join(lines[1:])
                kg.add(run, "condition", cond)
                feed_setup(kg, run, cond)
                if "no CaCO3" in cond:
                    kg.add(run, "buffered", False)
                elif flask:
                    kg.add(run, "buffered", True)


def parse_calibration(kg: KG, corpus: Path, ctx: dict):
    """Pump calibration (reactor set-up) and feed densities (fed-batch SOP): the numbers behind the unit conversions."""
    rel = "equipment/reactor_setup.md"
    text = read_text(corpus / rel)
    doc = document(kg, rel, "equipment_setup", text)
    for pump, ml in re.findall(r"- (Acid|Base|Feed) pump: ([\d.]+) mL per count", text):
        ctx["ml_per_count"][pump.lower()] = float(ml)
        cal = kg.typed("calibration", f"BioCtrl-{pump.lower()}-pump", "PumpCalibration", f"BioCtrl {pump.lower()} pump")
        kg.add(cal, "pump", pump.lower())
        kg.add(cal, "controller", "BioCtrl B1")
        kg.add(cal, "mlPerCount", float(ml))
        kg.add(cal, "reportedIn", doc)
    for vessels, controller in re.findall(r"^\| (R\d(?:, R\d)*) \| ([^|]+?) \|", text, re.M):
        for v in vessels.split(", "):
            kg.add(kg.vessel(v), "controller", controller.strip())
    rel = "sops/SOP-FB-001_v3.md"
    text = read_text(corpus / rel)
    doc = kg.source(rel)
    for feed, dens in re.findall(r"(Feed [AB]|C5 feeds) ([\d.]+) g/mL", text):
        code = feed.replace("Feed ", "").replace(" feeds", "")
        ctx["density"][code] = float(dens)
        f = kg.typed("feed", code, "FeedType", f"Feed {code}")
        kg.add(f, "feedCode", code)
        kg.add(f, "densityGPerML", float(dens))
        kg.add(f, "reportedIn", doc)


def inoculation(ctx: dict, run_id: str) -> dtm.datetime | None:
    return ctx["inoc"].get(run_id) or ctx["inoc_lims"].get(run_id)


def parse_online(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "online").iterdir()):
        rel = f"online/{path.name}"
        lines = read_text(path).splitlines()
        if path.suffix == ".txt":                                  # BioCtrl, tab separated, pump counts
            doc = document(kg, rel, "online_log_bioctrl")
            head = "\n".join(lines[:6])
            run_id = norm_run(re.search(r"Batch name\s*:\s*(\S+)", head).group(1))
            run = kg.run(run_id)
            kg.add(run, "inVessel", kg.vessel(re.search(r"Vessel\s*:\s*(R\d)", head).group(1)))
            kg.add(run, "operatedBy", kg.person(re.search(r"Operator\s*:\s*(\w+)", head).group(1)))
            inoc = to_datetime(re.search(r"Inoculation\s*:\s*(.+)", head).group(1))
            start = to_datetime(re.search(r"Log start\s*:\s*(.+)", head).group(1))
            kg.add(run, "inoculatedAt", inoc)
            ctx["inoc"][run_id] = inoc
            cols = lines[5].split("\t")
            for i, ln in enumerate(lines[6:]):
                c = dict(zip(cols, ln.split("\t")))
                ts = to_datetime(f"{c['Date']} {c['Time']}")
                rec = kg.typed("online", f"{run_id}/{i:05d}", "OnlineRecord")
                kg.add(rec, "ofRun", run)
                kg.add(rec, "timestamp", ts)
                kg.add(rec, "processTimeH", (ts - inoc).total_seconds() / 3600)
                kg.add(rec, "loggerElapsedH", (ts - start).total_seconds() / 3600)
                kg.add(rec, "stirrerRpm", float(c["Stirrer[rpm]"]))
                kg.add(rec, "temperatureC", float(c["Temp[C]"]))
                kg.add(rec, "pH", float(c["pH"]))
                kg.add(rec, "pO2Pct", float(c["pO2[%]"]))
                for pump, name in (("acid", "Acid_Pump[cnt]"), ("base", "Base_Pump[cnt]"), ("feed", "Feed_Pump[cnt]")):
                    kg.add(rec, f"{pump}PumpCounts", int(c[name]))
                    kg.add(rec, {"acid": "acidML", "base": "baseML", "feed": "feedPumpML"}[pump],
                           int(c[name]) * ctx["ml_per_count"][pump])
        else:                                                       # VesselView, csv, mL and feed balance
            doc = document(kg, rel, "online_log_vesselview")
            head = lines[1]
            run_id = norm_run(re.search(r"recipe (\S+)", head).group(1))
            run = kg.run(run_id)
            kg.add(run, "inVessel", kg.vessel(re.search(r"unit (R\d)", head).group(1)))
            inoc = inoculation(ctx, run_id)
            if inoc is None:
                kg.notes.append(f"{rel}: no inoculation time known for {run_id}; processTimeH not set")
            for i, r in enumerate(csv.DictReader(lines[2:])):
                ts = to_datetime(r["Timestamp"])
                rec = kg.typed("online", f"{run_id}/{i:05d}", "OnlineRecord")
                kg.add(rec, "ofRun", run)
                kg.add(rec, "timestamp", ts)
                if inoc:
                    kg.add(rec, "processTimeH", (ts - inoc).total_seconds() / 3600)
                kg.add(rec, "loggerElapsedH", float(r["Elapsed_h"]))
                kg.add(rec, "stirrerRpm", float(r["AGIT_rpm"]))
                kg.add(rec, "temperatureC", float(r["TEMP_C"]))
                kg.add(rec, "pH", float(r["PH"]))
                kg.add(rec, "pO2Pct", float(r["DO_pct"]))
                kg.add(rec, "feedBalanceG", float(r["FEED_BAL_g"]))
                kg.add(rec, "baseML", float(r["BASE_mL"]))
                kg.add(rec, "acidML", float(r["ACID_mL"]))
        kg.add(run, "describedIn", doc)
        kg.add(run, "hasOnlineLog", doc)


def parse_offgas(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "offgas").glob("*.xlsx")):
        rel = f"offgas/{path.name}"
        doc = document(kg, rel, "offgas_export")
        rows = sheet_rows(path)["GA4"]
        run_id = norm_run(rows[1][1])
        run = kg.run(run_id)
        kg.add(run, "describedIn", doc)
        kg.add(run, "hasOffgasData", doc)
        ga = kg.typed("equipment", "GA-4", "Equipment", "GasLab GA-4 off-gas analyser")
        kg.add(ga, "airflowLPerMin", to_number(rows[1][4])[0])
        inlet = re.search(r"O2 ([\d.]+) %, CO2 ([\d.]+) %", str(rows[2][4]))
        if inlet:
            kg.add(ga, "inletO2Pct", float(inlet.group(1)))
            kg.add(ga, "inletCO2Pct", float(inlet.group(2)))
        inoc = inoculation(ctx, run_id)
        head = rows[4]
        for i, r in enumerate(rows[5:]):
            if not r or r[0] is None:
                continue
            c = dict(zip(head, r))
            ts = to_datetime(f"{c['Date']} {c['Time']}")
            rec = kg.typed("offgas", f"{run_id}/{i:05d}", "OffgasRecord")
            kg.add(rec, "ofRun", run)
            kg.add(rec, "timestamp", ts)
            if inoc:
                kg.add(rec, "processTimeH", (ts - inoc).total_seconds() / 3600)
            kg.add(rec, "analyserRelTimeS", float(c["RelTime[s]"]))
            kg.add(rec, "o2Pct", float(c["O2 [%]"]))
            kg.add(rec, "co2Pct", float(c["CO2 [%]"]))
            kg.add(rec, "n2Pct", float(c["N2 [%]"]))
            kg.add(rec, "arPct", float(c["Ar [%]"]))


def parse_strains(kg: KG, corpus: Path, ctx: dict):
    rel = "strains/strain_registry.xlsx"
    doc = document(kg, rel, "strain_registry")
    rows = sheet_rows(corpus / rel)["Strains"]
    for r in rows[1:]:
        if not r or r[0] is None:
            continue
        c = dict(zip(rows[0], r))
        s = kg.strain(c["Strain ID"])
        kg.add(s, "alias", c["Alias"])
        kg.add(s, SKOS_ALT, c["Alias"])
        if c["Parent"] and c["Parent"] != "-":
            kg.strain_ref(s, c["Parent"], strain_p="parentStrain", clone_p="builtOnClone")
        kg.add(s, "genotype", c["Genotype / change"])
        d = to_datetime(c["Constructed"])
        kg.add(s, "constructedOn", d.date() if d else None)
        kg.add(s, "constructedBy", kg.person(c["By"]))
        kg.add(s, "clonesTested", int(c["Clones tested"]) if c["Clones tested"] else None)
        if c["Selected clone"]:
            kg.add(s, "selectedClone", kg.clone(c["Strain ID"], c["Selected clone"]))
        kg.add(s, "cryoPosition", c["Cryo position"])
        kg.add(s, "status", c["Status"])
        kg.add(s, "describedIn", doc)
    for path in sorted((corpus / "strains").glob("construction_*.md")):
        rel = f"strains/{path.name}"
        text = read_text(path)
        doc = document(kg, rel, "construction_record", text)
        f = md_fields(text)
        s = kg.strain(re.search(r"YX-\d{3}", path.name).group(0))
        kg.add(s, "describedIn", doc)
        kg.add(s, "constructedBy", kg.person(f["Author"]) if f.get("Author") else None)
        if f.get("Parent") and "YX-" in f["Parent"]:
            kg.strain_ref(s, f["Parent"], strain_p="parentStrain", clone_p="builtOnClone")
        kg.add(s, "change", f.get("Change"))
        for part in re.findall(r"^- (.+)$", text, re.M):
            kg.add(s, "geneticPart", part.strip())
        sel = re.search(r"Clone (c\d) selected after .*? on (\d{4}-\d{2}-\d{2})", text)
        if sel:
            kg.add(s, "selectedClone", kg.clone(s.value.rsplit("/", 1)[-1], sel.group(1)))
            kg.add(s, "cloneSelectedOn", to_datetime(sel.group(2)).date())
        body = text.split("\n\n")
        kg.add(s, "constructionNote", " ".join(p.strip() for p in body if p.startswith(("Transformed", "Note:"))))


def parse_materials(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "materials").glob("CoA_*.md")):
        rel = f"materials/{path.name}"
        text = read_text(path)
        doc = document(kg, rel, "certificate_of_analysis", text)
        f = md_fields(text)
        lot_id = re.search(r"C5-\d{4}", path.name).group(0)
        lot = kg.lot(lot_id)
        kg.add(lot, "supplier", f.get("Supplier"))
        d = to_datetime(f.get("Date shipped", ""))
        kg.add(lot, "shippedOn", d.date() if d else None)
        kg.add(lot, "describedIn", doc)
        for comp, val in re.findall(r"^\| ([^|]+?) \| ([\d.]+) \|$", text, re.M):
            if comp == "pH":
                kg.add(lot, "pH", float(val))
                continue
            c = kg.typed("coa", f"{lot_id}/{re.sub(r'[^A-Za-z0-9]+', '_', comp)}", "CoAValue", f"{lot_id} {comp}")
            kg.add(c, "ofLot", lot)
            kg.add(c, "componentLabel", comp)
            an = norm_analyte(comp)
            kg.add(c, "analyte", kg.analyte(an) if an else None)
            kg.add(c, "value", float(val))
            kg.add(c, "unit", "g/L")
            kg.add(c, "reportedIn", doc)


def parse_deviations(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "deviations").glob("DEV-*.md")):
        rel = f"deviations/{path.name}"
        text = read_text(path)
        doc = document(kg, rel, "deviation_report", text)
        f, sec = md_fields(text), md_sections(text)
        dev_id = path.stem
        dev = kg.typed("deviation", dev_id, "Deviation", dev_id)
        kg.add(dev, "deviationId", dev_id)
        kg.add(dev, "title", re.search(r"^# [^:]+:\s*(.+)$", text, re.M).group(1))
        kg.add(dev, "raisedBy", kg.person(f["Author"]) if f.get("Author") else None)
        d = to_datetime(f.get("Date", ""))
        kg.add(dev, "date", d.date() if d else None)
        affected = f.get("Affected", "")
        kg.add(dev, "affected", affected)
        if norm_run(affected):
            kg.add(dev, "affectsRun", kg.node("run", norm_run(affected)))
        for seq in re.findall(r"SEQ-\d{8}", affected):
            kg.add(dev, "affectsSequence", kg.node("sequence", seq))
        for v in re.findall(r"\bR[1-4]\b", affected):
            kg.add(dev, "affectsVessel", kg.vessel(v))
        for label, p in (("Description", "description"), ("Root cause", "rootCause"),
                         ("Impact on data", "impactOnData"), ("Corrective action", "correctiveAction")):
            kg.add(dev, p, sec.get(label))
        kg.add(dev, "reportedIn", doc)


def parse_maintenance(kg: KG, corpus: Path, ctx: dict):
    rel = "equipment/maintenance_log.csv"
    text = read_text(corpus / rel)
    doc = document(kg, rel, "maintenance_log", text)
    _, rows = csv_rows(corpus / rel)
    for i, r in enumerate(rows):
        ev = kg.typed("maintenance", f"{r['date']}-{i:02d}", "MaintenanceEvent", f"{r['date']} {r['equipment']}")
        kg.add(ev, "date", to_datetime(r["date"]).date())
        kg.add(ev, "equipment", r["equipment"])
        kg.add(ev, "description", r["event"])
        kg.add(ev, "performedBy", kg.person(r["by"]))
        m = re.match(r"R(\d)(?:-R(\d))?", r["equipment"])
        if m:
            for v in range(int(m.group(1)), int(m.group(2) or m.group(1)) + 1):
                kg.add(ev, "concernsVessel", kg.vessel(f"R{v}"))
        kg.add(ev, "reportedIn", doc)


def parse_sops(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "sops").glob("SOP-*.md")):
        rel = f"sops/{path.name}"
        text = read_text(path)
        doc = document(kg, rel, "sop", text)
        f, sec = md_fields(text), md_sections(text)
        sop_id, ver = re.match(r"(SOP-[A-Z]+-\d+)_v(\d+)", path.stem).groups()
        sop = kg.typed("sop", f"{sop_id}_v{ver}", "SOP", f"{sop_id} v{ver}")
        kg.add(sop, "sopId", sop_id)
        kg.add(sop, "version", f"v{ver}")
        kg.add(sop, "title", re.search(r"^# (.+)$", text, re.M).group(1))
        kg.add(sop, "author", kg.person(f["Author"]) if f.get("Author") else None)
        d = to_datetime(f.get("Effective from", "") or f.get("Date", ""))
        kg.add(sop, "effectiveFrom", d.date() if d else None)
        scope = re.search(r"Quantification of (.+?) in culture", sec.get("Scope", ""))
        if scope:
            for name in scope.group(1).split(","):
                an = norm_analyte(name.strip())
                kg.add(sop, "coversAnalyte", kg.analyte(an) if an else None)
        lod = re.search(r"LOD ([\d.]+)-([\d.]+) g/L in the vial", text)
        if lod:
            kg.add(sop, "lodVialLowGPerL", float(lod.group(1)))
            kg.add(sop, "lodVialHighGPerL", float(lod.group(2)))
        changes = next((v for k, v in sec.items() if k.startswith("Changes")), None)
        kg.add(sop, "changeNote", changes)
        kg.add(sop, "reportedIn", doc)


def parse_notes(kg: KG, corpus: Path, ctx: dict):
    for path in sorted((corpus / "notes").glob("*.md")):
        rel = f"notes/{path.name}"
        text = read_text(path)
        doc = document(kg, rel, "lab_note", text)
        f = md_fields(text)
        run_id = norm_run(path.name) if path.name.startswith("FB") else None
        author_ini = path.stem.split("_")[-1]                       # file naming: <run>_<initials>.md
        kg.add(doc, "author", kg.person(f.get("Author") or author_ini))
        if run_id is None:
            continue
        run = kg.run(run_id)
        kg.add(run, "describedIn", doc)
        kg.add(run, "labNote", doc)
        kg.add(run, "operatedBy", kg.person(f.get("Author") or author_ini))
        if f.get("Vessel"):
            kg.add(run, "inVessel", kg.vessel(re.match(r"R\d", f["Vessel"]).group(0)))
        if f.get("Strain"):
            kg.strain_ref(run, f["Strain"])
        if f.get("Inoculated"):
            inoc = to_datetime(f["Inoculated"])
            kg.add(run, "inoculatedAt", inoc)
            ctx["inoc"].setdefault(run_id, inoc)
        title = re.search(r"^# (.+)$", text, re.M).group(1)
        m = re.fullmatch(r"FB\d+ (R\d) (YX-\d{3})", title)          # terse style: '# FB16 R1 YX-005'
        if m:
            kg.add(run, "inVessel", kg.vessel(m.group(1)))
            kg.strain_ref(run, m.group(2))
        m = re.search(r"^- (\d{2}/\d{2}/\d{2}) inoc (\d{2}:\d{2})", text, re.M)
        if m:
            inoc = to_datetime(f"{m.group(1)} {m.group(2)}", day_first_short=True)
            kg.add(run, "inoculatedAt", inoc)
            ctx["inoc"].setdefault(run_id, inoc)
        m = re.search(r"^- feed (.+?) on at ([\d,\.]+)\s*h(?:, lot (C5-\d{4}))?", text, re.M)
        if m:
            feed_setup(kg, run, m.group(1))
            kg.add(run, "feedStartH", to_number(m.group(2))[0])
            if m.group(3):
                kg.add(run, "statedLot", kg.lot(m.group(3)))
        m = re.search(r"^- end (\d+)\s*h: (.+)$", text, re.M)       # 'end 168h: Xol 62,2 / Xyl 11,2 / DW 3,3'
        if m:
            for item in m.group(2).split("/"):
                parts = item.strip().split()
                if len(parts) == 2 and norm_analyte(parts[0]):
                    reported_value(kg, run_id, doc, f"end {parts[0]}", parts[1], at_h=float(m.group(1)),
                                   analyte=norm_analyte(parts[0]))


def parse_narrative(kg: KG, corpus: Path, ctx: dict):
    for folder, doc_type in (("meetings", "meeting_minutes"), ("reports", "weekly_report")):
        for path in sorted((corpus / folder).glob("*.md")):
            document(kg, f"{folder}/{path.name}", doc_type, read_text(path))


def add_vocabulary(kg: KG):
    for pref, alts in ANALYTES.items():
        n = kg.typed("analyte", pref.replace(" ", "_"), "Analyte", pref)
        kg.add(n, "analyteName", pref)
        for a in alts:
            kg.add(n, SKOS_ALT, a)


def add_inoculation(kg: KG, ctx: dict):
    """Runs without an explicit inoculation record get the time implied by the LIMS register."""
    for run_id, when in ctx["inoc_lims"].items():
        run = kg.node("run", run_id)
        explicit = ctx["inoc"].get(run_id)
        if explicit is None:
            kg.add(run, "inoculatedAt", when)
            kg.add(run, "inoculationDerivedFromLIMS", True)
        elif abs((explicit - when).total_seconds()) > 600:
            kg.notes.append(f"{run_id}: stated inoculation {explicit} differs from the LIMS-implied {when}")
    for run_id in ctx["inoc_lims"]:
        when = inoculation(ctx, run_id)
        kg.add(kg.node("run", run_id), "startDate", when.date())


# ---------------------------------------------------------------------------
def build(corpus: Path = CORPUS, store_path: Path | None = STORE) -> tuple[Store, KG]:
    if "answer_key" in corpus.parts:
        raise ValueError("The answer key must never be indexed.")
    kg = KG()
    ctx = dict(sample_at={}, inoc={}, inoc_lims={}, sequences={}, current=[], ml_per_count={}, density={})
    add_vocabulary(kg)
    parse_lims(kg, corpus, ctx)                 # runs, samples, sample times
    parse_calibration(kg, corpus, ctx)          # needed before the online logs
    parse_ev_workbooks(kg, corpus, ctx)         # explicit inoculation times
    parse_notes(kg, corpus, ctx)
    parse_md_master(kg, corpus, ctx)
    parse_pn_workbooks(kg, corpus, ctx)
    parse_hplc_raw(kg, corpus, ctx)
    parse_hplc_processed(kg, corpus, ctx)
    parse_online(kg, corpus, ctx)
    parse_offgas(kg, corpus, ctx)
    parse_strains(kg, corpus, ctx)
    parse_materials(kg, corpus, ctx)
    parse_deviations(kg, corpus, ctx)
    parse_maintenance(kg, corpus, ctx)
    parse_sops(kg, corpus, ctx)
    parse_narrative(kg, corpus, ctx)
    add_inoculation(kg, ctx)

    if store_path is not None:
        shutil.rmtree(store_path, ignore_errors=True)
        store_path.parent.mkdir(parents=True, exist_ok=True)
        store = Store(str(store_path))
    else:
        store = Store()
    store.bulk_extend(set(kg.quads))
    store.load(path=str(ONTOLOGY), format=RdfFormat.TURTLE, to_graph=NamedNode(GRAPH + "ontology"))
    store.flush()
    return store, kg


def rebuild(corpus: Path = CORPUS, store_path: Path = STORE) -> int:
    """Build the on-disk store and release it, so that readers can open it. Returns the triple count.
    Callers rebuild at start-up (a few seconds) rather than trusting a store left by older code."""
    store, _ = build(corpus, store_path)
    n = len(store)
    del store
    gc.collect()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--corpus", type=Path, default=CORPUS)
    ap.add_argument("--store", type=Path, default=STORE)
    args = ap.parse_args()
    store, kg = build(args.corpus, args.store)
    print(f"{len(store)} triples in {args.store}")
    q = PREFIXES + "SELECT ?c (COUNT(DISTINCT ?s) AS ?n) WHERE { ?s a ?c . FILTER(STRSTARTS(STR(?c), STR(bp:))) } " \
                   "GROUP BY ?c ORDER BY DESC(?n)"
    for row in store.query(q):
        print(f"  {row['c'].value.split('#')[-1]:<18} {row['n'].value}")
    if kg.notes:
        print("Build notes:")
        for n in kg.notes:
            print("  -", n)


if __name__ == "__main__":
    main()
