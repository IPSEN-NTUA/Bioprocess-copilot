"""
The three systems under comparison. They are one agent loop with different tools, so the only thing that
varies between arms is the access path to the record:

  rag     search_documents  - passages of the raw files (text, tables as text)
  kg      sparql            - the knowledge graph built by build_kg.py
  hybrid  both

Each arm sees the same question, the same model, the same step budget and the same answer format.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import re
from pathlib import Path

from pyoxigraph import NamedNode, RdfFormat, Store

from copilot import llm
from copilot.build_kg import BP, ONTOLOGY, PREFIXES, RES, STORE
from copilot.rag import Index

ARMS = ("rag", "kg", "hybrid")
MAX_STEPS = 10
MAX_ROWS = 60
SPARQL_TIMEOUT_S = 60

COMMON = """You answer questions about the records of a fermentation development lab (a xylitol strain-engineering \
campaign: bioreactor runs FB-001.., shake flasks SF-001.., strains YX-001.., C5 hydrolysate lots, HPLC analytics).

Rules:
- Use only what the tools return. Do not use outside knowledge for facts about this lab; general science and \
arithmetic are allowed and should be shown.
- The records come from several people and instruments and can disagree or have been corrected later. When they \
do, say which value is valid and why.
- If the records cannot answer the question, or the question assumes something the records contradict, say so and \
explain what is missing or wrong. Do not guess.
- Work step by step with the tools, then stop calling tools and write the answer.
- Finish with one line starting with 'ANSWER:' that states the final answer concisely, with units. Write dates as \
YYYY-MM-DD and identifiers in full (FB-002, YX-003, C5-2501, clone c3).
"""

RAG_NOTE = """
You have a search tool over the lab's files (notes, spreadsheets, instrument exports, SOPs, minutes). It returns \
passages; tables come as text rows with their original headers. Different files use different identifier styles, \
date formats, abbreviations and units. Search again with other wording or identifiers when the first passages do \
not contain what you need.
"""

KG_NOTE = """
You have a SPARQL tool over a knowledge graph built from the lab's files. The schema follows. All prefixes are \
already declared; write SELECT queries only.

Conventions:
- Identifiers are canonical: runId 'FB-002', sampleId 'FB006-S08', strainId 'YX-003', cloneId 'YX-002-c1', \
lotId 'C5-2501', sequenceId 'SEQ-20250321', deviationId 'DEV-2025-02'.
- Compare numbers with FILTER (e.g. FILTER(?t = 96), FILTER(ABS(?t - 47.5) < 0.1)), never by writing a number in \
a triple pattern. Dates: FILTER(?d >= "2025-05-01"^^xsd:date).
- A sample usually has several Measurement nodes for the same analyte, one per source. Always select \
bp:sourceType, and for 'hplc_processed' also bp:isCurrent.
- bp:value is absent when bp:belowLOD is true, so use OPTIONAL { ?m bp:value ?v } when below-LOD results matter.
- Aggregates (MIN, MAX, AVG, COUNT, GROUP BY) and arithmetic are supported. Use LIMIT when exploring.

Two patterns:
  # all statements about one analyte in one sample
  SELECT ?src ?current ?v ?below WHERE {
    ?m bp:ofSample/bp:sampleId "FB001-S04" ; bp:analyte/bp:analyteName "Xylose" ; bp:sourceType ?src ; bp:belowLOD ?below .
    OPTIONAL { ?m bp:value ?v } OPTIONAL { ?m bp:isCurrent ?current } }
  # a statistic of an online signal in a process-time window
  SELECT (AVG(?x) AS ?mean) WHERE {
    ?o a bp:OnlineRecord ; bp:ofRun/bp:runId "FB-001" ; bp:processTimeH ?t ; bp:temperatureC ?x . FILTER(?t >= 10 && ?t <= 20) }

SCHEMA
"""

SEARCH_TOOL = {
    "name": "search_documents",
    "description": "Search the lab's files and return the most relevant passages.",
    "parameters": {"type": "object", "properties": {"query": {"type": "string", "description": "What to look for."}},
                   "required": ["query"]},
}
SPARQL_TOOL = {
    "name": "sparql",
    "description": "Run a SPARQL SELECT query on the knowledge graph. Prefixes bp:, res:, rdfs:, skos:, xsd: are declared.",
    "parameters": {"type": "object", "properties": {"query": {"type": "string"}}, "required": ["query"]},
}


# ---------------------------------------------------------------------------
# Schema text for the prompt, generated from the ontology file
# ---------------------------------------------------------------------------
def schema_text(ontology: Path = ONTOLOGY) -> str:
    st = Store()
    st.load(path=str(ontology), format=RdfFormat.TURTLE)
    q = """PREFIX owl: <http://www.w3.org/2002/07/owl#> PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
           SELECT ?x ?kind ?dom ?rng ?sup ?c WHERE { ?x a ?kind . FILTER(?kind IN (owl:Class, owl:ObjectProperty, owl:DatatypeProperty))
           OPTIONAL { ?x rdfs:domain ?dom } OPTIONAL { ?x rdfs:range ?rng } OPTIONAL { ?x rdfs:subClassOf ?sup }
           OPTIONAL { ?x rdfs:comment ?c } }"""
    short = lambda t: None if t is None else "bp:" + t.value.split("#")[-1] if t.value.startswith(BP) else "xsd:" + t.value.split("#")[-1]
    classes, props = [], []
    for r in st.query(q):
        name, comment = short(r["x"]), (r["c"].value if r["c"] else "")
        if r["kind"].value.endswith("Class"):
            sup = f" (a {short(r['sup'])})" if r["sup"] else ""
            classes.append(f"{name}{sup} - {comment}")
        else:
            dom = short(r["dom"]) or "any"
            rng = short(r["rng"]) or ("node" if r["kind"].value.endswith("ObjectProperty") else "literal")
            props.append((dom, f"  {name} -> {rng}" + (f"   # {comment}" if comment else "")))
    out = ["Classes:"] + sorted(classes) + ["", "Properties (grouped by the class they belong to; 'any' = several classes):"]
    for dom in sorted({d for d, _ in props}):
        out.append(f"{dom}:")
        out += sorted(line for d, line in props if d == dom)
    return "\n".join(out)


# ---------------------------------------------------------------------------
# SPARQL execution in a worker process, so that a runaway query can be stopped
# ---------------------------------------------------------------------------
def _fmt(term) -> str:
    if term is None:
        return ""
    v = term.value
    if isinstance(term, NamedNode):
        return "res:" + v[len(RES):] if v.startswith(RES) else "bp:" + v[len(BP):] if v.startswith(BP) else f"<{v}>"
    dt = getattr(term, "datatype", None)
    if dt is not None and dt.value.endswith("double"):
        return f"{float(v):.6g}"
    return v


def _sparql_worker(store_path: str, conn):
    store = Store.read_only(store_path)
    while True:
        query = conn.recv()
        if query is None:
            return
        try:
            if not re.match(r"\s*(#.*\n\s*)*(PREFIX\s+\S+\s*<[^>]*>\s*)*(SELECT|ASK)\b", query, re.I):
                raise ValueError("only SELECT or ASK queries are allowed")
            res = store.query(PREFIXES + query)
            if isinstance(res, bool) or not hasattr(res, "variables"):
                conn.send(str(bool(res)))
                continue
            names = [v.value for v in res.variables]
            rows = []
            for i, sol in enumerate(res):
                if i >= MAX_ROWS:
                    rows.append(f"... more rows not shown (limit {MAX_ROWS}); aggregate or filter in the query")
                    break
                rows.append("\t".join(_fmt(sol[n]) for n in names))
            conn.send("\t".join(names) + "\n" + ("\n".join(rows) if rows else "(no rows)"))
        except Exception as e:                                     # query errors go back to the model
            conn.send(f"SPARQL error: {str(e)[:600]}")


class Sparql:
    def __init__(self, store_path: Path = STORE):
        if not Path(store_path).exists():
            raise FileNotFoundError(f"No KG at {store_path}; run: python -m copilot.build_kg")
        self.store_path, self.proc = str(store_path), None

    def _start(self):
        ctx = mp.get_context("spawn")
        self.conn, child = ctx.Pipe()
        self.proc = ctx.Process(target=_sparql_worker, args=(self.store_path, child), daemon=True)
        self.proc.start()

    def run(self, query: str) -> str:
        if self.proc is None or not self.proc.is_alive():
            self._start()
        self.conn.send(query)
        try:
            if self.conn.poll(SPARQL_TIMEOUT_S):
                return self.conn.recv()
        except (EOFError, OSError):
            self.proc = None
            return "SPARQL error: the query engine stopped unexpectedly; try a simpler query"
        self.proc.terminate()
        self.proc = None
        return f"SPARQL error: query stopped after {SPARQL_TIMEOUT_S} s; make it more selective"

    def close(self):
        if self.proc is not None and self.proc.is_alive():
            self.conn.send(None)
            self.proc.join(2)


# ---------------------------------------------------------------------------
# Agent loop
# ---------------------------------------------------------------------------
class Copilot:
    def __init__(self, index: Index | None, sparql: Sparql | None, model: str | None = None,
                 retrieval: str = "hybrid", top_k: int = 8):
        self.index, self.sparql, self.model = index, sparql, model or llm.MODEL
        self.retrieval, self.top_k = retrieval, top_k
        self._schema = schema_text() if sparql else ""

    def _search(self, query: str) -> str:
        hits = self.index.search(query, k=self.top_k, mode=self.retrieval)
        return "\n\n---\n\n".join(c.text for c in hits)

    def _call(self, name: str, args: dict) -> str:
        if name == "search_documents" and self.index is not None:
            return self._search(str(args.get("query", "")))
        if name == "sparql" and self.sparql is not None:
            return self.sparql.run(str(args.get("query", "")))
        return f"Unknown tool {name!r}."

    def ask(self, question: str, arm: str) -> dict:
        """Returns {'answer', 'final', 'steps', 'trace'}; 'final' is the text of the ANSWER: line."""
        if arm not in ARMS:
            raise ValueError(arm)
        tools = ([SEARCH_TOOL] if arm in ("rag", "hybrid") else []) + ([SPARQL_TOOL] if arm in ("kg", "hybrid") else [])
        system = COMMON + (RAG_NOTE if arm in ("rag", "hybrid") else "") + \
            (KG_NOTE + self._schema if arm in ("kg", "hybrid") else "")
        if arm == "hybrid":
            system += ("\n\nThe graph holds the structured data with identifiers, units and time bases already reconciled; "
                       "the files hold the narrative (reasons, decisions, procedures). Use both.")
        messages = [{"role": "system", "content": system}, {"role": "user", "content": question}]
        trace = []
        if arm == "rag":                                            # classic RAG: first retrieval is on the question
            passages = self._search(question)
            trace.append({"tool": "search_documents", "args": {"query": question}, "auto": True, "result": passages})
            messages[1]["content"] = f"{question}\n\nPassages retrieved for the question:\n\n{passages}"
        answer = ""
        for step in range(MAX_STEPS):
            last = step == MAX_STEPS - 1
            if last:
                messages.append({"role": "user", "content": "Step budget used up. Give your best answer now, ending with the ANSWER: line."})
            msg = llm.chat(messages, tools=None if last else tools, model=self.model)
            messages.append(msg)
            if not msg.get("tool_calls"):
                answer = msg["content"]
                break
            for call in msg["tool_calls"]:
                result = self._call(call["name"], call["arguments"])
                trace.append({"tool": call["name"], "args": call["arguments"], "result": result})
                messages.append({"role": "tool", "tool_call_id": call["id"], "name": call["name"], "content": result})
        final = re.findall(r"ANSWER:\s*(.+)", answer)
        return {"answer": answer, "final": final[-1].strip() if final else answer.strip(), "steps": len(trace), "trace": trace}


def trace_json(trace: list[dict], max_chars: int = 4000) -> str:
    """Trace with long tool results shortened, for the results file."""
    return json.dumps([{**t, "result": t["result"][:max_chars]} for t in trace])
