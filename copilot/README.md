# Copilot: the three systems under test

RAG, knowledge-graph and hybrid question answering over `benchmark/pilot/corpus`, with a runner and grader for
the 76 benchmark questions. Nothing here reads `benchmark/pilot/answer_key/`; both the graph builder and the
indexer refuse that folder.

```
copilot/
  build_kg.py            corpus -> knowledge graph (Oxigraph), deterministic, no LLM
  ontology/bioprocess.ttl  schema; also the source of the schema text in the KG prompt
  rag.py                 corpus -> text chunks -> dense + BM25 index
  arms.py                the three arms: one agent loop, different tools
  llm.py                 chat and embeddings (Ollama or any OpenAI-compatible endpoint)
  grade.py               grading rules of benchmark/README.md
  run_eval.py            run, grade, tabulate
  ask.py                 one question, three answers side by side (demo)
  tests/                 graph structure, competency queries against the gold answers, grader, index
```

## Set-up

```
pip install -r requirements.txt -r copilot/requirements.txt
ollama pull nomic-embed-text            # embeddings for the RAG index
python -m pytest copilot/tests          # 42 tests, no LLM needed
```

## Running

```
python -m copilot.build_kg                                   # ~4 s, 525k triples, prints class counts
python -m copilot.rag                                        # builds and caches the index (~5 min the first time)
python -m copilot.ask "Which C5 lot was actually used in FB-024?" --show-steps
python -m copilot.run_eval --name pilot-<model>              # all questions x three arms -> copilot/results/
```

The model is set with environment variables (see `llm.py`):

```
# local Ollama (default)
export COPILOT_LLM_MODEL=qwen2.5:14b
# any OpenAI-compatible endpoint
export COPILOT_LLM_API=openai COPILOT_LLM_BASE_URL=https://.../v1 COPILOT_LLM_API_KEY=... COPILOT_LLM_MODEL=...
export COPILOT_JUDGE_MODEL=...          # optional: a different model for grading
```

`run_eval` rebuilds the graph at start, writes one JSON line per question and arm (answer, tool trace, grade) and
can be interrupted and resumed. `summary.md` breaks the scores down by question set, difficulty tier, grading type
and interoperability tag.

## Design

**One loop, three toolsets.** The arms differ only in how they reach the record: `search_documents` (RAG),
`sparql` (KG) or both (hybrid). Model, prompt rules, step budget (10) and answer format are the same. The RAG arm
gets the passages for the question up front, as in classic RAG, and may search again.

**The graph reconciles form, not content.** `build_kg.py` normalises identifiers, analyte names, number and date
formats, time bases and units, and links re-run results to the values they supersede. It does not decide which of
two disagreeing sources is right: every source's value is kept as its own `bp:Measurement` with `bp:sourceType`
and `bp:reportedIn`. So the planted faults (decimal slip, stale values in the master log, lot mix-up, sample swap)
are still there to be found.

**The RAG index normalises nothing.** Chunks keep each file's own headers, identifiers and formats; tables become
row windows with the header repeated. Retrieval interleaves a dense and a BM25 ranking, at most three chunks per
file.

**Narrative stays out of the graph.** Notes, minutes, reports and SOP text are only `bp:Document` nodes with
`bp:mentions` links to the runs, strains, lots and deviations they name. Reasons and decisions are therefore
reachable by the RAG and hybrid arms only. Header fields and the structured sections of deviation reports are in
the graph.

## What the tests establish

`tests/test_kg.py` checks the graph against the corpus (counts, every sample resolves to a LIMS sample, every
predicate is declared in the ontology) and runs reference SPARQL for about 35 benchmark questions, comparing with
the gold answers in `questions.json`: lookups, aggregates, final titres, identifier and vocabulary resolution,
below-LOD handling, superseded values, source conflicts, time alignment, unit conversion and strain lineage.
These reference queries are in the tests only. The KG arm never sees them; it writes its own from the ontology.

## Limits to keep in mind for the paper

- The graph builder was written by someone who had read the questions. The parsers are per document type, not
  per question, but this is a real threat to validity; a held-out question set written after the builder was
  frozen would remove it.
- 40 of the 76 questions (key_points, rubric) and the 7 unanswerable ones are graded by an LLM judge. The judge
  has not been validated against human grading yet.
- `numeric` grading passes if any number in the `ANSWER:` line is within tolerance; the numbers found are stored
  in each record so that lenient passes can be audited.
- The KG prompt contains two generic example queries (not benchmark questions).
- Only smoke runs with `llama3.2` (3B) have been made, to test the plumbing. It cannot write SPARQL reliably and
  its scores mean nothing.
