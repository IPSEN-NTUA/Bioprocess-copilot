"""
Run the benchmark questions through the RAG, KG and hybrid arms, grade them and tabulate the results.

  python -m copilot.run_eval --name pilot                       all 76 questions, three arms
  python -m copilot.run_eval --name smoke --ids L1-01 IO-04     selected questions
  python -m copilot.run_eval --name pilot --report-only         re-tabulate an existing run

Results go to copilot/results/<name>/: answers.jsonl (one line per question and arm, with the tool trace),
summary.md and summary.csv. A run can be interrupted and resumed; finished question/arm pairs are skipped.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from collections import defaultdict
from pathlib import Path

from copilot import llm
from copilot.arms import ARMS, Copilot, Sparql, trace_json
from copilot.build_kg import ROOT, STORE, rebuild
from copilot.grade import grade, load_questions
from copilot.rag import load_index

QUESTIONS = ROOT / "benchmark" / "pilot" / "questions" / "questions.json"
RESULTS = ROOT / "copilot" / "results"


def run(name: str, arms: list[str], ids: list[str] | None, sets: list[str] | None, model: str, judge: str,
        retrieval: str, top_k: int):
    out_dir = RESULTS / name
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "answers.jsonl"
    done = set()
    if path.exists():
        done = {(r["id"], r["arm"]) for r in map(json.loads, path.read_text().splitlines()) if r.get("error") is None}
    questions = [q for q in load_questions(QUESTIONS)
                 if (not ids or q["id"] in ids) and (not sets or q["set"] in sets)]
    need_rag = any(a in ("rag", "hybrid") for a in arms)
    need_kg = any(a in ("kg", "hybrid") for a in arms)
    index = load_index(embed=retrieval != "lexical") if need_rag else None
    if need_kg:
        print(f"knowledge graph rebuilt: {rebuild()} triples")
    sparql = Sparql(STORE) if need_kg else None
    bot = Copilot(index, sparql, model=model, retrieval=retrieval, top_k=top_k)
    (out_dir / "config.json").write_text(json.dumps(dict(
        model=model, judge=judge, api=llm.API, retrieval=retrieval, top_k=top_k, embed_model=llm.EMBED_MODEL,
        num_ctx=llm.NUM_CTX, arms=arms), indent=1))
    todo = [(q, a) for q in questions for a in arms if (q["id"], a) not in done]
    print(f"{len(todo)} question/arm pairs to run ({len(done)} already done), model {model}, judge {judge}")
    with path.open("a") as f:
        for i, (q, arm) in enumerate(todo, 1):
            t0 = time.time()
            rec = {"id": q["id"], "arm": arm, "set": q["set"], "tier": q["tier"], "grading": q["grading"],
                   "interop": q["interop"] or [], "model": model, "error": None}
            try:
                res = bot.ask(q["question"], arm)
                rec.update(answer=res["answer"], final=res["final"], steps=res["steps"], trace=trace_json(res["trace"]))
                rec.update(grade(q, res, judge))
            except llm.LLMError as e:
                rec.update(error=str(e), correct=None, score=None)
            rec["seconds"] = round(time.time() - t0, 1)
            f.write(json.dumps(rec) + "\n")
            f.flush()
            mark = {True: "ok", False: "--", None: "??"}[rec.get("correct")]
            print(f"[{i}/{len(todo)}] {q['id']:<6} {arm:<6} {mark}  {rec['seconds']:>6.1f}s  {str(rec.get('final', rec['error']))[:90]}")
    if sparql:
        sparql.close()
    report(name)


def _table(records: list[dict], key, arms: list[str]) -> list[list]:
    groups = defaultdict(lambda: defaultdict(list))
    for r in records:
        for k in key(r):
            groups[k][r["arm"]].append(r)
    rows = []
    for k in sorted(groups):
        row = [k, max(len(v) for v in groups[k].values())]
        for a in arms:
            rs = [r for r in groups[k].get(a, []) if r.get("correct") is not None]
            row.append(f"{sum(r['correct'] for r in rs)}/{len(rs)}" if rs else "-")
        rows.append(row)
    return rows


def report(name: str):
    out_dir = RESULTS / name
    latest = {}
    for r in map(json.loads, (out_dir / "answers.jsonl").read_text().splitlines()):
        latest[(r["id"], r["arm"])] = r                       # a re-run of a failed pair replaces the old line
    records = list(latest.values())
    arms = [a for a in ARMS if any(r["arm"] == a for r in records)]
    cfg = json.loads((out_dir / "config.json").read_text()) if (out_dir / "config.json").exists() else {}
    sections = [
        ("All questions", lambda r: ["all"]),
        ("By question set", lambda r: [r["set"]]),
        ("By difficulty tier", lambda r: [r["tier"]]),
        ("By grading type", lambda r: [r["grading"]]),
        ("By interoperability problem", lambda r: r["interop"] or ["(none tagged)"]),
    ]
    md = [f"# Results: {name}", "",
          f"Model `{cfg.get('model')}`, judge `{cfg.get('judge')}`, retrieval `{cfg.get('retrieval')}` (top {cfg.get('top_k')}). "
          "Cells are correct / graded.", ""]
    flat = []
    for title, key in sections:
        rows = _table(records, key, arms)
        md += [f"## {title}", "", "| | n | " + " | ".join(arms) + " |", "|---|---|" + "---|" * len(arms)]
        md += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows] + [""]
        flat += [[title, *row] for row in rows]
    errors = [r for r in records if r.get("error")]
    ungraded = [r for r in records if r.get("correct") is None and not r.get("error")]
    judged = [r for r in records if r.get("judge")]
    md += ["## Notes", "",
           f"- {len(judged)} of {len(records)} answers were graded by the LLM judge (key_points, unanswerable, rubric); "
           "spot-check them before reporting.",
           f"- {len(errors)} runs failed with an LLM error; {len(ungraded)} could not be graded.",
           "- Mean tool calls per question: " + ", ".join(
               f"{a} {sum(r.get('steps', 0) for r in records if r['arm'] == a) / max(sum(r['arm'] == a for r in records), 1):.1f}"
               for a in arms), ""]
    (out_dir / "summary.md").write_text("\n".join(md))
    with (out_dir / "summary.csv").open("w", newline="") as f:
        csv.writer(f).writerows([["breakdown", "group", "n", *arms], *flat])
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True, help="results folder name")
    ap.add_argument("--arms", nargs="+", default=list(ARMS), choices=ARMS)
    ap.add_argument("--ids", nargs="+", help="question ids to run (default: all)")
    ap.add_argument("--sets", nargs="+", choices=["core", "IO", "CA", "XD"])
    ap.add_argument("--model", default=llm.MODEL)
    ap.add_argument("--judge", default=None, help="judge model (default: COPILOT_JUDGE_MODEL or --model)")
    ap.add_argument("--retrieval", default="hybrid", choices=["dense", "lexical", "hybrid"])
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    if a.report_only:
        report(a.name)
    else:
        judge = a.judge or os.environ.get("COPILOT_JUDGE_MODEL") or a.model
        run(a.name, a.arms, a.ids, a.sets, a.model, judge, a.retrieval, a.top_k)


if __name__ == "__main__":
    main()
