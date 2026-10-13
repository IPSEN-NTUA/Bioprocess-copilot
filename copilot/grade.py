"""
Grading, following the rules in benchmark/README.md.

  exact, numeric, ranking   deterministic, from the ANSWER: line
  key_points, unanswerable, rubric   an LLM judge compares the full answer with the gold answer

Deterministic rules:
  numeric   correct if a number in the ANSWER: line is within the question's tolerance of the gold value
            (the count of numbers in that line is recorded, so lenient passes can be audited)
  exact     the gold value appears in the ANSWER: line (numbers as numbers, text case-insensitively)
  ranking   the run identifiers in the ANSWER: line, in order, equal the gold order

Judge-graded items carry the judge's reasons. They should be spot-checked by a person before the numbers
are reported; `judge` in each record names the model that graded it.
"""
from __future__ import annotations

import json
import re

from copilot import llm
from copilot.build_kg import norm_run

_NUM = r"(?<![\w.])-?\d+(?:[.,]\d+)?"


def numbers(text: str) -> list[float]:
    """Numbers in a text. Identifier digits (FB-002, C5-2501, SEQ-20250321, c3, 2025-06-01) are skipped."""
    cleaned = re.sub(r"\b(?:FB|SF|YX|DEV|SEQ|SOP|C5|R|S|c|v)[-_]?[\w-]*\d[\w-]*|\b\d{4}-\d{2}-\d{2}\b", " ", text)
    return [float(n.replace(",", ".")) for n in re.findall(_NUM, cleaned)]


def _norm(s) -> str:
    return re.sub(r"\s+", " ", str(s)).strip().lower()


def grade_numeric(q: dict, final: str) -> dict:
    found = numbers(final)
    tol = (q["tolerance"] or 0) + 1e-9
    hit = any(abs(n - float(q["value"])) <= tol for n in found)
    return {"correct": hit, "score": float(hit), "detail": {"numbers_in_answer": found}}


def grade_exact(q: dict, final: str) -> dict:
    v = q["value"]
    if isinstance(v, (int, float)):
        hit = any(abs(n - float(v)) < 1e-9 for n in numbers(final))
    else:
        hit = re.search(rf"(?<![\w-]){re.escape(_norm(v))}(?![\w-])", _norm(final)) is not None
    return {"correct": hit, "score": float(hit), "detail": {}}


def grade_ranking(q: dict, final: str) -> dict:
    order = []
    for m in re.finditer(r"\b(?:FB|SF)[-_ ]?\d{1,3}\b", final):
        r = norm_run(m.group(0))
        if r not in order:
            order.append(r)
    hit = order == list(q["value"])
    return {"correct": hit, "score": float(hit), "detail": {"order": order}}


JUDGE_HEAD = """You are grading an answer to a question about a laboratory's records. Judge only against the \
reference; do not use your own knowledge of the subject. Be strict about numbers, identifiers and dates, and \
lenient about wording, rounding to the reference's precision, and extra correct detail.

Question:
{question}

Reference answer:
{gold}

Answer to grade:
{answer}
"""

JUDGE_KEY_POINTS = JUDGE_HEAD + """
The answer is correct only if it conveys every one of these key points:
{points}

Reply with JSON: {{"points": [{{"point": "...", "present": true or false}}], "contradicts_reference": true or false, \
"reason": "one sentence"}}"""

JUDGE_UNANSWERABLE = JUDGE_HEAD + """
The records cannot answer this question as asked, or its premise is false. The answer is correct only if it says \
so and gives a reason consistent with the reference. An answer that supplies a value, or accepts the false premise, \
is wrong.

Reply with JSON: {{"declines_or_corrects_premise": true or false, "reason_consistent": true or false, \
"reason": "one sentence"}}"""

JUDGE_RUBRIC = JUDGE_HEAD + """
Rubric items (an item counts only if the answer proposes it AND supports it with evidence from the records):
{rubric}

Penalised claims (any one of these in the answer fails it):
{penalties}

Reply with JSON: {{"items": [{{"id": "...", "proposed_with_evidence": true or false}}], \
"penalised_claims_made": ["..."], "reason": "one sentence"}}"""


def grade_judged(q: dict, answer: str, model: str | None = None) -> dict:
    model = model or llm.JUDGE_MODEL
    base = dict(question=q["question"], gold=q["gold_answer"], answer=answer or "(empty)")
    if q["grading"] == "unanswerable":
        j = llm.chat_json(JUDGE_UNANSWERABLE.format(**base), model=model)
        hit = bool(j.get("declines_or_corrects_premise")) and bool(j.get("reason_consistent"))
        return {"correct": hit, "score": float(hit), "detail": j, "judge": model}
    if q["grading"] == "rubric":
        rubric = "\n".join(f"- {r['id']}: {r['suggestion']} (evidence: {r['evidence']})" for r in q["rubric"])
        penalties = "\n".join(f"- {p}" for p in (q["penalties"] or [])) or "- none"
        j = llm.chat_json(JUDGE_RUBRIC.format(rubric=rubric, penalties=penalties, **base), model=model)
        n = sum(1 for it in j.get("items", []) if it.get("proposed_with_evidence"))
        hit = n >= 3 and not j.get("penalised_claims_made")
        return {"correct": hit, "score": n / max(len(q["rubric"]), 1), "detail": j, "judge": model}
    points = q["key_points"] or [q["gold_answer"]]
    j = llm.chat_json(JUDGE_KEY_POINTS.format(points="\n".join(f"- {p}" for p in points), **base), model=model)
    present = [bool(p.get("present")) for p in j.get("points", [])]
    hit = len(present) >= len(points) and all(present) and not j.get("contradicts_reference")
    return {"correct": hit, "score": sum(present) / max(len(points), 1), "detail": j, "judge": model}


def grade(q: dict, result: dict, judge_model: str | None = None) -> dict:
    """result is the dict returned by Copilot.ask (needs 'answer' and 'final')."""
    kind = q["grading"]
    if kind == "numeric":
        return grade_numeric(q, result["final"])
    if kind == "exact":
        return grade_exact(q, result["final"])
    if kind == "ranking":
        return grade_ranking(q, result["final"])
    try:
        return grade_judged(q, result["answer"], judge_model)
    except llm.LLMError as e:
        return {"correct": None, "score": None, "detail": {"judge_error": str(e)}, "judge": judge_model or llm.JUDGE_MODEL}


def load_questions(path) -> list[dict]:
    return json.loads(open(path).read())
