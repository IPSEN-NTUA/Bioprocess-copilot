"""
Ask one question and see the three arms side by side (the comparison demo).

  python -m copilot.ask "What xylitol concentration was measured in FB-006 at 96 h?"
  python -m copilot.ask --arms kg hybrid --show-steps "Which C5 lot was actually used in FB-024?"
"""
from __future__ import annotations

import argparse
import textwrap

from copilot import llm
from copilot.arms import ARMS, Copilot, Sparql
from copilot.build_kg import rebuild
from copilot.rag import load_index

TITLES = {"rag": "[1] RAG - passages of the raw files", "kg": "[2] Knowledge graph - SPARQL",
          "hybrid": "[3] Hybrid - graph + passages"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("question")
    ap.add_argument("--arms", nargs="+", default=list(ARMS), choices=ARMS)
    ap.add_argument("--model", default=llm.MODEL)
    ap.add_argument("--retrieval", default="hybrid", choices=["dense", "lexical", "hybrid"])
    ap.add_argument("--show-steps", action="store_true", help="print the searches and SPARQL queries each arm made")
    a = ap.parse_args()
    index = load_index(embed=a.retrieval != "lexical") if any(x in ("rag", "hybrid") for x in a.arms) else None
    sparql = None
    if any(x in ("kg", "hybrid") for x in a.arms):
        rebuild()
        sparql = Sparql()
    bot = Copilot(index, sparql, model=a.model, retrieval=a.retrieval)
    for arm in a.arms:
        res = bot.ask(a.question, arm)
        print("\n" + "=" * 100 + f"\n{TITLES[arm]}   ({res['steps']} tool calls, model {a.model})\n" + "=" * 100)
        if a.show_steps:
            for t in res["trace"]:
                arg = t["args"].get("query", t["args"])
                print(f"  > {t['tool']}: {' '.join(str(arg).split())[:300]}")
                print(textwrap.indent(t["result"][:500], "      "))
            print("-" * 100)
        print(res["answer"].strip())
    if sparql:
        sparql.close()


if __name__ == "__main__":
    main()
