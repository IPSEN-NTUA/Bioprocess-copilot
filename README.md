# Bioprocess-copilot: RAG vs knowledge graph vs hybrid on bioprocess records

Benchmark and simulator for a paper comparing retrieval-augmented generation (RAG), knowledge-graph (KG) and
hybrid question answering over heterogeneous bioprocess data, prepared for Frontiers in Chemical Engineering
(Research Topic "Semantic Interoperability in Industrial Operations").

The benchmark is a synthetic but simulator-grounded record of a xylitol strain-engineering campaign in a fictional
lab: lab notes in four personal styles, spreadsheets, controller and off-gas exports, raw and processed HPLC data,
LIMS, SOPs, deviations, minutes and a strain registry, with deliberate faults. It comes with 76 questions whose gold
answers are computed from the documents themselves.

## Layout

```
simulator/      fed-batch xylitol process simulator (extended kinetic model 31) and strain-design layer
  simulator.py          process model: growth, xylitol, overflow, O2/CO2, pH, hydrolysate inhibitors
  strain_design.py      genetic modifications -> parameter changes (network response analysis)
  MODEL.md              model description
  validate_bc18.py, tune_bc18.py, compare_fit.py, predict_runs.py
                        calibration / validation against partner data (not included, see data/README.md)
benchmark/      corpus generator, the pilot corpus, answer key and questions (see benchmark/README.md)
  generator/            code
  pilot/corpus/         145 documents - the only thing a system under test may index
  pilot/answer_key/     simulator truth - never index
  pilot/questions/      76 questions (json / csv / md)
paper/          LaTeX source (Frontiers template)
data/           placeholder for the partner data (git-ignored)
```

## Quick start

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

cd benchmark
python generator/questions.py          # regenerate the questions from the committed corpus (seconds)
python generator/verify_questions.py   # independent check of every gold answer
```

Rebuild the whole pilot from the simulator (about 6 minutes):

```
cd benchmark
python generator/build_pilot.py && python generator/audit_pilot.py && python generator/design_checks.py \
  && python generator/questions.py && python generator/verify_questions.py
```

Expected: audit 151 PASS / 0 FAIL / 5 NOTE (the notes are the planted faults), question check 129 PASS / 0 FAIL.
The build is deterministic (seed 20250217).

## Questions

| Set | n | What it tests |
|---|---|---|
| core | 50 | difficulty tiers L1 (lookup) to L5 (synthesis), plus 7 unanswerable / false-premise questions |
| IO | 13 | semantic interoperability: identifiers, below-LOD encodings, date and number formats, time bases, units, superseded values |
| CA | 10 | calculations: unit conversion, RQ from off-gas, yield from a mass balance, acid-base, power analysis |
| XD | 3 | experiment design ("what next?"), graded with a rubric; each rubric item also has a simulator check |

Every question is tagged with the interoperability problems it involves, so results can be reported per problem.
Details and grading rules: `benchmark/README.md`.

## Notes for contributors

- Do not commit partner data or plots made from it (`.gitignore` covers `data/partner/` and the validation plots).
- The corpus and questions are generated: change the generator, rebuild, and commit the regenerated `pilot/`.
- Never put `benchmark/pilot/answer_key/` into a retrieval index.
- The simulator code and `MODEL.md` still name the partner and its run IDs; scrub them before the repository is made
  public.
