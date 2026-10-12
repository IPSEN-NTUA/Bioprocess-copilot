# Bioprocess document benchmark (pilot)

Synthetic, simulator-grounded record of a fictional xylitol strain-engineering campaign (Feb to Jul 2025, lab
"Larkspur Bio"), plus 76 question/answer pairs for comparing RAG, knowledge-graph and hybrid retrieval.
Every corpus file is marked as synthetic.

```
benchmark/
  generator/            code (needs ../simulator)
    world.py            cast, vessels, C5 lots, strain plan
    campaign.py         runs, flasks, faults, clone decisions (taken from measured data)
    measure.py          measurement layer: HPLC, dry weight, OD, online signals, off-gas
    data_files.py       file renderers (xlsx / csv / instrument exports)
    narrative.py        notes, minutes, reports, SOPs, CoAs, deviations, construction records
    build_pilot.py      builds pilot/corpus and pilot/answer_key        (seed 20250217)
    audit_pilot.py      physics + cross-document consistency audit
    design_checks.py    simulates candidate next experiments on the true strains (answer key only)
    questions.py        builds pilot/questions (gold answers read from the corpus only)
    questions_extra.py  interoperability (IO), calculation (CA) and experiment-design (XD) sets
    verify_questions.py re-derives the gold answers from raw files and checks the unanswerables
  pilot/
    corpus/             145 files - the only thing a system under test may index
    answer_key/         simulator truth, faults, noise settings - never index
    questions/          questions.json / .csv / .md
```

Rebuild everything (about 6 min; build_pilot clears answer_key, so design_checks must run before questions):

```
python generator/build_pilot.py && python generator/audit_pilot.py && python generator/design_checks.py \
  && python generator/questions.py && python generator/verify_questions.py
```

Current pilot: audit 151 PASS / 0 FAIL / 5 NOTE (the NOTEs are the deliberate faults); question verification
129 PASS / 0 FAIL.

## Corpus

| Folder | Files | Content |
|---|---|---|
| online | 25 | controller logs; R1/R2 tab-separated BioCtrl export, R3/R4 VesselView csv |
| offgas | 14 | gas analyser xlsx, R1/R2 only |
| hplc/raw, hplc/processed | 26 + 6 | one csv per HPLC sequence; monthly long-format results workbooks |
| lims | 1 | sample register |
| spreadsheets | 14 | EV per-run workbooks, MD master log, PN flask and lot-comparison workbooks |
| notes | 27 | lab notes per run (four personal styles) and two flask-screen notes |
| meetings, reports | 7 + 6 | minutes, weekly reports |
| sops, materials, strains, deviations, equipment | 3 + 3 + 7 + 4 + 2 | SOPs, CoAs, strain registry + construction records, deviations, maintenance log, reactor set-up (pump calibration, time bases) |

Deliberate faults (answer_key/faults.csv): old preculture (FB-002), feed-pump over-delivery (FB-005), pO2 probe
drift (FB-009), power cut (FB-010), antifoam (FB-014), wrong C5 lot (FB-024), biased HPLC xylitol standard
(SEQ-20250321, re-run SEQ-20250409), S07 sample swap between FB-012 and FB-013, and one decimal-point slip in the
MD master log (FB-015 S08).

## Questions

76 questions in four sets (pilot/questions/questions.md lists them all with evidence and grading).

| Set | n | Levels | What it tests |
|---|---|---|---|
| core | 50 | L1 9, L2 9, L3 10, L4 9, L5 6, U 7 | lookup to synthesis, plus unanswerable / false-premise questions |
| IO | 13 | L2-L4 | semantic interoperability: the answer needs heterogeneous sources reconciled |
| CA | 10 | L2-L5 | calculations: unit conversion, gas balance (RQ), mass balance (yield), acid-base, power analysis |
| XD | 3 | L4-L5 | experiment design: what to do next, graded with a rubric |

Difficulty levels: L1 one fact from one document; L2 one document plus a small calculation; L3 joins across 2-3
documents; L4 provenance, fault diagnosis or harder calculation; L5 campaign-level synthesis; U not answerable.

Every question carries `interop` tags (also on the core set) so results can be reported per interoperability
problem: `id_resolution`, `lod_encoding`, `number_date_format`, `time_alignment`, `vocabulary`, `units_conversion`,
`versioning`, `source_authority`, `implicit_context`. IO questions also record the `trap` - the answer a system gives
when it does not reconcile the sources (e.g. reading Elapsed_h as process time, 03/02/25 as 2 March).

Fields: `id, set, tier, category, question, gold_answer, value, tolerance, grading, key_points, evidence, hops,
interop, trap, rubric, penalties`. Grading: `exact`, `numeric` (|answer - value| <= tolerance), `key_points` (all
present), `ranking`, `unanswerable` (states that the corpus cannot answer and why), `rubric` (XD: count of rubric
items proposed with correct corpus evidence - at least three - and no penalised claim).

Gold answers are the values recorded in the corpus (measured, with noise), not the simulator truth. Where
measurement noise makes a conclusion unsafe (P3 strain ranking, P6 titre differences, enzyme effect), the gold
answer says so; the truth in answer_key/ confirms these caveats (true P3 means differ by < 4 g/L).

### Experiment design and the simulator

The XD rubric items are built from the corpus. Each also has a `simulator_check` from design_checks.py, which runs the
suggested experiment on the true strains (noise-free). This allows a second, stronger evaluation: not only whether a
system's suggestion is grounded, but whether it would have been informative. Main outcomes for this pilot:

- a higher feed rate (3 instead of 2 mL/h) raises titre (YX-006 145 -> 168 g/L) and enlarges strain differences
  (YX-003 vs YX-001 5.6 -> 8.9 g/L); at standard conditions ~10 runs per strain are needed to resolve them;
- with feed B, 700 rpm or a pO2 cascade avoids the oxygen limitation seen at 500 rpm;
- YX-006 shows no titre loss up to 100 % C5 fraction - titre is set by the xylose top-up, not by the inhibitors;
- the planned C5-2503 no-enzyme repeat would show only ~1.4 g/L at 5 U/L; a dose-response (20 U/L: +5 g/L) is more
  informative;
- PN's idea that the PPP strain needs more glucose is not supported (advantage +2.8 g/L on feed B vs +5.6 on feed A).

These checks are in answer_key/design_checks.json and must never be indexed.
