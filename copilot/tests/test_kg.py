"""
Is the knowledge graph built correctly?

Two kinds of checks:
  * structure - identifiers, counts and vocabulary agree with the corpus and the ontology;
  * competency - hand-written reference SPARQL reproduces the benchmark's gold answers.

The reference queries live only here. They are never shown to the systems under test: the KG arm has to
write its own queries from the ontology.
"""
import json
import statistics

import pytest

from copilot import build_kg as B

QUESTIONS = {q["id"]: q for q in json.loads((B.ROOT / "benchmark/pilot/questions/questions.json").read_text())}


@pytest.fixture(scope="session")
def built():
    return B.build(store_path=None)


@pytest.fixture(scope="session")
def store(built):
    return built[0]


def rows(store, query):
    out = []
    res = store.query(B.PREFIXES + query)
    names = [v.value for v in res.variables]
    for sol in res:
        out.append(tuple(_py(sol[n]) for n in names))
    return out


def _py(term):
    if term is None:
        return None
    dt = getattr(term, "datatype", None)
    if dt is not None and dt.value.endswith(("double", "integer", "decimal")):
        return float(term.value)
    return term.value.rsplit("/", 1)[-1] if term.value.startswith(B.RES) else term.value


def one(store, query):
    r = rows(store, query)
    assert len(r) == 1, r
    return r[0][0] if len(r[0]) == 1 else r[0]


def gold(qid):
    return QUESTIONS[qid]["value"]


def close(qid, value):
    q = QUESTIONS[qid]
    assert value == pytest.approx(q["value"], abs=(q["tolerance"] or 0) + 1e-9), (qid, value, q["gold_answer"])


FINAL = """
  ?s bp:ofRun ?run ; bp:sampleType "broth" ; bp:processTimeH ?t .
  FILTER NOT EXISTS { ?s2 bp:ofRun ?run ; bp:sampleType "broth" ; bp:processTimeH ?t2 . FILTER(?t2 > ?t) }
"""


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("raw,expected", [("FB02", "FB-002"), ("fb001", "FB-001"), ("FB-013", "FB-013"),
                                          ("FB003_R3.csv", "FB-003"), ("SF-7", "SF-007"), ("YX-002", None)])
def test_norm_run(raw, expected):
    assert B.norm_run(raw) == expected


def test_norm_sample():
    assert B.norm_sample("FB005-S07") == "FB005-S07"
    assert B.norm_sample("S08", run="FB15") == "FB015-S08"
    assert B.norm_sample("FB004-F") == "FB004-F"


@pytest.mark.parametrize("raw,expected", [("3,3", (3.3, False, None)), ("<0.05", (None, True, 0.05)),
                                          ("<LOD", (None, True, None)), ("n.d.", (None, True, None)),
                                          (None, (None, False, None)), (12, (12.0, False, None))])
def test_to_number(raw, expected):
    assert B.to_number(raw) == expected


def test_day_first_dates_only_when_declared():
    assert B.to_datetime("03/02/25", day_first_short=True).date().isoformat() == "2025-02-03"
    assert B.to_datetime("03/02/25") is None
    assert B.to_datetime("23 Jun 2025").date().isoformat() == "2025-06-23"


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------
def test_no_build_notes(built):
    assert built[1].notes == []


def test_every_predicate_and_class_is_in_the_ontology(store):
    used = rows(store, "SELECT DISTINCT ?p WHERE { ?s ?p ?o FILTER(STRSTARTS(STR(?p), STR(bp:))) }")
    declared = {r[0] for r in rows(store, f"SELECT ?p WHERE {{ GRAPH <{B.GRAPH}ontology> {{ ?p a ?t }} }}")}
    assert {p for (p,) in used} - declared == set()
    classes = rows(store, "SELECT DISTINCT ?c WHERE { ?s a ?c }")
    assert {c for (c,) in classes} - declared == set()


def test_entity_counts_match_the_corpus(store):
    lims = [ln for ln in (B.CORPUS / "lims/sample_register.csv").read_text().splitlines()[2:] if ln]
    assert one(store, "SELECT (COUNT(DISTINCT ?s) AS ?n) WHERE { ?s a bp:Sample }") == len(lims)
    assert one(store, "SELECT (COUNT(DISTINCT ?r) AS ?n) WHERE { ?r a bp:BioreactorRun }") == 25
    assert one(store, "SELECT (COUNT(DISTINCT ?r) AS ?n) WHERE { ?r a bp:FlaskRun }") == 35
    assert one(store, "SELECT (COUNT(DISTINCT ?d) AS ?n) WHERE { ?d a bp:Document }") == \
        sum(1 for p in B.CORPUS.rglob("*") if p.is_file())


def test_identifiers_resolve_to_single_nodes(store):
    # every sample that any source mentions is a sample of the LIMS register, and every run has its basics
    assert rows(store, "SELECT ?s WHERE { ?m bp:ofSample ?s . FILTER NOT EXISTS { ?s bp:sampledAt ?t } }") == []
    assert rows(store, """SELECT ?r WHERE { ?r a bp:Run . FILTER NOT EXISTS { ?r bp:inoculatedAt ?t ; bp:usesStrain ?x } }""") == []
    assert rows(store, """SELECT ?r WHERE { ?r a bp:BioreactorRun .
                          FILTER NOT EXISTS { ?r bp:inVessel ?v ; bp:operatedBy ?p ; bp:hasOnlineLog ?d } }""") == []
    for single in ("usesStrain", "inVessel", "operatedBy", "inoculatedAt"):
        assert rows(store, f"SELECT ?r WHERE {{ ?r a bp:Run ; bp:{single} ?a, ?b . FILTER(?a != ?b) }}") == [], single


def test_process_time_agrees_between_sources(store):
    # an online row logged at the wall-clock time of a sample must carry that sample's process time
    # (within one 10-min step); a wrong time base would be off by hours
    drift = rows(store, """SELECT ?sid ?t ?to WHERE { ?s bp:sampleId ?sid ; bp:ofRun ?run ; bp:processTimeH ?t ; bp:sampledAt ?when .
                           ?o a bp:OnlineRecord ; bp:ofRun ?run ; bp:processTimeH ?to ; bp:timestamp ?when .
                           FILTER(ABS(?to - ?t) > 0.11) }""")
    assert drift == [], drift[:5]


# ---------------------------------------------------------------------------
# Competency: lookups
# ---------------------------------------------------------------------------
def test_lookups(store):
    assert one(store, 'SELECT ?id WHERE { ?r bp:runId "FB-013" ; bp:usesStrain/bp:strainId ?id }') == gold("L1-01")
    assert one(store, 'SELECT ?id WHERE { ?r bp:runId "FB-022" ; bp:inVessel/bp:vesselId ?id }') == gold("L1-02")
    close("L1-03", one(store, 'SELECT ?v WHERE { ?c bp:ofLot/bp:lotId "C5-2503" ; bp:analyte/bp:analyteName "Formate" ; bp:value ?v }'))
    assert one(store, 'SELECT ?n WHERE { ?s bp:strainId "YX-004" ; bp:constructedBy/bp:name ?n }') == gold("L1-04")
    assert one(store, 'SELECT ?d WHERE { ?s bp:sopId "SOP-ANA-003" ; bp:version "v2" ; bp:effectiveFrom ?d }') == gold("L1-05")
    assert one(store, 'SELECT ?c WHERE { ?s bp:strainId "YX-006" ; bp:selectedClone/bp:cloneLabel ?c }') == gold("L1-07")
    close("L1-08", one(store, 'SELECT ?v WHERE { ?x bp:ofRun/bp:runId "FB-001" ; bp:valueLabel "Final xylitol (g/L)" ; bp:value ?v }'))
    assert one(store, 'SELECT DISTINCT ?d WHERE { ?m bp:ofSample/bp:sampleId "FB018-S06" ; bp:dilution ?d }') == gold("L1-09")
    assert one(store, 'SELECT ?t WHERE { ?r bp:runId "FB-002" ; bp:inoculatedAt ?t }') == "2025-02-03T10:15:00"


def test_aggregates(store):
    close("L2-04", one(store, """SELECT (MAX(?v) AS ?m) WHERE { ?x bp:ofSample/bp:ofRun/bp:runId "FB-018" ;
                                 bp:analyte/bp:analyteName "Dry weight" ; bp:value ?v }"""))
    close("L2-07", one(store, """SELECT (COUNT(?s) AS ?n) WHERE { ?s a bp:HPLCSequence ; bp:runDate ?d .
                                 FILTER(?d >= "2025-05-01"^^xsd:date && ?d < "2025-06-01"^^xsd:date) }"""))
    close("L2-09", one(store, """SELECT (COUNT(DISTINCT ?r) AS ?n) WHERE { ?r a bp:FlaskRun ; bp:describedIn/bp:path ?p .
                                 FILTER(CONTAINS(?p, "clone_screen_Mar2025")) }"""))
    runs = [r for (r,) in rows(store, 'SELECT ?id WHERE { ?r a bp:BioreactorRun ; bp:runId ?id ; bp:usesStrain/bp:strainId "YX-006" } ORDER BY ?id')]
    assert runs == [f"FB-{n:03d}" for n in (14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25)]
    top = rows(store, """SELECT ?n (COUNT(DISTINCT ?r) AS ?k) WHERE { ?r a bp:BioreactorRun ; bp:operatedBy/bp:name ?n }
                         GROUP BY ?n ORDER BY DESC(?k)""")
    assert top[0] == ("Mark Doherty", 13) and dict(top) == {"Mark Doherty": 13, "Elena Varga": 10, "Priya Nair": 2}
    offgas = rows(store, "SELECT DISTINCT ?v WHERE { ?r bp:hasOffgasData ?d ; bp:inVessel/bp:vesselId ?v } ORDER BY ?v")
    assert offgas == [("R1",), ("R2",)]
    assert one(store, "SELECT (COUNT(DISTINCT ?r) AS ?n) WHERE { ?r bp:hasOffgasData ?d }") == 14


def test_final_titres(store):
    q = f"""SELECT ?id ?v WHERE {{ ?run bp:runId ?id . {FINAL}
            ?m bp:ofSample ?s ; bp:sourceType "operator_log" ; bp:analyte/bp:analyteName "Xylitol" ; bp:value ?v .
            FILTER(?id >= "FB-012" && ?id <= "FB-017") }} ORDER BY DESC(?v)"""
    assert [r[0] for r in rows(store, q)] == gold("L3-02")
    q = f"""SELECT ?id ?v WHERE {{ ?run bp:runId ?id . {FINAL}
            ?m bp:ofSample ?s ; bp:sourceType "hplc_processed" ; bp:isCurrent true ;
               bp:analyte/bp:analyteName "Xylitol" ; bp:value ?v . FILTER(?id <= "FB-003") }}"""
    finals = dict(rows(store, q))
    close("L3-03", statistics.mean(round(v, 1) for v in finals.values()))
    q = f"""SELECT ?id ?v WHERE {{ VALUES ?id {{ "FB-022" "FB-023" }} ?run bp:runId ?id . {FINAL}
            ?m bp:ofSample ?s ; bp:sourceType "run_workbook" ; bp:analyte/bp:analyteName "Xylitol" ; bp:value ?v }}"""
    f = dict(rows(store, q))
    close("L3-07", round(f["FB-023"], 1) - round(f["FB-022"], 1))


# ---------------------------------------------------------------------------
# Competency: interoperability
# ---------------------------------------------------------------------------
def test_identifier_and_vocabulary_resolution(store):
    # 'FB06' / 'S08' / 'Xol' in the master log and 'FB006-S08' / 'Xylitol' in the HPLC results are one sample, one analyte
    got = rows(store, """SELECT ?src ?v WHERE { ?s bp:ofRun/bp:runId "FB-006" ; bp:processTimeH ?t . FILTER(?t = 96)
                         ?m bp:ofSample ?s ; bp:analyte/bp:analyteName "Xylitol" ; bp:sourceType ?src ; bp:value ?v }""")
    assert {s for s, _ in got} == {"operator_log", "hplc_processed", "hplc_raw"}
    for _, v in got:
        close("IO-01", v)
    close("IO-06", one(store, f"""SELECT ?v WHERE {{ ?run bp:runId "FB-016" . {FINAL}
                                 ?m bp:ofSample ?s ; bp:sourceType "hplc_processed" ; bp:analyte/bp:analyteName "Arabinose" ; bp:value ?v }}"""))
    seq = one(store, """SELECT ?when ?seq WHERE { ?s bp:sampleId "FB015-S08" ; bp:sampledAt ?when .
                        ?m bp:ofSample ?s ; bp:sourceType "hplc_processed" ; bp:fromSequence/bp:sequenceId ?seq } LIMIT 1""")
    assert seq == ("2025-05-09T11:15:00", "SEQ-20250513")


def test_below_lod_is_a_result_not_a_gap(store):
    got = rows(store, """SELECT ?src ?below ?v WHERE { ?s bp:ofRun/bp:runId "FB-002" ; bp:processTimeH ?t . FILTER(?t = 23.5)
                         ?m bp:ofSample ?s ; bp:analyte/bp:analyteName "Glycerol" ; bp:sourceType ?src ; bp:belowLOD ?below .
                         OPTIONAL { ?m bp:value ?v } }""")
    assert len(got) == 3 and all(below == "true" and v is None for _, below, v in got)
    first = one(store, """SELECT (MIN(?t) AS ?first) WHERE { ?s bp:ofRun/bp:runId "FB-003" ; bp:processTimeH ?t .
                          ?m bp:ofSample ?s ; bp:sourceType "hplc_processed" ; bp:belowLOD true ; bp:analyte/bp:analyteName "Glucose" }""")
    close("L2-01", first)


def test_versioning(store):
    cur = one(store, """SELECT ?v WHERE { ?m bp:ofSample/bp:sampleId "FB006-S05" ; bp:sourceType "hplc_processed" ;
                        bp:isCurrent true ; bp:analyte/bp:analyteName "Xylitol" ; bp:value ?v }""")
    close("IO-07", cur)
    old, new, seq = one(store, """SELECT ?old ?new ?seq WHERE { ?m bp:ofSample/bp:sampleId "FB005-S07" ; bp:isCurrent false ;
                                  bp:value ?old ; bp:supersededBy ?n . ?n bp:value ?new ; bp:fromSequence/bp:sequenceId ?seq }""")
    assert (round(old, 2), round(new, 2), seq) == (103.88, 97.82, "SEQ-20250409")
    assert one(store, "SELECT (COUNT(?m) AS ?n) WHERE { ?m bp:isCurrent false }") == 16
    # succinate was not part of method v1, so there is nothing for FB-006 - not even a below-LOD result
    assert rows(store, 'SELECT ?m WHERE { ?m bp:ofSample/bp:ofRun/bp:runId "FB-006" ; bp:analyte/bp:analyteName "Succinate" }') == []


def test_source_conflicts_are_kept(store):
    q = """SELECT ?sid ?log ?hplc WHERE {
             ?a bp:ofSample ?s ; bp:sourceType "operator_log" ; bp:analyte ?an ; bp:value ?log .
             ?b bp:ofSample ?s ; bp:sourceType "hplc_processed" ; bp:isCurrent true ; bp:analyte ?an ; bp:value ?hplc .
             ?an bp:analyteName "Xylitol" . ?s bp:sampleId ?sid . FILTER(ABS(?log - ?hplc) > 0.02) } ORDER BY ?sid"""
    diff = {sid: (log, hplc) for sid, log, hplc in rows(store, q)}
    assert diff["FB015-S08"] == pytest.approx((490.7, 49.07), abs=0.01)          # decimal-point slip
    assert len(diff) == 9                                                        # plus 8 values not updated after the re-run
    stated, used = one(store, """SELECT ?stated ?used WHERE { ?r bp:runId "FB-024" ; bp:statedLot/bp:lotId ?stated .
                                 ?fb bp:forRun ?r ; bp:lotUsed/bp:lotId ?used }""")
    assert (stated, used) == ("C5-2503", "C5-2502")


def test_time_alignment(store):
    v = one(store, """SELECT ?v WHERE { ?o a bp:OnlineRecord ; bp:ofRun/bp:runId "FB-021" ; bp:processTimeH ?t ; bp:pO2Pct ?v .
                      FILTER(ABS(?t - 47.5) < 0.01) }""")
    close("IO-04", v)
    lo = rows(store, """SELECT ?t ?v WHERE { ?o a bp:OnlineRecord ; bp:ofRun/bp:runId "FB-005" ; bp:processTimeH ?t ; bp:pO2Pct ?v .
                        FILTER(?t >= 2 && ?t <= 40) } ORDER BY ?v LIMIT 1""")[0]
    assert lo == pytest.approx((15.0, 8.7), abs=0.5)
    peak = rows(store, """SELECT ?t ?v WHERE { ?o a bp:OffgasRecord ; bp:ofRun/bp:runId "FB-001" ; bp:processTimeH ?t ; bp:co2Pct ?v }
                          ORDER BY DESC(?v) LIMIT 1""")[0]
    close("IO-05", peak[0])
    times = [t for (t,) in rows(store, """SELECT ?t WHERE { ?o a bp:OnlineRecord ; bp:ofRun/bp:runId "FB-010" ; bp:processTimeH ?t } ORDER BY ?t""")]
    close("L4-02", max(b - a for a, b in zip(times, times[1:])))


def test_units(store):
    base = {rid: v for rid, v in rows(store, """SELECT ?id (MAX(?b) AS ?v) WHERE { VALUES ?id { "FB-002" "FB-006" }
                                                ?o a bp:OnlineRecord ; bp:ofRun/bp:runId ?id ; bp:baseML ?b } GROUP BY ?id""")}
    assert [base["FB-002"], base["FB-006"]] == pytest.approx(gold("IO-10"), abs=0.01)
    close("CA-02", one(store, """SELECT ?v WHERE { ?m bp:ofSample/bp:sampleId "FB018-S06" ; bp:sourceType "hplc_raw" ;
                                 bp:analyte/bp:analyteName "Xylitol" ; bp:value ?v }"""))
    feed = one(store, """SELECT ((MAX(?g) - MIN(?g)) / ?rho AS ?ml) WHERE { ?r bp:runId "FB-006" ; bp:feedType ?code .
                         ?o a bp:OnlineRecord ; bp:ofRun ?r ; bp:feedBalanceG ?g . ?f bp:feedCode ?code ; bp:densityGPerML ?rho }
                         GROUP BY ?rho""")
    close("CA-03", feed)
    pump, start, end = one(store, """SELECT (MAX(?p) AS ?pump) ?start ?end WHERE { ?r bp:runId "FB-005" .
                                     ?o a bp:OnlineRecord ; bp:ofRun ?r ; bp:feedPumpML ?p .
                                     ?a bp:ofRun ?r ; bp:valueLabel "Feed bottle start (g)" ; bp:value ?start .
                                     ?b bp:ofRun ?r ; bp:valueLabel "Feed bottle end (g)" ; bp:value ?end } GROUP BY ?start ?end""")
    close("CA-04", ((start - end) / 1.15 / pump - 1) * 100)


def test_lineage_and_clones(store):
    chain = rows(store, """SELECT ?child ?parent ?clone WHERE { ?s bp:strainId "YX-006" ; bp:parentStrain* ?c .
                           ?c bp:strainId ?child ; bp:parentStrain/bp:strainId ?parent . OPTIONAL { ?c bp:builtOnClone/bp:cloneId ?clone } }
                           ORDER BY DESC(?child)""")
    assert chain == [("YX-006", "YX-003", "YX-003-c3"), ("YX-003", "YX-002", "YX-002-c1"), ("YX-002", "YX-001", None)]
    used = {c for (c,) in rows(store, 'SELECT DISTINCT ?c WHERE { VALUES ?id { "FB-005" "FB-008" } ?r bp:runId ?id ; bp:usesClone/bp:cloneId ?c }')}
    assert used == {"YX-002-c3"}                                              # IO-08: not the clone YX-003 was built on
    assert one(store, 'SELECT ?d WHERE { ?s bp:strainId "YX-005" ; bp:cloneSelectedOn ?d }') == "2025-04-29"
    od = one(store, """SELECT ?id ?v WHERE { ?r a bp:FlaskRun ; bp:runId ?id ; bp:usesClone/bp:cloneId "YX-006-c3" ; bp:c5FractionPct ?c5 .
                       ?s bp:ofRun ?r ; bp:processTimeH ?t . FILTER(?c5 = 70 && ?t = 72) ?m bp:ofSample ?s ; bp:analyte/bp:analyteName "OD600" ; bp:value ?v }""")
    assert od[0] == "SF-033"
    close("IO-12", od[1])
