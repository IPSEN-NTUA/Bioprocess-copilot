import pytest

from copilot import arms, build_kg, rag


@pytest.fixture(scope="module")
def index():
    return rag.load_index(embed=False)          # lexical only: no model server needed for the tests


def test_every_corpus_file_is_chunked(index):
    files = {p.relative_to(build_kg.CORPUS).as_posix() for p in build_kg.CORPUS.rglob("*") if p.is_file()}
    chunked = {f for f in files if any(c.source.startswith(f) for c in index.chunks)}
    assert chunked == files
    assert max(len(c.text) for c in index.chunks) < 4000


def test_answer_key_is_refused():
    key = build_kg.ROOT / "benchmark/pilot/answer_key"
    with pytest.raises(ValueError):
        rag.chunk_corpus(key)
    with pytest.raises(ValueError):
        build_kg.build(key, None)


def test_chunks_keep_the_raw_formats(index):
    text = "\n".join(c.text for c in index.chunks if "MD_fermentation_master" in c.source)
    assert "FB02 | S01 | 03/02/25" in text and "<0.05" in text      # nothing is normalised on the RAG side


def test_lexical_search_and_per_file_cap(index):
    hits = index.search("FB018-S06 Xylitol dilution", k=8, mode="lexical")
    assert any("FB018-S06" in c.text for c in hits[:3])
    files = [c.source.split(" ")[0] for c in hits]
    assert max(files.count(f) for f in files) <= 3


def test_schema_prompt_covers_the_ontology():
    text = arms.schema_text()
    for name in ("bp:Measurement", "bp:OnlineRecord", "bp:processTimeH", "bp:isCurrent", "bp:supersededBy", "bp:lotUsed"):
        assert name in text


def test_sparql_tool(tmp_path):
    build_kg.rebuild(store_path=tmp_path / "kg")
    sp = arms.Sparql(tmp_path / "kg")
    try:
        assert sp.run('SELECT ?id WHERE { ?r bp:runId "FB-013" ; bp:usesStrain/bp:strainId ?id }').split() == ["id", "YX-005"]
        assert sp.run("SELECT ?x WHERE { ?x a bp:Nothing }").endswith("(no rows)")
        assert sp.run("DELETE WHERE { ?s ?p ?o }").startswith("SPARQL error")
        assert sp.run("SELECT ?x WHERE { ?x a").startswith("SPARQL error")
        assert "more rows not shown" in sp.run("SELECT ?s WHERE { ?s a bp:Sample }")
    finally:
        sp.close()
