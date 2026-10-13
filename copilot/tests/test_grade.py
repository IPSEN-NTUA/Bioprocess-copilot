import pytest

from copilot import grade as G


def q(**kw):
    return {"tolerance": None, "value": None, **kw}


def test_numbers_skip_identifiers_and_dates():
    assert G.numbers("FB-006 S08 reached 112.63 g/L on 2025-03-21 (SEQ-20250321, clone c3, R3, SOP v1)") == [112.63]
    assert G.numbers("3,3 g/L at 70 % C5") == [3.3, 70.0]
    assert G.numbers("1:10") == [1.0, 10.0]


@pytest.mark.parametrize("final,ok", [("112.63 g/L", True), ("about 112.6 g/L", False), ("FB-112 gave 41.63", False)])
def test_numeric(final, ok):
    assert G.grade_numeric(q(value=112.63, tolerance=0.02), final)["correct"] is ok


def test_numeric_zero_tolerance():
    assert G.grade_numeric(q(value=6, tolerance=0), "6 sequences")["correct"]
    assert not G.grade_numeric(q(value=6, tolerance=0), "5 sequences in May, 26 in total")["correct"]


def test_exact():
    assert G.grade_exact(q(value="YX-005"), "Strain YX-005.")["correct"]
    assert not G.grade_exact(q(value="YX-005"), "YX-0055")["correct"]
    assert G.grade_exact(q(value="Mark Doherty"), "mark doherty (MD)")["correct"]
    assert G.grade_exact(q(value=10), "1:10")["correct"]
    assert G.grade_exact(q(value="2025-06-01"), "effective from 2025-06-01")["correct"]
    assert not G.grade_exact(q(value="c3"), "clone c1")["correct"]


def test_ranking():
    gold = ["FB-013", "FB-014", "FB-012"]
    assert G.grade_ranking(q(value=gold), "FB-013 (67.6) > FB14 > FB-012")["correct"]
    assert not G.grade_ranking(q(value=gold), "FB-014 > FB-013 > FB-012")["correct"]
