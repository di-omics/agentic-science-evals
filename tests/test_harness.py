"""The shared harness on a toy task: grader shape, both graders, controls, eval.json checks."""
from science_evals.controls import Control, coverage_problems, run_controls
from science_evals.evalfile import make_eval, validate
from science_evals.grading import grader_config, own_check
from science_evals.spec import Field, Task, canary

REF = {"n_cells": 100, "frac": 0.25, "markers": ["A", "B", "C", "D"]}
FILES = [{"name": "x.csv", "source": "https://example.org/x.csv", "sha256": "0" * 64, "bytes": 1}]
PAPER = {"citation": "Toy paper", "short": "Toy"}


def _task() -> Task:
    return Task(
        id="toy_task",
        title="Toy task",
        paper_anchor="Fig. 1",
        prompt="Count the cells and name the markers.",
        fields=[
            Field("n_cells", "int", "cells kept", tolerance={"value": 5}),
            Field("frac", "float", "a fraction", tolerance={"value": 0.02}),
            Field("markers", "labels", "marker genes", pass_threshold=0.5),
            Field("decision", "choice"),
        ],
        data=["x.csv"],
        failure_modes=["qc_control_artifact.context_inappropriate_qc_threshold"],
        category="toy",
        kit="scrna",
        choice_question="Which reading is right?",
        choices=["A", "B"],
        correct_choice=2,
        intent="i",
        trap="t",
        tolerance_rationale="r",
    )


def test_grader_shape_and_own_check():
    t = _task()
    assert t.check() == []
    g = grader_config(t, REF)
    assert g["type"] == "all_of" and len(g["config"]["children"]) == 3
    good = {"n_cells": 103, "frac": 0.26, "markers": ["A", "B", "C"], "decision": 2}
    assert own_check(good, g)["passed"]
    assert not own_check({**good, "decision": 1}, g)["passed"]
    assert not own_check({**good, "n_cells": 110}, g)["passed"]
    assert not own_check({**good, "markers": ["X", "Y"]}, g)["passed"]
    assert not own_check({**good, "frac": True}, g)["passed"]
    assert not own_check({k: v for k, v in good.items() if k != "frac"}, g)["passed"]


def test_controls_eval_file_and_latch_agreement():
    t = _task()
    controls = [
        Control("ref", "reference", {"n_cells": 100, "frac": 0.25, "markers": ["A", "B", "C", "D"], "decision": 2}, "the analysis"),
        Control("var", "variant", {"n_cells": 96, "frac": 0.24, "markers": ["A", "B", "C"], "decision": 2}, "lighter filter"),
        Control("trap", "trap", {"n_cells": 150, "frac": 0.4, "markers": ["A", "B", "C", "D"], "decision": 1}, "no filter"),
    ]
    assert coverage_problems(t.id, controls) == []
    e = make_eval(t, REF, FILES, PAPER, ["x"])
    assert validate(e) == []
    assert e["canary"] == canary("toy_task")
    rows = run_controls(e, controls)  # also compares with latch-eval-tools when installed
    assert all(r["as_expected"] for r in rows), rows


def test_validate_catches_problems():
    e = make_eval(_task(), REF, FILES, PAPER, ["x"])
    e["data_files"][0]["source"] = "ftp://example.org/x.csv"
    assert any("source" in p for p in validate(e))
    e = make_eval(_task(), REF, FILES, PAPER, ["x"])
    e["notes"] += " " + chr(8212)
    assert any("non-ASCII" in p for p in validate(e))
