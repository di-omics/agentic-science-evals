"""Recorded agent runs: every verdict still follows from the answer and the current grader."""
import json
from pathlib import Path

import pytest

from science_evals.grading import own_check

ROOT = Path(__file__).resolve().parents[1]
RUNS = sorted((ROOT / "results" / "agent_runs").glob("*/runs.json"))


@pytest.mark.parametrize("path", RUNS, ids=[p.parent.name for p in RUNS])
def test_recorded_verdicts_regrade(path):
    runs = json.loads(path.read_text())["runs"]
    assert runs
    for r in runs:
        e = json.loads((ROOT / "evals" / r["paper"] / r["task"] / "eval.json").read_text())
        assert own_check(r["answer"] or {}, e["grader"])["passed"] == r["passed"], (r["model"], r["task"])
        assert r["latch_passed"] in (None, r["passed"]), (r["model"], r["task"])


@pytest.mark.parametrize("path", RUNS, ids=[p.parent.name for p in RUNS])
def test_no_local_paths_recorded(path):
    text = path.read_text()
    assert "/Users/" not in text and "/private/tmp" not in text
