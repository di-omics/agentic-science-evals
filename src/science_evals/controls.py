"""Scripted controls: analyses that must pass or must fail a task's grader.

reference  the analysis the task is built around. Must pass.
variant    a defensible alternative (another convention, a lighter filter). Must pass.
trap       a plausible wrong analysis, named in the task notes. Must fail.

Each control is an answer produced by running code on the task's data, never a
number typed in to sit just inside or outside a window.
"""
from __future__ import annotations

from dataclasses import dataclass

from .grading import latch_check, own_check

KINDS = ("reference", "variant", "trap")


@dataclass
class Control:
    name: str
    kind: str  # reference | variant | trap
    answer: dict
    what: str  # one line: what this analysis does differently


def run_controls(eval_json: dict, controls: list[Control]) -> list[dict]:
    rows = []
    for c in controls:
        if c.kind not in KINDS:
            raise ValueError(f"{c.name}: kind must be one of {KINDS}")
        own = own_check(c.answer, eval_json["grader"])
        latch = latch_check(c.answer, eval_json["grader"])
        expected = c.kind != "trap"
        agree = latch is None or latch["passed"] == own["passed"]
        rows.append({
            "task": eval_json["id"],
            "control": c.name,
            "kind": c.kind,
            "what": c.what,
            "expected_pass": expected,
            "own_pass": own["passed"],
            "latch_pass": None if latch is None else latch["passed"],
            "graders_agree": agree,
            "as_expected": own["passed"] == expected and agree,
            "answer": c.answer,
            "fields": own["fields"],
        })
    return rows


def coverage_problems(task_id: str, controls: list[Control]) -> list[str]:
    kinds = {c.kind for c in controls}
    missing = [k for k in KINDS if k not in kinds]
    return [f"{task_id}: no {k} control" for k in missing]
