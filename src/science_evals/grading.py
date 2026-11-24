"""Grader configs in LatchBio's public eval format, and two ways to apply them.

own_check     a few lines of arithmetic written here, with no dependencies.
latch_check   the published latch-eval-tools package, when it is installed.

Every scripted control is graded both ways and the two verdicts must agree.
Supported grader types: numeric_tolerance, label_set_jaccard, and all_of over them.
A multiple-choice decision is a numeric_tolerance field with tolerance 0.
"""
from __future__ import annotations

from .spec import Field, Task


def _round(value, f: Field):
    if f.kind in ("int", "choice") or f.digits == 0:
        return int(round(value))
    return round(float(value), f.digits)


def _tolerance(f: Field) -> dict:
    tol = {"type": "absolute"}
    tol.update(f.tolerance or {"value": 0})
    return tol


def grader_config(task: Task, reference: dict) -> dict:
    """Build the eval.json grader block from a task and its reference answer."""
    children = []
    numbers = [f for f in task.fields if f.kind in ("int", "float")]
    if numbers:
        children.append({
            "type": "numeric_tolerance",
            "config": {
                "ground_truth": {f.name: _round(reference[f.name], f) for f in numbers},
                "tolerances": {f.name: _tolerance(f) for f in numbers},
            },
        })
    for f in task.fields:
        if f.kind == "labels":
            labels = sorted({str(x) for x in reference[f.name]})
            children.append({
                "type": "label_set_jaccard",
                "config": {
                    "ground_truth_labels": labels,
                    "scoring": {"pass_threshold": f.pass_threshold},
                    "answer_field": f.name,
                },
            })
    choices = [f for f in task.fields if f.kind == "choice"]
    if choices:
        children.append({
            "type": "numeric_tolerance",
            "config": {
                "ground_truth": {f.name: task.correct_choice for f in choices},
                "tolerances": {f.name: {"type": "absolute", "value": 0} for f in choices},
            },
        })
    if len(children) == 1:
        return children[0]
    return {"type": "all_of", "config": {"children": children}}


def _leaves(grader: dict) -> list[dict]:
    if grader["type"] == "all_of":
        out = []
        for child in grader["config"]["children"]:
            out.extend(_leaves(child))
        return out
    return [grader]


def windows(grader: dict) -> dict[str, tuple[float, float]]:
    """Accepted closed interval per numeric answer field."""
    out = {}
    for leaf in _leaves(grader):
        if leaf["type"] != "numeric_tolerance":
            continue
        cfg = leaf["config"]
        for name, truth in cfg["ground_truth"].items():
            tol = cfg["tolerances"][name]
            value = tol.get("value", 0)
            out[name] = (truth - tol.get("lower", value), truth + tol.get("upper", value))
    return out


def own_check(answer: dict, grader: dict) -> dict:
    """Independent verdict: every number inside its window, every label set close enough."""
    detail = {}
    for name, (lo, hi) in windows(grader).items():
        value = answer.get(name)
        ok = isinstance(value, (int, float)) and not isinstance(value, bool) and lo <= value <= hi
        gap = 0.0
        if isinstance(value, (int, float)) and not isinstance(value, bool) and not ok:
            gap = lo - value if value < lo else value - hi
        detail[name] = {"value": value, "window": [lo, hi], "ok": ok, "gap": gap}
    for leaf in _leaves(grader):
        if leaf["type"] != "label_set_jaccard":
            continue
        cfg = leaf["config"]
        name = cfg.get("answer_field", "cell_types_predicted")
        truth = {str(x) for x in cfg["ground_truth_labels"]}
        given = answer.get(name)
        if not isinstance(given, list):
            detail[name] = {"value": given, "jaccard": 0.0, "ok": False}
            continue
        pred = {str(x) for x in given}
        union = truth | pred
        jac = len(truth & pred) / len(union) if union else 0.0
        threshold = cfg.get("scoring", {}).get("pass_threshold", 0.9)
        detail[name] = {"jaccard": round(jac, 4), "threshold": threshold, "ok": jac >= threshold}
    known = {"numeric_tolerance", "label_set_jaccard"}
    unknown = [leaf["type"] for leaf in _leaves(grader) if leaf["type"] not in known]
    if unknown:
        raise ValueError(f"own_check does not support grader types {unknown}")
    return {"passed": all(d["ok"] for d in detail.values()), "fields": detail}


def latch_available() -> bool:
    try:
        import latch_eval_tools.graders  # noqa: F401
    except Exception:
        return False
    return True


def latch_version() -> str | None:
    try:
        from importlib.metadata import version

        return version("latch-eval-tools")
    except Exception:
        return None


def latch_check(answer: dict, grader: dict) -> dict | None:
    """Verdict from LatchBio's own grader. None when the package is absent."""
    if not latch_available():
        return None
    from latch_eval_tools.graders import get_grader

    result = get_grader(grader["type"]).evaluate(answer, grader["config"])
    return {"passed": bool(result.passed), "score": float(result.score)}
