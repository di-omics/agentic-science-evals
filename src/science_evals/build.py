"""Build a paper's eval.json files from its data, run its controls, write the results.

    python -m science_evals.build ave_migration      one paper
    python -m science_evals.build --index            the top-level manifest and summary

Outputs per paper (all generated, never edited by hand):
    evals/<paper>/<task>/eval.json, evals/<paper>/manifest.json
    results/<paper>/controls.json, results/<paper>/CONTROLS.md
    results/<paper>/paper_anchors.json, results/<paper>/PAPER_ANCHORS.md (if the paper has anchors)
    docs/<paper>/TASKS.md
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

from .controls import coverage_problems, run_controls
from .data import DATA_ROOT, ROOT
from .evalfile import make_eval, validate
from .grading import latch_version
from .papers import PAPERS


def _dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=True) + "\n")


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, list):
        return f"{len(v)} labels"
    return str(v)


def control_line(row: dict) -> str:
    verdict = "passes" if row["own_pass"] else "fails"
    shown = ", ".join(f"{k}={_fmt(v)}" for k, v in list(row["answer"].items())[:4])
    return f"{row['kind']} {row['control']}: {row['what']} -> {shown} ({verdict})"


def build_paper(slug: str) -> dict:
    mod = importlib.import_module(f"science_evals.papers.{slug}")
    data_dir = DATA_ROOT / slug
    built = mod.build(data_dir)
    problems, all_rows, manifest, task_docs = [], [], [], []
    for b in built:
        problems += b.task.check() + coverage_problems(b.task.id, b.controls)
        probe = make_eval(b.task, b.reference, b.data_files, mod.PAPER, [])
        rows = run_controls(probe, b.controls)
        e = make_eval(b.task, b.reference, b.data_files, mod.PAPER, [control_line(r) for r in rows])
        problems += [f"{b.task.id}: {p}" for p in validate(e)]
        rel = Path("evals") / slug / b.task.id / "eval.json"
        _dump(ROOT / rel, e)
        manifest.append({"id": b.task.id, "path": str(rel)})
        all_rows += rows
        task_docs.append(f"## `{b.task.id}`: {b.task.title}\n\nFrom the paper: {b.task.paper_anchor}\n\n"
                         f"```\n{e['task']}\n```\n\nGrader:\n\n```json\n{json.dumps(e['grader'], indent=1)}\n```\n\n"
                         f"Notes:\n\n```\n{e['notes']}\n```\n")
    _dump(ROOT / "evals" / slug / "manifest.json", {"paper": mod.PAPER["citation"], "evals": manifest})
    ok = sum(r["as_expected"] for r in all_rows)
    summary = {
        "paper": mod.PAPER["citation"],
        "tasks": len(built),
        "controls": len(all_rows),
        "as_expected": ok,
        "latch_eval_tools": latch_version(),
        "rows": all_rows,
    }
    _dump(ROOT / "results" / slug / "controls.json", summary)
    lines = [f"# Controls: {mod.PAPER['short']}", "",
             f"{ok} of {len(all_rows)} scripted analyses behave as expected "
             f"(graded with this repository's own check"
             + (f" and latch-eval-tools {latch_version()}" if latch_version() else "") + ").", "",
             "| Task | Control | Kind | Should | Own grader | Latch grader | What it does |", "|---|---|---|---|---|---|---|"]
    for r in all_rows:
        lp = "n/a" if r["latch_pass"] is None else ("pass" if r["latch_pass"] else "fail")
        lines.append(f"| `{r['task']}` | {r['control']} | {r['kind']} | {'pass' if r['expected_pass'] else 'fail'} | "
                     f"{'pass' if r['own_pass'] else 'fail'} | {lp} | {r['what']} |")
    (ROOT / "results" / slug / "CONTROLS.md").write_text("\n".join(lines) + "\n")
    if hasattr(mod, "anchors"):
        anchors = mod.anchors(data_dir)
        _dump(ROOT / "results" / slug / "paper_anchors.json", anchors)
        alines = [f"# Paper anchors: {mod.PAPER['short']}", "", "| The paper says | Where | Recomputed here |", "|---|---|---|"]
        alines += [f"| {a['paper_says']} | {a['where']} | {a['recomputed']} |" for a in anchors]
        (ROOT / "results" / slug / "PAPER_ANCHORS.md").write_text("\n".join(alines) + "\n")
    docs = [f"# Tasks: {mod.PAPER['short']}", "", mod.PAPER["citation"], ""] + task_docs
    (ROOT / "docs" / slug).mkdir(parents=True, exist_ok=True)
    (ROOT / "docs" / slug / "TASKS.md").write_text("\n".join(docs))
    summary["problems"] = problems
    return summary


def build_index() -> dict:
    papers = []
    for slug in PAPERS:
        m = ROOT / "evals" / slug / "manifest.json"
        if m.exists():
            papers.append({"paper": slug, **json.loads(m.read_text())})
    _dump(ROOT / "evals" / "manifest.json", {"papers": papers})
    return {"papers": len(papers), "tasks": sum(len(p["evals"]) for p in papers)}


def main(argv: list[str]) -> int:
    if argv == ["--index"]:
        print(json.dumps(build_index()))
        return 0
    bad = 0
    for slug in argv or PAPERS:
        s = build_paper(slug)
        print(f"{slug}: {s['tasks']} tasks, {s['as_expected']}/{s['controls']} controls as expected")
        for p in s["problems"]:
            print(f"  PROBLEM {p}")
        bad += len(s["problems"]) + (s["controls"] - s["as_expected"])
    build_index()
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
