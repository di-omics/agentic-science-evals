"""Checks on what is committed: recorded controls, eval files, house rules. Needs no data."""
import json
import subprocess
from pathlib import Path

import pytest

from science_evals.evalfile import validate
from science_evals.grading import own_check
from science_evals.papers import PAPERS

ROOT = Path(__file__).resolve().parents[1]
TEXT = {".py", ".md", ".json", ".toml", ".yml", ".yaml", ".txt", ".csv", ".tsv", ".svg", ".html", ".cfg", ".sha256", ""}


def _files() -> list[Path]:
    try:
        out = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                             cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
        return [ROOT / p for p in out if (ROOT / p).is_file()]
    except Exception:
        skip = {"data", ".venv", ".git"}
        return [p for p in ROOT.rglob("*") if p.is_file() and not skip & set(p.relative_to(ROOT).parts)]


@pytest.mark.parametrize("slug", PAPERS)
def test_recorded_controls_hold(slug):
    path = ROOT / "results" / slug / "controls.json"
    if not path.exists():
        pytest.skip(f"{slug} not built yet")
    s = json.loads(path.read_text())
    assert s["controls"] > 0 and s["as_expected"] == s["controls"]
    kinds = {}
    for r in s["rows"]:
        kinds.setdefault(r["task"], set()).add(r["kind"])
        e = json.loads((ROOT / "evals" / slug / r["task"] / "eval.json").read_text())
        assert own_check(r["answer"], e["grader"])["passed"] == r["own_pass"], (r["task"], r["control"])
    for task, ks in kinds.items():
        assert ks == {"reference", "variant", "trap"}, (task, ks)


def test_eval_files_are_valid():
    for f in sorted((ROOT / "evals").glob("*/*/eval.json")):
        assert validate(json.loads(f.read_text())) == [], f


def test_tracked_text_is_ascii():
    bad = []
    for p in _files():
        if p.suffix.lower() in TEXT or p.name in ("LICENSE", ".gitignore"):
            data = p.read_bytes()
            if any(b > 126 or (b < 32 and b not in (9, 10, 13)) for b in data):
                bad.append(str(p.relative_to(ROOT)))
    assert bad == [], bad


def test_no_forbidden_names():
    assert [str(p) for p in _files() if "cuttag" in str(p.relative_to(ROOT)).lower()] == []
