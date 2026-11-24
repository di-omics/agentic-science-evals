"""Run AI agents on the tasks: stage one clean folder per (model, task), then grade and audit.

    python -m science_evals.runs stage <run_root> <model> [<model> ...]
    python -m science_evals.runs collect <run_root> <transcripts_dir> <out_dir>

stage    copies each task's input files (APFS clones where possible) into
         <run_root>/<model>/<task>/ and writes TASK.md there: the eval's task text, unchanged.
collect  reads every agent transcript (JSONL) in <transcripts_dir>, finds the run folder named
         in its prompt, takes the answer (the last <EVAL_ANSWER> block of the final message,
         else answer.json in the folder), grades it with both graders, and audits every tool
         call: web use, and any command or path that does not stay inside the run folder.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from .data import DATA_ROOT, ROOT
from .grading import latch_check, latch_version, own_check

ANSWER_RE = re.compile(r"<EVAL_ANSWER>\s*(\{.*?\})\s*</EVAL_ANSWER>", re.S)
FOLDER_RE = re.compile(r"Your working folder: (\S+)")
PATH_RE = re.compile(r"(?<![\w.])(~/[^\s'\"`;|&<>)]*|/(?:Users|private|tmp|var|etc|opt|Library|Volumes|home|root)\b[^\s'\"`;|&<>)]*)")
ALLOWED_PREFIXES = ("/dev/null",)


def _evals():
    manifest = json.loads((ROOT / "evals" / "manifest.json").read_text())
    for paper in manifest["papers"]:
        for entry in paper["evals"]:
            yield paper["slug"], json.loads((ROOT / entry["path"]).read_text())


def stage(run_root: Path, models: list[str]) -> dict:
    index = {}
    for slug, e in _evals():
        for model in models:
            ws = run_root / model / e["id"]
            ws.mkdir(parents=True, exist_ok=False)
            for f in e["data_files"]:
                src = DATA_ROOT / slug / f["name"]
                if not src.exists():
                    raise FileNotFoundError(src)
                dst = ws / f["name"]
                if subprocess.run(["cp", "-c", str(src), str(dst)], capture_output=True).returncode != 0:
                    shutil.copy2(src, dst)
            (ws / "TASK.md").write_text(e["task"] + "\n")
            index[str(ws)] = {"model": model, "paper": slug, "task": e["id"]}
    (run_root / "index.json").write_text(json.dumps(index, indent=1))
    return index


def _text_blocks(msg) -> list[str]:
    content = msg.get("content")
    if isinstance(content, str):
        return [content]
    return [b.get("text", "") for b in content or [] if isinstance(b, dict) and b.get("type") == "text"]


def _read_transcript(path: Path) -> dict:
    prompt, final, tools, models, stamps = None, "", [], set(), []
    usage = {"input_tokens": 0, "output_tokens": 0, "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0}
    with open(path) as fh:
        for line in fh:
            try:
                d = json.loads(line)
            except Exception:
                continue
            if not isinstance(d, dict):
                continue
            if d.get("timestamp"):
                stamps.append(d["timestamp"])
            msg = d.get("message") if isinstance(d.get("message"), dict) else {}
            if d.get("type") == "user" and prompt is None and isinstance(msg.get("content"), str):
                prompt = msg["content"]
            if d.get("type") == "assistant":
                if msg.get("model"):
                    models.add(msg["model"])
                for k in usage:
                    usage[k] += int((msg.get("usage") or {}).get(k) or 0)
                texts = _text_blocks(msg)
                if any(t.strip() for t in texts):
                    final = "\n".join(texts)
                for b in msg.get("content") or []:
                    if isinstance(b, dict) and b.get("type") == "tool_use":
                        tools.append({"name": b.get("name"), "input": b.get("input")})
    dur = None
    if len(stamps) >= 2:
        t = [datetime.fromisoformat(s.replace("Z", "+00:00")) for s in (min(stamps), max(stamps))]
        dur = round((t[1] - t[0]).total_seconds())
    return {"prompt": prompt or "", "final": final, "tools": tools, "models": sorted(models), "usage": usage, "seconds": dur}


def audit(tools: list[dict], ws: str, python: str, harness_dir: str = "") -> dict:
    """Web use, and paths outside the run folder. Paths in the harness's own task folder (where
    background shell output lands) are listed apart: they are not inputs, keys or other runs'
    folders, but an agent that lists that folder can see other runs' file names."""
    allowed = (ws.rstrip("/"), str(Path(python).parent.parent)) + ALLOWED_PREFIXES
    web, outside, no_folder = 0, set(), []
    for t in tools:
        name = t["name"] or ""
        if name in ("WebFetch", "WebSearch") or (name.startswith("mcp__") and not name.startswith("mcp__ccd_turn__")):
            web += 1
        text = json.dumps(t["input"])
        for p in PATH_RE.findall(text.replace("\\/", "/")):
            p = p.rstrip(".,:")
            if not p.startswith(allowed):
                outside.add(p)
        if t["name"] == "Bash":
            cmd = (t["input"] or {}).get("command", "")
            if ws not in cmd and not cmd.strip().startswith(str(Path(python))):
                no_folder.append(cmd[:200])
    harness = sorted(p for p in outside if harness_dir and p.startswith(harness_dir))
    other = sorted(p for p in outside if p not in harness)
    return {"web_calls": web, "paths_outside": other, "harness_paths": harness, "commands_without_folder": no_folder}


def collect(run_root: Path, transcripts: Path, out_dir: Path, python: str) -> list[dict]:
    index = json.loads((run_root / "index.json").read_text())
    evals = {e["id"]: e for _, e in _evals()}
    rows = []
    for path in sorted(transcripts.glob("a*.output")):
        try:
            tr = _read_transcript(path)
        except (FileNotFoundError, UnicodeDecodeError):
            continue
        m = FOLDER_RE.search(tr["prompt"])
        if not m or m.group(1).rstrip("/") not in index:
            continue
        ws = m.group(1).rstrip("/")
        meta = index[ws]
        answer, source = None, None
        found = ANSWER_RE.findall(tr["final"])
        if found:
            try:
                answer, source = json.loads(found[-1]), "final message"
            except Exception:
                answer = None
        if answer is None and (Path(ws) / "answer.json").exists():
            try:
                answer, source = json.loads((Path(ws) / "answer.json").read_text()), "answer.json"
            except Exception:
                answer = None
        e = evals[meta["task"]]
        own = own_check(answer or {}, e["grader"])
        latch = latch_check(answer or {}, e["grader"])
        rows.append({
            **meta,
            "model_ids": tr["models"],
            "passed": own["passed"],
            "latch_passed": None if latch is None else latch["passed"],
            "answer_source": source,
            "answer": answer,
            "fields": own["fields"],
            "tool_calls": len(tr["tools"]),
            "seconds": tr["seconds"],
            "usage": tr["usage"],
            "audit": audit(tr["tools"], ws, python, str(transcripts)),
        })
    out_dir.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (r["paper"], r["task"], r["model"]))
    summary = {"latch_eval_tools": latch_version(), "runs": rows}
    text = json.dumps(summary, indent=1, ensure_ascii=True)
    # Local paths say nothing about the result: replace them with neutral names.
    for local, neutral in ((str(run_root), "<run>"), (str(transcripts), "<harness>"),
                           (str(Path(python).parent.parent), "<agent-env>"), (str(Path.home()), "~")):
        text = text.replace(local.replace("/", "\\/"), neutral).replace(local, neutral)
    text = re.sub(r"/private/tmp/claude-\d+/[^\s\"'\\]*", "<tmp>", text)
    text = re.sub(r"/Users/[A-Za-z0-9._-]+", "~", text)
    rows = json.loads(text)["runs"]
    (out_dir / "runs.json").write_text(text + "\n")
    (out_dir / "RESULTS.md").write_text(report(rows))
    return rows


def _miss(row: dict) -> str:
    """The fields that failed, and by how much."""
    bad = []
    for name, d in row["fields"].items():
        if d.get("ok"):
            continue
        if "jaccard" in d:
            bad.append(f"{name} (Jaccard {d['jaccard']})")
        elif isinstance(d.get("value"), (int, float)):
            lo, hi = d["window"]
            bad.append(f"{name} {d['value']:.4g} (window {lo:.4g} to {hi:.4g})")
        else:
            bad.append(f"{name} missing")
    return "; ".join(bad)


def report(rows: list[dict]) -> str:
    models = sorted({r["model"] for r in rows}, key=["haiku", "sonnet", "opus"].index)
    tasks = sorted({(r["paper"], r["task"]) for r in rows})
    cell = {(r["task"], r["model"]): r for r in rows}
    lines = ["# Agent runs", "",
             "One fresh agent per model and task, on a folder holding only TASK.md and the task's input files. "
             "Graded by this repository's own check and by latch-eval-tools; an audit of every tool call looks "
             "for web use and for paths outside the run folder.", "",
             "| Task | " + " | ".join(models) + " |", "|---|" + "---|" * len(models)]
    for _, task in tasks:
        marks = []
        for m in models:
            r = cell.get((task, m))
            marks.append("not run" if r is None else ("pass" if r["passed"] else "fail"))
        lines.append(f"| `{task}` | " + " | ".join(marks) + " |")
    totals = [f"{m}: {sum(1 for r in rows if r['model'] == m and r['passed'])} of "
              f"{sum(1 for r in rows if r['model'] == m)}" for m in models]
    lines += ["", "Passed: " + ", ".join(totals) + ".", "", "## Each run", "",
              "| Task | Model | Verdict | Graders agree | Tool calls | Seconds | Audit | What missed |",
              "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        a = r["audit"]
        flags = []
        if a["web_calls"]:
            flags.append(f"{a['web_calls']} web calls")
        if a["paths_outside"]:
            flags.append(f"{len(a['paths_outside'])} paths outside")
        if a.get("harness_paths"):
            flags.append("touched the harness's output folder")
        agree = "n/a" if r["latch_passed"] is None else ("yes" if r["latch_passed"] == r["passed"] else "NO")
        lines.append(f"| `{r['task']}` | {r['model']} | {'pass' if r['passed'] else 'fail'} | {agree} | "
                     f"{r['tool_calls']} | {r['seconds']} | {', '.join(flags) or 'clean'} | {_miss(r) or ''} |")
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    if argv[:1] == ["stage"]:
        idx = stage(Path(argv[1]), argv[2:])
        print(f"staged {len(idx)} run folders")
        return 0
    if argv[:1] == ["collect"]:
        rows = collect(Path(argv[1]), Path(argv[2]), Path(argv[3]), argv[4])
        print(f"collected {len(rows)} runs, {sum(r['passed'] for r in rows)} passed")
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
