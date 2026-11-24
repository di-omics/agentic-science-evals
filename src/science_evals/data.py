"""Fetch public input files and check them against their recorded sha256.

Every file a task hands to an agent is listed in its paper's DATA_FILES as
{"name", "source", "sha256", "bytes"}. "source" is a public https URL, or
"derived:<function>" for a file this repository makes from public files (the
paper module then says how, and the sha256 still pins the result).
"""
from __future__ import annotations

import hashlib
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "data"
USER_AGENT = "agentic-science-evals/0.1 (public data fetch)"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(entry: dict, dest_dir: Path) -> Path:
    """Download one https entry if it is missing, then verify its hash."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / entry["name"]
    if not path.exists():
        source = entry["source"]
        if not source.startswith("https://"):
            raise ValueError(f"{entry['name']}: only https sources are fetched, got {source}")
        tmp = path.with_suffix(path.suffix + ".part")
        req = urllib.request.Request(source, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
            shutil.copyfileobj(resp, out, 1 << 20)
        tmp.rename(path)
    want = entry.get("sha256")
    if want:
        got = sha256_file(path)
        if got != want:
            raise ValueError(f"{entry['name']}: sha256 {got} does not match {want}")
    return path


def describe(path: Path) -> dict:
    return {"sha256": sha256_file(path), "bytes": path.stat().st_size}
