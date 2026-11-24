"""AVE migration paper: the package contract, the decision rule, and (when the data are present)
the numbers the tasks rest on. Tests that need data/ave_migration skip without it, as in CI."""
import json
import re
from pathlib import Path

import pytest

from science_evals.papers import ave_migration as paper
from science_evals.papers.ave_migration import inputs as io
from science_evals.papers.ave_migration.choices import AVE04_CHOICES, decide_dkk1

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "ave_migration"
PUBLIC_NAMES = [e["name"] for e in io.PUBLIC]
HAVE_DATA = all((DATA / n).exists() for n in PUBLIC_NAMES)
needs_data = pytest.mark.skipif(not HAVE_DATA, reason="data/ave_migration not downloaded")


def test_paper_contract():
    for key in ("short", "citation", "doi", "author_role", "data", "code", "license"):
        assert paper.PAPER[key]
    assert "D.H. and S.R. did sample collections for the transcriptomics and proteomics, respectively" \
        in paper.PAPER["author_role"]
    assert io.SCAVE_COMMIT in paper.PAPER["code"]


def test_data_files_are_pinned():
    names = [f["name"] for f in paper.DATA_FILES]
    assert len(names) == len(set(names))
    for f in paper.DATA_FILES:
        assert f["source"].startswith("https://") or f["source"].startswith("derived:"), f
        assert re.fullmatch(r"[0-9a-f]{64}", f["sha256"]), f
        assert f["bytes"] > 0, f
        assert "cuttag" not in f["name"].lower()
    for f in paper.DATA_FILES:
        if f["source"].startswith("https://raw.githubusercontent.com/"):
            assert io.SCAVE_COMMIT in f["source"]


def test_ave04_options_and_rule():
    lengths = [len(c) for c in AVE04_CHOICES]
    assert len(AVE04_CHOICES) == 6
    assert all(c.isascii() for c in AVE04_CHOICES)
    key = decide_dkk1({"dkk1_auc_e55": 0.52, "dkk1_auc_e625": 0.78}, detected_e55=0.8)
    assert key == 2 and lengths[key - 1] not in (max(lengths), min(lengths))
    assert decide_dkk1({"dkk1_auc_e55": 0.80, "dkk1_auc_e625": 0.85}, 0.8) == 3
    assert decide_dkk1({"dkk1_auc_e55": 0.20, "dkk1_auc_e625": 0.25}, 0.8) == 5
    assert decide_dkk1({"dkk1_auc_e55": 0.52, "dkk1_auc_e625": 0.78}, 0.1) == 6


def test_generated_evals_if_built():
    manifest = ROOT / "evals" / "ave_migration" / "manifest.json"
    if not manifest.exists():
        pytest.skip("ave_migration not built yet")
    ids = [e["id"] for e in json.loads(manifest.read_text())["evals"]]
    assert ids == ["ave01_plate_qc", "ave02_pseudotime_root", "ave03_ephrin_bootstrap_signs", "ave04_dkk1_timing"]
    for i in ids:
        e = json.loads((ROOT / "evals" / "ave_migration" / i / "eval.json").read_text())
        names = {f["name"] for f in e["data_files"]}
        for name in names:
            assert f"`{name}`" in e["task"], (i, name)


@needs_data
def test_derived_files_match_their_pins():
    for name in io.DERIVED_ORDER:
        if name == io.ENA_COUNTS and not (DATA / name).exists():
            pytest.skip("ENA run report not derived yet (needs network)")
        got = io.describe(io.derive(DATA, name))
        assert (got["sha256"], got["bytes"]) == io.DERIVED_PINS[name], name


@needs_data
def test_qc_reference_reproduces_the_authors_wells():
    pytest.importorskip("scipy")
    from science_evals.papers.ave_migration import reference as ref

    answers = ref.ave01(DATA)
    sdrf = io.read_sdrf(DATA)
    authors = sorted(sdrf.index[sdrf["authors_pass"] & (sdrf["stage"] == "5.5 dpc")])
    assert answers["reference"]["wells_kept"] == authors
    assert answers["droplet_defaults"]["n_wells_kept"] > 350


@needs_data
def test_ephrin_numbers_match_the_paper():
    pytest.importorskip("scipy")
    from science_evals.papers.ave_migration import reference as ref

    got = ref._ave03_answer(str(DATA))
    assert (got["n_genes_kept_e55"], got["n_genes_kept_e625"]) == (14, 13)
    assert (got["n_significant_e55"], got["n_significant_e625"]) == (6, 8)
    assert ref._ave03_answer(str(DATA), reflect=False)["n_significant_e55"] == 0


@needs_data
def test_dkk1_is_enriched_only_after_migration():
    pytest.importorskip("anndata")
    from science_evals.papers.ave_migration import reference as ref

    aucs = ref.stage_aucs(str(DATA))
    assert abs(aucs["dkk1_auc_e55"] - 0.5) < 0.1 and aucs["dkk1_auc_e625"] > 0.7
    assert aucs["cer1_auc_e55"] > 0.75 and aucs["cer1_auc_e625"] > 0.75
