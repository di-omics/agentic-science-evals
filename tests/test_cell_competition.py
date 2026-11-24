"""Cell competition (Nat Metab 2021) package: metadata, pinned files, reference
functions on toy data, and, when data/cell_competition is present, a rerun of
the recorded controls."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from science_evals.papers import cell_competition as cc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cell_competition"
AUTHOR_ROLE = "D.H. performed embryo dissections, treatments and cell dissociation prior to scRNA-seq experiments."


def test_paper_metadata():
    for k in ("short", "citation", "doi", "author_role", "data", "code", "license"):
        assert cc.PAPER.get(k), k
    assert cc.PAPER["author_role"] == AUTHOR_ROLE
    assert cc.PAPER["doi"] in cc.PAPER["citation"]
    text = json.dumps(cc.PAPER)
    assert all(ord(ch) < 127 for ch in text)


def test_data_files_are_pinned():
    from science_evals.papers.cell_competition import prepare

    names = [f["name"] for f in cc.DATA_FILES]
    assert len(names) == len(set(names))
    for f in cc.DATA_FILES:
        assert set(f) == {"name", "source", "sha256", "bytes"}, f
        assert len(f["sha256"]) == 64 and f["bytes"] > 0, f["name"]
        if f["source"].startswith("derived:"):
            assert hasattr(prepare, f["source"].split(".", 1)[1]), f["source"]
        else:
            assert f["source"].startswith("https://"), f["source"]
            assert "/main/" not in f["source"] and "/master/" not in f["source"], "pin a commit"


def test_qc_outlier_rules_on_toy_table():
    pytest.importorskip("scipy")
    from science_evals.papers.cell_competition import reference as R

    rng = np.random.default_rng(0)
    n = 200
    qc = pd.DataFrame({
        "treatment": ["CI"] * 100 + ["DMSO"] * 100,
        "log10_total_reads": rng.normal(6, 0.1, n),
        "n_genes_above_10rpm": rng.normal(6000, 200, n),
        "frac_reads_ercc": rng.normal(0.02, 0.002, n),
        "frac_reads_mito": rng.normal(0.04, 0.004, n),
    })
    # a continuum of stressed cells making up just over half of one treatment, as
    # the intermediate and CI-only cells do in the real data
    qc.loc[:54, "frac_reads_ercc"] = rng.uniform(0.025, 0.08, 55)
    pooled = R.qc_outliers(qc, 3)
    within = R.qc_outliers_by(qc, 5, "treatment")
    assert pooled[:55].mean() > 0.8
    assert within[:55].sum() == 0
    assert R.qc_answer(qc, pooled)["n_ci_cells_kept"] <= 60
    assert R.qc_answer(qc, within) == {"n_cells_kept": 200, "n_ci_cells_kept": 100, "n_dmso_cells_kept": 100}


def test_heteroplasmy_rules_on_toy_counts():
    pytest.importorskip("scipy")
    from science_evals.papers.cell_competition import reference as R

    pos = np.arange(1, 401)
    counts = np.zeros((4, len(pos), 4), dtype=np.int64)
    counts[:, :, 0] = 100  # all A, deep coverage
    cols = [int(np.flatnonzero(pos == p)[0]) for p in R.RNR1_POSITIONS]
    counts[1, cols, 1] = 5  # cell 1: about 5% C at the six positions
    counts[2, cols[:3], :] = 0  # cell 2: three positions uncovered
    counts[2, cols[3:], 1] = 5
    depth, h = R.heteroplasmy(counts)
    cluster = np.array(["winner", "loser", "loser", "intermediate"])
    full = R.rnr1_fractions(depth, h, pos, cluster, min_positions=None)
    assert full["n_cells_used"] == 3 and full["frac_heteroplasmic_loser"] == 1.0
    zero = R.rnr1_fractions(depth, h, pos, cluster, min_positions=None, mode="zero_fill")
    assert zero["n_cells_used"] == 4 and zero["frac_heteroplasmic_loser"] == 1.0
    avail = R.rnr1_fractions(depth, h, pos, cluster, min_positions=None, mode="available")
    assert avail["n_cells_used"] == 4
    nan = R.rnr1_fractions(depth, h, pos, cluster, min_positions=None, mode="nan_false")
    assert nan["n_cells_used"] == 4 and nan["frac_heteroplasmic_loser"] == 0.5


def test_diffusion_pseudotime_orders_a_curve():
    pytest.importorskip("scipy")
    from science_evals.papers.cell_competition.trajectory import diffusion_map, spearman

    t = np.linspace(0, 1, 120)
    x = np.c_[np.cos(2.5 * t), np.sin(2.5 * t), 0.01 * np.sin(40 * t)]
    dm = diffusion_map(x, n_eigs=10)
    tips = dm.tips(0)
    root = tips[int(np.argmin(t[tips]))]
    assert spearman(dm.dpt(root), t) > 0.95


def _have_data() -> bool:
    needed = [f["name"] for f in cc.DATA_FILES if not f["source"].startswith("derived:")]
    return all((DATA / n).exists() for n in needed)


@pytest.mark.skipif(not _have_data(), reason="data/cell_competition not downloaded")
def test_recorded_controls_reproduce_from_data():
    pytest.importorskip("rdata")
    pytest.importorskip("scipy")
    recorded = ROOT / "results" / "cell_competition" / "controls.json"
    if not recorded.exists():
        pytest.skip("cell_competition not built yet")
    rows = json.loads(recorded.read_text())["rows"]
    built = {b.task.id: b for b in cc.build(DATA)}
    for r in rows:
        ctrl = next(c for c in built[r["task"]].controls if c.name == r["control"])
        assert json.loads(json.dumps(ctrl.answer)) == r["answer"], (r["task"], r["control"])
