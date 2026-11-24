"""APOE microglia (Nat Commun 2025) package: metadata, pinned files, the DESeq2
engine and the reference functions on toy data, and, when data/apoe_microglia is
present, a recomputation of one paper anchor."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from science_evals.papers import apoe_microglia as am

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "apoe_microglia"
AUTHOR_ROLE = "D.H., L.W., S.R., G.L.F., and R.M. performed experiments."


def test_paper_metadata():
    for k in ("short", "citation", "doi", "author_role", "data", "code", "license"):
        assert am.PAPER.get(k), k
    assert am.PAPER["author_role"] == AUTHOR_ROLE
    assert am.PAPER["doi"] in am.PAPER["citation"]
    assert all(ord(ch) < 127 for ch in json.dumps(am.PAPER))


def test_data_files_are_pinned():
    from science_evals.papers.apoe_microglia import inputs

    names = [f["name"] for f in am.DATA_FILES]
    assert len(names) == len(set(names))
    for f in am.DATA_FILES:
        assert set(f) == {"name", "source", "sha256", "bytes"}, f
        assert len(f["sha256"]) == 64 and f["bytes"] > 0, f["name"]
        assert ("cut" + "tag") not in f["name"].lower()
        if f["source"].startswith("derived:"):
            assert hasattr(inputs, f["source"].split(".", 1)[1]), f["source"]
        else:
            assert f["source"].startswith("https://"), f["source"]
            assert "/main/" not in f["source"] and "/master/" not in f["source"], "pin a commit"


def _toy_counts(seed=0, n_genes=400, per_group=4, lfc_genes=40):
    rng = np.random.default_rng(seed)
    base = rng.lognormal(4, 1, n_genes)
    lfc = np.zeros(n_genes)
    lfc[:lfc_genes] = rng.choice([-2.0, 2.0], lfc_genes)
    sf = rng.uniform(0.6, 1.6, 2 * per_group)
    groups = ["A"] * per_group + ["B"] * per_group
    mu = np.stack([base * sf[j] * (2 ** lfc if g == "B" else 1) for j, g in enumerate(groups)], 1)
    disp = 0.05
    counts = rng.negative_binomial(1 / disp, 1 / (1 + disp * mu))
    return counts, groups, lfc


def test_deseq_engine_recovers_toy_effects():
    pytest.importorskip("scipy")
    from science_evals.papers.apoe_microglia import deseq as D

    counts, groups, lfc = _toy_counts()
    fit = D.fit(counts, groups)
    r = fit.contrast("B", "A")
    q = D.bh(r["pvalue"])
    called = q < 0.05
    assert called[:40].mean() > 0.8  # most true effects found
    assert called[40:].mean() < 0.05  # few false calls
    assert np.all(np.sign(r["log2FoldChange"][:40][called[:40]]) == np.sign(lfc[:40][called[:40]]))
    flipped = fit.contrast("A", "B")
    assert np.allclose(flipped["log2FoldChange"], -r["log2FoldChange"], equal_nan=True)
    # normalisation factors equal to the size factors give the same answer
    nf = np.broadcast_to(D.size_factors(counts)[None, :], counts.shape)
    same = D.fit(counts, groups, norm_factors=nf).contrast("B", "A")
    assert np.allclose(same["pvalue"], r["pvalue"], equal_nan=True)


def test_bh_and_filter_by_expr():
    from science_evals.papers.apoe_microglia import deseq as D

    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5])
    q = D.bh(p)
    assert np.isnan(q[3])
    assert np.allclose(q[[0, 1, 2, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.5])  # R: p.adjust(p, "BH")
    counts = np.array([[0, 0, 0, 0], [50, 60, 0, 0], [100, 90, 80, 70], [1, 0, 0, 2]])
    keep = D.filter_by_expr(counts, ["a", "a", "b", "b"])
    assert keep.tolist() == [False, True, True, False]


def _toy_rna():
    from science_evals.papers.apoe_microglia import reference as R

    samples = [f"{g}-{i}" for g in ("E2", "E3", "E4", "KO") for i in range(1, 4)]
    rng = np.random.default_rng(1)
    counts = pd.DataFrame(rng.integers(50, 150, (5, len(samples))), index=[348, 1, 2, 3, 4], columns=samples)
    counts.loc[348] = [1000] * 9 + [900, 60, 50]  # KO-1 still expresses APOE
    lengths = pd.DataFrame(1000.0, index=counts.index, columns=samples)
    genotype = {"E2": "APOE2", "E3": "APOE3", "E4": "APOE4", "KO": "APOE-KO"}
    sheet = pd.DataFrame({"sample": samples, "genotype": [genotype[s.split("-")[0]] for s in samples]})
    return R.Rna(counts, lengths, sheet, pd.Series({348: "APOE"}))


def test_failed_knockout_rule_on_toy_samples():
    from science_evals.papers.apoe_microglia import reference as R

    rna = _toy_rna()
    assert R.failed_knockouts(rna) == ["KO-1"]
    assert rna.columns(["E3", "KO"], exclude=["KO-1"]) == ["E3-1", "E3-2", "E3-3", "KO-2", "KO-3"]


def test_vdr_direction_on_toy_sheet():
    pytest.importorskip("scipy")
    from science_evals.papers.apoe_microglia import reference as R

    sheet = pd.DataFrame({
        "hgnc_symbol": [f"G{i}" for i in range(20)],
        "log2FoldChange": [-1.0] * 5 + [1.0] * 5 + [0.1] * 10,
        "fdr": [0.01] * 10 + [0.9] * 10,
    })
    targets = {"G0", "G1", "G2", "G7", "X"}
    up_in_apoe2 = R.vdr_test(sheet, targets)  # APOE4 over APOE2: negative is higher in APOE2
    assert up_in_apoe2["vdr_targets_higher_in_apoe2"] == ["G0", "G1", "G2"]
    assert up_in_apoe2["n_higher_in_apoe2"] == 5
    reversed_ = R.vdr_test(sheet, targets, higher_in_apoe2_sign=1)
    assert reversed_["vdr_targets_higher_in_apoe2"] == ["G7"]


def test_concordance_matched_differs_from_pooled_on_toy_tables():
    pytest.importorskip("scipy")
    from science_evals.papers.apoe_microglia import reference as R

    rng = np.random.default_rng(2)
    genes = np.arange(1, 41)
    degs, rows = {}, []
    for k, comp in enumerate(R.SHEETS):
        lfc = rng.normal(0, 2, genes.size)
        degs[comp] = pd.DataFrame({"entrezgene_id": genes, "hgnc_symbol": [f"G{g}" for g in genes],
                                   "log2FoldChange": lfc, "fdr": 0.01})
        acc = lfc + rng.normal(0, 0.3, genes.size)  # accessibility follows expression within a comparison
        rows += [{"comparison": comp, "geneId": g, "SYMBOL": f"G{g}", "log2FoldChange": a} for g, a in zip(genes, acc)]
    peaks = pd.DataFrame(rows)
    matched = R.concordance(degs, peaks)
    pooled = R.concordance(degs, peaks, pooled=True)
    assert min(matched[k] for k in matched if k.startswith("r_")) > 0.9
    assert max(pooled[k] for k in pooled if k.startswith("r_")) < 0.9
    assert matched["n_e2_vs_e3"] == pooled["n_e2_vs_e3"] == 40


@pytest.mark.skipif(not (DATA / "GSE271384_APOE_microglia_ATACseq_counts.csv.gz").exists(),
                    reason="data/apoe_microglia not downloaded (CI)")
def test_filter_by_expr_reproduces_the_paper_peak_count():
    from science_evals.papers.apoe_microglia import deseq as D

    counts = pd.read_csv(DATA / "GSE271384_APOE_microglia_ATACseq_counts.csv.gz", index_col=0).iloc[:, 5:]
    cols = [c for c in counts.columns if c not in ("KO_1", "KO_2")]
    keep = D.filter_by_expr(counts[cols].to_numpy(dtype=float), [c.split("_")[0] for c in cols])
    assert int(keep.sum()) == 167951
