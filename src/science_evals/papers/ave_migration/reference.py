"""The analyses behind each task: the reference, defensible variants and plausible wrong ones.

Every function here runs on the task's input files and returns answers in the shape of the
task's answer block. Nothing is typed in by hand.

ave01  which wells of a Smart-seq2 plate hold one good cell (QC)
ave02  which end of the embryonic VE trajectory is the AVE (pseudotime root and direction)
ave03  which ephrin/Eph genes drive the pathway's principal components (bootstrap with axis signs)
ave04  is Dkk1 an AVE signal before or after migration (stage-resolved enrichment)
"""
from __future__ import annotations

import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as ss

from . import inputs as io
from .choices import AVE04_CHOICES, ENRICHED, NULL_BAND, decide_dkk1  # noqa: F401

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def _mad_low(x, n=3.0):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return x < med - n * mad


def _mad_high(x, n=3.0):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return x > med + n * mad


# ============================================================================
# ave01: plate QC
# ============================================================================

PAPER_QC = {"min_genes_10rpm": 3000, "min_log10_reads": 4.0, "min_frac_mapped": 0.5,
            "max_frac_mito": 0.1, "max_frac_ercc": 0.3}


@lru_cache(maxsize=2)
def qc_table(data_dir: str) -> pd.DataFrame:
    """Per-well metrics from the counts, the well sheet and the gene table."""
    d = Path(data_dir)
    counts = io.read_counts(d, "e55")
    wells = pd.read_csv(io.derive(d, io.WELLS_E55), sep="\t").set_index("well_id")
    genes = pd.read_csv(io.derive(d, io.GENES), sep="\t").set_index("gene_id")
    ercc = counts.index.str.startswith("ERCC-")
    g = counts.loc[~ercc]
    mito = genes.loc[g.index, "symbol"].str.startswith("mt-").to_numpy()
    gs = g.sum(0).to_numpy()
    es = counts.loc[ercc].sum(0).to_numpy()
    ms = g.loc[mito].sum(0).to_numpy()
    w = wells.loc[counts.columns]
    rp = w["read_pairs"].to_numpy().astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rpm = g.to_numpy() / gs * 1e6
        t = pd.DataFrame(index=counts.columns)
        t["well_type"] = w["well_type"].to_numpy()
        t["read_pairs"] = rp
        t["counted"] = gs + es
        t["n_genes"] = (g.to_numpy() > 0).sum(0)
        t["n_genes_10rpm"] = np.nan_to_num(rpm > 10).sum(0)
        t["frac_mapped"] = (gs + es) / rp
        t["frac_ercc"] = es / rp
        t["frac_mito"] = ms / rp
        t["frac_ercc_counted"] = np.where(gs + es > 0, es / (gs + es), 1.0)
        t["frac_mito_counted"] = np.where(gs + es > 0, ms / (gs + es), 0.0)
    return t


def _paper_thresholds(t: pd.DataFrame, min_genes: int = 3000) -> pd.Series:
    return ((t.n_genes_10rpm > min_genes) & (np.log10(t.read_pairs) > PAPER_QC["min_log10_reads"])
            & (t.frac_mapped > PAPER_QC["min_frac_mapped"]) & (t.frac_mito < PAPER_QC["max_frac_mito"])
            & (t.frac_ercc < PAPER_QC["max_frac_ercc"]))


def _qc_answer(keep: pd.Series) -> dict:
    kept = sorted(keep.index[keep.to_numpy()])
    return {"n_wells_kept": len(kept), "wells_kept": kept}


def ave01(data_dir: Path) -> dict[str, dict]:
    t = qc_table(str(data_dir))
    single = t.well_type == "single cell"
    out = {
        "reference": _qc_answer(_paper_thresholds(t) & single),
        "thresholds_all_wells": _qc_answer(_paper_thresholds(t)),
        "counted_reads_only": _qc_answer(
            (t.n_genes_10rpm > 3000) & (np.log10(t.counted.clip(lower=1)) > 4)
            & (t.frac_ercc_counted < 0.3) & (t.frac_mito_counted < 0.1) & single),
        "gene_floor_2000": _qc_answer(_paper_thresholds(t, min_genes=2000) & single),
    }
    s = t[single]
    out_mad = (_mad_low(np.log10(s.read_pairs)) | _mad_low(np.log10(s.n_genes.clip(lower=1)))
               | _mad_high(s.frac_ercc) | _mad_high(s.frac_mito))
    keep = pd.Series(False, index=t.index)
    keep.loc[s.index[~out_mad.to_numpy()]] = True
    out["mad_outliers"] = _qc_answer(keep)
    out["droplet_defaults"] = _qc_answer((t.n_genes >= 200) & (t.frac_mito_counted < 0.2))
    out["genes_and_mito_only"] = _qc_answer((t.n_genes > 1000) & (t.frac_mito_counted < 0.1) & single)
    return out


def ave01_facts(data_dir: Path) -> dict:
    t = qc_table(str(data_dir))
    single = t.well_type == "single cell"
    ref = _paper_thresholds(t) & single
    removed = t[single & ~ref]
    drop_kept = (t.n_genes >= 200) & (t.frac_mito_counted < 0.2)
    controls_pass = t[(t.well_type != "single cell") & _paper_thresholds(t)]
    return {
        "n_wells": len(t),
        "n_single": int(single.sum()),
        "n_no_cell": int((t.well_type == "no cell").sum()),
        "n_multi": int((t.well_type == "multiple cells").sum()),
        "n_ref": int(ref.sum()),
        "n_removed_single": len(removed),
        "removed_median_ercc": float(removed.frac_ercc.median()),
        "removed_median_mapped": float(removed.frac_mapped.median()),
        "removed_median_genes": float(removed.n_genes.median()),
        "kept_median_ercc": float(t[ref].frac_ercc.median()),
        "kept_median_mapped": float(t[ref].frac_mapped.median()),
        "kept_median_genes": float(t[ref].n_genes.median()),
        "removed_passing_droplet": int((drop_kept & single & ~ref).sum()),
        "n_droplet": int(drop_kept.sum()),
        "droplet_annotated": int((drop_kept & ~single).sum()),
        "fail_ercc": int((removed.frac_ercc >= 0.3).sum()),
        "fail_mapped": int((removed.frac_mapped <= 0.5).sum()),
        "fail_genes": int((removed.n_genes_10rpm <= 3000).sum()),
        "controls_passing": len(controls_pass),
        "controls_passing_types": controls_pass.well_type.value_counts().to_dict(),
        "max_mito_all": float(t.frac_mito.max()),
        "strict_floor_no_spikein": int(((t.n_genes_10rpm > 3000) & (np.log10(t.read_pairs) > 4)
                                        & (t.frac_mito < 0.1) & single).sum()),
    }


# ============================================================================
# ave02: pseudotime root and direction
# ============================================================================

AVE_MARKERS = ("Cer1", "Hhex", "Lefty1")
PANEL_SIZE = 5


@lru_cache(maxsize=2)
def _embryonic_ve(data_dir: str):
    import scanpy as sc

    d = Path(data_dir)
    a = io.read_norm(d, "e55")
    cells = io.derive(d, io.EMB_VE_E55).read_text().split()
    a = a[cells].copy()
    sc.pp.log1p(a)
    return a


def _spearman_rep(pcs: np.ndarray) -> np.ndarray:
    """Rows scaled so that Euclidean distance is 2 * sqrt((1 - rho) / 2), rho the Spearman
    correlation between two cells' principal-component vectors (the lab's distance, times 2;
    the kNN graph's weights do not change under a constant factor)."""
    r = np.apply_along_axis(ss.rankdata, 1, pcs)
    r = r - r.mean(1, keepdims=True)
    return r / np.linalg.norm(r, axis=1, keepdims=True)


@lru_cache(maxsize=4)
def _graph(data_dir: str, metric: str = "spearman", scaled: bool = False):
    import scanpy as sc

    a = _embryonic_ve(data_dir)
    hv = a.copy()
    sc.pp.highly_variable_genes(hv, n_top_genes=3000)
    hv = hv[:, hv.var["highly_variable"]].copy()
    if scaled:
        sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, svd_solver="arpack", random_state=0)
    pcs = np.asarray(hv.obsm["X_pca"][:, :10], dtype=np.float64)
    hv.obsm["X_rep"] = _spearman_rep(pcs) if metric == "spearman" else pcs
    sc.pp.neighbors(hv, n_neighbors=15, use_rep="X_rep", random_state=0)
    sc.tl.diffmap(hv, n_comps=20, random_state=0)
    return hv


def _ave_score(a) -> np.ndarray:
    x = np.column_stack([np.asarray(a[:, g].X).ravel() for g in AVE_MARKERS])
    z = (x - x.mean(0)) / x.std(0)
    return z.mean(1)


def _dc1_ends(hv, a) -> tuple[int, int]:
    """(AVE end, emVE end) of the first non-trivial diffusion component."""
    dc1 = np.asarray(hv.obsm["X_diffmap"][:, 1])
    s = np.sign(ss.spearmanr(dc1, _ave_score(a))[0]) or 1.0
    return int(np.argmax(s * dc1)), int(np.argmin(s * dc1))


def _pseudotime(hv, root: int) -> np.ndarray:
    import scanpy as sc

    hv.uns["iroot"] = int(root)
    sc.tl.dpt(hv, n_dcs=10)
    return hv.obs["dpt_pseudotime"].to_numpy(dtype=float).copy()


def _ns_basis(x: np.ndarray, df: int) -> np.ndarray:
    """Natural cubic spline basis with an intercept, df - 1 interior knots at quantiles."""
    knots = np.quantile(x, np.linspace(0, 1, df + 1))
    k = len(knots)

    def dk(j):
        return (np.clip(x - knots[j], 0, None) ** 3 - np.clip(x - knots[-1], 0, None) ** 3) / (knots[-1] - knots[j])

    cols = [np.ones_like(x), x] + [dk(j) - dk(k - 2) for j in range(k - 2)]
    return np.column_stack(cols)


def _ftest(p: np.ndarray, y: np.ndarray, df: int) -> np.ndarray:
    x1 = _ns_basis(p, df) if df > 1 else np.column_stack([np.ones_like(p), p])
    beta, *_ = np.linalg.lstsq(x1, y, rcond=None)
    rss1 = ((y - x1 @ beta) ** 2).sum(0)
    rss0 = ((y - y.mean(0)) ** 2).sum(0)
    n, q = x1.shape
    f = ((rss0 - rss1) / (q - 1)) / (rss1 / (n - q))
    return ss.f.sf(f, q - 1, n - q)


def _spearman_all(p: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rp = ss.rankdata(p)
    ry = np.apply_along_axis(ss.rankdata, 0, y)
    rp = (rp - rp.mean()) / rp.std()
    sd = ry.std(0)
    ry = (ry - ry.mean(0)) / np.where(sd == 0, 1, sd)
    rho = (rp[:, None] * ry).mean(0)
    n = len(p)
    t = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
    return rho, 2 * ss.t.sf(np.abs(t), n - 2)


@lru_cache(maxsize=1)
def _tested_genes(data_dir: str):
    a = _embryonic_ve(data_dir)
    x = np.asarray(a.X, dtype=np.float64)
    keep = (x > 0).sum(0) >= 10
    return a.var_names[keep].to_numpy(), x[:, keep]


def _dpt_result(data_dir: str, root: int, metric: str = "spearman", test: str = "spline4",
                scaled: bool = False) -> dict:
    a = _embryonic_ve(data_dir)
    hv = _graph(data_dir, metric, scaled)
    p = _pseudotime(hv, root)
    genes, y = _tested_genes(data_dir)
    rho, p_rank = _spearman_all(p, y)
    if test == "spearman":
        pv = p_rank
    else:
        pv = _ftest(p, y, {"spline4": 4, "spline3": 3, "linear": 1}[test])
    sig = _bh(pv) < 0.01
    high = sig & (rho < 0)
    cer1 = np.asarray(a[:, "Cer1"].X).ravel()
    return {"pseudotime": p, "genes": genes, "rho": rho, "sig": sig, "high": high,
            "rho_cer1": float(ss.spearmanr(p, cer1)[0]), "n_sig": int(sig.sum()),
            "frac_high": float(high.sum() / max(sig.sum(), 1))}


@lru_cache(maxsize=1)
def ave02_panel(data_dir: str) -> tuple[list[str], list[str]]:
    """Panel genes: the strongest high-in-AVE and low-in-AVE genes of the reference, the
    markers named in the prompt left out. Returns (panel in alphabetical order, keyed subset)."""
    ref = _dpt_result(data_dir, _dc1_ends(_graph(data_dir), _embryonic_ve(data_dir))[0])
    genes, rho, sig = ref["genes"], ref["rho"], ref["sig"]
    ok = sig & ~np.isin(genes, AVE_MARKERS)
    neg = np.flatnonzero(ok & (rho < 0))
    pos = np.flatnonzero(ok & (rho > 0))
    high = list(genes[neg[np.argsort(rho[neg])][:PANEL_SIZE]])
    low = list(genes[pos[np.argsort(-rho[pos])][:PANEL_SIZE]])
    return sorted(high + low), sorted(high)


def _ave02_answer(data_dir: str, res: dict) -> dict:
    panel, _ = ave02_panel(data_dir)
    high = set(res["genes"][res["high"]])
    return {"rho_pseudotime_cer1": res["rho_cer1"], "frac_high_in_ave": res["frac_high"],
            "high_in_ave_panel": sorted(g for g in panel if g in high)}


def ave02(data_dir: Path) -> dict[str, dict]:
    d = str(data_dir)
    a = _embryonic_ve(d)
    hv = _graph(d)
    ave_end, emve_end = _dc1_ends(hv, a)
    marker_root = int(np.argmax(_ave_score(a)))
    hv_e = _graph(d, "euclidean")
    ave_end_e, _ = _dc1_ends(hv_e, a)
    ave_end_s, _ = _dc1_ends(_graph(d, "spearman", True), a)
    dc2 = np.asarray(hv.obsm["X_diffmap"][:, 2])
    score = _ave_score(a)
    dc2_ends = (int(np.argmin(dc2)), int(np.argmax(dc2)))
    dc2_root = max(dc2_ends, key=lambda i: score[i])
    authors = pd.read_csv(io.path(Path(d), io.AUTHORS_DC_E55), index_col=0).loc[a.obs_names]
    authors_root = int(np.argmin(authors["DC2"].to_numpy()))
    return {
        "reference": _ave02_answer(d, _dpt_result(d, ave_end)),
        "root_max_marker_score": _ave02_answer(d, _dpt_result(d, marker_root)),
        "euclidean_knn": _ave02_answer(d, _dpt_result(d, ave_end_e, metric="euclidean")),
        "scaled_pca": _ave02_answer(d, _dpt_result(d, ave_end_s, scaled=True)),
        "root_dc2_extreme_in_ave": _ave02_answer(d, _dpt_result(d, dc2_root)),
        "authors_root_cell": _ave02_answer(d, _dpt_result(d, authors_root)),
        "linear_trend_test": _ave02_answer(d, _dpt_result(d, ave_end, test="linear")),
        "spearman_test": _ave02_answer(d, _dpt_result(d, ave_end, test="spearman")),
        "spline_3df": _ave02_answer(d, _dpt_result(d, ave_end, test="spline3")),
        "root_at_emve_end": _ave02_answer(d, _dpt_result(d, emve_end)),
        "root_first_cell": _ave02_answer(d, _dpt_result(d, 0)),
    }


def ave02_facts(data_dir: Path) -> dict:
    d = str(data_dir)
    a = _embryonic_ve(d)
    hv = _graph(d)
    ave_end, emve_end = _dc1_ends(hv, a)
    ref = _dpt_result(d, ave_end)
    rev = _dpt_result(d, emve_end)
    lab = pd.read_csv(io.derive(Path(d), io.LABELS)).set_index("cell_id")
    is_ave = (lab.loc[a.obs_names, "cluster"] == "AVE").to_numpy()
    first = _dpt_result(d, 0)
    dm = np.asarray(hv.obsm["X_diffmap"][:, 1])
    score = _ave_score(a)
    return {
        "n_cells": a.n_obs,
        "n_tested": len(ref["genes"]),
        "n_sig": ref["n_sig"],
        "n_high": int(ref["high"].sum()),
        "frac_high": ref["frac_high"],
        "rev_frac_high": rev["frac_high"],
        "rho_cer1": ref["rho_cer1"],
        "rev_rho_cer1": rev["rho_cer1"],
        "first_cell": a.obs_names[0],
        "first_cell_cluster": lab.loc[a.obs_names[0], "cluster"],
        "first_cell_rho": first["rho_cer1"],
        "root_cluster": lab.loc[a.obs_names[ave_end], "cluster"],
        "emve_end_cluster": lab.loc[a.obs_names[emve_end], "cluster"],
        "n_ave": int(is_ave.sum()),
        "mean_pt_ave": float(ref["pseudotime"][is_ave].mean()),
        "mean_pt_emve": float(ref["pseudotime"][~is_ave].mean()),
        "dc1_marker_rho": float(ss.spearmanr(dm, score)[0]),
        "authors_root": str(pd.read_csv(io.path(Path(d), io.AUTHORS_DC_E55), index_col=0)["DC2"].idxmin()),
        "ave_root_rhos": _ave_root_rhos(d, is_ave),
    }


def _ave_root_rhos(data_dir: str, is_ave: np.ndarray) -> list[float]:
    """Cer1 correlation of the pseudotime rooted at each AVE cell in turn."""
    a = _embryonic_ve(data_dir)
    hv = _graph(data_dir)
    cer1 = np.asarray(a[:, "Cer1"].X).ravel()
    return [float(ss.spearmanr(_pseudotime(hv, int(r)), cer1)[0]) for r in np.flatnonzero(is_ave)]


# ============================================================================
# ave03: ephrin/Eph PCA loadings, bootstrap with axis signs
# ============================================================================

N_BOOT = 1000


@lru_cache(maxsize=1)
def _ephrin_table(data_dir: str) -> pd.DataFrame:
    return pd.read_csv(io.derive(Path(data_dir), io.EPHRIN))


def ephrin_matrix(data_dir: str, stage: str) -> tuple[list[str], np.ndarray]:
    t = _ephrin_table(data_dir)
    t = t[t["stage"] == stage]
    genes = [c for c in t.columns if c not in ("cell_id", "stage", "cell_type")]
    x = np.log1p(t[genes].to_numpy(dtype=np.float64))
    keep = x.mean(0) > 0.5
    return [g for g, k in zip(genes, keep) if k], x[:, keep]


def _loadings(x: np.ndarray, k: int = 3, solver: str = "numpy") -> np.ndarray:
    if solver == "sklearn":
        from sklearn.decomposition import PCA

        return PCA(n_components=k, svd_solver="full").fit(x).components_
    xc = x - x.mean(0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    return vt[:k]


def bootstrap_pvalues(x: np.ndarray, n_boot: int = N_BOOT, seed: int = 0, reflect: bool = True,
                      reorder: bool = False, solver: str = "numpy", k: int = 3) -> np.ndarray:
    """Fraction of bootstrap samples in which each loading has the other sign (k x genes)."""
    rng = np.random.default_rng(seed)
    l0 = _loadings(x, k, solver)
    flips = np.zeros_like(l0)
    n = len(x)
    for _ in range(n_boot):
        lb = _loadings(x[rng.integers(0, n, n)], k, solver)
        if reorder:
            c = np.abs(l0 @ lb.T)
            order, used = [], set()
            for i in range(k):
                j = max((j for j in range(k) if j not in used), key=lambda j: c[i, j])
                used.add(j)
                order.append(j)
            lb = lb[order]
        if reflect:
            s = np.sign(np.sum(l0 * lb, axis=1))
            s[s == 0] = 1
            lb = lb * s[:, None]
        flips += np.sign(lb) != np.sign(l0)
    return flips / n_boot


def _ave03_answer(data_dir: str, **kw) -> dict:
    out = {}
    for key, stage in (("e55", "5.5 dpc"), ("e625", "6.25 dpc")):
        genes, x = ephrin_matrix(data_dir, stage)
        p = bootstrap_pvalues(x, **kw)
        sig = sorted(g for j, g in enumerate(genes) if (p[:, j] < 0.01).any())
        out[f"n_genes_kept_{key}"] = len(genes)
        out[f"n_significant_{key}"] = len(sig)
        out[f"significant_genes_{key}"] = sig
    return out


def ave03(data_dir: Path) -> dict[str, dict]:
    d = str(data_dir)
    return {
        "reference": _ave03_answer(d),
        "axis_reordering_too": _ave03_answer(d, reorder=True),
        "other_seed": _ave03_answer(d, seed=1),
        "sklearn_pca": _ave03_answer(d, solver="sklearn"),
        "two_thousand_samples": _ave03_answer(d, n_boot=2000),
        "signs_not_aligned": _ave03_answer(d, reflect=False),
        "sklearn_signs_not_aligned": _ave03_answer(d, reflect=False, solver="sklearn"),
    }


def ave03_facts(data_dir: Path) -> dict:
    d = str(data_dir)
    out = {}
    for key, stage in (("e55", "5.5 dpc"), ("e625", "6.25 dpc")):
        genes, x = ephrin_matrix(d, stage)
        p_raw = bootstrap_pvalues(x, reflect=False)
        p_ref = bootstrap_pvalues(x)
        out[f"median_p_raw_{key}"] = float(np.median(p_raw))
        out[f"min_p_raw_{key}"] = float(p_raw.min())
        out[f"median_p_ref_{key}"] = float(np.median(p_ref))
        out[f"n_cells_{key}"] = int(len(x))
        ev = np.linalg.svd(x - x.mean(0), compute_uv=False) ** 2
        out[f"var_top3_{key}"] = [float(v) for v in ev[:3] / ev.sum()]
    return out


# ============================================================================
# ave04: Dkk1 timing
# ============================================================================

@lru_cache(maxsize=2)
def _stage_cells(data_dir: str, stage: str):
    import scanpy as sc

    d = Path(data_dir)
    a = io.read_norm(d, stage)
    lab = pd.read_csv(io.derive(d, io.LABELS)).set_index("cell_id")
    a = a[a.obs_names.isin(lab.index)].copy()
    a.obs["cluster"] = lab.loc[a.obs_names, "cluster"].to_numpy()
    a = a[a.obs["cluster"].isin(["AVE", "emVE"])].copy()
    sc.pp.log1p(a)
    return a


def _auc(x: np.ndarray, is_ave: np.ndarray) -> float:
    u = ss.mannwhitneyu(x[is_ave], x[~is_ave], alternative="two-sided").statistic
    return float(u / (is_ave.sum() * (~is_ave).sum()))


def _gene_by_stage(data_dir: str, gene: str, log2: bool = False):
    for stage in ("e55", "e625"):
        a = _stage_cells(data_dir, stage)
        x = np.asarray(a[:, gene].X).ravel()
        if log2:
            x = x / np.log(2)
        yield stage, x, (a.obs["cluster"] == "AVE").to_numpy()


def stage_aucs(data_dir: str, log2: bool = False, pooled: bool = False, auc_fn=_auc) -> dict:
    out = {}
    for gene, key in (("Dkk1", "dkk1"), ("Cer1", "cer1")):
        parts = list(_gene_by_stage(data_dir, gene, log2))
        if pooled:
            x = np.concatenate([p[1] for p in parts])
            y = np.concatenate([p[2] for p in parts])
            v = auc_fn(x, y)
            out[f"{key}_auc_e55"] = v
            out[f"{key}_auc_e625"] = v
        else:
            for stage, x, y in parts:
                out[f"{key}_auc_{stage}"] = auc_fn(x, y)
    return out


def _sk_auc(x, y) -> float:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(y.astype(int), x))


def ave04_facts(data_dir: Path) -> dict:
    d = str(data_dir)
    out = {}
    for stage in ("e55", "e625"):
        a = _stage_cells(d, stage)
        y = (a.obs["cluster"] == "AVE").to_numpy()
        for gene in ("Dkk1", "Cer1"):
            x = np.asarray(a[:, gene].X).ravel()
            out[f"{gene}_{stage}_mean_ave"] = float(x[y].mean())
            out[f"{gene}_{stage}_mean_emve"] = float(x[~y].mean())
            out[f"{gene}_{stage}_detected"] = float((x > 0).mean())
            out[f"{gene}_{stage}_p"] = float(ss.mannwhitneyu(x[y], x[~y]).pvalue)
        out[f"n_ave_{stage}"] = int(y.sum())
        out[f"n_emve_{stage}"] = int((~y).sum())
    pooled = stage_aucs(d, pooled=True)
    out["pooled_dkk1_auc"] = pooled["dkk1_auc_e55"]
    xs, ys = zip(*[(x, y) for _, x, y in _gene_by_stage(d, "Dkk1")])
    x, y = np.concatenate(xs), np.concatenate(ys)
    out["pooled_dkk1_p"] = float(ss.mannwhitneyu(x[y], x[~y]).pvalue)
    return out


def ave04(data_dir: Path) -> tuple[dict[str, dict], int]:
    d = str(data_dir)
    facts = ave04_facts(data_dir)
    ref = stage_aucs(d)
    key = decide_dkk1(ref, facts["Dkk1_e55_detected"])
    with_choice = lambda ans, c: {**ans, "conclusion": c}  # noqa: E731
    controls = {
        "reference": with_choice(ref, key),
        "log2_scale": with_choice(stage_aucs(d, log2=True), key),
        "roc_auc_score": with_choice(stage_aucs(d, auc_fn=_sk_auc), key),
        "stages_pooled": with_choice(stage_aucs(d, pooled=True), 1),
        "prior_driven_reading": with_choice(ref, 4),
    }
    return controls, key


# ============================================================================
# anchors
# ============================================================================

def _pam(dist: np.ndarray, k: int) -> np.ndarray:
    """Partitioning around medoids (BUILD then SWAP), deterministic."""
    n = len(dist)
    med = [int(np.argmin(dist.sum(1)))]
    while len(med) < k:
        costs = [np.minimum(dist[:, med].min(1), dist[:, c]).sum() if c not in med else np.inf for c in range(n)]
        med.append(int(np.argmin(costs)))
    best = dist[:, med].min(1).sum()
    improved = True
    while improved:
        improved = False
        for i in range(k):
            for c in range(n):
                if c in med:
                    continue
                trial = med.copy()
                trial[i] = c
                cost = dist[:, trial].min(1).sum()
                if cost < best - 1e-12:
                    med, best, improved = trial, cost, True
    return dist[:, med].argmin(1)


def pam_split(data_dir: Path) -> tuple[int, int, float]:
    """Split the 5.5 dpc embryonic VE in 2 with pam on the Spearman distance of its HVGs."""
    import scanpy as sc

    a = _embryonic_ve(str(data_dir))
    hv = a.copy()
    sc.pp.highly_variable_genes(hv, n_top_genes=3000)
    x = np.asarray(hv[:, hv.var["highly_variable"]].X, dtype=np.float64)
    r = np.apply_along_axis(ss.rankdata, 1, x)
    dist = np.sqrt(np.clip((1 - np.corrcoef(r)) / 2, 0, None))
    cl = _pam(dist, 2)
    cer1 = np.asarray(a[:, "Cer1"].X).ravel()
    ave = 0 if cer1[cl == 0].mean() > cer1[cl == 1].mean() else 1
    lab = pd.read_csv(io.derive(Path(data_dir), io.LABELS)).set_index("cell_id")
    truth = (lab.loc[a.obs_names, "cluster"] == "AVE").to_numpy()
    mine = cl == ave
    agree = float((mine == truth).mean())
    return int(mine.sum()), int((~mine).sum()), agree


def qc_both_stages(data_dir: Path) -> dict:
    """The paper's five thresholds at both stages, from the counts and the ENA read counts."""
    d = Path(data_dir)
    sdrf = io.read_sdrf(d)
    ena = pd.read_csv(io.derive(d, io.ENA_COUNTS), sep="\t").set_index("run_accession")["read_count"]
    genes = pd.read_csv(io.derive(d, io.GENES), sep="\t").set_index("gene_id")
    out = {}
    for stage in ("e55", "e625"):
        c = io.read_counts(d, stage)
        ercc = c.index.str.startswith("ERCC-")
        g = c.loc[~ercc]
        mito = genes.loc[g.index, "symbol"].str.startswith("mt-").to_numpy()
        gs, es, ms = g.sum(0).to_numpy(), c.loc[ercc].sum(0).to_numpy(), g.loc[mito].sum(0).to_numpy()
        meta = sdrf.loc[c.columns]
        rp = ena.loc[meta["Comment[ENA_RUN]"]].to_numpy() / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            n10 = np.nan_to_num(g.to_numpy() / gs * 1e6 > 10).sum(0)
        keep = ((n10 > 3000) & (np.log10(rp) > 4) & ((gs + es) / rp > 0.5) & (ms / rp < 0.1) & (es / rp < 0.3))
        single = (meta["well_type"] == "single cell").to_numpy()
        authors = meta["authors_pass"].to_numpy()
        out[stage] = {"thresholds": int(keep.sum()), "single_cell_wells": int((keep & single).sum()),
                      "agree_with_authors": int(((keep & single) == authors).sum()), "wells": len(keep),
                      "labelled": int(meta["cluster"].notna().sum())}
    return out


def dpt_gene_counts(data_dir: Path) -> dict:
    """Genes changing along pseudotime (natural spline F-test, FDR < 0.01) at both stages,
    rooted at the AVE end of DC1, on the authors' AVE and emVE cells."""
    import scanpy as sc

    d = Path(data_dir)
    lab = pd.read_csv(io.derive(d, io.LABELS)).set_index("cell_id")
    out = {}
    for stage in ("e55", "e625"):
        a = io.read_norm(d, stage)
        a = a[a.obs_names.isin(lab.index)].copy()
        a = a[lab.loc[a.obs_names, "cluster"].isin(["AVE", "emVE"]).to_numpy()].copy()
        sc.pp.log1p(a)
        hv = a.copy()
        sc.pp.highly_variable_genes(hv, n_top_genes=3000)
        hv = hv[:, hv.var["highly_variable"]].copy()
        sc.tl.pca(hv, svd_solver="arpack", random_state=0)
        hv.obsm["X_rep"] = _spearman_rep(np.asarray(hv.obsm["X_pca"][:, :10], dtype=np.float64))
        sc.pp.neighbors(hv, n_neighbors=15, use_rep="X_rep", random_state=0)
        sc.tl.diffmap(hv, n_comps=20, random_state=0)
        root, _ = _dc1_ends(hv, a)
        p = _pseudotime(hv, root)
        x = np.asarray(a.X, dtype=np.float64)
        keep = (x > 0).sum(0) >= 10
        sig = _bh(_ftest(p, x[:, keep], 4)) < 0.01
        rho, _ = _spearman_all(p, x[:, keep])
        out[stage] = {"n_sig": int(sig.sum()), "n_high": int((sig & (rho < 0)).sum()), "n_tested": int(keep.sum())}
    return out


def sema6d_plxna1(data_dir: Path) -> dict:
    import scanpy as sc

    d = Path(data_dir)
    lab = pd.read_csv(io.derive(d, io.LABELS)).set_index("cell_id")
    a = io.read_norm(d, "e55")
    a = a[a.obs_names.isin(lab.index)].copy()
    sc.pp.log1p(a)
    cl = lab.loc[a.obs_names, "cluster"].to_numpy()
    out = {}
    for g in ("Sema6d", "Plxna1"):
        x = np.asarray(a[:, g].X).ravel()
        out[g] = {c: float(x[cl == c].mean()) for c in ("AVE", "emVE", "Epi", "exVE", "ExE")}
    return out
