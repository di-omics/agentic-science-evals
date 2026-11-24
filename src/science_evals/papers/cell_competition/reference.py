"""Reference analyses, variants and traps for the cell competition tasks.

Every function here takes data loaded from the public files (see inputs.py and
prepare.py) and returns an answer dict, so that each control in a task is code
run on the data. Nothing here types in a number.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .trajectory import diffusion_map, spearman, spearman_dissimilarity

# --- cc01: which cells to keep --------------------------------------------------


def _mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x))) * 1.4826)


def _out_hi(x: np.ndarray, n: float) -> np.ndarray:
    return x > np.median(x) + n * _mad(x)


def _out_lo(x: np.ndarray, n: float) -> np.ndarray:
    return x < np.median(x) - n * _mad(x)


def qc_outliers(qc: pd.DataFrame, nmads: float, metrics=("lib", "genes", "ercc", "mito")) -> np.ndarray:
    """Outlier flags in the style of scater's perCellQCFilters: low library size
    and low detected genes on the log scale, high ERCC and mitochondrial
    fractions on the raw scale, each nmads median absolute deviations."""
    flag = np.zeros(len(qc), bool)
    if "lib" in metrics:
        flag |= _out_lo(qc["log10_total_reads"].to_numpy(float), nmads)
    if "genes" in metrics:
        flag |= _out_lo(np.log10(qc["n_genes_above_10rpm"].to_numpy(float)), nmads)
    if "ercc" in metrics:
        flag |= _out_hi(qc["frac_reads_ercc"].to_numpy(float), nmads)
    if "mito" in metrics:
        flag |= _out_hi(qc["frac_reads_mito"].to_numpy(float), nmads)
    return flag


def qc_outliers_by(qc: pd.DataFrame, nmads: float, by: str) -> np.ndarray:
    flag = np.zeros(len(qc), bool)
    for _, idx in qc.groupby(by).indices.items():
        flag[idx] = qc_outliers(qc.iloc[idx], nmads)
    return flag


def qc_answer(qc: pd.DataFrame, drop: np.ndarray) -> dict:
    kept = qc[~np.asarray(drop, bool)]
    return {
        "n_cells_kept": int(len(kept)),
        "n_ci_cells_kept": int((kept["treatment"] == "CI").sum()),
        "n_dmso_cells_kept": int((kept["treatment"] == "DMSO").sum()),
    }


# --- cc02: which way the trajectory runs ------------------------------------------


def log_tpm(tpm: pd.DataFrame) -> pd.DataFrame:
    return np.log10(tpm + 1.0)


def gene_correlations(score: np.ndarray, log_expr: pd.DataFrame, genes: dict[str, str]) -> dict:
    """Spearman rho between a per-cell score and log expression of named genes.
    genes: answer field -> gene id."""
    return {field: spearman(score, log_expr.loc[gid].to_numpy()) for field, gid in genes.items()}


# --- cc03: heteroplasmy in mt-Rnr1 by cluster -------------------------------------

RNR1_POSITIONS = (300, 303, 304, 305, 326, 327)
CLUSTERS = ("winner", "intermediate", "loser")


def heteroplasmy(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """counts: cells x positions x 4 alleles. Returns (depth, H) with
    H = 1 - frequency of the most common allele (NaN where depth is 0)."""
    depth = counts.sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        h = 1.0 - counts.max(axis=2) / depth
    return depth, h


def rnr1_fractions(depth: np.ndarray, h: np.ndarray, positions: np.ndarray, cluster: np.ndarray,
                   min_reads: int = 50, min_positions: int | None = 2000, mode: str = "complete",
                   cut: float = 0.01, strict: bool = True) -> dict:
    """Fraction of cells per cluster whose mean heteroplasmy over the six
    mt-Rnr1 positions is above `cut`.

    mode="complete": a cell counts only if all six positions have more than
                     min_reads reads (the paper's Fig. 6h rule).
    mode="available": average over the positions that are covered (at least 1).
    mode="zero_fill": every cell counts; positions without coverage count as no
                      heteroplasmy.
    mode="nan_false": every cell counts; cells without full coverage have a
                      missing mean, which never exceeds the cut (no data read
                      as no variant).
    """
    valid = depth > min_reads if strict else depth >= min_reads
    keep = np.ones(len(cluster), bool)
    if min_positions is not None:
        keep &= valid.sum(axis=1) > min_positions
    cols = [int(np.flatnonzero(positions == p)[0]) for p in RNR1_POSITIONS]
    v, hh = valid[:, cols], h[:, cols]
    masked = np.where(v, hh, np.nan)
    n_cov = v.sum(axis=1)
    mean_cov = np.divide(np.nansum(masked, axis=1), n_cov, out=np.full(len(n_cov), np.nan), where=n_cov > 0)
    if mode == "complete":
        use = keep & v.all(axis=1)
        mean_h = np.where(use, mean_cov, np.nan)
    elif mode == "available":
        use = keep & v.any(axis=1)
        mean_h = mean_cov
    elif mode == "zero_fill":
        use = keep
        mean_h = np.where(v, hh, 0.0).mean(axis=1)
    elif mode == "nan_false":
        # positions without coverage are missing; a missing mean compares as
        # False against the cut, and the cell stays in the denominator
        use = keep
        mean_h = np.where(v.all(axis=1), mean_cov, np.nan)
    else:
        raise ValueError(mode)
    out = {}
    for name in CLUSTERS:
        s = use & (cluster == name)
        with np.errstate(invalid="ignore"):
            above = np.nan_to_num(mean_h[s], nan=-1.0) > cut
        out[f"frac_heteroplasmic_{name}"] = float(np.mean(above)) if s.any() else 0.0
    out["n_cells_used"] = int(use.sum())
    return out


# --- cc04: over-representation of respiratory chain genes ------------------------

RESP_CHAIN_PATTERN = r"^(Nduf[abcsv][0-9]|Sdh[abcd]$|Uqcr|Cyc1$|Cox[4-8]|Atp5|mt-(Nd|Co|Cytb|Atp))"


def trend_test(log_expr: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares slope of each row of log_expr on score, with a two-sided
    t-test p-value: the linear part of the authors' GAM fit."""
    from scipy.stats import t as tdist

    y = np.asarray(log_expr, dtype=float)
    x = np.asarray(score, dtype=float) - np.mean(score)
    n = len(x)
    sxx = float((x ** 2).sum())
    yc = y - y.mean(axis=1, keepdims=True)
    slope = yc @ x / sxx
    resid = yc - np.outer(slope, x)
    s2 = (resid ** 2).sum(axis=1) / (n - 2)
    se = np.sqrt(s2 / sxx)
    with np.errstate(divide="ignore", invalid="ignore"):
        tval = np.where(se > 0, slope / se, 0.0)
    p = 2 * tdist.sf(np.abs(tval), n - 2)
    return slope, p


def spearman_test(log_expr: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from scipy.stats import rankdata, t as tdist

    ry = np.apply_along_axis(rankdata, 1, np.asarray(log_expr, dtype=float))
    rx = rankdata(score)
    ryc = ry - ry.mean(axis=1, keepdims=True)
    rxc = rx - rx.mean()
    denom = np.sqrt((ryc ** 2).sum(axis=1) * (rxc ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(denom > 0, ryc @ rxc / denom, 0.0)
    n = len(rx)
    tval = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
    p = 2 * tdist.sf(np.abs(tval), n - 2)
    return rho, p


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def over_representation(in_set: np.ndarray, called: np.ndarray) -> dict:
    """One-sided Fisher exact test of `in_set` among `called`, both boolean
    over the same universe of genes."""
    from scipy.stats import fisher_exact

    a = int((in_set & called).sum())
    b = int((in_set & ~called).sum())
    c = int((~in_set & called).sum())
    dd = int((~in_set & ~called).sum())
    odds, p = fisher_exact([[a, b], [c, dd]], alternative="greater")
    return {
        "frac_genes_decreasing": float(called.mean()),
        "frac_resp_chain_decreasing": float(a / (a + b)) if a + b else 0.0,
        "odds_ratio": float(odds),
        "n_universe": int(len(called)),
        "n_resp_chain": int(a + b),
        "neg_log10_p": float(-np.log10(max(p, 1e-300))),
    }


def heteroplasmy_scan(counts: np.ndarray, positions: np.ndarray, cluster: np.ndarray, score: np.ndarray,
                      min_reads: int = 50, min_positions: int = 2000, min_het: float = 0.01,
                      min_cells: int = 10, fdr: float = 1e-3) -> dict:
    """The paper's heteroplasmy filters, then Spearman's test of heteroplasmy
    against the losing score at each remaining position (BH-adjusted)."""
    from scipy.stats import spearmanr

    depth, h = heteroplasmy(counts)
    valid = depth > min_reads
    keep = valid.sum(axis=1) > min_positions
    ok = np.ones(depth.shape[1], bool)
    for name in CLUSTERS:
        s = keep & (cluster == name)
        ok &= valid[s].sum(axis=0) > 0.5 * s.sum()
    n_het = (np.where(valid & keep[:, None], h, 0.0) > min_het).sum(axis=0)
    relevant = np.flatnonzero(ok & (n_het > min_cells))
    rho, p = [], []
    for j in relevant:
        s = keep & valid[:, j]
        r = spearmanr(score[s], h[s, j])
        rho.append(float(r.statistic))
        p.append(float(r.pvalue))
    q = bh(np.array(p))
    rho = np.array(rho)
    inc = sorted(int(positions[j]) for j, qq, rr in zip(relevant, q, rho) if qq < fdr and rr > 0)
    return {"n_cells": int(keep.sum()), "n_positions": int(ok.sum()), "n_tested": int(len(relevant)),
            "increasing": inc}


# --- the epiblast trajectory, shared by cc02, cc04 and the anchors ----------------------

N_DMSO_NEIGHBOURS = 20


def epiblast_trajectory(tpm: pd.DataFrame, hvg: list[str], epi: pd.DataFrame) -> dict:
    """The authors' diffusion map of CI epiblast cells, rebuilt, with the root
    chosen without cluster labels: the tip whose nearest neighbours among all
    epiblast cells include the most cells from DMSO-treated embryos."""
    cells = list(epi.index)
    ci_mask = (epi["treatment"] == "CI").to_numpy()
    lxh = log_tpm(tpm.loc[[g for g in hvg if g in tpm.index], cells])
    diss = spearman_dissimilarity(lxh.to_numpy())
    ci_idx = np.flatnonzero(ci_mask)
    dm = diffusion_map(diss[np.ix_(ci_idx, ci_idx)])
    tips = dm.tips(0)

    def dmso_nbrs(i_all: int) -> int:
        order = np.argsort(diss[i_all], kind="stable")[1:N_DMSO_NEIGHBOURS + 1]
        return int((~ci_mask[order]).sum())

    nbrs = [dmso_nbrs(ci_idx[t]) for t in tips]
    root, far = tips[int(np.argmax(nbrs))], tips[int(np.argmin(nbrs))]
    # the same map built on all epiblast cells, evaluated on the CI cells
    dm_all = diffusion_map(diss)
    tips_all = dm_all.tips(0)
    nbrs_all = [dmso_nbrs(t) for t in tips_all]
    root_all = tips_all[int(np.argmax(nbrs_all))]
    return {
        "ci_cells": [cells[i] for i in ci_idx],
        "ci_idx": ci_idx,
        "ci_mask": ci_mask,
        "lx_hvg": lxh,
        "tips": tips,
        "tip_dmso_neighbours": nbrs,
        "root": root,
        "score": dm.dpt(root),
        "reversed": dm.dpt(far),
        "score_all_epiblast_map": dm_all.dpt(root_all)[ci_idx],
        "dmso_nbrs": dmso_nbrs,
    }


def scanpy_dpt(x_cells_by_genes: np.ndarray, end_dmso_nbrs, n_pcs: int = 30, n_neighbors: int = 15) -> np.ndarray:
    """Diffusion pseudotime the scanpy way (PCA, kNN graph, diffmap, dpt), rooted
    at the end of the first diffusion component closer to DMSO cells."""
    import anndata as ad
    import scanpy as sc

    a = ad.AnnData(X=np.asarray(x_cells_by_genes, dtype=np.float32))
    sc.pp.pca(a, n_comps=n_pcs, random_state=0)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, random_state=0)
    sc.tl.diffmap(a, n_comps=15)
    dc1 = a.obsm["X_diffmap"][:, 1]
    ends = [int(np.argmin(dc1)), int(np.argmax(dc1))]
    a.uns["iroot"] = ends[int(np.argmax([end_dmso_nbrs(e) for e in ends]))]
    sc.tl.dpt(a)
    return a.obs["dpt_pseudotime"].to_numpy(dtype=float)


def scanpy_own_hvg_dpt(tpm_ci: pd.DataFrame, end_dmso_nbrs, n_top: int = 2000) -> np.ndarray:
    import anndata as ad
    import scanpy as sc

    b = ad.AnnData(X=np.log1p(tpm_ci.to_numpy().T).astype(np.float32))
    sc.pp.highly_variable_genes(b, n_top_genes=n_top, flavor="seurat")
    return scanpy_dpt(b.X[:, b.var["highly_variable"].to_numpy()], end_dmso_nbrs)


# --- genes along the trajectory (cc04 and the anchors) ---------------------------------


def trend_genes(tpm_ci: pd.DataFrame, score: np.ndarray, keep: np.ndarray, test: str = "linear",
                fdr: float = 0.01) -> dict:
    sub = tpm_ci[np.asarray(keep, bool)]
    lx = np.log10(sub.to_numpy(dtype=float) + 1.0)
    stat, p = trend_test(lx, score) if test == "linear" else spearman_test(lx, score)
    q = bh(p)
    ids = sub.index.to_numpy()
    return {"ids": ids, "down": ids[(stat < 0) & (q < fdr)], "up": ids[(stat > 0) & (q < fdr)]}


def expressed(tpm_ci: pd.DataFrame, rule: str) -> np.ndarray:
    t = tpm_ci.to_numpy(dtype=float)
    if rule == "gt15_in_gt10":  # the paper: more than 15 TPM in more than 10 cells
        return (t > 15).sum(axis=1) > 10
    if rule == "gt1_in_10pct":
        return (t > 1).mean(axis=1) >= 0.10
    if rule == "mean_ge10":
        return t.mean(axis=1) >= 10
    if rule == "detected_in_3":
        return (t > 0).sum(axis=1) >= 3
    raise ValueError(rule)


def enrichment_answer(universe: np.ndarray, down: np.ndarray, in_set: np.ndarray) -> dict:
    called = np.isin(universe, down)
    r = over_representation(np.asarray(in_set, bool), called)
    return {
        "frac_genes_decreasing": r["frac_genes_decreasing"],
        "frac_resp_chain_decreasing": r["frac_resp_chain_decreasing"],
        "fold_enrichment": r["frac_resp_chain_decreasing"] / r["frac_genes_decreasing"],
        "neg_log10_p": r["neg_log10_p"],
    }


def marker_agreement(a: np.ndarray, b: np.ndarray, ref: np.ndarray, positions: np.ndarray,
                     markers=(12217, 13536), min_reads: int = 20) -> dict:
    """Do two allele tables of the same cells agree on which cells carry an
    inherited non-reference allele? A check that rows were matched to the
    right cells."""
    out = {"cells": 0, "carriers": 0, "agree": 0}
    for p in markers:
        j = int(np.flatnonzero(positions == p)[0])
        r = "ACGT".index(str(ref[j]))
        da, db = a[:, j].sum(axis=1), b[:, j].sum(axis=1)
        ok = (da > min_reads) & (db > min_reads)
        ca = a[ok, j, r] / da[ok] < 0.5
        cb = b[ok, j, r] / db[ok] < 0.5
        out["cells"] += int(ok.sum())
        out["carriers"] += int(ca.sum())
        out["agree"] += int((ca == cb).sum())
    return out
