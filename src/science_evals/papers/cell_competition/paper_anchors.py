"""Statements from the paper, recomputed from the public files."""
from __future__ import annotations

import numpy as np

from . import reference as R
from .inputs import fetch_public, read_rdata
from .prepare import Data, ci_epiblast_cells


def _row(paper_says: str, where: str, recomputed: str, agrees: bool) -> dict:
    verdict = "Reproduces." if agrees else "Does not fully reproduce."
    return {"paper_says": paper_says, "where": where, "recomputed": f"{verdict} {recomputed}",
            "agrees": bool(agrees)}


def anchor_qc(d: Data) -> dict:
    meta = read_rdata(fetch_public(d.data_dir, "input_for_figures.Rdata"))["meta_dati"]
    n_all = len(meta)
    n_good = int((meta["Quality"].astype(str) == "good").sum())
    qc = d.qc
    passes = ((qc["log10_total_reads"] >= 4) & (qc["frac_mapped_reads"] >= 0.8)
              & (qc["n_genes_above_10rpm"] >= 3000) & (qc["frac_reads_endogenous"] >= 0.5))
    return _row(
        "1,495 cells were sequenced and 723 passed the quality control thresholds.",
        "Methods, data processing, quality control and normalization",
        f"{n_all:,} cells in the authors' cell table, {n_good} marked good; all {int(passes.sum())} of the "
        f"{len(qc)} cells in their QC table meet the 4 stated thresholds.",
        n_all == 1495 and n_good == 723 and int(passes.sum()) == len(qc) == 723,
    )


def anchor_ci_only_cluster(d: Data) -> dict:
    c = d.cells
    los = c[c["cluster"] == "loser"]
    n_ci, n_dmso = int((los["treatment"] == "CI").sum()), int((los["treatment"] == "DMSO").sum())
    return _row(
        "One of the three epiblast clusters (cluster 4) is only composed of cells from CI-treated embryos.",
        "Results and Fig. 1c-e",
        f"The authors' cluster 4 holds {n_ci} cells from CI-treated embryos and {n_dmso} from DMSO-treated "
        f"embryos.",
        n_dmso == 0 and n_ci > 0,
    )


def anchor_rnr1_fractions(d: Data) -> dict:
    e = ci_epiblast_cells(d)
    clu = e["cluster"].to_numpy()

    def fractions(allele: dict) -> dict:
        idx = np.array([allele["cells"].index(c) for c in e.index])
        depth, h = R.heteroplasmy(allele["counts"][idx])
        return R.rnr1_fractions(depth, h, allele["positions"], clu)

    r, cap = fractions(d.allele), fractions(d.allele_capped)
    w, i, lo = (r[f"frac_heteroplasmic_{k}"] for k in R.CLUSTERS)
    cw, ci_, clo = (cap[f"frac_heteroplasmic_{k}"] for k in R.CLUSTERS)
    return _row(
        "About 5% of winner, 55% of intermediate and 87% of loser epiblast cells have a mean mt-Rnr1 "
        "heteroplasmy above 0.01 (215 cells).",
        "Fig. 6h and legend; Methods, analysis of heteroplasmy",
        f"{w:.0%}, {i:.0%} and {lo:.0%} of {r['n_cells_used']} cells in the E-MTAB-8640 allele table with the "
        f"paper's coverage rules. The re-count shipped with MitoHEAR, made with a pileup that keeps about "
        f"250 reads per position, gives {cw:.0%}, {ci_:.0%} and {clo:.0%} of {cap['n_cells_used']} cells.",
        abs(w - 0.05) <= 0.03 and abs(i - 0.55) <= 0.08 and abs(lo - 0.87) <= 0.08,
    )


PAPER_POSITIONS = {300, 303, 304, 305, 326, 327, 2031, 2074, 2077, 2079, 2081}


def anchor_positions(d: Data) -> dict:
    """Coverage filters and positions whose heteroplasmy rises with the losing score."""
    e = ci_epiblast_cells(d)
    idx = np.array([d.allele["cells"].index(c) for c in e.index])
    res = R.heteroplasmy_scan(d.allele["counts"][idx], d.allele["positions"], e["cluster"].to_numpy(),
                              e["losing_score"].to_numpy())
    found = set(res["increasing"])
    shared = found & PAPER_POSITIONS
    return _row(
        "After coverage filters 259 cells and 5,192 mtDNA positions remain; 11 positions, 6 in mt-Rnr1 and 5 "
        "in mt-Rnr2, have heteroplasmy that rises with the losing score (FDR < 0.001); Spearman's test gives "
        "substantially the same result.",
        "Methods, analysis of heteroplasmy; Fig. 6 and Extended Data Fig. 7",
        f"{res['n_cells']} cells and {res['n_positions']:,} positions pass the same filters; Spearman's test "
        f"finds {len(found)} positions rising with the authors' losing score at FDR < 0.001, {len(shared)} of "
        f"them among the paper's 11 ({', '.join(str(p) for p in sorted(found))}).",
        len(shared) >= 9 and abs(res["n_cells"] - 259) <= 5,
    )


def anchor_trajectory(d: Data, traj: dict) -> dict:
    e = ci_epiblast_cells(d)
    med = e.groupby("cluster")["losing_score"].median()
    return _row(
        "The winner and loser clusters sit at the two ends of the diffusion trajectory of CI epiblast "
        "cells, so the diffusion pseudotime from the winner end reads as a losing score.",
        "Results and Fig. 2a, 2d-e",
        f"Median losing score in the authors' table: winner {med['winner']:.2f}, intermediate "
        f"{med['intermediate']:.2f}, loser {med['loser']:.2f}. A Python rebuild of their diffusion map "
        f"(Spearman distance over their 1,921 variable genes, destiny-style kernel) rooted at the winner "
        f"tip gives Spearman rho {traj['rho_with_authors']:.4f} with the authors' score.",
        med["winner"] < med["intermediate"] < med["loser"] and traj["rho_with_authors"] > 0.9,
    )


def anchor_trend_genes(d: Data, trend: dict) -> dict:
    return _row(
        "5,311 genes are differentially expressed along the trajectory (FDR < 0.01), out of the genes with "
        "more than 15 TPM in more than 10 cells.",
        "Methods, differential gene expression analysis along the trajectory",
        f"{trend['n_tested']:,} genes pass the expression filter; a linear trend test of log10(TPM + 1) on the "
        f"authors' losing score (the parametric part of their GAM) calls {trend['n_de']:,} at FDR < 0.01, "
        f"{trend['n_down']:,} decreasing and {trend['n_up']:,} increasing toward the loser end.",
        abs(trend["n_de"] - 5311) / 5311 <= 0.15,
    )


def anchors(d: Data, traj: dict | None = None, trend: dict | None = None) -> list[dict]:
    out = [anchor_qc(d), anchor_ci_only_cluster(d)]
    if traj is not None:
        out.append(anchor_trajectory(d, traj))
    if trend is not None:
        out.append(anchor_trend_genes(d, trend))
    out += [anchor_positions(d), anchor_rnr1_fractions(d)]
    return out
