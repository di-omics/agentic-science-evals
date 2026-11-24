"""Statements from the paper, recomputed from the public files."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import deseq as D
from . import inputs as I
from . import reference as R

PAIRS = (("E2", "E3"), ("E4", "E3"), ("E4", "E2"), ("KO", "E3"))


def _atac_fit(counts: pd.DataFrame, drop=("KO_1", "KO_2")):
    cols = [c for c in counts.columns if c not in drop]
    groups = [c.split("_")[0] for c in cols]
    keep = D.filter_by_expr(counts[cols].to_numpy(dtype=float), groups)
    fit = D.fit(counts[cols].to_numpy(dtype=float)[keep], groups)
    return keep, fit, counts.index[keep]


def _calls(fit, a, b):
    r = fit.contrast(a, b)
    sig = D.bh(r["pvalue"]) < 0.05
    return int(sig.sum()), int((sig & (r["log2FoldChange"] > 0)).sum()), int((sig & (r["log2FoldChange"] < 0)).sum())


def _fmt3(t):
    return f"{t[0]} ({t[1]} up, {t[2]} down)"


def anchors(data_dir: Path) -> list[dict]:
    data_dir = Path(data_dir)
    out = []

    # 1. peaks kept by filterByExpr, on the GEO matrix
    geo = pd.read_csv(data_dir / I.ATAC_COUNTS_GEO, index_col=0).iloc[:, 5:]
    repo = pd.read_csv(data_dir / I.ATAC_COUNTS_REPO, index_col=0).iloc[:, 5:]
    keep_geo, fit_geo, _ = _atac_fit(geo)
    keep_repo, fit_repo, peaks_repo = _atac_fit(repo)
    out.append({
        "paper_says": "Peaks filtered with filterByExpr; this left 167,951 peaks for downstream analyses.",
        "where": "Methods, ATAC-seq QC and processing",
        "recomputed": (f"{int(keep_geo.sum()):,} peaks of the {len(geo):,} in the GEO matrix (GSE271384) pass edgeR's "
                       f"filterByExpr rule with the 16 samples the paper kept. The same rule on the count matrix in "
                       f"the authors' repository, a different consensus peak set of {len(repo):,} peaks, keeps "
                       f"{int(keep_repo.sum()):,}, the set behind the published DAR tables."),
        "agrees": bool(int(keep_geo.sum()) == 167951),
    })

    # 2. differentially accessible regions
    paper_dars = {("E2", "E3"): (40, 24, 16), ("E4", "E3"): (50, 38, 12), ("E4", "E2"): (72, 52, 20),
                  ("KO", "E3"): (22, 14, 8)}
    got = {p: _calls(fit_repo, *p) for p in PAIRS}
    got_geo = {p: _calls(fit_geo, *p) for p in PAIRS}
    out.append({
        "paper_says": "DARs at FDR < 0.05: APOE2 vs APOE3 40 (24 up, 16 down), APOE4 vs APOE3 50 (38 up, 12 down), "
                      "APOE4 vs APOE2 72 (52 up, 20 down), APOE-KO vs APOE3 14 up and 8 down.",
        "where": "Results, chromatin accessibility; Fig. 3a-b, Supplementary Fig. 4",
        "recomputed": ("From the authors' repository counts with DESeq2's default steps (this package's deseq.py) "
                       "and BH over all filtered peaks: "
                       + "; ".join(f"{a} vs {b} {_fmt3(got[(a, b)])}" for a, b in PAIRS)
                       + ". The GEO matrix gives "
                       + "; ".join(f"{a} vs {b} {_fmt3(got_geo[(a, b)])}" for a, b in PAIRS) + "."),
        "agrees": all(got[p] == paper_dars[p] for p in PAIRS),
    })

    # 3. differentially expressed genes, with the tximport lengths
    rna = R.load_rna(data_dir)
    failed = R.failed_knockouts(rna)
    cols = [c for c in rna.counts.columns if c not in failed]
    degs = R.load_degs(data_dir)
    res, tables = {}, {}
    for a, b in PAIRS:
        # BH over the genes the authors tested after their symbol merge: the Entrez IDs of their table
        universe = degs[f"{a}vs{b}"]["entrezgene_id"]
        r = R.de(rna, cols, a, b, "offsets")
        r = r.loc[r.index.isin(universe) & r["pvalue"].notna()].copy()
        r["fdr"] = D.bh(r["pvalue"])
        tables[(a, b)] = r
        sig = r["fdr"] < 0.05
        res[(a, b)] = (int(sig.sum()), int((sig & (r["log2FoldChange"] > 0)).sum()),
                       int((sig & (r["log2FoldChange"] < 0)).sum()))
    plain = {p: int((R.de(rna, cols, *p, "counts")["fdr"] < 0.05).sum()) for p in PAIRS}
    paper_degs = {("E2", "E3"): 470, ("E4", "E3"): 332, ("E4", "E2"): 1639, ("KO", "E3"): 332}
    close = all(abs(res[p][0] - paper_degs[p]) / paper_degs[p] < 0.02 for p in PAIRS)
    within = max(abs(res[p][0] - paper_degs[p]) / paper_degs[p] for p in PAIRS)
    out.append({
        "paper_says": "DEGs at FDR < 0.05: APOE2 vs APOE3 470 (282 up, 188 down), APOE4 vs APOE3 332 (170 up, 158 "
                      "down), APOE4 vs APOE2 1639 (751 up, 888 down); APOE-KO vs APOE3 126 up and 206 down.",
        "where": "Results, first section; Fig. 1d, Fig. 2a-c, Supplementary Fig. 3",
        "recomputed": ("GEO counts with the tximport lengths from the authors' DESeqDataSet as offsets, the 17 "
                       f"samples kept ({', '.join(failed)} excluded), BH over the genes in the authors' table: "
                       + "; ".join(f"{a} vs {b} {_fmt3(res[(a, b)])}" for a, b in PAIRS)
                       + f". Largest difference from the text {100 * within:.1f}%; Supplementary Data 2 lists "
                       + "; ".join(f"{a} vs {b} {int((degs[f'{a}vs{b}']['fdr'] < 0.05).sum())}" for a, b in PAIRS)
                       + ". Without the lengths the counts give "
                       + "; ".join(f"{a} vs {b} {plain[(a, b)]}" for a, b in PAIRS)
                       + ". The bar labels of Fig. 1d (374, 217, 1424, 461) match neither."),
        "agrees": bool(close),
    })

    # 4. CHCHD2: expression and promoter accessibility
    r43, r42 = tables[("E4", "E3")], tables[("E4", "E2")]
    chchd2 = 51142
    peaks = R.load_promoter_peaks(data_dir)
    pk = peaks[(peaks["SYMBOL"] == "CHCHD2") & (peaks["distanceToTSS"] == 0)]["Peak"].iloc[0]
    i = list(peaks_repo).index(pk)
    a43 = fit_repo.contrast("E4", "E3")
    a42 = fit_repo.contrast("E4", "E2")
    fdr43 = D.bh(a43["pvalue"])[i]
    out.append({
        "paper_says": "CHCHD2 lower in APOE4 than APOE3 (logFC -7.8, p 1.3e-11) and APOE2 (logFC -7.4, p 7.8e-11); "
                      "accessibility at its TSS lower than APOE3 (logFC -7.8, p 1.2e-7) and APOE2 (logFC -6.5, "
                      "p 8.3e-9).",
        "where": "Results, chromatin accessibility; Fig. 3c-d",
        "recomputed": (f"Expression: log2FC {r43.loc[chchd2, 'log2FoldChange']:.2f} (p {r43.loc[chchd2, 'pvalue']:.1e}, "
                       f"FDR {r43.loc[chchd2, 'fdr']:.1e}) against APOE3 and {r42.loc[chchd2, 'log2FoldChange']:.2f} "
                       f"(p {r42.loc[chchd2, 'pvalue']:.1e}, FDR {r42.loc[chchd2, 'fdr']:.1e}) against APOE2. Promoter "
                       f"peak {pk.rsplit('_', 1)[1]}: log2FC {a43['log2FoldChange'][i]:.2f} (p {a43['pvalue'][i]:.1e}, "
                       f"FDR {fdr43:.1e}) against APOE3 and {a42['log2FoldChange'][i]:.2f} (p {a42['pvalue'][i]:.1e}) "
                       f"against APOE2. The text's expression p-values and its accessibility p-value against APOE3 "
                       f"are the FDRs. CHCHD2 counts are near zero in all 5 APOE4 and all 5 APOE-KO grafts."),
        "agrees": bool(abs(r43.loc[chchd2, "log2FoldChange"] + 7.8) < 0.15 and abs(a43["log2FoldChange"][i] + 7.8) < 0.1
                       and abs(a42["log2FoldChange"][i] + 6.5) < 0.1),
    })

    # 5. knockout confirmed in 3 of 5 grafts
    cpm = R.apoe_cpm(rna)
    e3 = cpm[rna.columns(["E3"])].median()
    rel = {c: cpm[c] / e3 for c in rna.columns(["KO"])}
    out.append({
        "paper_says": "APOE expression was significantly lower in three out of five knockout samples, and only these "
                      "were used in downstream analyses.",
        "where": "Results, first section; Fig. 1b, Supplementary Fig. 2",
        "recomputed": ("APOE counts per million relative to the APOE3 median: "
                       + ", ".join(f"{c} {100 * v:.0f}%" for c, v in rel.items())
                       + f". The rule used here (at least 25% of the APOE3 median) flags {', '.join(failed)}, the two "
                         "samples the authors' code removes."),
        "agrees": bool(sorted(failed) == ["KO-1", "KO-2"]),
    })

    # 6. concordance of expression and promoter accessibility
    matched = R.concordance(degs, peaks)
    pooled = R.concordance(degs, peaks, pooled=True)
    rank = R.spearman_concordance(degs, peaks)
    tags = (("e2_vs_e3", "E2vsE3"), ("e4_vs_e3", "E4vsE3"), ("e4_vs_e2", "E4vsE2"), ("ko_vs_e3", "KOvsE3"))
    out.append({
        "paper_says": "We correlated the logFC values between DEGs and DARs at the corresponding promoter peaks and "
                      "observed a strong correlation across all APOE comparisons.",
        "where": "Results, chromatin accessibility; Supplementary Fig. 6",
        "recomputed": ("Pearson r over DEGs with a promoter peak, same comparison in both assays: "
                       + "; ".join(f"{c} {matched['r_' + t]:.2f} (n {matched['n_' + t]})" for t, c in tags)
                       + ". With each gene's accessibility pooled over the four comparisons, as the code in the "
                         "authors' repository computes it: "
                       + "; ".join(f"{c} {pooled['r_' + t]:.2f}" for t, c in tags)
                       + ". Spearman correlations (same comparison) are lower: "
                       + "; ".join(f"{c} {rank[c]:.2f}" for t, c in tags)
                       + ", so genes with large changes in both assays carry much of the Pearson correlation."),
        "agrees": bool(min(matched["r_" + t] for t, _ in tags) > 0.4),
    })
    return out
