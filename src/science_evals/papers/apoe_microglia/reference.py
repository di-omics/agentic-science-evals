"""Reference analyses, defensible variants and traps, as functions of the task files.

Each task's controls are computed here from the files the agent receives, so
every number in an eval.json comes from running code on the data.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import deseq as D
from . import inputs as I

APOE = 348  # Entrez gene ID of APOE
SHORT = {"APOE2": "E2", "APOE3": "E3", "APOE4": "E4", "APOE-KO": "KO"}
SHEETS = {"E2vsE3": "E2vsE3_DEGs", "E4vsE3": "E4vsE3_DEGs", "E4vsE2": "E4vE2_DEGs", "KOvsE3": "KOvsE3_DEGs"}


# --- RNA-seq: counts, tximport lengths, samples --------------------------------


@dataclass
class Rna:
    counts: pd.DataFrame  # genes (Entrez ID) x samples
    lengths: pd.DataFrame  # same shape: tximport average transcript lengths
    samples: pd.DataFrame  # one row per sample
    symbols: pd.Series  # Entrez ID -> symbol

    def group(self, cols) -> list[str]:
        g = self.samples.set_index("sample")["genotype"].map(SHORT)
        return [g[c] for c in cols]

    def columns(self, groups, exclude=()) -> list[str]:
        return [c for c, g in zip(self.counts.columns, self.group(self.counts.columns))
                if g in groups and c not in exclude]


def load_rna(data_dir: Path) -> Rna:
    counts = pd.read_csv(data_dir / I.RNA_COUNTS, index_col=0)
    lengths = pd.read_csv(data_dir / I.LENGTHS, index_col=0)
    if list(lengths.columns) != list(counts.columns) or not lengths.index.equals(counts.index):
        raise ValueError("length matrix and count matrix do not line up")
    samples = pd.read_csv(data_dir / I.SAMPLES)
    genes = pd.read_csv(data_dir / I.GENES)
    return Rna(counts, lengths, samples, genes.set_index("entrez_id")["symbol"])


def apoe_cpm(rna: Rna) -> pd.Series:
    return rna.counts.loc[APOE] / rna.counts.sum() * 1e6


def failed_knockouts(rna: Rna, frac: float = 0.25) -> list[str]:
    """KO samples whose APOE level is at least `frac` of the APOE3 median (counts per million).

    A graft of the knockout line should carry little APOE mRNA; the three
    samples the paper kept sit near 5% of the APOE3 level, the other two near 50%.
    """
    cpm = apoe_cpm(rna)
    e3 = cpm[rna.columns(["E3"])].median()
    return [c for c in rna.columns(["KO"]) if cpm[c] >= frac * e3]


def de(rna: Rna, cols, num: str, den: str, mode: str = "offsets", adjust: str = "bh") -> pd.DataFrame:
    """Negative binomial GLM with genotype as the only factor (DESeq2 defaults), one Wald contrast.

    mode: "offsets" (tximport lengths as DESeq2 normalisation factors), "counts"
    (plain counts, size factors only) or "length_scaled" (lengths divided out
    of the counts first). adjust: "bh" over every gene with a p-value, or
    "independent_filtering" (DESeq2's results() padj, alpha 0.1).
    """
    y = rna.counts[cols].to_numpy(dtype=float)
    lengths = rna.lengths[cols].to_numpy(dtype=float)
    nf = None
    if mode == "offsets":
        nf = D.norm_factors_from_lengths(y, lengths)
    elif mode == "length_scaled":
        y = D.length_scaled_counts(y, lengths)
    elif mode != "counts":
        raise ValueError(mode)
    fit = D.fit(y, rna.group(cols), norm_factors=nf)
    r = pd.DataFrame(fit.contrast(num, den), index=rna.counts.index)
    if adjust == "bh":
        r["fdr"] = D.bh(r["pvalue"])
    else:
        r["fdr"] = D.padj_independent_filtering(r["baseMean"], r["pvalue"])
    return r


def pydeseq2_de(rna: Rna, cols, num: str, den: str, mode: str, cache_dir: Path) -> pd.DataFrame:
    """The same test run with the pydeseq2 package (BH over every gene with a p-value).

    pydeseq2 takes size factors only, so lengths can enter only as length-scaled
    counts. Results are cached under data/ keyed by the inputs; the first run
    writes them with fixed precision and every run reads them back the same way.
    """
    if mode == "counts":
        y = rna.counts[cols].to_numpy(dtype=float)
    elif mode == "length_scaled":
        y = D.length_scaled_counts(rna.counts[cols].to_numpy(dtype=float), rna.lengths[cols].to_numpy(dtype=float))
    else:
        raise ValueError(mode)
    import pydeseq2

    key = hashlib.sha256(y.tobytes() + repr((list(cols), rna.group(cols), num, den, pydeseq2.__version__)).encode()).hexdigest()
    path = Path(cache_dir) / f"pydeseq2_{key[:20]}.csv"
    if not path.exists():
        import warnings

        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        warnings.filterwarnings("ignore")
        counts = pd.DataFrame(y.T.astype(int), index=list(cols), columns=[str(g) for g in rna.counts.index])
        meta = pd.DataFrame({"group": rna.group(cols)}, index=list(cols))
        dds = DeseqDataSet(counts=counts, metadata=meta, design="~group", quiet=True, n_cpus=4)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["group", num, den], quiet=True, n_cpus=4)
        st.summary()
        res = st.results_df[["baseMean", "log2FoldChange", "pvalue", "padj"]]
        path.parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(path, float_format="%.10g", lineterminator="\n")
    r = pd.read_csv(path, index_col=0)
    r.index = r.index.astype(int)
    r = r.reindex(rna.counts.index)
    r["fdr"] = D.bh(r["pvalue"])
    return r


def de_summary(r: pd.DataFrame, up_name: str) -> dict:
    sig = r["fdr"] < 0.05
    return {"n_deg": int(sig.sum()), up_name: int((sig & (r["log2FoldChange"] > 0)).sum()),
            "apoe_log2fc": float(r.loc[APOE, "log2FoldChange"])}


# --- VDR targets among genes higher in APOE2 -----------------------------------


def load_degs(data_dir: Path) -> dict:
    sheets = pd.read_excel(data_dir / I.DEGS, sheet_name=list(SHEETS.values()))
    return {k: sheets[v] for k, v in SHEETS.items()}


def load_vdr(data_dir: Path) -> pd.DataFrame:
    return pd.read_csv(data_dir / I.VDR, encoding="utf-8-sig")


def vdr_induced(vdr: pd.DataFrame) -> set:
    t = vdr.dropna(subset=["GeneSymbol", "log2FoldChange", "padj"])
    return set(t.loc[(t["log2FoldChange"] > 0) & (t["padj"] < 0.05), "GeneSymbol"])


def vdr_test(sheet: pd.DataFrame, targets: set, higher_in_apoe2_sign: int = -1, unit: str = "symbol") -> dict:
    """One-sided hypergeometric test of VDR targets among genes higher in APOE2.

    The sheet is APOE4 against APOE2, so genes higher in APOE2 have a negative
    log2 fold change (higher_in_apoe2_sign=-1). unit "symbol": unique symbols,
    universe = symbols in the sheet, targets counted only inside it.
    unit "authors": rows, universe = all rows, every target counted (the
    authors' code).
    """
    from scipy.stats import hypergeom

    sig = sheet[(sheet["fdr"] < 0.05) & (np.sign(sheet["log2FoldChange"]) == higher_in_apoe2_sign)]
    if unit == "symbol":
        universe = set(sheet["hgnc_symbol"].dropna())
        chosen = set(sig["hgnc_symbol"].dropna())
        m, total, n_chosen = len(targets & universe), len(universe), len(chosen)
    elif unit == "authors":
        chosen = set(sig["hgnc_symbol"].dropna())
        m, total, n_chosen = len(targets), len(sheet), len(sig)
    else:
        raise ValueError(unit)
    hits = sorted(chosen & targets)
    k = len(hits)
    p = float(hypergeom.sf(k - 1, total, m, n_chosen)) if k > 0 else 1.0
    return {"n_higher_in_apoe2": n_chosen, "n_vdr_targets": k, "vdr_targets_higher_in_apoe2": hits,
            "neglog10_p": -math.log10(max(p, 1e-300))}


# --- promoter accessibility against expression ---------------------------------


def load_promoter_peaks(data_dir: Path) -> pd.DataFrame:
    p = pd.read_csv(data_dir / I.PROMOTER_PEAKS, na_values=["NA", ""])
    p["geneId"] = pd.to_numeric(p["geneId"], errors="coerce")
    return p


def concordance(degs: dict, peaks: pd.DataFrame, pooled: bool = False, agg: str = "median",
                key: str = "entrez") -> dict:
    """Pearson r between expression and median promoter accessibility log2 fold changes, per comparison.

    pooled=True takes each gene's accessibility summary over the rows of all
    four comparisons at once, as a filter on the comparison column that matches
    nothing would.
    """
    from scipy.stats import pearsonr

    out, sizes = {}, {}
    for comp, tag in (("E2vsE3", "e2_vs_e3"), ("E4vsE3", "e4_vs_e3"), ("E4vsE2", "e4_vs_e2"), ("KOvsE3", "ko_vs_e3")):
        d = degs[comp]
        sig = d[d["fdr"] < 0.05]
        p = peaks if pooled else peaks[peaks["comparison"] == comp]
        if key == "entrez":
            acc = p.dropna(subset=["geneId"]).groupby("geneId")["log2FoldChange"].agg(agg)
            expr = sig.set_index("entrezgene_id")["log2FoldChange"]
        else:
            acc = p.dropna(subset=["SYMBOL"]).groupby("SYMBOL")["log2FoldChange"].agg(agg)
            expr = sig.dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")["log2FoldChange"]
        shared = expr.index.intersection(acc.index)
        out[f"r_{tag}"] = float(pearsonr(expr.loc[shared], acc.loc[shared])[0])
        sizes[f"n_{tag}"] = int(len(shared))
    return {**out, **sizes}


def spearman_concordance(degs: dict, peaks: pd.DataFrame) -> dict:
    from scipy.stats import spearmanr

    out = {}
    for comp in SHEETS:
        d = degs[comp]
        expr = d[d["fdr"] < 0.05].set_index("entrezgene_id")["log2FoldChange"]
        acc = peaks[peaks["comparison"] == comp].dropna(subset=["geneId"]).groupby("geneId")["log2FoldChange"].median()
        shared = expr.index.intersection(acc.index)
        out[comp] = float(spearmanr(expr.loc[shared], acc.loc[shared])[0])
    return out
