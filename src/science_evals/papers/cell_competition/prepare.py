"""Load the public inputs once and write the files handed to agents.

Everything an agent receives is written by a function in this module (named
in DATA_FILES as derived:<function>) from the public files listed in
inputs.PUBLIC. The outputs are plain CSV with fixed row and column order and
fixed float formatting, so two builds give identical bytes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .inputs import cell_id, fetch_public, read_rdata, write_csv

TREATMENT = {"CI-Treated": "CI", "DMSO": "DMSO"}
CLUSTER_NAME = {1: "winner", 3: "intermediate", 4: "loser", 2: "extraembryonic_ectoderm", 5: "visceral_endoderm"}
EPIBLAST = ("winner", "intermediate", "loser")


@dataclass
class Data:
    data_dir: Path
    cells: pd.DataFrame  # 723 cells: treatment, batch, cluster label, losing score
    qc: pd.DataFrame  # agent-facing QC table (cc01)
    genes: pd.DataFrame  # gene_id, gene_name, chromosome, start, end
    hvg: list[str]
    tpm: pd.DataFrame | None = None  # genes x 723 cells (TPM)
    allele: dict = field(default_factory=dict)  # E-MTAB-8640 counts, CI epiblast cells x positions x 4
    allele_capped: dict = field(default_factory=dict)  # MitoHEAR re-count of the same cells (cross-check)


def load(data_dir: Path, with_tpm: bool = True, with_allele: bool = True) -> Data:
    figs = read_rdata(fetch_public(data_dir, "input_for_figures.Rdata"))
    meta = figs["meta_dati"]
    qc_start = figs["QC_start"]
    after = read_rdata(fetch_public(data_dir, "after_qc.Rda"))["after_qc"]
    meta_all = pd.read_csv(fetch_public(data_dir, "meta_all.txt"), sep=" ")

    # one row per QC-passing cell, in the authors' order (QC_start)
    names = [str(x) for x in qc_start.index]
    cells = pd.DataFrame(index=[cell_id(n) for n in names])
    cells.index.name = "cell_id"
    cells["treatment"] = [TREATMENT[str(meta.loc[n, "medium"])] for n in names]
    cells["batch"] = [int(meta.loc[n, "Batch"]) for n in names]
    ma = meta_all.set_index(meta_all["cell_name"].map(cell_id))
    cells["cluster"] = [CLUSTER_NAME[int(ma.loc[c, "cluster"])] for c in cells.index]
    aq = after.set_index(after["new_name"].astype(str))
    cells["losing_score"] = [float(aq.loc[c, "pseudo_time"]) for c in cells.index]
    cells["after_qc_cluster"] = [CLUSTER_NAME[int(aq.loc[c, "cluster"])] for c in cells.index]
    if not (cells["cluster"] == cells["after_qc_cluster"]).all():
        raise ValueError("cluster labels in meta_all.txt and after_qc.Rda disagree")
    cells = cells.drop(columns="after_qc_cluster")

    qc = pd.DataFrame({
        "cell_id": cells.index,
        "treatment": cells["treatment"].to_numpy(),
        "batch": cells["batch"].to_numpy(),
        "log10_total_reads": qc_start["totReadsLog10"].to_numpy(float),
        "frac_mapped_reads": qc_start["fracMappedReads"].to_numpy(float) / 100.0,
        "frac_reads_endogenous": qc_start["fracGenes"].to_numpy(float),
        "frac_reads_mito": qc_start["fracMTs"].to_numpy(float),
        "frac_reads_ercc": qc_start["fracSpikes"].to_numpy(float),
        "n_genes_above_10rpm": qc_start["numAbove10rpm"].to_numpy(int),
    })

    ann = pd.read_csv(fetch_public(data_dir, "useful_mt_names.txt"), sep="\t", dtype=str)
    ann.columns = ["gene_id", "gene_name", "start", "end", "chromosome"]
    genes = ann[["gene_id", "gene_name", "chromosome", "start", "end"]].copy()
    hvg = [str(g) for g in np.asarray(figs["HVGC"])]

    d = Data(data_dir=data_dir, cells=cells, qc=qc, genes=genes, hvg=hvg)
    if with_tpm:
        d.tpm = load_tpm(data_dir, list(cells.index))
    if with_allele:
        ci = cells[(cells["treatment"] == "CI") & cells["cluster"].isin(EPIBLAST)]
        d.allele = load_allele_ebi(data_dir, list(ci.index))
        d.allele_capped = load_allele(data_dir, list(ci.index))
    return d


def load_tpm(data_dir: Path, cell_order: list[str]) -> pd.DataFrame:
    obj = read_rdata(fetch_public(data_dir, "data_tpm.Rda"))
    tpm = next(iter(obj.values()))
    if not isinstance(tpm, pd.DataFrame):  # an R matrix arrives as an xarray DataArray
        tpm = pd.DataFrame(np.asarray(tpm.values, dtype=float),
                           index=[str(x) for x in tpm.coords[tpm.dims[0]].values],
                           columns=[str(x) for x in tpm.coords[tpm.dims[1]].values])
    tpm.columns = [cell_id(c) for c in tpm.columns]
    return tpm[cell_order].astype(float)


def load_allele(data_dir: Path, cell_order: list[str]) -> dict:
    obj = read_rdata(fetch_public(data_dir, "output_SNP_mt.Rda"))["output_SNP_mt"]
    arr = obj[0]
    names = [str(x) for x in arr.coords[arr.dims[0]].values]
    cols = [str(x) for x in arr.coords[arr.dims[1]].values]
    alleles = [c.split("_")[1] for c in cols[:4]]
    if alleles != ["A", "C", "G", "T"]:
        raise ValueError(f"unexpected allele order {alleles}")
    counts = np.asarray(arr.values, dtype=np.int64).reshape(len(names), -1, 4)
    positions = np.array([int(c.split("_")[0]) for c in cols[::4]])
    ref = np.array([c.split("_")[2] for c in cols[::4]])
    order = pd.Index(names).get_indexer(cell_order)
    if (order < 0).any():
        raise ValueError("cells missing from the allele table")
    return {"cells": list(cell_order), "positions": positions, "ref": ref, "counts": counts[order]}


# --- files handed to agents ----------------------------------------------------


def epiblast_cells(d: Data) -> pd.DataFrame:
    return d.cells[d.cells["cluster"].isin(EPIBLAST)]


def ci_epiblast_cells(d: Data) -> pd.DataFrame:
    e = epiblast_cells(d)
    return e[e["treatment"] == "CI"]


def write_qc_metrics(d: Data, path: Path) -> dict:
    return write_csv(d.qc, path)


def write_epiblast_cells(d: Data, path: Path) -> dict:
    e = epiblast_cells(d)
    df = pd.DataFrame({"cell_id": e.index, "treatment": e["treatment"].to_numpy(), "batch": e["batch"].to_numpy()})
    return write_csv(df, path)


def write_epiblast_tpm(d: Data, path: Path) -> dict:
    e = epiblast_cells(d)
    tpm = d.tpm[list(e.index)]
    df = tpm.reset_index().rename(columns={"index": "gene_id"})
    df.columns = ["gene_id"] + list(e.index)
    return write_csv(df, path, float_format="%.5g")


def write_genes(d: Data, path: Path) -> dict:
    g = d.genes.copy()
    g["authors_hvg"] = g["gene_id"].isin(set(d.hvg)).astype(int)
    return write_csv(g, path)


def write_ci_epiblast_cells(d: Data, path: Path) -> dict:
    e = ci_epiblast_cells(d)
    df = pd.DataFrame({
        "cell_id": e.index,
        "cluster": e["cluster"].to_numpy(),
        "batch": e["batch"].to_numpy(),
        "losing_score": e["losing_score"].to_numpy(float),
    })
    return write_csv(df, path, float_format="%.6f")


def write_mt_allele_counts(d: Data, path: Path) -> dict:
    e = ci_epiblast_cells(d)
    idx = pd.Index(d.allele["cells"]).get_indexer(e.index)
    counts = d.allele["counts"][idx]
    pos = d.allele["positions"]
    cols = [f"{p}_{a}" for p in pos for a in "ACGT"]
    df = pd.DataFrame(counts.reshape(len(idx), -1), columns=cols)
    df.insert(0, "cell_id", list(e.index))
    return write_csv(df, path)


def old_to_new_names(data_dir: Path) -> dict:
    """Old BAM-file cell names -> ArrayExpress cell ids, by the positional join in
    the authors' figure code (name_full.Rda order matches provo.Rda rows)."""
    provo = read_rdata(fetch_public(data_dir, "provo.Rda"))["provo"]
    full = [str(x) for x in np.asarray(read_rdata(fetch_public(data_dir, "name_full.Rda"))["nomi_full"])]
    pw = read_rdata(fetch_public(data_dir, "provo_ward_no_un.Rda"))["provo_ward_no_un"]
    if sorted(full) != sorted(str(x) for x in pw["Cellname"]) or len(full) != len(provo):
        raise ValueError("name tables do not line up")
    return {o: cell_id(n) for o, n in zip(full, provo.index)}


def load_allele_ebi(data_dir: Path, cell_order: list[str]) -> dict:
    """raw_counts_allele_mt.csv from E-MTAB-8640: one row per pileup file (named in
    old_2_cell_names.csv), columns '<pos> <ref>><alt>' with alt in G, T, C, A order."""
    files = pd.read_csv(fetch_public(data_dir, "old_2_cell_names.csv"))
    pile = [str(x) for x in files[files.columns[1]]]
    prefix, suffix = "pile_files/", ".unique.bam.MT.bam.pileup"
    old = [p[len(prefix):-len(suffix)] if p.startswith(prefix) and p.endswith(suffix) else p for p in pile]
    tab = pd.read_csv(fetch_public(data_dir, "raw_counts_allele_mt.csv"), index_col=0, dtype=np.int32)
    if len(tab) != len(old):
        raise ValueError(f"{len(tab)} rows in the allele table, {len(old)} names")
    cols = list(tab.columns)
    pos = np.array([int(c.split(" ")[0]) for c in cols[::4]])
    alts = [c.split(">")[1] for c in cols[:4]]
    ref = np.array([c.split(" ")[1].split(">")[0] for c in cols[::4]])
    counts = tab.to_numpy(dtype=np.int32).reshape(len(tab), -1, 4)
    counts = counts[:, :, [alts.index(a) for a in "ACGT"]]
    mapping = old_to_new_names(data_dir)
    new = [mapping.get(o) for o in old]
    index = {n: i for i, n in enumerate(new) if n is not None}
    missing = [c for c in cell_order if c not in index]
    if missing:
        raise ValueError(f"{len(missing)} cells missing from the EBI allele table")
    order = np.array([index[c] for c in cell_order])
    return {"cells": list(cell_order), "positions": pos, "ref": ref, "counts": counts[order],
            "source": "raw_counts_allele_mt.csv"}
