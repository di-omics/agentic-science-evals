"""Public inputs for the AVE migration paper, and the small files derived from them.

Public files (fetched by URL and checked against their sha256):
    counts_e55.txt, counts_e625.txt   Salmon gene-level counts, 384 wells per stage (ArrayExpress)
    E-MTAB-9645.sdrf.txt              per-well annotation, including the authors' cluster labels
    norm_data_55.h5ad, norm_data_625.h5ad
                                      scran-normalized expression of the cells that passed QC
                                      (authors' repository, pinned commit)
    DC_55.csv                         the authors' first 2 diffusion components of the 5.5 dpc AVE
                                      and emVE cells (same repository; used for a control, not given
                                      to the agent)

Derived files (made here from public files; their sha256 pins the result):
    ena_read_counts_ERP124375.tsv     read count per sequencing run, from the ENA run report of
                                      the same accession, sorted by run so the file is stable
    wells_e55.tsv                     per-well plate design and read pairs, 5.5 dpc
    genes_grcm38.tsv                  Ensembl gene ID to symbol, for the rows of counts_e55.txt
    embryonic_ve_e55.txt              the 173 embryonic VE cells (AVE and emVE) at 5.5 dpc
    ephrin_family_expression.csv      22 ephrin and Eph genes in the labelled cells, both stages
    cell_labels.csv                   stage and cluster of every labelled cell, both stages
"""
from __future__ import annotations

import io
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from ...data import USER_AGENT, describe, fetch

EBI = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/645/E-MTAB-9645/Files/"
SCAVE_COMMIT = "6221759bfa0e103f4a037433185cb3384149a67e"
GITHUB = f"https://raw.githubusercontent.com/ScialdoneLab/scAVE/{SCAVE_COMMIT}/DATA/"
ENA_REPORT = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP124375"
              "&result=read_run&fields=run_accession,read_count&format=tsv")

COUNTS_E55 = "counts_e55.txt"
COUNTS_E625 = "counts_e625.txt"
SDRF = "E-MTAB-9645.sdrf.txt"
NORM_E55 = "norm_data_55.h5ad"
NORM_E625 = "norm_data_625.h5ad"
AUTHORS_DC_E55 = "DC_55.csv"  # the authors' diffusion coordinates of the 5.5 dpc AVE and emVE cells
ENA_COUNTS = "ena_read_counts_ERP124375.tsv"
WELLS_E55 = "wells_e55.tsv"
GENES = "genes_grcm38.tsv"
EMB_VE_E55 = "embryonic_ve_e55.txt"
EPHRIN = "ephrin_family_expression.csv"
LABELS = "cell_labels.csv"

PUBLIC = [
    {"name": COUNTS_E55, "source": EBI + COUNTS_E55,
     "sha256": "d6f234936da8c584b1139565ddd9d4d56e922b4104fe40e8ba179f63afdac6de", "bytes": 46123537},
    {"name": COUNTS_E625, "source": EBI + COUNTS_E625,
     "sha256": "578da6d1788d8554f3c4714ca038dce7ae15c288f6a65428280e1df3796283ed", "bytes": 45808002},
    {"name": SDRF, "source": EBI + SDRF,
     "sha256": "ba65f0c2eece22f20d35c19f594d47bdfc11bdc2de5a6095ad1316d844664aec", "bytes": 1121192},
    {"name": NORM_E55, "source": GITHUB + NORM_E55,
     "sha256": "ed98526f77daff92f72922c05dd4100e3c3807c535e3746a27e75dd1ec840c41", "bytes": 58126979},
    {"name": NORM_E625, "source": GITHUB + NORM_E625,
     "sha256": "9add07f74cbe863ede948c30bfc5e17b5c770128f4a60c4ce62041345c6c133d", "bytes": 54403350},
    {"name": AUTHORS_DC_E55, "source": GITHUB + AUTHORS_DC_E55,
     "sha256": "bf0b0542303dcee4f52ba10ab812b8488137b88e55c1c838d700811698f8dc9e", "bytes": 7190},
]

# sha256 and size of each derived file, recorded from a build; a rebuild must reproduce them.
DERIVED_PINS = {
    ENA_COUNTS: ("58f63ec818742e6ebcf85b1f25e76b46187ba975fb2a741c31b6eed765537158", 14241),
    WELLS_E55: ("39cf7ae09cdbc9de8c903ef79974821158c457d5a426cb8e14f9113ff1627462", 12311),
    GENES: ("28210d447b716485c792109227444b961691cb2782c4a19d42f47877d57bce4c", 1479722),
    EMB_VE_E55: ("0cb4a8cd421662d3a4d296843dd53161b9df7bb80e4c49c687c52525748a3f19", 2249),
    EPHRIN: ("d3950ddd09f1eb4000b28c9bf7c2ee62ed12a0b2a67e870d066517ffaa027428", 57915),
    LABELS: ("112c4fb625ab7225dc2b16e28604aec7cf5e670b702684bf905dba249f774720", 12360),
}

STAGES = {"e55": "5.5 dpc", "e625": "6.25 dpc"}
CLUSTER_NAMES = {
    "anterior visceral endoderm": "AVE",
    "emVE": "emVE",
    "exVE": "exVE",
    "epiblast": "Epi",
    "extraembryonic ectoderm": "ExE",
}
WELL_TYPES = {"OK": "single cell", "no cell": "no cell", "multiple cells": "multiple cells"}
FAMILY = re.compile(r"(Efn[ab]\d+|Eph[ab]\d+)")


def _public(name: str) -> dict:
    return next(e for e in PUBLIC if e["name"] == name)


def path(data_dir: Path, name: str) -> Path:
    """Path of a public file, fetched and checked if needed."""
    return fetch(_public(name), Path(data_dir))


# --- readers ------------------------------------------------------------------


def read_sdrf(data_dir: Path) -> pd.DataFrame:
    s = pd.read_csv(path(data_dir, SDRF), sep="\t", dtype=str)
    s = s.drop_duplicates("Source Name").copy()
    s["well_id"] = s["Factor Value[single cell identifier]"]
    s["stage"] = s["Characteristics[age]"].map({"embryonic day 5.5": "5.5 dpc", "embryonic day 6.25": "6.25 dpc"})
    s["cluster"] = s["Characteristics[cell type]"].map(CLUSTER_NAMES)
    s["well_type"] = s["Characteristics[single cell well quality]"].map(WELL_TYPES)
    s["authors_pass"] = s["Characteristics[post analysis well quality]"] == "pass"
    return s.set_index("well_id")


def read_counts(data_dir: Path, stage: str = "e55") -> pd.DataFrame:
    name = COUNTS_E55 if stage == "e55" else COUNTS_E625
    return pd.read_csv(path(data_dir, name), sep="\t", index_col=0)


def read_norm(data_dir: Path, stage: str = "e55"):
    import anndata as ad

    name = NORM_E55 if stage == "e55" else NORM_E625
    return ad.read_h5ad(path(data_dir, name))


def well_id_of(obs_name: str) -> str:
    """norm_data cell names use '.', the counts files and the SDRF use '_'."""
    return obs_name.replace(".", "_")


# --- derived files --------------------------------------------------------------


def _ena_read_counts(data_dir: Path) -> bytes:
    req = urllib.request.Request(ENA_REPORT, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as resp:
        raw = resp.read().decode()
    t = pd.read_csv(io.StringIO(raw), sep="\t", dtype={"run_accession": str, "read_count": "int64"})
    t = t[["run_accession", "read_count"]].sort_values("run_accession").reset_index(drop=True)
    return t.to_csv(sep="\t", index=False, lineterminator="\n").encode()


def _wells_e55(data_dir: Path) -> bytes:
    sdrf = read_sdrf(data_dir)
    ena = pd.read_csv(Path(data_dir) / ENA_COUNTS, sep="\t").set_index("run_accession")["read_count"]
    wells = pd.read_csv(path(data_dir, COUNTS_E55), sep="\t", index_col=0, nrows=1).columns
    meta = sdrf.loc[list(wells)]
    reads = ena.loc[meta["Comment[ENA_RUN]"].values].to_numpy()
    if (reads % 2).any():
        raise ValueError("ENA read counts are expected to count both mates of each pair")
    t = pd.DataFrame({"well_id": list(wells), "well_type": meta["well_type"].to_numpy(),
                      "read_pairs": reads // 2})
    if t["well_type"].isna().any():
        raise ValueError("unmapped well type in the SDRF")
    return t.to_csv(sep="\t", index=False, lineterminator="\n").encode()


def _genes(data_dir: Path) -> bytes:
    """Symbols for the Ensembl IDs of counts_e55.txt.

    The authors' normalized file names its genes by symbol (Ensembl ID where there is no
    symbol) in ascending Ensembl ID order. Pairing it with the sorted IDs of the count file
    gives the map; every gene that the normalized file names by its Ensembl ID must land on
    that same ID, which checks the pairing.
    """
    ids = pd.read_csv(path(data_dir, COUNTS_E55), sep="\t", index_col=0, usecols=[0]).index
    ids = [i for i in ids if not i.startswith("ERCC-")]
    norm = read_norm(data_dir, "e55")
    if len(ids) != norm.n_vars:
        raise ValueError("gene count mismatch between counts and normalized data")
    mapping = dict(zip(sorted(ids), norm.var_names))
    bad = [i for i, s in mapping.items() if s.startswith("ENSMUSG") and s != i]
    if bad:
        raise ValueError(f"symbol pairing fails for {len(bad)} genes")
    t = pd.DataFrame({"gene_id": ids, "symbol": [mapping[i] for i in ids]})
    return t.to_csv(sep="\t", index=False, lineterminator="\n").encode()


def _labelled(data_dir: Path, stage: str):
    """Normalized AnnData of one stage with the authors' cluster labels; unassigned cells dropped."""
    sdrf = read_sdrf(data_dir)
    a = read_norm(data_dir, stage)
    a.obs["cluster"] = sdrf.loc[[well_id_of(n) for n in a.obs_names], "cluster"].to_numpy()
    a.obs["stage"] = STAGES[stage]
    return a[a.obs["cluster"].notna()].copy()


def _embryonic_ve(data_dir: Path) -> bytes:
    a = _labelled(data_dir, "e55")
    cells = a.obs_names[a.obs["cluster"].isin(["AVE", "emVE"])]
    return ("\n".join(cells) + "\n").encode()


def _labels(data_dir: Path) -> bytes:
    rows = []
    for stage in ("e55", "e625"):
        a = _labelled(data_dir, stage)
        rows.append(pd.DataFrame({"cell_id": a.obs_names, "stage": a.obs["stage"].to_numpy(),
                                  "cluster": a.obs["cluster"].to_numpy()}))
    return pd.concat(rows).to_csv(index=False, lineterminator="\n").encode()


def _ephrin(data_dir: Path) -> bytes:
    parts = []
    for stage in ("e55", "e625"):
        a = _labelled(data_dir, stage)
        genes = sorted((g for g in a.var_names if FAMILY.fullmatch(g)), key=_family_key)
        x = np.asarray(a[:, genes].X, dtype=np.float64)
        t = pd.DataFrame(x, columns=genes)
        t.insert(0, "cell_type", a.obs["cluster"].to_numpy())
        t.insert(0, "stage", a.obs["stage"].to_numpy())
        t.insert(0, "cell_id", a.obs_names)
        parts.append(t)
    return pd.concat(parts).to_csv(index=False, float_format="%.6g", lineterminator="\n").encode()


def _family_key(g: str):
    m = re.fullmatch(r"(Efn|Eph)([ab])(\d+)", g)
    return (m.group(1), m.group(2), int(m.group(3)))


DERIVERS = {
    ENA_COUNTS: _ena_read_counts,
    WELLS_E55: _wells_e55,
    GENES: _genes,
    EMB_VE_E55: _embryonic_ve,
    EPHRIN: _ephrin,
    LABELS: _labels,
}
DERIVED_ORDER = [ENA_COUNTS, WELLS_E55, GENES, EMB_VE_E55, EPHRIN, LABELS]


def derive(data_dir: Path, name: str) -> Path:
    """Make a derived file if it is missing, then check it against its pin."""
    data_dir = Path(data_dir)
    out = data_dir / name
    if not out.exists():
        content = DERIVERS[name](data_dir)
        tmp = out.with_suffix(out.suffix + ".part")
        tmp.write_bytes(content)
        tmp.rename(out)
    want, size = DERIVED_PINS[name]
    if want:
        got = describe(out)
        if got["sha256"] != want or got["bytes"] != size:
            raise ValueError(f"{name}: derived file {got} does not match its pin {want} ({size} bytes)")
    return out


def prepare(data_dir: Path) -> None:
    """Fetch every public file and make every derived file."""
    for e in PUBLIC:
        fetch(e, Path(data_dir))
    for name in DERIVED_ORDER:
        derive(data_dir, name)


def entry(data_dir: Path, name: str) -> dict:
    """The DATA_FILES entry of one file."""
    for e in PUBLIC:
        if e["name"] == name:
            return dict(e)
    p = derive(data_dir, name)
    return {"name": name, "source": f"derived:{DERIVERS[name].__name__.lstrip('_')}", **describe(p)}


def derived_entries() -> list[dict]:
    """DATA_FILES entries of the derived files, from their pins."""
    return [{"name": n, "source": f"derived:{DERIVERS[n].__name__.lstrip('_')}",
             "sha256": DERIVED_PINS[n][0], "bytes": DERIVED_PINS[n][1]} for n in DERIVED_ORDER]
