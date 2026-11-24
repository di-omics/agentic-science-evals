"""Public inputs for the cell competition paper, and the files derived from them.

Public files come from ArrayExpress E-MTAB-8640 (through the BioStudies file
service) and from the authors' figure repository at a fixed commit. The files
handed to agents are derived from them here, deterministically, and written as
plain CSV (no compression, so the bytes do not depend on the zlib build).

The .Rda / .Rdata files are read with the pure-Python `rdata` package, imported
only when a derivation needs it.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ...data import fetch, sha256_file

COMMIT = "6f5ae9cc7a525c1f8b58e1a3cb2947df4e47aef7"
REPO = "https://github.com/ScialdoneLab/Cell-Competition-Paper-Figures"
GITHUB_RAW = (
    "https://raw.githubusercontent.com/ScialdoneLab/Cell-Competition-Paper-Figures/"
    f"{COMMIT}/Input_data_figures/"
)
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-8640/"
# The same lab's heteroplasmy package; its vignette reruns this paper's Fig. 6
# analysis and ships the per-cell losing score and the allele counts it used.
MITOHEAR_COMMIT = "dc38f39d518ae6b7b15eb1067fc77e06819c27cd"
MITOHEAR_REPO = "https://github.com/ScialdoneLab/MitoHEAR"
MITOHEAR_RAW = (
    f"https://raw.githubusercontent.com/ScialdoneLab/MitoHEAR/{MITOHEAR_COMMIT}/inst/extdata/"
)

# Public files, pinned by sha256. Sizes and hashes were measured on the files
# as served on 2026-11-24.
PUBLIC = {
    "data_tpm.Rda": {
        "source": BIOSTUDIES + "data_tpm.Rda",
        "sha256": "719bb38160f3b720d505a000689780f35c671a9a6a6e8bc9c414bb10d10ebc77",
        "bytes": 45696709,
    },
    "raw_counts_allele_mt.csv": {
        "source": BIOSTUDIES + "raw_counts_allele_mt.csv",
        "sha256": "61eca63920a399ed4f0dbf3ec60e68a3c610f6ddda8336c8cb38863f6c12ea88",
        "bytes": 77475975,
    },
    "input_for_figures.Rdata": {
        "source": GITHUB_RAW + "input_for_figures.Rdata",
        "sha256": "a3e80a6a17491b58b3d4d07c2a186a6514717b32c9e11607077a267048fd3431",
        "bytes": 692611,
    },
    "meta_all.txt": {
        "source": GITHUB_RAW + "meta_all.txt",
        "sha256": "3c09e76e390a01edb35863c41ea33bc3a9c4cfdaabd148558be45f87765bb441",
        "bytes": 37625,
    },
    "old_2_cell_names.csv": {
        "source": GITHUB_RAW + "old_2_cell_names.csv",
        "sha256": "a588ce921d23fa4e44ee36471cb1bcc1ec70c0c8169bab2360c3a08e86fba49e",
        "bytes": 27033,
    },
    "name_full.Rda": {
        "source": GITHUB_RAW + "name_full.Rda",
        "sha256": "ddc9e6b4a6aca251c6fd233f3368db286347183486b90fa303f031423bda7a16",
        "bytes": 1159,
    },
    "provo.Rda": {
        "source": GITHUB_RAW + "provo.Rda",
        "sha256": "b9ecf5b8793371a5b024ed098889b59d564817ff296b58434d083ba438903eae",
        "bytes": 2259,
    },
    "provo_ward_no_un.Rda": {
        "source": GITHUB_RAW + "provo_ward_no_un.Rda",
        "sha256": "4a1e0fdbacd1e844adf4f10db30b163453784b1ae2719099cf37e62784b96669",
        "bytes": 3295,
    },
    "useful_mt_names.txt": {
        "source": GITHUB_RAW + "useful_mt_names.txt",
        "sha256": "fe7806473e69f9a1b8c93e0ef10f20c8e0d6f799656cec956203fee93f37de13",
        "bytes": 2562551,
    },
    "after_qc.Rda": {
        "source": MITOHEAR_RAW + "after_qc.Rda",
        "sha256": "12da5b80fb60ea326554b544573cd56510a8634d1a314e09c5105703d88cf6bf",
        "bytes": 13292,
    },
    "output_SNP_mt.Rda": {
        "source": MITOHEAR_RAW + "output_SNP_mt.Rda",
        "sha256": "986ec19695fd16477eb51fb9c8ea6c2ac7c37638db56beb7f09514d3716dbf7e",
        "bytes": 17536734,
    },
}


def public_entry(name: str) -> dict:
    return {"name": name, **PUBLIC[name]}


def fetch_public(data_dir: Path, name: str) -> Path:
    return fetch(public_entry(name), data_dir)


# --- readers -----------------------------------------------------------------


def read_rdata(path: Path) -> dict:
    import rdata  # optional dependency, see README

    return rdata.conversion.convert(rdata.parser.parse_file(str(path)))


def cell_id(name: str) -> str:
    """The ArrayExpress single cell identifier: '24538_8#6' -> '24538_8_6'."""
    return str(name).replace("#", "_")


# --- deterministic writer ----------------------------------------------------


def write_csv(df: pd.DataFrame, path: Path, float_format: str = "%.6g") -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = df.to_csv(index=False, float_format=float_format, lineterminator="\n")
    path.write_bytes(text.encode("ascii"))
    return {"sha256": sha256_file(path), "bytes": path.stat().st_size}

