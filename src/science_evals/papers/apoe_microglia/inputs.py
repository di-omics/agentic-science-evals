"""Public files for the APOE microglia paper, and the files derived from them.

Every public file is pinned by URL, sha256 and size. Derived files are written
by the functions below from the pinned public files only; their sha256 pins the
bytes, and build() stops if a derived file comes out different.
"""
from __future__ import annotations

import hashlib
import io
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn"
REPO = "https://github.com/Marzi-Lab/APOE_microglia"
COMMIT = "26a1dc41eb0f489e98bc4ab2876fab53e575a08c"  # tag v1.0.0, archived as Zenodo 10.5281/zenodo.15176132
RAW = f"https://raw.githubusercontent.com/Marzi-Lab/APOE_microglia/{COMMIT}/data"
ZENODO = "https://zenodo.org/api/records/15202660/files"

RNA_COUNTS = "GSE271385_APOE_microglia_RNAseq_counts.csv.gz"
ATAC_COUNTS_GEO = "GSE271384_APOE_microglia_ATACseq_counts.csv.gz"
ATAC_COUNTS_REPO = "APOE_microglia_ATACseq_counts.csv"
PROMOTER_PEAKS_REPO = "APOE_microglia_ATACseq_promoter_peaks.csv"
RNA_RDS = "APOE_microglia_RNAseq_DESeqDataset_obj.rds"
VDR = "VDR_target_genes_24h_monocytes_Warwick2021.csv"
DEGS = "Supplementary_Table2_DEGs.xlsx"

LENGTHS = "rnaseq_avg_tx_length.csv"
SAMPLES = "rnaseq_samples.csv"
GENES = "rnaseq_genes.csv"
PROMOTER_PEAKS = "atac_promoter_peaks.csv"

PUBLIC = [
    {"name": RNA_COUNTS, "source": f"{GEO}/GSE271385/suppl/{RNA_COUNTS}",
     "sha256": "55daea740870257d41a10c1cb50c2f6333da033252485e1963596573229c3996", "bytes": 602839},
    {"name": ATAC_COUNTS_GEO, "source": f"{GEO}/GSE271384/suppl/{ATAC_COUNTS_GEO}",
     "sha256": "4c158d010927be365ac3fcc8886b1155f747ca00d3123f39500fa9383bce5562", "bytes": 7012528},
    {"name": ATAC_COUNTS_REPO, "source": f"{RAW}/{ATAC_COUNTS_REPO}",
     "sha256": "d716336b10bb94fe9acfbab1ca83bf9bc165ab0c3ce47f210f5cc0af48bc5d59", "bytes": 35883323},
    {"name": PROMOTER_PEAKS_REPO, "source": f"{RAW}/{PROMOTER_PEAKS_REPO}",
     "sha256": "bb4341916f58ec467f983d432fa597b1357bcfe5094f3d0e0ea2a8e7251c2c42", "bytes": 82949462},
    {"name": RNA_RDS, "source": f"{RAW}/{RNA_RDS}",
     "sha256": "1572d05d18edd9ea671b129ffdfdac48684e4e3edd39022fb89e0a11ac9bb492", "bytes": 4567960},
    {"name": VDR, "source": f"{RAW}/{VDR}",
     "sha256": "b499b1122466292cd362c5b1534fdf080176d538e81e1cabb1feefde1b6c5ec2", "bytes": 1552751},
    {"name": DEGS, "source": f"{ZENODO}/{DEGS}/content",
     "sha256": "3c225676ed29c5e9ad24db33991dfaf35417925dfb8026ad03169cdecc9fbd79", "bytes": 6198093},
]

# Pinned after the first build; build() refuses a derived file whose bytes change.
DERIVED = {
    LENGTHS: {"writer": "write_lengths",
              "sha256": "43cc606ad1190a0352952fe6ef2c0d4f720887e8c3cde09738daa2e7a0e07b2d", "bytes": 6083582},
    SAMPLES: {"writer": "write_samples",
              "sha256": "5b006f3f7ccb11c3f31195e6da5172cd8663a6585f2fd0c6ee40aac807286874", "bytes": 1591},
    GENES: {"writer": "write_genes",
            "sha256": "1ad0153ff11906a1267aaf0c77dbbd79b33033e062f74e2f675abd13b57818b3", "bytes": 220831},
    PROMOTER_PEAKS: {"writer": "write_promoter_peaks",
                     "sha256": "a60121ae7ab8652dd76f911b61b46136236c291286c79b59df126358f0916ce2", "bytes": 39522715},
}

# GEO sample records of GSE271385 (RNA-seq); iPSC lines and their APOE alleles from the paper's Methods.
LINES = {
    "E2": ("APOE2", "BIONi010-C-6", "e2/KO"),
    "E3": ("APOE3", "BIONi010-C-2", "e3/KO"),
    "E4": ("APOE4", "BIONi010-C-4", "e4/KO"),
    "KO": ("APOE-KO", "BIONi010-C-3", "KO/KO"),
}
RNA_GSM = {f"{g}-{i}": None for g in ("E2", "E3", "E4", "KO") for i in range(1, 6) if not (g == "E3" and i == 5)}
for _n, _s in enumerate(RNA_GSM):
    RNA_GSM[_s] = f"GSM{8375887 + _n}"
HOST = "Rag2-/- Il2rg-/- hCSF1KI AppNL-G-F"


def entry(name: str) -> dict:
    for f in PUBLIC:
        if f["name"] == name:
            return dict(f)
    info = DERIVED[name]
    return {"name": name, "source": f"derived:inputs.{info['writer']}", "sha256": info["sha256"],
            "bytes": info["bytes"]}


def _describe(path: Path) -> dict:
    data = path.read_bytes()
    return {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}


def fetch_public(data_dir: Path) -> None:
    from ...data import fetch

    for f in PUBLIC:
        fetch(f, data_dir)


# --- derived files -----------------------------------------------------------


def _rds_assays(path: Path) -> dict:
    """The assay matrices of a serialised DESeqDataSet, by name (needs the rdata package)."""
    import sys

    import numpy as np
    import rdata

    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))
    parsed = rdata.parser.parse_file(path)

    def resolve(o):
        while o is not None and getattr(o, "info", None) is not None and o.info.type.name == "REF":
            o = o.referenced_object
        return o

    def tag(t):
        t = resolve(t)
        if t is None:
            return None
        v = t.value
        v = resolve(v) if hasattr(v, "info") else v
        v = getattr(v, "value", v)
        return v.decode() if isinstance(v, bytes) else None

    def pairs(o):
        out, o = [], resolve(o)
        while o is not None and getattr(o, "info", None) is not None and o.info.type.name == "LIST":
            out.append((tag(o.tag), resolve(o.value[0])))
            o = resolve(o.value[1])
        return out

    def slot(o, name):
        for k, v in pairs(o.attributes):
            if k == name:
                return v
        raise KeyError(name)

    assays = slot(slot(slot(resolve(parsed.object), "assays"), "data"), "listData")
    names = [resolve(x).value.decode() for x in slot(assays, "names").value]
    out = {}
    for name, obj in zip(names, assays.value):
        obj = resolve(obj)
        dim = tuple(int(x) for x in np.asarray(slot(obj, "dim").value))
        out[name] = np.asarray(obj.value).reshape(dim, order="F")
    return out


def write_lengths(data_dir: Path, out: Path) -> dict:
    """tximport's avgTxLength from the authors' DESeqDataSet, labelled like the GEO counts.

    The object's counts are checked to be the GEO count matrix, column for column
    (GEO renumbered the APOE3 samples 2-5 as 1-4), so the lengths line up with it.
    """
    import numpy as np
    import pandas as pd

    counts = pd.read_csv(data_dir / RNA_COUNTS, index_col=0)
    assays = _rds_assays(data_dir / RNA_RDS)
    if not np.array_equal(assays["counts"], counts.to_numpy()):
        raise ValueError("the DESeqDataSet counts are not the GEO count matrix")
    lengths = pd.DataFrame(assays["avgTxLength"], index=counts.index, columns=counts.columns)
    lengths.index.name = "entrez_id"
    lengths.to_csv(out, float_format="%.4f", lineterminator="\n")
    return _describe(out)


def write_samples(data_dir: Path, out: Path) -> dict:
    import pandas as pd

    counts = pd.read_csv(data_dir / RNA_COUNTS, index_col=0, nrows=1)
    rows = []
    for s in counts.columns:
        grp = s.split("-")[0]
        genotype, line, alleles = LINES[grp]
        rows.append({"sample": s, "geo_accession": RNA_GSM[s], "genotype": genotype, "ipsc_line": line,
                     "line_apoe_alleles": alleles, "host_mouse": HOST, "age_months": 12})
    pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")
    return _describe(out)


def write_genes(data_dir: Path, out: Path) -> dict:
    import pandas as pd

    t = pd.read_excel(data_dir / DEGS, sheet_name="E2vsE3_DEGs", usecols=["entrezgene_id", "hgnc_symbol"])
    t = t.dropna().drop_duplicates("entrezgene_id").sort_values("entrezgene_id", kind="stable")
    t.columns = ["entrez_id", "symbol"]
    t.to_csv(out, index=False, lineterminator="\n")
    return _describe(out)


PROMOTER_COLUMNS = ["comparison", "Peak", "Chr", "Start", "End", "log2FoldChange", "pvalue", "fdr",
                    "annotation", "distanceToTSS", "geneId", "SYMBOL"]


def write_promoter_peaks(data_dir: Path, out: Path) -> dict:
    """The authors' promoter-peak table with the columns the task needs.

    Some GENENAME values in the source contain unquoted commas, which shifts
    their rows; GENENAME (the last free-text field before `comparison`) is
    dropped, so those rows parse. Values are copied as text, unchanged.
    """
    import csv

    src = data_dir / PROMOTER_PEAKS_REPO
    with open(src, newline="") as fh:
        header = fh.readline().rstrip("\n").split(",")
        n_fixed = header.index("GENENAME")
        keep = [header.index(c) if c != "comparison" else None for c in PROMOTER_COLUMNS]
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(PROMOTER_COLUMNS)
        for line in fh:
            f = line.rstrip("\n").split(",")
            fixed, comparison = f[:n_fixed], f[-1]
            w.writerow([comparison if i is None else fixed[i] for i in keep])
    out.write_text(buf.getvalue())
    return _describe(out)


def write_derived(data_dir: Path, check: bool = True) -> dict:
    """Write every derived file; return name -> DATA_FILES entry."""
    from . import inputs as me

    out = {}
    for name, info in DERIVED.items():
        got = getattr(me, info["writer"])(data_dir, data_dir / name)
        if check and info["sha256"] and got["sha256"] != info["sha256"]:
            raise ValueError(f"{name}: derived file changed (sha256 {got['sha256']}, pinned {info['sha256']})")
        out[name] = {"name": name, "source": f"derived:inputs.{info['writer']}", **got}
    return out
