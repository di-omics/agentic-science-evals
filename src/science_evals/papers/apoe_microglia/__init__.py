"""Murphy, Hu et al. Nat Commun 2025: APOE isoforms shape the transcriptome and
chromatin of human microglia xenografted into an Alzheimer's disease model.
Tasks built from GEO GSE271384/GSE271385, the authors' Zenodo supplementary data
and their code repository. See README.md in this folder."""
from __future__ import annotations

from pathlib import Path

from . import inputs as I

PAPER = {
    "short": "APOE microglia, Nat Commun 2025",
    "citation": (
        "Murphy KB, Hu D, Wolfs L, Rohde SK, Leguia Fauro G, Geric I, Mancuso R, De Strooper B, Marzi SJ. "
        "The APOE isoforms differentially shape the transcriptomic and epigenomic landscapes of human microglia "
        "xenografted into a mouse model of Alzheimer's disease. Nat Commun 16, 4883 (2025). "
        "doi:10.1038/s41467-025-60099-4"
    ),
    "doi": "10.1038/s41467-025-60099-4",
    "author_role": "D.H., L.W., S.R., G.L.F., and R.M. performed experiments.",
    "data": (
        "GEO GSE271384 (ATAC-seq) and GSE271385 (RNA-seq), processed count matrices; Supplementary Data 2 on Zenodo "
        "(doi:10.5281/zenodo.15202660); the count matrices, promoter-peak table, saved DESeqDataSet and VDR target "
        "list in the authors' repository. No reads."
    ),
    "code": f"{I.REPO} at commit {I.COMMIT} (tag v1.0.0, archived as doi:10.5281/zenodo.15176132)",
    "license": (
        "Paper: Creative Commons Attribution 4.0 International. Supplementary data on Zenodo: CC BY 4.0. GEO "
        "records: NCBI public data. Authors' repository: GPL-3.0."
    ),
}

DATA_FILES = [I.entry(f["name"]) for f in I.PUBLIC] + [I.entry(n) for n in I.DERIVED]


def build(data_dir: Path):
    from . import reference as R
    from . import tasks as T

    data_dir = Path(data_dir)
    I.fetch_public(data_dir)
    files = {f["name"]: dict(f) for f in I.PUBLIC}
    files.update(I.write_derived(data_dir))
    cache = data_dir / "cache"
    rna = R.load_rna(data_dir)
    degs = R.load_degs(data_dir)
    return [
        T.t1_failed_knockouts(rna, files, cache),
        T.t2_length_offsets(rna, files, cache),
        T.t3_vdr_direction(degs, R.load_vdr(data_dir), files),
        T.t4_matched_concordance(degs, R.load_promoter_peaks(data_dir), files),
    ]


def anchors(data_dir: Path) -> list[dict]:
    from .paper_anchors import anchors as compute

    return compute(Path(data_dir))
