"""Lima, Lubatti et al. Nat Metab 2021: cell competition eliminates epiblast cells
with mitochondrial defects. Tasks built from E-MTAB-8640 and the authors'
processed inputs. See README.md in this folder."""
from __future__ import annotations

from pathlib import Path

from . import prepare
from .inputs import COMMIT, MITOHEAR_COMMIT, MITOHEAR_REPO, PUBLIC, REPO, public_entry

PAPER = {
    "short": "Cell competition, Nat Metab 2021",
    "citation": (
        "Lima A, Lubatti G, Burgstaller J, Hu D, et al. Cell competition acts as a purifying selection "
        "to eliminate cells with mitochondrial defects during early mouse development. Nat Metab 3, "
        "1091-1108 (2021). doi:10.1038/s42255-021-00422-7"
    ),
    "doi": "10.1038/s42255-021-00422-7",
    "author_role": (
        "D.H. performed embryo dissections, treatments and cell dissociation prior to scRNA-seq experiments."
    ),
    "data": (
        "ArrayExpress E-MTAB-8640 (single-cell RNA-seq, Smart-seq2) and E-MTAB-8692 (bulk RNA-seq). "
        "These tasks use E-MTAB-8640 processed data and the authors' processed inputs in their two "
        "repositories; no reads."
    ),
    "code": (
        f"{REPO} at commit {COMMIT} (figure code and inputs); {MITOHEAR_REPO} at commit "
        f"{MITOHEAR_COMMIT} (the same lab's heteroplasmy package, whose vignette reruns Fig. 6)"
    ),
    "license": (
        "Paper: author manuscript in Europe PMC (PMC7611553) under the publisher's terms for academic "
        "research use, not an open licence. Data: public ArrayExpress records under the EMBL-EBI terms of "
        "use. Code: the figure repository has no licence file; MitoHEAR is GPL-3.0."
    ),
}

# Files written by prepare.py from the public files above. Hashes pin the bytes;
# build() rewrites each file and stops if a hash changes.
DERIVED = {
    "cc_qc_metrics.csv": {
        "writer": "write_qc_metrics",
        "sha256": "a8e307deacbc4c289a91467a66d7a80f50fb25ebd3d2b7e1591c3686c946afb4",
        "bytes": 49476,
    },
    "cc_epiblast_cells.csv": {
        "writer": "write_epiblast_cells",
        "sha256": "7128f9bbaab7eae4968e8c4ba2f85ef64118db541f7ec4e8acccaa2391972e80",
        "bytes": 9082,
    },
    "cc_epiblast_tpm.csv": {
        "writer": "write_epiblast_tpm",
        "sha256": "d8704452c9a8b549555870731d7cc8f2cb2691ab00c13120da0997080992f334",
        "bytes": 60011363,
    },
    "cc_genes.csv": {
        "writer": "write_genes",
        "sha256": "188d3d197944358f1c76eed3ad0505f4a35fd86802ed202c1f581265d777f1ff",
        "bytes": 2669922,
    },
    "cc_ci_epiblast_cells.csv": {
        "writer": "write_ci_epiblast_cells",
        "sha256": "55b0788d773bbf4ff6b984640d832987bae1a5a96cf99a2b2cc7162338aa4772",
        "bytes": 8871,
    },
    "cc_mt_allele_counts.csv": {
        "writer": "write_mt_allele_counts",
        "sha256": "4db20c4e74d0bb83c6d6f530cb13986731c7180f481029e817405e91b616c22e",
        "bytes": 40442819,
    },
}


def _derived_entry(name: str) -> dict:
    info = DERIVED[name]
    return {"name": name, "source": f"derived:prepare.{info['writer']}", "sha256": info["sha256"],
            "bytes": info["bytes"]}


DATA_FILES = [public_entry(n) for n in PUBLIC] + [_derived_entry(n) for n in DERIVED]


def write_derived(d: prepare.Data, data_dir: Path, check: bool = True) -> dict:
    """Write every derived file and return its DATA_FILES entry, name -> entry."""
    out = {}
    for name, info in DERIVED.items():
        got = getattr(prepare, info["writer"])(d, data_dir / name)
        if check and info["sha256"] and got["sha256"] != info["sha256"]:
            raise ValueError(f"{name}: derived file changed (sha256 {got['sha256']}, pinned {info['sha256']})")
        out[name] = {"name": name, "source": f"derived:prepare.{info['writer']}", **got}
    return out


_CACHE: dict = {}


def _loaded(data_dir: Path):
    key = str(Path(data_dir).resolve())
    if key not in _CACHE:
        from . import tasks

        d = prepare.load(Path(data_dir))
        _CACHE.clear()
        _CACHE[key] = (d, tasks.context(d))
    return _CACHE[key]


def build(data_dir: Path):
    from . import tasks

    d, ctx = _loaded(data_dir)
    files = write_derived(d, Path(data_dir))
    return [
        tasks.make_cc01(d, files),
        tasks.make_cc02(d, files, ctx),
        tasks.make_cc03(d, files),
        tasks.make_cc04(d, files, ctx),
    ]


def anchors(data_dir: Path) -> list[dict]:
    from . import paper_anchors as A

    d, ctx = _loaded(data_dir)
    return A.anchors(d, traj=ctx["traj"], trend=ctx["trend"])
