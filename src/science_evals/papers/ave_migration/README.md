# AVE migration, Dev Cell 2024

Thowfeequ S, Fiorentino J, Hu D, Solovey M, Ruane S, Whitehead M, Zhou F, Godwin J,
Mateo-Otero Y, Vanhaesebroeck B, Scialdone A, Srinivas S. An integrated approach identifies the
molecular underpinnings of murine anterior visceral endoderm migration. Developmental Cell 59,
2347-2363.e9 (2024). doi:10.1016/j.devcel.2024.05.014. Open access, CC BY 4.0.

**The repository owner's part**, from the paper's author contributions: "D.H. and S.R. did
sample collections for the transcriptomics and proteomics, respectively". D.H. (Di Hu) did the
sample collections for the transcriptomics, the single-cell RNA-seq that these tasks use.

**The paper in one line.** Smart-seq2 profiles of visceral endoderm (VE) cells from 5.5 and 6.25
dpc mouse embryos, with phosphoproteomics, imaging and lineage labelling, describe how the
anterior visceral endoderm (AVE) emerges, migrates and returns to an emVE-like state, and point
to semaphorin 6D signalling in its migration.

Every key, window and control below was checked against the data by script; none has been
reviewed by a person.

## Tasks

| Task | The decision | Reference finds | Trap, and where it lands |
|---|---|---|---|
| `ave01_plate_qc` | Which wells of a Smart-seq2 plate (384 wells, ERCC spike-ins) hold one good cell | 255 wells, the authors' set well for well | Droplet defaults (200 genes, mitochondrial cap) keep 378 wells: 105 nearly empty wells, mostly spike-in and unmapped reads, get through. Window 235 to 275 |
| `ave02_pseudotime_root` | Where to root the diffusion pseudotime so that it starts in the AVE | Cer1 correlation with pseudotime -0.79; 37.7% of 1,064 changing genes high-in-AVE | Root at an unchecked diffusion-component extreme or at the first cell (both emVE): correlation +0.81 or +0.72, 60 to 62% called high-in-AVE, every panel gene in the wrong group |
| `ave03_ephrin_bootstrap_signs` | Align each bootstrap PCA axis with the full-data axis before comparing loading signs | 14 and 13 genes kept; 6 and 8 significant, as in the paper | Signs compared as the solver returns them: no gene significant at 5.5 dpc (median p 0.50) |
| `ave04_dkk1_timing` | Whether Dkk1 is an AVE signal while the AVE migrates | Dkk1 AUC (AVE vs emVE) 0.52 at 5.5 dpc and 0.78 at 6.25 dpc; keyed option: the data argue against AVE-derived DKK1 starting the migration | Prior-driven reading of the 6.25 dpc (or pooled, AUC 0.62) enrichment as support for DKK1 guidance |

`results/ave_migration/CONTROLS.md` lists every scripted control (30), its answer and both
graders' verdicts. `docs/ave_migration/TASKS.md` has the full prompts, graders and notes.

## Paper anchors

Recomputed from the public files (`results/ave_migration/PAPER_ANCHORS.md`):

- The paper's five QC thresholds, with read pairs from the ENA run report, keep 263 wells at
  5.5 dpc and 238 at 6.25 dpc; setting aside the 8 wells annotated as no cell or multiple cells
  leaves 255, and both stages then match the authors' pass or fail in the SDRF for all 384 wells.
- pam (k = 2) on the Spearman distance splits the 5.5 dpc embryonic VE into the authors' 95 AVE
  and 78 emVE cells exactly.
- Genes changing along pseudotime: 1,064 and 789 with a natural spline F-test, against 952 and
  915 with the paper's R gam loess fit (within 15%; the test differs, so exact counts are not
  expected).
- Ephrin/Eph genes: 14 and 13 kept, 6 and 8 significant, as in the paper.
- Dkk1 enriched in the AVE only at 6.25 dpc; Sema6d highest in the AVE and Plxna1 higher in emVE
  and Epi than in the AVE at 5.5 dpc, as in the paper.

## Data

Public files (pinned by sha256 in `inputs.py`; about 205 MB in all):

| File | Source | Bytes | Used by |
|---|---|---|---|
| `counts_e55.txt` | ArrayExpress E-MTAB-9645, processed data | 46,123,537 | ave01, derived files |
| `counts_e625.txt` | ArrayExpress E-MTAB-9645, processed data | 45,808,002 | anchors |
| `E-MTAB-9645.sdrf.txt` | ArrayExpress E-MTAB-9645 | 1,121,192 | derived files |
| `norm_data_55.h5ad` | github.com/ScialdoneLab/scAVE, commit 6221759 | 58,126,979 | ave02, ave04 |
| `norm_data_625.h5ad` | same repository and commit | 54,403,350 | ave04 |
| `DC_55.csv` | same repository and commit | 7,190 | a control of ave02 (the authors' root cell) |

Derived files (made by `inputs.py`, pinned by sha256):

- `ena_read_counts_ERP124375.tsv`: read count per run from the ENA run report of the same
  accession (metadata only, no reads), sorted by run so the file is stable. Read pairs per well
  are half of it; this matches the authors' total reads for every cell in their normalized files.
- `wells_e55.tsv`: per-well annotation (`single cell`, `no cell`, `multiple cells`, from the
  SDRF's single cell well quality column) and read pairs.
- `genes_grcm38.tsv`: Ensembl ID to symbol. The authors' normalized file lists symbols in
  ascending Ensembl ID order; the pairing is checked against the 726 genes it names by ID.
- `embryonic_ve_e55.txt`: the 173 cells of the authors' AVE and emVE clusters at 5.5 dpc, not
  labelled.
- `ephrin_family_expression.csv`: the 22 Efn and Eph genes in the authors' labelled cells.
- `cell_labels.csv`: stage and cluster (AVE, emVE, exVE, Epi, ExE) of the labelled cells.

No raw reads were downloaded.

## Not built, and why

- **Ligand-receptor screen (Fig. 6B).** It starts from a curated list of 2,548 pairs from
  another publication; no copy at a pinned public URL was found, so the 413, 471 and 43 pair
  counts could not be recomputed.
- **AVE sub-clusters in the reused 10x dataset (Fig. 5).** The paper reuses the 10x data of its
  reference 17, a 2019 atlas of mouse gut endoderm (Human Cell Atlas project
  4e6f083b-5b9a-4393-9890-2a83da8188f1). Its
  contributor matrix `sc_endoderm_all_cells.h5ad` is 4,059,291,777 bytes, over this paper's
  3 GB download budget, so it was not downloaded.
- **Proteomics and phosphoproteomics.** No deposit was found in PRIDE; the peptide tables are
  in the journal's supplement (Tables S14 and S15), which was not fetched.
- **RNA velocity.** The loom files are in the code repository, but a velocity readout stable
  enough to grade was not attempted.
- **AVE versus emVE at 5.5 dpc as a task.** The split is a continuum: reasonable clustering
  choices recover the authors' AVE cells with Jaccard indices from about 0.3 to 1.0, too
  spread for a window that separates them from a careless analysis. It is kept as an anchor.

## Rebuild

```
.venv/bin/python -m science_evals.build ave_migration
.venv/bin/python -m pytest -q tests/test_ave_migration.py
```

The first build downloads the public files (EBI and GitHub) and the ENA run report; later builds
reuse `data/ave_migration/`. Two builds give identical files (fixed seeds; the ENA report is
sorted before it is pinned).

## Requirements

Only the repository's `analysis` and `grade` extras (numpy, pandas, scipy, scikit-learn, anndata,
scanpy, latch-eval-tools). No other package was installed. Network access is needed for the first
build.
