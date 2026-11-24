# Cell competition eliminates epiblast cells with mitochondrial defects (Nat Metab 2021)

Lima A, Lubatti G, Burgstaller J, Hu D, et al. Cell competition acts as a purifying
selection to eliminate cells with mitochondrial defects during early mouse development.
Nat Metab 3, 1091-1108 (2021). doi:10.1038/s42255-021-00422-7

The text used here is the author manuscript in Europe PMC (PMC7611553), available under
the publisher's terms for academic research use; it is not an open licence.

## The author's role

From the paper's author contributions, verbatim: "D.H. performed embryo dissections,
treatments and cell dissociation prior to scRNA-seq experiments."

## The experiment

Mouse embryos were dissected at E5.5 and cultured overnight with a pan-caspase inhibitor
(CI) or vehicle (DMSO), so that epiblast cells that would normally die are kept and can
be profiled. 1,495 single cells were sequenced (Smart-seq2, ERCC spike-ins) and 723
passed quality control. Epiblast cells from CI-treated embryos form a trajectory from
"winner" cells, common in both treatments, to "loser" cells found only when apoptosis is
blocked; the diffusion pseudotime along it is the paper's "losing score". Genes falling
along that score are enriched for mitochondrial functions, and loser cells carry
heteroplasmic sequence changes in mt-Rnr1 and mt-Rnr2.

## Data

Public inputs, all pinned by sha256 in `inputs.py` (about 144 MB in all, no reads):

| File | Source | What it gives |
|---|---|---|
| `data_tpm.Rda` | ArrayExpress E-MTAB-8640, through the BioStudies file service | TPM for the 723 cells that passed QC (37,351 genes) |
| `raw_counts_allele_mt.csv` | same accession | A, C, G and T read counts at every mtDNA position for 524 cells, the table behind Fig. 6 (77 MB) |
| `input_for_figures.Rdata` | the authors' figure repository, [ScialdoneLab/Cell-Competition-Paper-Figures](https://github.com/ScialdoneLab/Cell-Competition-Paper-Figures) at commit 6f5ae9cc | per-cell QC metrics (723 cells), the 1,495-cell table, the 1,921 variable genes |
| `meta_all.txt` | same repository | the authors' clusters, treatment and batch per cell |
| `useful_mt_names.txt` | same repository | Ensembl gene annotation (53,700 genes, with mtDNA coordinates) |
| `old_2_cell_names.csv`, `name_full.Rda`, `provo.Rda`, `provo_ward_no_un.Rda` | same repository | the names that join the allele table's rows (older sequencing-file names) to cells, as the authors' figure code does |
| `after_qc.Rda` | the same lab's heteroplasmy package, [ScialdoneLab/MitoHEAR](https://github.com/ScialdoneLab/MitoHEAR) at commit dc38f39d | the authors' losing score per cell (the package's vignette reruns Fig. 6) |
| `output_SNP_mt.Rda` | same package | the package's own re-count of the allele table; used only to check the name join and in the Fig. 6h anchor |

The files handed to agents are written from these by `prepare.py` (plain CSV, fixed
order and formatting, hashes pinned in `__init__.py`):

| File | Rows x columns | Used by |
|---|---|---|
| `cc_qc_metrics.csv` | 723 cells x 9 | cc01 |
| `cc_epiblast_tpm.csv` | 37,351 genes x 512 epiblast cells (60 MB) | cc02, cc04 |
| `cc_epiblast_cells.csv` | 512 cells: treatment, batch | cc02 |
| `cc_genes.csv` | 53,700 genes: id, symbol, chromosome, start, end, `authors_hvg` | cc02, cc04 |
| `cc_ci_epiblast_cells.csv` | 274 CI epiblast cells: cluster, batch, losing score | cc03, cc04 |
| `cc_mt_allele_counts.csv` | 274 cells x 65,196 allele counts from the E-MTAB-8640 table (40 MB) | cc03 |

The join between the allele table and the cells was checked by script: the cells that
carry inherited non-reference alleles at mtDNA positions 12217 and 13536 are the same in
the accession's table and in MitoHEAR's independent re-count (333 of 333 calls agree).

Not used, and why:

- Raw reads (ENA ERP119309): never downloaded.
- E-MTAB-8692 (bulk RNA-seq of the BG and HB stem cell lines): the accession holds no
  processed files. The authors' repository has the count matrix (`GeneCounts.csv`, 16
  samples); a comparator task on it was tried and dropped, because the co-culture effect
  on respiratory chain genes is small (median log2 fold change about -0.03) next to the
  line difference (about -0.18), too close to method noise for a clean window.
- `counts.txt` (124 MB, all 1,495 cells) and `data_raw.Rda` (24 MB) from E-MTAB-8640: not
  needed by any task. The EBI file service delivered 10 to 40 KB/s during this build.
- The journal's supplementary tables: the publisher's article page answers scripted
  requests with a cookie check, so it was skipped.
- The ArrayExpress SDRF (per-cell metadata, including the authors' inferred cell type)
  was read to check cell counts; each task gives agents only the per-cell columns it
  needs, and no task gives the winner, intermediate and loser labels to an agent asked
  to orient the trajectory.

## Tasks

| Task | The decision | Reference | Trap, and where it lands |
|---|---|---|---|
| `cc01_cell_qc` | Which of the 723 first-pass cells to keep | Keep all 723 (367 CI, 356 DMSO), as the authors did | Pooled 3-MAD adaptive QC keeps 543 (241 CI): it removes 52 of the 71 cells of the CI-only cluster, which carry more ERCC reads and fewer detected genes. 106 CI cells outside a +/-20 window |
| `cc02_losing_score_direction` | Which end of the diffusion trajectory is its start | DPT from the tip next to DMSO cells: Samm50 -0.71, mt-Atp6 -0.64, Atf4 +0.65, Sesn2 +0.53 | Root at the CI-only end: every sign flips, at least 0.76 outside +/-0.3 windows |
| `cc03_rnr1_heteroplasmy` | Which cells a per-cell heteroplasmy fraction is based on | Cells covered by more than 50 reads at all six mt-Rnr1 positions: 6%, 54% and 87% of 215 cells, the paper's Fig. 6h | Cells without full coverage get a missing mean, which compares as False, and stay in the denominator: intermediate 0.40 and loser 0.68, 0.08 and 0.12 outside the windows, over 274 cells, 39 outside |
| `cc04_resp_chain_enrichment` | The background of an over-representation test | Tested genes (9,916): 41% decreasing, 95% of respiratory chain genes decreasing, fold enrichment 2.3, -log10 p 30 | Whole annotation as background: 7% decreasing, fold enrichment 10.6, -log10 p 81; outside every window |

Each task's notes give the controls, the windows and what the prompt leaves unsaid; the
controls are listed in `results/cell_competition/CONTROLS.md`.

## Paper anchors

`results/cell_competition/PAPER_ANCHORS.md` recomputes 6 statements from the public files,
and all 6 reproduce:

- 1,495 cells sequenced and 723 kept; all 723 meet the 4 stated thresholds.
- The CI-only cluster: 71 cells from CI-treated embryos, none from DMSO-treated ones.
- The trajectory: a Python rebuild of the authors' destiny diffusion map gives Spearman
  0.9998 with their published losing score.
- 5,278 genes change along it at FDR < 0.01 against 5,311 reported (a linear trend test,
  the parametric part of their GAM, instead of the GAM itself).
- Heteroplasmy filters leave 259 cells and 5,192 positions, as reported; Spearman's test
  finds the paper's 11 positions and 1 more (1862).
- Fig. 6h: 6%, 54% and 87% of 215 cells against about 5%, 55% and 87% of 215. MitoHEAR's
  re-count of the same cells gives 6%, 51% and 69%, from fewer reads per position (its
  depth histogram drops at 240 to 250 reads; the paper reports about 700). The tasks use
  the accession's table.

## Other notes

- The authors' figure repository has no licence file. MitoHEAR is GPL-3.0. Only data
  files are taken from either, and both are the authors' own deposits.
- The losing score given in cc03 and cc04 is the authors' published one. cc02 asks the
  agent to build it.
- A high mitochondrial read fraction does not mark the CI-only cells: their median is
  0.015 against 0.040 in winner cells. What adaptive QC catches is their ERCC fraction
  and detected-gene count.

## Requirements

On top of the repository's `analysis` extra (numpy, pandas, scipy, scikit-learn,
anndata, scanpy): `rdata` (pure-Python reader for .Rda and .Rdata files; 1.1.0 used),
installed with

```
uv pip install --python .venv/bin/python rdata
```

## Rebuild

```
.venv/bin/python -m science_evals.build cell_competition
```

The build fetches the public inputs into `data/cell_competition/` (git-ignored), writes
the derived files there and stops if any of their hashes change, runs every control
under both graders, and writes `evals/cell_competition/`, `results/cell_competition/`
and `docs/cell_competition/TASKS.md`. It takes about 90 seconds once the inputs are
local; two builds give identical files.

Keys, windows and controls were checked against the data by script; not reviewed by a
person.
