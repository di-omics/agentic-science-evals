# APOE microglia, Nat Commun 2025

Murphy KB, Hu D, Wolfs L, Rohde SK, Leguia Fauro G, Geric I, Mancuso R, De Strooper B,
Marzi SJ. The APOE isoforms differentially shape the transcriptomic and epigenomic
landscapes of human microglia xenografted into a mouse model of Alzheimer's disease.
Nat Commun 16, 4883 (2025). doi:10.1038/s41467-025-60099-4

**Author's role.** From the paper's author contributions: "D.H., L.W., S.R., G.L.F., and R.M.
performed experiments."

**Licenses.** Paper: Creative Commons Attribution 4.0 International (stated in the article).
Supplementary data on Zenodo (doi:10.5281/zenodo.15202660): CC BY 4.0. GEO records: NCBI
public data. Authors' repository (Marzi-Lab/APOE_microglia): GPL-3.0.

Keys, windows and controls were checked against the data by script; not reviewed by a person.

## The study in one paragraph

Isogenic human iPSC lines carrying one APOE2, APOE3 or APOE4 allele, or no APOE (KO), were
differentiated into microglial precursors and grafted into AppNL-G-F mice. After 12 months
the human microglia were sorted and profiled by bulk ATAC-seq and RNA-seq. APOE3 is the
baseline for most comparisons, and APOE4 is also compared directly with APOE2.

## Tasks

| Task | The decision | Reference | Trap |
|---|---|---|---|
| `apoe01_failed_knockouts` | Confirm that the knockout grafts lack APOE before using them (Fig. 1b, Supplementary Fig. 2) | KO-1 and KO-2 carry APOE at 55% and 48% of the APOE3 level and are left out; APOE-KO against APOE3 gives 338 DEGs (133 up) and an APOE log2 fold change of -4.08 | all five KO samples: 458 DEGs (231 up), APOE log2 fold change -2.13 |
| `apoe02_transcript_length_offsets` | Use tximport's per-sample transcript lengths with Salmon gene counts, as the authors' DESeqDataSetFromTximport does | APOE4 against APOE3 with the lengths as offsets: 270 DEGs (153 up), APOE log2 fold change -0.463 | plain counts: 1,220 DEGs (628 up), -0.726; 87% of the 993 extra genes move with the APOE4 libraries' length differences |
| `apoe03_vdr_targets_direction` | Read which way a stored comparison runs (Fig. 7b) | genes higher in APOE2 are the negative log2 fold changes of the APOE4 against APOE2 sheet: 890 genes, 82 VDR targets, -log10 p 10.44 | reversed sign: 735 genes, 30 VDR targets, none shared, -log10 p 0.20 |
| `apoe04_matched_contrast_concordance` | Take promoter accessibility from the same comparison as expression (Supplementary Fig. 6) | Pearson r 0.44, 0.64, 0.45 and 0.63 for E2vsE3, E4vsE3, E4vsE2 and KOvsE3 | accessibility pooled over all four comparisons: r 0.14, 0.61, 0.39 and 0.51 |

Prompts, keys, windows and notes: `docs/apoe_microglia/TASKS.md`. Every control and its
verdict under both graders: `results/apoe_microglia/CONTROLS.md`. Paper statements recomputed
from the public files: `results/apoe_microglia/PAPER_ANCHORS.md`.

Candidate decisions that were tried and dropped, because on these data they did not
separate a careless analysis from a careful one (checked by script during development, not
part of the build): ATAC normalisation (median of ratios, reads in peaks, total reads after
filtering and TMM move each comparison's region count by at most 9, E4 vs E2 63 to 80
against 72); keeping the two failed knockout grafts in the ATAC comparison (28 regions
against 22); and dropping the positive-marker filter in the microglial-state enrichment of
Fig. 4a. Counts of differentially accessible regions also differ by 14 to 28% between DESeq2
and pydeseq2, which keeps them out of the answer fields. All mice are male and every line
comes from one donor, so there is no sex or donor covariate to adjust for.

## Data

Public files (pinned by sha256 in `inputs.py`; about 139 MB in all; no reads):

| File | Source | Used for |
|---|---|---|
| `GSE271385_APOE_microglia_RNAseq_counts.csv.gz` | GEO GSE271385 supplementary file | apoe01, apoe02, anchors |
| `GSE271384_APOE_microglia_ATACseq_counts.csv.gz` | GEO GSE271384 supplementary file | anchors |
| `APOE_microglia_ATACseq_counts.csv` | authors' repository, commit 26a1dc4 | anchors |
| `APOE_microglia_ATACseq_promoter_peaks.csv` | authors' repository, commit 26a1dc4 | source of `atac_promoter_peaks.csv` |
| `APOE_microglia_RNAseq_DESeqDataset_obj.rds` | authors' repository, commit 26a1dc4 | source of `rnaseq_avg_tx_length.csv` |
| `VDR_target_genes_24h_monocytes_Warwick2021.csv` | authors' repository, commit 26a1dc4 | apoe03 |
| `Supplementary_Table2_DEGs.xlsx` | Zenodo 15202660 | apoe03, apoe04, source of `rnaseq_genes.csv` |

Derived files (written by `inputs.py`, pinned by sha256):

- `rnaseq_avg_tx_length.csv`: the avgTxLength assay of the authors' saved DESeqDataSet. Its
  counts assay is checked to equal the GEO matrix column for column (GEO renumbered the
  APOE3 samples 2-5 as 1-4), and the lengths take the GEO sample names.
- `rnaseq_samples.csv`: sample name, GEO accession, genotype, iPSC line and its APOE alleles
  (from the paper's Methods), host strain, age.
- `rnaseq_genes.csv`: Entrez ID to symbol, from the first sheet of Supplementary Data 2.
- `atac_promoter_peaks.csv`: the authors' promoter-peak table reduced to 12 columns. Some
  rows of the source have unquoted commas in the free-text gene name; that column is dropped.

## The DESeq2 engine

`deseq.py` repeats DESeq2's default steps for a one-factor design in numpy (size factors or
tximport normalisation factors, Cox-Reid gene-wise dispersions with DESeq2's line search,
parametric trend, MAP dispersions, IRLS with DESeq2's ridge, Wald tests, Cook's distances).
Checked against the data by script: on the authors' ATAC count matrix it returns the published
differentially accessible regions exactly (40, 50, 72 and 22, with the same up and down
splits) and, on the RNA-seq matrix with the tximport lengths, lands within 1% of the gene
counts in Supplementary Data 2; anchors 2 and 3 recompute both in every build. A one-off
development script also compared it peak by peak with the authors' Supplementary Data 4:
median relative difference of the standard errors below 1e-6.
One detail matters for low counts: DESeq2 floors the fitted means at 0.5 when it computes the
coefficient covariance; pydeseq2 does not, and calls 14 to 28% fewer regions on the ATAC
data. The RNA tasks run pydeseq2 as a variant control, and their windows accept it.

## Observations from the public files

Checked by script while building these tasks:

- The GEO ATAC matrix (228,041 peaks) and the repository's matrix (301,631 peaks) are two
  different consensus peak sets. The Methods count of 167,951 filtered peaks matches the GEO
  matrix; the published region counts match the repository matrix (173,284 peaks after the
  same filter).
- GEO holds only the RNA-seq counts. Without the tximport lengths, APOE4 against APOE3 gives
  about 4 times as many genes as the paper (task apoe02).

## Requirements

The repository's `analysis` extra (scipy, statsmodels, openpyxl and pydeseq2 are used here)
plus `rdata`, which reads the authors' saved DESeqDataSet (installed with
`uv pip install --python .venv/bin/python rdata`; version 1.1.0 was used). pydeseq2 results
are cached in `data/apoe_microglia/cache/`.

## Rebuild

```
.venv/bin/python -m science_evals.build apoe_microglia
.venv/bin/python -m pytest -q tests/test_apoe_microglia.py
```

The first build downloads the public files above into `data/apoe_microglia/` and writes the
derived files; a later build reuses them and the pydeseq2 cache.
