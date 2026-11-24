# Tasks: Cell competition, Nat Metab 2021

Lima A, Lubatti G, Burgstaller J, Hu D, et al. Cell competition acts as a purifying selection to eliminate cells with mitochondrial defects during early mouse development. Nat Metab 3, 1091-1108 (2021). doi:10.1038/s42255-021-00422-7

## `cc01_cell_qc`: Which cells to keep after first-pass QC when the treatment keeps dying cells alive

From the paper: Methods, data processing, quality control and normalization; Extended Data Fig. 1a

```
Mouse embryos were dissected at embryonic day 5.5 and cultured for 16 hours with a pan-caspase inhibitor (CI) or with vehicle (DMSO). Blocking caspases keeps alive the epiblast cells that would otherwise be eliminated by apoptosis during this window, so that they can be profiled. The embryos were dissociated, pooled by treatment within each of 5 processing batches, and single cells were sequenced with Smart-seq2 with ERCC spike-ins in every well.

Of 1,495 sequenced cells, 723 passed the lab's first-pass filter: log10 total reads of at least 4, at least 80% of reads mapped, at least 3,000 genes above 10 reads per million, and at least half of all reads on endogenous genes. `cc_qc_metrics.csv` lists these 723 cells, one row per cell: the treatment of the embryos the cell came from (`treatment`, CI or DMSO), the processing batch (`batch`, 1 to 5), log10 total reads, the fraction of reads mapped, the fractions of reads on endogenous genes, on mitochondrial genes and on ERCC spike-ins, and the number of genes above 10 reads per million.

The next steps are to cluster these cells and to characterise the epiblast cells that survive only when apoptosis is blocked. Decide which of the 723 cells should go into that analysis, applying any further quality filtering you judge appropriate. Report how many cells you keep, and how many of the kept cells come from CI-treated and from DMSO-treated embryos.

Return EXACTLY:

<EVAL_ANSWER>
{"n_cells_kept": <integer>, "n_ci_cells_kept": <integer>, "n_dmso_cells_kept": <integer>}
</EVAL_ANSWER>

`n_cells_kept`: number of cells you carry into clustering
`n_ci_cells_kept`: how many of the kept cells come from CI-treated embryos
`n_dmso_cells_kept`: how many of the kept cells come from DMSO-treated embryos
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "n_cells_kept": 723,
   "n_ci_cells_kept": 367,
   "n_dmso_cells_kept": 356
  },
  "tolerances": {
   "n_cells_kept": {
    "type": "absolute",
    "value": 45
   },
   "n_ci_cells_kept": {
    "type": "absolute",
    "value": 20
   },
   "n_dmso_cells_kept": {
    "type": "absolute",
    "value": 35
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent checks what a standard adaptive QC rule would remove before applying it, in an experiment built to keep stressed, dying cells alive, and whether it notices that the outliers are concentrated in one treatment.

2) Primary trap: Adaptive outlier rules at the scater default of 3 MADs are applied to all cells pooled. They flag 180 cells, 126 of them from CI embryos, and remove 52 of the 71 cells that the authors' clustering later places in the cluster found only in CI embryos, and 56 of the 125 CI cells of the intermediate cluster. Those cells passed the fixed thresholds and are the population the experiment was built to capture; they carry more ERCC reads (median 0.060 against 0.019 in winner epiblast cells) and fewer genes above 10 reads per million (4872 against 6099), that is, less endogenous RNA per cell. A high mitochondrial fraction is not their signature (median 0.015 against 0.040), so a mitochondrial filter does not remove them. At 5 MADs on pooled cells the rules still remove 30 of the 71.

3) Scripted controls, run on the task data:
   reference authors_fixed_thresholds_only: Keep every cell that passed the fixed thresholds, as the authors did. -> n_cells_kept=723, n_ci_cells_kept=367, n_dmso_cells_kept=356 (passes)
   variant high_mito_outliers_removed: Also remove cells more than 3 MADs above the median mitochondrial fraction. -> n_cells_kept=705, n_ci_cells_kept=355, n_dmso_cells_kept=350 (passes)
   variant fixed_mito_cut_10pct: Also remove cells with more than 10% of reads on mitochondrial genes. -> n_cells_kept=711, n_ci_cells_kept=359, n_dmso_cells_kept=352 (passes)
   variant adaptive_qc_within_treatment_5mad: Adaptive outlier rules (5 MADs) computed separately within each treatment group. -> n_cells_kept=686, n_ci_cells_kept=359, n_dmso_cells_kept=327 (passes)
   trap adaptive_qc_pooled_3mad: scater-style outlier rules at the default 3 MADs (library size, detected genes, ERCC and mitochondrial fractions) on all cells pooled. -> n_cells_kept=543, n_ci_cells_kept=241, n_dmso_cells_kept=302 (fails)

4) Tolerance: The windows accept keeping all 723 cells, removing high-mitochondrial outliers (3 MADs, or a fixed 10% cut), and adaptive rules computed within each treatment group at 5 MADs. They do not accept adaptive rules that pool the treatments: at 5 MADs on pooled cells 303 CI cells are kept, within each processing batch at 3 MADs 248, and within each treatment at 3 MADs 308, against 367 in the key.

5) Not said in the prompt: The prompt does not say that the QC metrics differ between treatments, or that the cells of interest hold less endogenous RNA.

6) Provenance: The authors kept all 723 cells that passed the fixed thresholds (Methods). The QC table is the authors' QC_start object, with treatment and batch from meta_dati, both in input_for_figures.Rdata of their figure repository. Cluster labels used in these notes come from meta_all.txt in the same repository.
```

## `cc02_losing_score_direction`: Which way the winner-to-loser trajectory runs

From the paper: Fig. 2a, 2d-e, Fig. 3c and Extended Data Fig. 4d; Methods, identification of a single-cell trajectory in the epiblast

```
Mouse embryos were dissected at embryonic day 5.5 and cultured for 16 hours with a pan-caspase inhibitor (CI) or with vehicle (DMSO). Blocking caspases keeps alive the epiblast cells that would otherwise be eliminated by apoptosis during this window, so that they can be profiled. The embryos were dissociated, pooled by treatment within each of 5 processing batches, and single cells were sequenced with Smart-seq2 with ERCC spike-ins in every well.

`cc_epiblast_tpm.csv` holds TPM values (one row per gene, Ensembl gene id in `gene_id`, one column per cell) for the 512 epiblast cells of the experiment, and `cc_epiblast_cells.csv` gives each cell's treatment (`CI` or `DMSO`) and batch. `cc_genes.csv` gives gene symbols; its `authors_hvg` column marks the 1,921 highly variable genes the lab selected on all cells.

In CI-treated embryos the epiblast cells form one continuum, from a state that is also common in DMSO-treated embryos to states seen only when apoptosis is blocked. Using the 274 CI-treated epiblast cells, order the cells along this continuum by diffusion pseudotime and express it as a losing score: 0 at the end where cells survive without the inhibitor, increasing toward the cells that would otherwise have been eliminated. Report the Spearman correlation, across the 274 CI-treated epiblast cells, between your losing score and log10(TPM + 1) of each of Samm50, mt-Atp6, Atf4 and Sesn2.

Return EXACTLY:

<EVAL_ANSWER>
{"rho_samm50": <float>, "rho_mt_atp6": <float>, "rho_atf4": <float>, "rho_sesn2": <float>}
</EVAL_ANSWER>

`rho_samm50`: Spearman correlation of the losing score with Samm50
`rho_mt_atp6`: Spearman correlation of the losing score with mt-Atp6
`rho_atf4`: Spearman correlation of the losing score with Atf4
`rho_sesn2`: Spearman correlation of the losing score with Sesn2
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "rho_samm50": -0.713,
   "rho_mt_atp6": -0.64,
   "rho_atf4": 0.647,
   "rho_sesn2": 0.533
  },
  "tolerances": {
   "rho_samm50": {
    "type": "absolute",
    "value": 0.3
   },
   "rho_mt_atp6": {
    "type": "absolute",
    "value": 0.3
   },
   "rho_atf4": {
    "type": "absolute",
    "value": 0.3
   },
   "rho_sesn2": {
    "type": "absolute",
    "value": 0.3
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent decides which end of a pseudotime is its start from the experimental design (where cells from untreated embryos fall), instead of leaving it to the arbitrary sign of a diffusion component.

2) Primary trap: Diffusion pseudotime is rooted wherever the first diffusion component happens to put its extreme. The tips of the CI map have 0, 17 and 12 DMSO cells among their 20 nearest neighbours; rooting at the tip with none reverses the score (Spearman -0.999 with the authors' score) and every sign: Samm50 +0.71 instead of -0.71, Atf4 -0.65 instead of +0.65. Read that way, mitochondrial genes would rise and stress genes would fall as cells approach elimination.

3) Scripted controls, run on the task data:
   reference destiny_style_dpt_dmso_oriented: Rebuild the authors' diffusion map of CI epiblast cells (Spearman distance over the lab's variable genes) and root the pseudotime at the tip whose nearest neighbours include the most cells from DMSO-treated embryos. -> rho_samm50=-0.7132, rho_mt_atp6=-0.6401, rho_atf4=0.6475, rho_sesn2=0.5328 (passes)
   variant authors_published_losing_score: Use the authors' own losing score for these cells. -> rho_samm50=-0.7135, rho_mt_atp6=-0.6398, rho_atf4=0.6472, rho_sesn2=0.5325 (passes)
   variant map_on_all_epiblast_cells: Build the diffusion map on all 512 epiblast cells, root it the same way, keep the CI cells. -> rho_samm50=-0.7198, rho_mt_atp6=-0.6466, rho_atf4=0.6563, rho_sesn2=0.5461 (passes)
   variant scanpy_pca_knn_dpt: scanpy: PCA on the lab's variable genes, 15-nearest-neighbour graph, diffmap and dpt, rooted at the end of the first diffusion component next to DMSO cells. -> rho_samm50=-0.724, rho_mt_atp6=-0.6111, rho_atf4=0.6734, rho_sesn2=0.5658 (passes)
   variant scanpy_own_variable_genes: The same scanpy pipeline on its own 2,000 variable genes from ln(TPM + 1) of the CI cells. -> rho_samm50=-0.5795, rho_mt_atp6=-0.7097, rho_atf4=0.7968, rho_sesn2=0.6751 (passes)
   variant scanpy_no_gene_selection: The same scanpy pipeline on ln(TPM + 1) of every gene detected in at least 10 CI cells, with no variable-gene selection. -> rho_samm50=-0.6731, rho_mt_atp6=-0.4257, rho_atf4=0.4005, rho_sesn2=0.2955 (passes)
   trap root_at_ci_only_end: Root the pseudotime at the other end of the trajectory, the tip with no DMSO cell among its nearest neighbours. -> rho_samm50=0.7136, rho_mt_atp6=0.6387, rho_atf4=-0.6457, rho_sesn2=-0.5297 (fails)

4) Tolerance: Windows of +/-0.3 accept the authors' own score, a map built on all epiblast cells, and scanpy diffusion pseudotime on the lab's variable genes, on its own, or with no gene selection (largest shift 0.25): the point is the direction, not the strength of each correlation. A reversed score lands at least 0.76 outside every window.

5) Not said in the prompt: The prompt does not say how to tell the two ends of the trajectory apart.

6) Provenance: The paper roots destiny's DPT in the winner cluster of CI epiblast cells and calls it the losing score. The rebuild here (destiny-style kernel on Spearman distances over the lab's 1,921 variable genes) gives Spearman 0.9998 with the authors' published score (MitoHEAR after_qc.Rda). Samm50 and mt-Atp6 fall along the trajectory in Fig. 3c; Atf4 and Sesn2 rise in Extended Data Fig. 4d.
```

## `cc03_rnr1_heteroplasmy`: How common is mt-Rnr1 heteroplasmy in winner, intermediate and loser epiblast cells

From the paper: Fig. 6h and legend; Methods, analysis of heteroplasmy in single-cell RNA-seq

```
Mouse embryos were dissected at embryonic day 5.5 and cultured for 16 hours with a pan-caspase inhibitor (CI) or with vehicle (DMSO). Blocking caspases keeps alive the epiblast cells that would otherwise be eliminated by apoptosis during this window, so that they can be profiled. The embryos were dissociated, pooled by treatment within each of 5 processing batches, and single cells were sequenced with Smart-seq2 with ERCC spike-ins in every well.

The epiblast cells from CI-treated embryos fall into three clusters along a trajectory: winner cells, which are also common in DMSO-treated embryos, loser cells, found only when apoptosis is blocked, and intermediate cells between them. The lab asks whether mitochondrial DNA variants are more common toward the loser end.

`cc_mt_allele_counts.csv` gives, for each of the 274 epiblast cells from CI-treated embryos, the number of reads supporting A, C, G and T at every mtDNA position (columns `<position>_<base>`, positions 1 to 16299 of the mm10 mtDNA; reads uniquely mapped to the mtDNA, base quality above 33). `cc_ci_epiblast_cells.csv` gives each cell's cluster (`winner`, `intermediate` or `loser`), batch and losing score.

Heteroplasmy at a position is 1 minus the frequency of the most common base at that position in that cell, and it can only be estimated where the position is covered by more than 50 reads in that cell. Take the six mt-Rnr1 positions 300, 303, 304, 305, 326 and 327. For each cluster, report the fraction of cells whose mean heteroplasmy over these six positions is above 0.01, and report the total number of cells that your three fractions are based on.

Return EXACTLY:

<EVAL_ANSWER>
{"frac_heteroplasmic_winner": <float>, "frac_heteroplasmic_intermediate": <float>, "frac_heteroplasmic_loser": <float>, "n_cells_used": <integer>}
</EVAL_ANSWER>

`frac_heteroplasmic_winner`: share of winner cells with mean heteroplasmy over the six positions above 0.01, a fraction between 0 and 1
`frac_heteroplasmic_intermediate`: the same share among intermediate cells, a fraction between 0 and 1
`frac_heteroplasmic_loser`: the same share among loser cells, a fraction between 0 and 1
`n_cells_used`: number of cells, all three clusters together, behind the three fractions
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "frac_heteroplasmic_winner": 0.059,
   "frac_heteroplasmic_intermediate": 0.543,
   "frac_heteroplasmic_loser": 0.873,
   "n_cells_used": 215
  },
  "tolerances": {
   "frac_heteroplasmic_winner": {
    "type": "absolute",
    "value": 0.05
   },
   "frac_heteroplasmic_intermediate": {
    "type": "absolute",
    "value": 0.06
   },
   "frac_heteroplasmic_loser": {
    "type": "absolute",
    "value": 0.08
   },
   "n_cells_used": {
    "type": "absolute",
    "value": 20
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent keeps cells without enough reads out of a per-cell variant summary instead of counting them as cells without the variant.

2) Primary trap: Heteroplasmy is masked where coverage is too low, as the prompt says, but the per-cell mean is then missing for cells without full coverage, a missing value compares as False against 0.01 without any error, and those cells stay in their cluster's denominator: no data is read as no variant. Coverage is incomplete for 10 of 78 winner, 33 of 125 intermediate and 16 of 71 loser cells, so the dilution falls where heteroplasmy is common: the loser fraction drops from 0.87 to 0.68, the intermediate fraction from 0.54 to 0.40, and the cell count from 215 to 274.

3) Scripted controls, run on the task data:
   reference all_six_covered: Keep cells covered by more than 50 reads at all six positions (after the paper's filter of more than 2,000 covered mtDNA positions per cell), average heteroplasmy over the six, and divide by the cells kept in each cluster. -> frac_heteroplasmic_winner=0.05882, frac_heteroplasmic_intermediate=0.5435, frac_heteroplasmic_loser=0.8727, n_cells_used=215 (passes)
   variant all_six_covered_no_cell_filter: The same without the per-cell filter on covered mtDNA positions. -> frac_heteroplasmic_winner=0.05882, frac_heteroplasmic_intermediate=0.5435, frac_heteroplasmic_loser=0.8727, n_cells_used=215 (passes)
   variant at_least_50_reads: Count a position as covered at 50 reads or more instead of more than 50. -> frac_heteroplasmic_winner=0.05882, frac_heteroplasmic_intermediate=0.5484, frac_heteroplasmic_loser=0.8727, n_cells_used=216 (passes)
   variant mean_over_covered_positions: Average over whichever of the six positions are covered in each cell. -> frac_heteroplasmic_winner=0.05882, frac_heteroplasmic_intermediate=0.5758, frac_heteroplasmic_loser=0.8793, n_cells_used=225 (passes)
   trap uncovered_cells_kept_in_denominator: Positions with 50 reads or fewer are set to missing; a cell's mean over the six is then missing, compares as False against 0.01, and the cell stays in the denominator. -> frac_heteroplasmic_winner=0.05128, frac_heteroplasmic_intermediate=0.4, frac_heteroplasmic_loser=0.6761, n_cells_used=274 (fails)

4) Tolerance: The windows accept the paper's per-cell filter on covered mtDNA positions or none, a 50-read cut taken as at least 50, and averaging over the covered subset of the six positions. They do not accept keeping cells without full coverage in the denominator, whether their mean is missing or their uncovered positions are set to 0 (intermediate 0.44, loser 0.70, 274 cells).

5) Not said in the prompt: The prompt does not say what to do with cells that lack coverage at some of the six positions, or that coverage differs between clusters.

6) Provenance: The paper reports about 5%, 55% and 87% (Fig. 6h, 215 cells). The allele counts are raw_counts_allele_mt.csv from E-MTAB-8640, the table behind Fig. 6. Its rows are pileup files named by older sequencing-file names, joined to cells through the authors' own name tables (name_full.Rda and provo.Rda), as their figure code does. The join was checked by script: the cells carrying inherited non-reference alleles at mtDNA positions 12217 and 13536 are the same in this table and in the independent re-count shipped with the authors' MitoHEAR package (333 of 333 cell-position calls agree, 93 carriers). Cluster labels and losing scores are the authors' (MitoHEAR after_qc.Rda).
```

## `cc04_resp_chain_enrichment`: Are respiratory chain genes over-represented among genes that fall along the losing score

From the paper: Fig. 3a-b; Methods, differential gene expression analysis along the trajectory

```
Mouse embryos were dissected at embryonic day 5.5 and cultured for 16 hours with a pan-caspase inhibitor (CI) or with vehicle (DMSO). Blocking caspases keeps alive the epiblast cells that would otherwise be eliminated by apoptosis during this window, so that they can be profiled. The embryos were dissociated, pooled by treatment within each of 5 processing batches, and single cells were sequenced with Smart-seq2 with ERCC spike-ins in every well.

`cc_epiblast_tpm.csv` holds TPM values (one row per gene, Ensembl gene id in `gene_id`, one column per cell) for the epiblast cells. `cc_ci_epiblast_cells.csv` lists the 274 epiblast cells from CI-treated embryos with the lab's losing score, a diffusion pseudotime that is 0 in cells common to both treatments and highest in cells found only when apoptosis is blocked. `cc_genes.csv` is the Ensembl gene annotation (gene id, symbol, chromosome, start, end).

The lab wants to know whether cells switch off the respiratory chain as they approach elimination. Using the 274 CI-treated epiblast cells, call a gene decreasing if a test of trend of log10(TPM + 1) against the losing score gives a negative slope (or correlation) at a Benjamini-Hochberg FDR below 0.01. Then test, with a one-sided Fisher exact test, whether respiratory chain and ATP synthase genes are over-represented among the decreasing genes. Define that gene set as the gene symbols matching the regular expression `^(Nduf[abcsv][0-9]|Sdh[abcd]$|Uqcr|Cyc1$|Cox[4-8]|Atp5|mt-(Nd|Co|Cytb|Atp))`.

Report: the share of the genes in your test's background that you call decreasing; the share of the respiratory chain genes in that background that you call decreasing; the fold enrichment (the second share divided by the first); and the one-sided Fisher p-value as -log10(p).

Return EXACTLY:

<EVAL_ANSWER>
{"frac_genes_decreasing": <float>, "frac_resp_chain_decreasing": <float>, "fold_enrichment": <float>, "neg_log10_p": <float>}
</EVAL_ANSWER>

`frac_genes_decreasing`: share of all background genes called decreasing, a fraction between 0 and 1
`frac_resp_chain_decreasing`: share of respiratory chain genes in the background called decreasing, a fraction between 0 and 1
`fold_enrichment`: the second share divided by the first
`neg_log10_p`: -log10 of the one-sided Fisher exact test p-value
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "frac_genes_decreasing": 0.406,
   "frac_resp_chain_decreasing": 0.949,
   "fold_enrichment": 2.34,
   "neg_log10_p": 29.9
  },
  "tolerances": {
   "frac_genes_decreasing": {
    "type": "absolute",
    "value": 0.15
   },
   "frac_resp_chain_decreasing": {
    "type": "absolute",
    "value": 0.08
   },
   "fold_enrichment": {
    "type": "absolute",
    "value": 1.2
   },
   "neg_log10_p": {
    "type": "absolute",
    "value": 15
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent takes the background of an over-representation test from the genes that could have been called, not from the whole genome.

2) Primary trap: The decreasing genes are compared with every gene in the annotation (53,700 genes) instead of the 9,916 genes that were expressed enough to be tested. Most of the genome could never have been called, so the share of background genes called decreasing falls from 0.406 to 0.074, the fold enrichment rises from 2.34 to 10.61, and -log10 p from 30 to 81. The respiratory chain result holds either way; its strength is what the background inflates.

3) Scripted controls, run on the task data:
   reference tested_genes_background: Test genes with more than 15 TPM in more than 10 cells (the paper's filter) for a linear trend on the losing score; use those tested genes as the background. -> frac_genes_decreasing=0.4063, frac_resp_chain_decreasing=0.949, fold_enrichment=2.336, neg_log10_p=29.92 (passes)
   variant spearman_trend_test: Same genes and background, Spearman correlation test instead of the linear trend. -> frac_genes_decreasing=0.3937, frac_resp_chain_decreasing=0.9082, fold_enrichment=2.307, neg_log10_p=25.99 (passes)
   variant lighter_expression_filter: Test genes above 1 TPM in at least 10% of cells; background = those genes. -> frac_genes_decreasing=0.363, frac_resp_chain_decreasing=0.9485, fold_enrichment=2.613, neg_log10_p=33.88 (passes)
   variant stricter_expression_filter: Test genes with mean TPM of at least 10; background = those genes. -> frac_genes_decreasing=0.5187, frac_resp_chain_decreasing=0.9574, fold_enrichment=1.846, neg_log10_p=20.7 (passes)
   trap whole_genome_background: Same decreasing genes, but the background is every gene in the annotation, as web enrichment tools assume by default. -> frac_genes_decreasing=0.07428, frac_resp_chain_decreasing=0.7881, fold_enrichment=10.61, neg_log10_p=80.84 (fails)

4) Tolerance: The windows accept a Spearman trend test and expression filters from about 1 TPM in 10% of cells to a mean of 10 TPM, with the tested genes as the background. They do not accept a background of the whole annotation, nor testing every gene detected in at least 3 cells and using those as the background (share decreasing 0.172, fold enrichment 4.79, -log10 p 50).

5) Not said in the prompt: The prompt does not say which genes to test or which genes form the background.

6) Provenance: The paper ran IPA, a commercial tool, on the genes changing along the trajectory with all tested genes (more than 15 TPM in more than 10 cells) as the background, and reported Mitochondrial Dysfunction (-log10 p = 21.1) and Oxidative Phosphorylation (-log10 p = 18.6) as the top pathways. This task uses an open, explicitly defined gene set (98 genes among the tested genes, 118 in the annotation) and Fisher's exact test; it does not reproduce IPA. With the paper's filter, 9,916 genes are tested and 5,278 change at FDR < 0.01 (the paper reports 5,311). The losing score is the authors' (MitoHEAR after_qc.Rda).
```
