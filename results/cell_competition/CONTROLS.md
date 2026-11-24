# Controls: Cell competition, Nat Metab 2021

22 of 22 scripted analyses behave as expected (graded with this repository's own check and latch-eval-tools 0.4.54).

| Task | Control | Kind | Should | Own grader | Latch grader | What it does |
|---|---|---|---|---|---|---|
| `cc01_cell_qc` | authors_fixed_thresholds_only | reference | pass | pass | pass | Keep every cell that passed the fixed thresholds, as the authors did. |
| `cc01_cell_qc` | high_mito_outliers_removed | variant | pass | pass | pass | Also remove cells more than 3 MADs above the median mitochondrial fraction. |
| `cc01_cell_qc` | fixed_mito_cut_10pct | variant | pass | pass | pass | Also remove cells with more than 10% of reads on mitochondrial genes. |
| `cc01_cell_qc` | adaptive_qc_within_treatment_5mad | variant | pass | pass | pass | Adaptive outlier rules (5 MADs) computed separately within each treatment group. |
| `cc01_cell_qc` | adaptive_qc_pooled_3mad | trap | fail | fail | fail | scater-style outlier rules at the default 3 MADs (library size, detected genes, ERCC and mitochondrial fractions) on all cells pooled. |
| `cc02_losing_score_direction` | destiny_style_dpt_dmso_oriented | reference | pass | pass | pass | Rebuild the authors' diffusion map of CI epiblast cells (Spearman distance over the lab's variable genes) and root the pseudotime at the tip whose nearest neighbours include the most cells from DMSO-treated embryos. |
| `cc02_losing_score_direction` | authors_published_losing_score | variant | pass | pass | pass | Use the authors' own losing score for these cells. |
| `cc02_losing_score_direction` | map_on_all_epiblast_cells | variant | pass | pass | pass | Build the diffusion map on all 512 epiblast cells, root it the same way, keep the CI cells. |
| `cc02_losing_score_direction` | scanpy_pca_knn_dpt | variant | pass | pass | pass | scanpy: PCA on the lab's variable genes, 15-nearest-neighbour graph, diffmap and dpt, rooted at the end of the first diffusion component next to DMSO cells. |
| `cc02_losing_score_direction` | scanpy_own_variable_genes | variant | pass | pass | pass | The same scanpy pipeline on its own 2,000 variable genes from ln(TPM + 1) of the CI cells. |
| `cc02_losing_score_direction` | scanpy_no_gene_selection | variant | pass | pass | pass | The same scanpy pipeline on ln(TPM + 1) of every gene detected in at least 10 CI cells, with no variable-gene selection. |
| `cc02_losing_score_direction` | root_at_ci_only_end | trap | fail | fail | fail | Root the pseudotime at the other end of the trajectory, the tip with no DMSO cell among its nearest neighbours. |
| `cc03_rnr1_heteroplasmy` | all_six_covered | reference | pass | pass | pass | Keep cells covered by more than 50 reads at all six positions (after the paper's filter of more than 2,000 covered mtDNA positions per cell), average heteroplasmy over the six, and divide by the cells kept in each cluster. |
| `cc03_rnr1_heteroplasmy` | all_six_covered_no_cell_filter | variant | pass | pass | pass | The same without the per-cell filter on covered mtDNA positions. |
| `cc03_rnr1_heteroplasmy` | at_least_50_reads | variant | pass | pass | pass | Count a position as covered at 50 reads or more instead of more than 50. |
| `cc03_rnr1_heteroplasmy` | mean_over_covered_positions | variant | pass | pass | pass | Average over whichever of the six positions are covered in each cell. |
| `cc03_rnr1_heteroplasmy` | uncovered_cells_kept_in_denominator | trap | fail | fail | fail | Positions with 50 reads or fewer are set to missing; a cell's mean over the six is then missing, compares as False against 0.01, and the cell stays in the denominator. |
| `cc04_resp_chain_enrichment` | tested_genes_background | reference | pass | pass | pass | Test genes with more than 15 TPM in more than 10 cells (the paper's filter) for a linear trend on the losing score; use those tested genes as the background. |
| `cc04_resp_chain_enrichment` | spearman_trend_test | variant | pass | pass | pass | Same genes and background, Spearman correlation test instead of the linear trend. |
| `cc04_resp_chain_enrichment` | lighter_expression_filter | variant | pass | pass | pass | Test genes above 1 TPM in at least 10% of cells; background = those genes. |
| `cc04_resp_chain_enrichment` | stricter_expression_filter | variant | pass | pass | pass | Test genes with mean TPM of at least 10; background = those genes. |
| `cc04_resp_chain_enrichment` | whole_genome_background | trap | fail | fail | fail | Same decreasing genes, but the background is every gene in the annotation, as web enrichment tools assume by default. |
