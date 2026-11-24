# Controls: AVE migration, Dev Cell 2024

30 of 30 scripted analyses behave as expected (graded with this repository's own check and latch-eval-tools 0.4.54).

| Task | Control | Kind | Should | Own grader | Latch grader | What it does |
|---|---|---|---|---|---|---|
| `ave01_plate_qc` | reference | reference | pass | pass | pass | The paper's five thresholds (genes above 10 rpm > 3,000, log10 read pairs > 4, mapped > 0.5, mitochondrial < 0.1, ERCC < 0.3, as fractions of read pairs) on the wells annotated as single cells. |
| `ave01_plate_qc` | thresholds_all_wells | variant | pass | pass | pass | The same thresholds, wells annotated as no cell or multiple cells not set aside. |
| `ave01_plate_qc` | counted_reads_only | variant | pass | pass | pass | The thresholds on shares of counted reads, without read pairs or the mapping filter. |
| `ave01_plate_qc` | gene_floor_2000 | variant | pass | pass | pass | The thresholds with a gene floor of 2,000 instead of 3,000. |
| `ave01_plate_qc` | mad_outliers | variant | pass | pass | pass | 3-MAD outliers on log read pairs, log genes detected, ERCC share and mitochondrial share, among single-cell wells. |
| `ave01_plate_qc` | droplet_defaults | trap | fail | fail | fail | Droplet defaults: at least 200 genes detected and a mitochondrial share of counts under 20%. |
| `ave01_plate_qc` | genes_and_mito_only | trap | fail | fail | fail | More than 1,000 genes detected and a mitochondrial share under 10%, on single-cell wells; spike-in share and mapping rate not used. |
| `ave02_pseudotime_root` | reference | reference | pass | pass | pass | Root at the end of the first diffusion component where the Cer1, Hhex and Lefty1 score is high; natural spline F-test with 4 df. |
| `ave02_pseudotime_root` | root_max_marker_score | variant | pass | pass | pass | Root at the cell with the highest Cer1, Hhex and Lefty1 score. |
| `ave02_pseudotime_root` | euclidean_knn | variant | pass | pass | pass | Euclidean distance on the 10 principal components instead of the Spearman distance. |
| `ave02_pseudotime_root` | scaled_pca | variant | pass | pass | pass | Genes scaled to unit variance before the PCA. |
| `ave02_pseudotime_root` | root_dc2_extreme_in_ave | variant | pass | pass | pass | The paper's recipe, an extreme of the second diffusion component, taken at the end with the higher AVE marker score. |
| `ave02_pseudotime_root` | authors_root_cell | variant | pass | pass | pass | The cell the authors rooted at (minimum of their second component, from DC_55.csv), placed in this graph. |
| `ave02_pseudotime_root` | linear_trend_test | variant | pass | pass | pass | Linear regression on pseudotime instead of the spline. |
| `ave02_pseudotime_root` | spearman_test | variant | pass | pass | pass | Spearman correlation test instead of the spline. |
| `ave02_pseudotime_root` | spline_3df | variant | pass | pass | pass | Natural spline with 3 degrees of freedom. |
| `ave02_pseudotime_root` | root_at_emve_end | trap | fail | fail | fail | Root at the extreme of the first diffusion component that lies in the emVE, as an unchecked argmin or argmax gives when the sign falls that way. |
| `ave02_pseudotime_root` | root_first_cell | trap | fail | fail | fail | Root at the first cell of the file (iroot = 0), an emVE cell. |
| `ave03_ephrin_bootstrap_signs` | reference | reference | pass | pass | pass | Each bootstrap axis reflected to match the full-data axis before signs are compared (numpy SVD, 1,000 samples, seed 0). |
| `ave03_ephrin_bootstrap_signs` | axis_reordering_too | variant | pass | pass | pass | Bootstrap axes also matched to the full-data axes by their loadings (Peres-Neto reordering) before reflection. |
| `ave03_ephrin_bootstrap_signs` | other_seed | variant | pass | pass | pass | Another random seed. |
| `ave03_ephrin_bootstrap_signs` | sklearn_pca | variant | pass | pass | pass | scikit-learn PCA (full SVD) with reflection. |
| `ave03_ephrin_bootstrap_signs` | two_thousand_samples | variant | pass | pass | pass | 2,000 bootstrap samples. |
| `ave03_ephrin_bootstrap_signs` | signs_not_aligned | trap | fail | fail | fail | Loading signs compared as the SVD returns them, axes not aligned. |
| `ave03_ephrin_bootstrap_signs` | sklearn_signs_not_aligned | trap | fail | fail | fail | Loading signs compared under scikit-learn's own sign convention, axes not aligned. |
| `ave04_dkk1_timing` | reference | reference | pass | pass | pass | AUC of ln(x + 1), AVE against emVE, at each stage; conclusion chosen from the measured pattern. |
| `ave04_dkk1_timing` | log2_scale | variant | pass | pass | pass | log2(x + 1) instead of ln(x + 1). |
| `ave04_dkk1_timing` | roc_auc_score | variant | pass | pass | pass | AUC from scikit-learn's roc_auc_score. |
| `ave04_dkk1_timing` | stages_pooled | trap | fail | fail | fail | AVE and emVE cells of both stages pooled, and the pooled enrichment read as support for guidance. |
| `ave04_dkk1_timing` | prior_driven_reading | trap | fail | fail | fail | Correct AUCs, but the enrichment at 6.25 dpc and DKK1's published role taken as support. |
