# Controls: APOE microglia, Nat Commun 2025

20 of 20 scripted analyses behave as expected (graded with this repository's own check and latch-eval-tools 0.4.54).

| Task | Control | Kind | Should | Own grader | Latch grader | What it does |
|---|---|---|---|---|---|---|
| `apoe01_failed_knockouts` | reference | reference | pass | pass | pass | Drop the knockout grafts whose APOE level is at least a quarter of the APOE3 median, then DESeq2 with the tximport lengths as offsets, BH over all genes. |
| `apoe01_failed_knockouts` | plain_counts | variant | pass | pass | pass | Same samples; counts without the transcript lengths. |
| `apoe01_failed_knockouts` | length_scaled_counts | variant | pass | pass | pass | Same samples; lengths divided out of the counts first. |
| `apoe01_failed_knockouts` | independent_filtering | variant | pass | pass | pass | Same samples and offsets; DESeq2's padj with independent filtering instead of BH over all genes. |
| `apoe01_failed_knockouts` | pydeseq2 | variant | pass | pass | pass | Same samples; the pydeseq2 package on the counts. |
| `apoe01_failed_knockouts` | all_knockout_grafts | trap | fail | fail | fail | Keep all five APOE-KO samples, two of which express APOE; otherwise the reference analysis. |
| `apoe02_transcript_length_offsets` | reference | reference | pass | pass | pass | APOE2, APOE3 and APOE4 samples; DESeq2 with the tximport lengths as offsets; BH over all genes. |
| `apoe02_transcript_length_offsets` | length_scaled_counts | variant | pass | pass | pass | Lengths divided out of the counts, then DESeq2 with size factors. |
| `apoe02_transcript_length_offsets` | pydeseq2_length_scaled | variant | pass | pass | pass | Length-scaled counts with the pydeseq2 package. |
| `apoe02_transcript_length_offsets` | independent_filtering | variant | pass | pass | pass | Offsets; DESeq2's padj with independent filtering. |
| `apoe02_transcript_length_offsets` | paper_model | variant | pass | pass | pass | The paper's model: offsets, all samples except the two knockout grafts that express APOE. |
| `apoe02_transcript_length_offsets` | counts_without_lengths | trap | fail | fail | fail | The gene-level counts tested as plain read counts, size factors only. |
| `apoe03_vdr_targets_direction` | reference | reference | pass | pass | pass | Genes with fdr < 0.05 and a negative log2 fold change in the APOE4 against APOE2 sheet; unique symbols; hypergeometric test inside the sheet's symbols. |
| `apoe03_vdr_targets_direction` | rows_not_symbols | variant | pass | pass | pass | Count table rows instead of unique symbols. |
| `apoe03_vdr_targets_direction` | authors_universe | variant | pass | pass | pass | The authors' code: all rows as the universe and every VDR target counted, inside the sheet or not. |
| `apoe03_vdr_targets_direction` | reversed_direction | trap | fail | fail | fail | Take a positive log2 fold change in the APOE4 against APOE2 sheet as higher in APOE2. |
| `apoe04_matched_contrast_concordance` | reference | reference | pass | pass | pass | For each comparison, median promoter log2 fold change per gene from the same comparison; Pearson r over DEGs with a promoter peak; genes matched by Entrez ID. |
| `apoe04_matched_contrast_concordance` | mean_per_gene | variant | pass | pass | pass | Mean instead of median of a gene's promoter peaks. |
| `apoe04_matched_contrast_concordance` | match_by_symbol | variant | pass | pass | pass | Genes matched by symbol instead of Entrez ID. |
| `apoe04_matched_contrast_concordance` | pooled_comparisons | trap | fail | fail | fail | Each gene's promoter summary taken over the rows of all four comparisons at once. |
