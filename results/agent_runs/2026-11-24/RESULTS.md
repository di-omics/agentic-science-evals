# Agent runs

One fresh agent per model and task, on a folder holding only TASK.md and the task's input files. Graded by this repository's own check and by latch-eval-tools; an audit of every tool call looks for web use and for paths outside the run folder.

| Task | haiku | sonnet | opus |
|---|---|---|---|
| `apoe01_failed_knockouts` | fail | fail | pass |
| `apoe02_transcript_length_offsets` | fail | pass | pass |
| `apoe03_vdr_targets_direction` | pass | pass | pass |
| `apoe04_matched_contrast_concordance` | pass | pass | pass |
| `ave01_plate_qc` | fail | pass | pass |
| `ave02_pseudotime_root` | pass | pass | pass |
| `ave03_ephrin_bootstrap_signs` | fail | pass | pass |
| `ave04_dkk1_timing` | pass | pass | pass |
| `cc01_cell_qc` | fail | pass | pass |
| `cc02_losing_score_direction` | pass | pass | pass |
| `cc03_rnr1_heteroplasmy` | pass | pass | pass |
| `cc04_resp_chain_enrichment` | fail | pass | pass |

Passed: haiku: 6 of 12, sonnet: 11 of 12, opus: 12 of 12.

## Each run

| Task | Model | Verdict | Graders agree | Tool calls | Seconds | Audit | What missed |
|---|---|---|---|---|---|---|---|
| `apoe01_failed_knockouts` | haiku | fail | yes | 10 | 358 | clean | n_deg 9 (window 293 to 413); n_higher_in_ko 4 (window 108 to 168); apoe_log2fc -2.134 (window -4.477 to -3.677) |
| `apoe01_failed_knockouts` | opus | pass | yes | 14 | 1121 | clean |  |
| `apoe01_failed_knockouts` | sonnet | fail | yes | 4 | 109 | clean | n_deg 455 (window 293 to 413); n_higher_in_ko 212 (window 108 to 168); apoe_log2fc -2.13 (window -4.477 to -3.677) |
| `apoe02_transcript_length_offsets` | haiku | fail | yes | 21 | 948 | touched the harness's output folder | n_deg 521 (window 220 to 360); n_higher_in_apoe4 281 (window 123 to 203) |
| `apoe02_transcript_length_offsets` | opus | pass | yes | 15 | 1033 | clean |  |
| `apoe02_transcript_length_offsets` | sonnet | pass | yes | 7 | 221 | touched the harness's output folder |  |
| `apoe03_vdr_targets_direction` | haiku | pass | yes | 5 | 49 | clean |  |
| `apoe03_vdr_targets_direction` | opus | pass | yes | 5 | 132 | clean |  |
| `apoe03_vdr_targets_direction` | sonnet | pass | yes | 3 | 64 | clean |  |
| `apoe04_matched_contrast_concordance` | haiku | pass | yes | 3 | 49 | clean |  |
| `apoe04_matched_contrast_concordance` | opus | pass | yes | 10 | 194 | clean |  |
| `apoe04_matched_contrast_concordance` | sonnet | pass | yes | 4 | 98 | clean |  |
| `ave01_plate_qc` | haiku | fail | yes | 9 | 93 | clean | wells_kept (Jaccard 0.4792) |
| `ave01_plate_qc` | opus | pass | yes | 20 | 514 | clean |  |
| `ave01_plate_qc` | sonnet | pass | yes | 9 | 285 | clean |  |
| `ave02_pseudotime_root` | haiku | pass | yes | 7 | 272 | touched the harness's output folder |  |
| `ave02_pseudotime_root` | opus | pass | yes | 13 | 533 | touched the harness's output folder |  |
| `ave02_pseudotime_root` | sonnet | pass | yes | 11 | 223 | touched the harness's output folder |  |
| `ave03_ephrin_bootstrap_signs` | haiku | fail | yes | 3 | 22 | clean | n_significant_e55 2 (window 4 to 8); significant_genes_e55 (Jaccard 0.3333) |
| `ave03_ephrin_bootstrap_signs` | opus | pass | yes | 12 | 392 | clean |  |
| `ave03_ephrin_bootstrap_signs` | sonnet | pass | yes | 5 | 31 | clean |  |
| `ave04_dkk1_timing` | haiku | pass | yes | 5 | 18 | clean |  |
| `ave04_dkk1_timing` | opus | pass | yes | 7 | 131 | clean |  |
| `ave04_dkk1_timing` | sonnet | pass | yes | 3 | 26 | clean |  |
| `cc01_cell_qc` | haiku | fail | yes | 5 | 24 | clean | n_cells_kept 548 (window 678 to 768); n_ci_cells_kept 247 (window 347 to 387); n_dmso_cells_kept 301 (window 321 to 391) |
| `cc01_cell_qc` | opus | pass | yes | 7 | 247 | clean |  |
| `cc01_cell_qc` | sonnet | pass | yes | 5 | 38 | clean |  |
| `cc02_losing_score_direction` | haiku | pass | yes | 18 | 405 | touched the harness's output folder |  |
| `cc02_losing_score_direction` | opus | pass | yes | 18 | 841 | touched the harness's output folder |  |
| `cc02_losing_score_direction` | sonnet | pass | yes | 19 | 319 | touched the harness's output folder |  |
| `cc03_rnr1_heteroplasmy` | haiku | pass | yes | 4 | 21 | clean |  |
| `cc03_rnr1_heteroplasmy` | opus | pass | yes | 8 | 262 | clean |  |
| `cc03_rnr1_heteroplasmy` | sonnet | pass | yes | 3 | 34 | clean |  |
| `cc04_resp_chain_enrichment` | haiku | fail | yes | 7 | 47 | clean | frac_genes_decreasing 0.1933 (window 0.256 to 0.556); frac_resp_chain_decreasing 0.8667 (window 0.869 to 1.029); fold_enrichment 4.483 (window 1.14 to 3.54); neg_log10_p 49.57 (window 14.9 to 44.9) |
| `cc04_resp_chain_enrichment` | opus | pass | yes | 10 | 392 | clean |  |
| `cc04_resp_chain_enrichment` | sonnet | pass | yes | 5 | 50 | clean |  |
