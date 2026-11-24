# Building a paper's tasks

How a paper is turned into tasks here: what its package provides, the rules every
task follows, and the commands.

## Layout per paper

```
src/science_evals/papers/<slug>/
  __init__.py    PAPER, DATA_FILES, build(), anchors()
  reference.py   the analyses (reference, variants, traps) as functions of the data
  README.md      the paper, the author's role, the data, the tasks, how to rebuild
tests/test_<slug>.py
data/<slug>/     downloaded inputs (git-ignored)
```

`python -m science_evals.build <slug>` writes, and nothing else edits:
`evals/<slug>/`, `results/<slug>/` and `docs/<slug>/TASKS.md`.

## The contract

**`PAPER`**: a dict with `short` (for example "AVE migration, Dev Cell 2024"),
`citation` (authors, title, journal, year, volume, pages, doi), `doi`,
`author_role` (the repository owner's part in the paper, quoted from the
paper's author contributions), `data` (accessions), `code` (the authors'
repository, with the commit used) and `license` (of the paper and the data).

**`DATA_FILES`**: a list of `{"name", "source", "sha256", "bytes"}`. `source`
is a public https URL that does not move (an accession file, or a raw URL at a
fixed commit), or `derived:<function>` for a file this package makes from
public files. A derived file must be deterministic; its sha256 pins it.

**`build(data_dir) -> list[Built]`**: fetches (`science_evals.data.fetch`) or
derives the inputs, runs the reference analysis and every control, and returns
one `Built(task, reference, controls, data_files)` per task. Fixed seeds only:
two builds give identical files.

**`anchors(data_dir) -> list[dict]`** (recommended): 3 to 6 statements from
the paper recomputed from the public files, each
`{"paper_says", "where", "recomputed", "agrees"}`. When a number does not
reproduce, say so and say why if known.

## What makes a task

- **One real decision from the paper**, where a competent-looking analysis can
  go wrong without any error message: which cells to keep, which comparison is
  the right one, which direction a trajectory runs, what a missing signal means.
- **The files a scientist had at that point**, and a question whose answer can
  be checked: numbers (absolute tolerance), a set of labels (Jaccard index), and
  at most 1 multiple-choice decision.
- **The key comes from the reference analysis**, run on the data. No number in
  an eval.json is typed by hand.
- **Controls**: at least 1 `reference` and 1 `variant` that must pass, and 1
  `trap` that must fail. Each is code run on the data, not a number placed near
  a window edge. Tolerances are set so the variants pass and the traps fail;
  `tolerance_rationale` says why.
- **One trap per task.** Differences that are not the point of a task go into
  its tolerance, as variants.
- **Multiple choice**: options of similar length; the keyed option is neither
  the longest nor the shortest; every wrong option is something a real analyst
  might do. Its `key_review` reads "Decision key checked against the data by
  script. Not yet reviewed by a person."
- **The prompt does not give the trap away.** What it leaves out on purpose
  goes in `not_said`.
- **Failure modes** come from taxonomy v1.1.0. Ids already in use:
  `qc_control_artifact.context_inappropriate_qc_threshold`,
  `qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology`,
  `metric_statistic_mismatch.wrong_metric_summary_rank_or_association`,
  `metric_statistic_mismatch.null_background_or_control_calibration_error`,
  `subset_scope_compartment.wrong_comparator_contrast_or_baseline`,
  `marker_reference_label_scope.dropout_panel_limit_or_absence_overclaim`,
  `biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error`,
  `biological_interpretation.prior_driven_or_overclaimed_biological_interpretation`,
  `assay_intervention_feature_linkage.barcode_feature_observation_or_metadata_linkage_error`.
- **Size**: a capable agent can finish a task in under 15 minutes on a laptop;
  inputs are ideally under 200 MB per task.

## Rules

- Public data only, from the accession, the journal's supplementary files or
  the authors' repository. No raw reads. A source that blocks scripts (a 403 or
  a bot check) is skipped, never worked around.
- ASCII only in every tracked file: no em or en dashes, no curly quotes, no
  Greek letters or symbols (write `+/-`, `um`, `deg`). Numbers as digits in prose.
- No file, folder or id named "cuttag". Write the assay as "CUT&Tag" in prose.
- No names of companies or people the work may be shown to.
- Say what checked what: "checked against the data by script". Never imply a
  person reviewed a key, a window or a control.
- The author's role is stated exactly as the paper's author contributions say it.

## Commands

```
.venv/bin/python -m science_evals.build <slug>
.venv/bin/python -m science_evals.build --index
.venv/bin/python -m pytest -q
```
