# agentic-science-evals

**Evaluation tasks for AI agents, built from published papers by one of their
authors.**

Each task takes one real analysis decision from a paper, gives an agent the
public files a scientist had at that point, and asks a question whose answer a
script can check. Every task ships with the scripted controls its grader was
tested against: the reference analysis and a defensible variant must pass, and
a plausible wrong analysis must fail.

- Tasks are written in LatchBio's public eval format (`eval.json`). Independent
  work, not affiliated with LatchBio.
- Written and checked by AI models under the author's direction. No person has
  reviewed the answer keys, the windows or the controls yet.

## Papers

| Paper | Author's part | Public data | Tasks |
|---|---|---|---|
| CUT&Tag recovers up to half of ENCODE ChIP-seq histone acetylation peaks. Nat Commun 2025 | Joint first author; performed the laboratory experiments with a co-author | GEO GSE286492 | 6 tasks and a pilot, in [epigenomics-evals](https://github.com/di-omics/epigenomics-evals) |
| An integrated approach identifies the molecular underpinnings of murine anterior visceral endoderm migration. Dev Cell 2024 | Sample collections for the transcriptomics | ArrayExpress E-MTAB-9645 | 4 tasks: [ave_migration](src/science_evals/papers/ave_migration/README.md) |
| Cell competition acts as a purifying selection to eliminate cells with mitochondrial defects during early mouse development. Nat Metab 2021 | Embryo dissections, treatments and cell dissociation before single-cell RNA-seq | ArrayExpress E-MTAB-8640 | 4 tasks: [cell_competition](src/science_evals/papers/cell_competition/README.md) |
| The APOE isoforms differentially shape the transcriptomic and epigenomic landscapes of human microglia xenografted into a mouse model of Alzheimer's disease. Nat Commun 2025 | Second author; performed experiments with 4 co-authors | GEO GSE271384, GSE271385 | 4 tasks: [apoe_microglia](src/science_evals/papers/apoe_microglia/README.md) |

## The tasks

| Task | The decision it tests |
|---|---|
| `ave01_plate_qc` | Which wells of a 384-well Smart-seq2 plate hold one good cell |
| `ave02_pseudotime_root` | Where to root a pseudotime so that it starts in the AVE |
| `ave03_ephrin_bootstrap_signs` | Aligning bootstrap PCA axes before comparing loading signs |
| `ave04_dkk1_timing` | Whether Dkk1 marks the AVE while it migrates |
| `cc01_cell_qc` | Whether extra QC filtering removes the very cells the study is about |
| `cc02_losing_score_direction` | Which way a winner-to-loser trajectory runs |
| `cc03_rnr1_heteroplasmy` | Which cells belong in the denominator of a heteroplasmy fraction |
| `cc04_resp_chain_enrichment` | The right background for a gene set enrichment test |
| `apoe01_failed_knockouts` | Confirming a knockout before comparing against it |
| `apoe02_transcript_length_offsets` | Using tximport transcript lengths with gene-level counts |
| `apoe03_vdr_targets_direction` | Reading which way a stored comparison runs |
| `apoe04_matched_contrast_concordance` | Taking accessibility from the same comparison as expression |

12 tasks, 72 scripted controls (56 that must pass, 16 that must fail), all behaving as
expected under both graders, and 18 statements of the papers recomputed from their public
files. None of it has been reviewed by a person yet, and no model has been run on the tasks.

Each paper's folder has the full table: what the reference analysis finds, the trap,
and where the trap lands against the grader's window. Prompts, keys and notes are in
`docs/<paper>/TASKS.md`; every control and its two verdicts in `results/<paper>/CONTROLS.md`;
statements of each paper recomputed from its public files in `results/<paper>/PAPER_ANCHORS.md`.

## Why the controls

A grader accepts a range of answers. Whether that range means anything depends
on what falls inside it and what falls outside. So before a task is accepted,
3 kinds of analysis are scripted, run on the task's data and graded:

- **reference**: the analysis the task is built around. Must pass.
- **variant**: a defensible alternative. Must pass. The tolerance is set so that it does.
- **trap**: a plausible wrong analysis, named in the task notes. Must fail.

This is assay development applied to benchmark writing: a readout is not
trusted until its positive and negative controls separate. Every control is
graded twice, by this repository's own check and by LatchBio's published
`latch-eval-tools`, and the two verdicts must agree.

## Layout

```
src/science_evals/            shared: task format, graders, controls, build
src/science_evals/papers/     one package per paper: data, reference analyses, tasks
evals/<paper>/<task>/         eval.json files (generated)
results/<paper>/              controls and paper anchors (generated)
docs/<paper>/TASKS.md         prompts, keys and notes in readable form (generated)
docs/BUILDING.md              how a paper becomes tasks, and the rules
```

## Rebuild

```
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[analysis,grade,dev]"
.venv/bin/python -m science_evals.build          # downloads public inputs into data/, then builds
.venv/bin/python -m pytest -q                    # needs no data
```

## License

Code: MIT. Each paper's data keep their own licenses; see the paper's folder.
