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
| An integrated approach identifies the molecular underpinnings of murine anterior visceral endoderm migration. Dev Cell 2024 | Sample collections for the transcriptomics | ArrayExpress E-MTAB-9645 | [ave_migration](src/science_evals/papers/ave_migration/README.md), in progress |
| Cell competition acts as a purifying selection to eliminate cells with mitochondrial defects during early mouse development. Nat Metab 2021 | Embryo dissections, treatments and cell dissociation before single-cell RNA-seq | see the paper's folder | [cell_competition](src/science_evals/papers/cell_competition/README.md), in progress |
| The APOE isoforms differentially shape the transcriptomic and epigenomic landscapes of human microglia xenografted into a mouse model of Alzheimer's disease. Nat Commun 2025 | Second author | GEO GSE271384, GSE271385 | [apoe_microglia](src/science_evals/papers/apoe_microglia/README.md), in progress |

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
