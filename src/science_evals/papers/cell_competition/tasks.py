"""The cell competition tasks: prompt, answer fields, controls and notes.

Each make_* function runs the reference analysis, its variants and its trap on
the data, and returns a Built. Numbers quoted in the notes are computed here.
"""
from __future__ import annotations

import numpy as np

from ...controls import Control
from ...spec import Field, Task
from .. import Built
from . import reference as R
from .prepare import Data

FRACTION = "a fraction between 0 and 1"
NUMERIC_REVIEW = "Numeric key, computed from the data by script. Not yet reviewed by a person. No decision field."

DESIGN = (
    "Mouse embryos were dissected at embryonic day 5.5 and cultured for 16 hours with a pan-caspase "
    "inhibitor (CI) or with vehicle (DMSO). Blocking caspases keeps alive the epiblast cells that would "
    "otherwise be eliminated by apoptosis during this window, so that they can be profiled. The embryos "
    "were dissociated, pooled by treatment within each of 5 processing batches, and single cells were "
    "sequenced with Smart-seq2 with ERCC spike-ins in every well."
)


# --- cc01 -------------------------------------------------------------------------


def make_cc01(d: Data, files: dict) -> Built:
    qc = d.qc
    clusters = d.cells["cluster"].to_numpy()
    ref = R.qc_answer(qc, np.zeros(len(qc), bool))
    trap_flag = R.qc_outliers(qc, 3)
    pooled5 = R.qc_outliers(qc, 5)
    within5 = R.qc_outliers_by(qc, 5, "treatment")
    within3 = R.qc_outliers_by(qc, 3, "treatment")
    by_batch3 = R.qc_outliers_by(qc, 3, "batch")
    mt3 = R.qc_outliers(qc, 3, metrics=("mito",))
    mt10 = qc["frac_reads_mito"].to_numpy(float) > 0.10
    controls = [
        Control("authors_fixed_thresholds_only", "reference",
                what="Keep every cell that passed the fixed thresholds, as the authors did.",
                answer=ref),
        Control("high_mito_outliers_removed", "variant",
                what="Also remove cells more than 3 MADs above the median mitochondrial fraction.",
                answer=R.qc_answer(qc, mt3)),
        Control("fixed_mito_cut_10pct", "variant",
                what="Also remove cells with more than 10% of reads on mitochondrial genes.",
                answer=R.qc_answer(qc, mt10)),
        Control("adaptive_qc_within_treatment_5mad", "variant",
                what="Adaptive outlier rules (5 MADs) computed separately within each treatment group.",
                answer=R.qc_answer(qc, within5)),
        Control("adaptive_qc_pooled_3mad", "trap",
                what="scater-style outlier rules at the default 3 MADs (library size, detected genes, ERCC "
                "and mitochondrial fractions) on all cells pooled.",
                answer=R.qc_answer(qc, trap_flag)),
    ]
    losers = clusters == "loser"
    ci_int = (clusters == "intermediate") & (qc["treatment"].to_numpy() == "CI")
    med = qc.groupby(clusters)[["frac_reads_ercc", "frac_reads_mito", "n_genes_above_10rpm"]].median()
    task = Task(
        id="cc01_cell_qc",
        title="Which cells to keep after first-pass QC when the treatment keeps dying cells alive",
        paper_anchor="Methods, data processing, quality control and normalization; Extended Data Fig. 1a",
        prompt=f"""{DESIGN}

Of 1,495 sequenced cells, 723 passed the lab's first-pass filter: log10 total reads of at least 4, at least 80% of reads mapped, at least 3,000 genes above 10 reads per million, and at least half of all reads on endogenous genes. `cc_qc_metrics.csv` lists these 723 cells, one row per cell: the treatment of the embryos the cell came from (`treatment`, CI or DMSO), the processing batch (`batch`, 1 to 5), log10 total reads, the fraction of reads mapped, the fractions of reads on endogenous genes, on mitochondrial genes and on ERCC spike-ins, and the number of genes above 10 reads per million.

The next steps are to cluster these cells and to characterise the epiblast cells that survive only when apoptosis is blocked. Decide which of the 723 cells should go into that analysis, applying any further quality filtering you judge appropriate. Report how many cells you keep, and how many of the kept cells come from CI-treated and from DMSO-treated embryos.""",
        fields=[
            Field("n_cells_kept", "int", "number of cells you carry into clustering", tolerance={"value": 45}),
            Field("n_ci_cells_kept", "int", "how many of the kept cells come from CI-treated embryos",
                  tolerance={"value": 20}),
            Field("n_dmso_cells_kept", "int", "how many of the kept cells come from DMSO-treated embryos",
                  tolerance={"value": 35}),
        ],
        data=["cc_qc_metrics.csv"],
        failure_modes=["qc_control_artifact.context_inappropriate_qc_threshold"],
        category="qc",
        kit="scrna",
        intent="Test whether the agent checks what a standard adaptive QC rule would remove before applying "
               "it, in an experiment built to keep stressed, dying cells alive, and whether it notices that "
               "the outliers are concentrated in one treatment.",
        trap=(
            f"Adaptive outlier rules at the scater default of 3 MADs are applied to all cells pooled. They "
            f"flag {int(trap_flag.sum())} cells, {int((trap_flag & (qc['treatment'] == 'CI').to_numpy()).sum())} "
            f"of them from CI embryos, and remove {int((trap_flag & losers).sum())} of the {int(losers.sum())} "
            f"cells that the authors' clustering later places in the cluster found only in CI embryos, and "
            f"{int((trap_flag & ci_int).sum())} of the {int(ci_int.sum())} CI cells of the intermediate "
            f"cluster. Those cells passed the fixed thresholds and are the population the experiment was "
            f"built to capture; they carry more ERCC reads (median "
            f"{med.loc['loser', 'frac_reads_ercc']:.3f} against {med.loc['winner', 'frac_reads_ercc']:.3f} "
            f"in winner epiblast cells) and fewer genes above 10 reads per million "
            f"({med.loc['loser', 'n_genes_above_10rpm']:.0f} against {med.loc['winner', 'n_genes_above_10rpm']:.0f}), "
            f"that is, less endogenous RNA per cell. A high mitochondrial fraction is not their signature "
            f"(median {med.loc['loser', 'frac_reads_mito']:.3f} against {med.loc['winner', 'frac_reads_mito']:.3f}), "
            f"so a mitochondrial filter does not remove them. At 5 MADs on pooled cells the rules still remove "
            f"{int((pooled5 & losers).sum())} of the {int(losers.sum())}."
        ),
        tolerance_rationale=(
            f"The windows accept keeping all 723 cells, removing high-mitochondrial outliers (3 MADs, or a "
            f"fixed 10% cut), and adaptive rules computed within each treatment group at 5 MADs. They do not "
            f"accept adaptive rules that pool the treatments: at 5 MADs on pooled cells "
            f"{R.qc_answer(qc, pooled5)['n_ci_cells_kept']} CI cells are kept, within each processing batch at "
            f"3 MADs {R.qc_answer(qc, by_batch3)['n_ci_cells_kept']}, and within each treatment at 3 MADs "
            f"{R.qc_answer(qc, within3)['n_ci_cells_kept']}, against {ref['n_ci_cells_kept']} in the key."
        ),
        not_said="The prompt does not say that the QC metrics differ between treatments, or that the cells "
                 "of interest hold less endogenous RNA.",
        provenance="The authors kept all 723 cells that passed the fixed thresholds (Methods). The QC table is "
                   "the authors' QC_start object, with treatment and batch from meta_dati, both in "
                   "input_for_figures.Rdata of their figure repository. Cluster labels used in these notes "
                   "come from meta_all.txt in the same repository.",
        key_review=NUMERIC_REVIEW,
    )
    return Built(task, ref, controls, [files["cc_qc_metrics.csv"]])


# --- cc03 -------------------------------------------------------------------------


def _rnr1_inputs(d: Data):
    from .prepare import ci_epiblast_cells

    e = ci_epiblast_cells(d)
    idx = np.array([d.allele["cells"].index(c) for c in e.index])
    depth, h = R.heteroplasmy(d.allele["counts"][idx])
    return e, depth, h, d.allele["positions"]


def make_cc03(d: Data, files: dict) -> Built:
    e, depth, h, pos = _rnr1_inputs(d)
    clu = e["cluster"].to_numpy()
    ref = R.rnr1_fractions(depth, h, pos, clu)
    trap = R.rnr1_fractions(depth, h, pos, clu, mode="nan_false", min_positions=None)
    chk = R.marker_agreement(d.allele["counts"], d.allele_capped["counts"], d.allele["ref"], pos)
    zero = R.rnr1_fractions(depth, h, pos, clu, mode="zero_fill", min_positions=None)
    controls = [
        Control("all_six_covered", "reference",
                what="Keep cells covered by more than 50 reads at all six positions (after the paper's "
                     "filter of more than 2,000 covered mtDNA positions per cell), average heteroplasmy "
                     "over the six, and divide by the cells kept in each cluster.",
                answer=ref),
        Control("all_six_covered_no_cell_filter", "variant",
                what="The same without the per-cell filter on covered mtDNA positions.",
                answer=R.rnr1_fractions(depth, h, pos, clu, min_positions=None)),
        Control("at_least_50_reads", "variant",
                what="Count a position as covered at 50 reads or more instead of more than 50.",
                answer=R.rnr1_fractions(depth, h, pos, clu, strict=False)),
        Control("mean_over_covered_positions", "variant",
                what="Average over whichever of the six positions are covered in each cell.",
                answer=R.rnr1_fractions(depth, h, pos, clu, mode="available")),
        Control("uncovered_cells_kept_in_denominator", "trap",
                what="Positions with 50 reads or fewer are set to missing; a cell's mean over the six is "
                     "then missing, compares as False against 0.01, and the cell stays in the denominator.",
                answer=trap),
    ]
    cols = [int(np.flatnonzero(pos == p)[0]) for p in R.RNR1_POSITIONS]
    covered = (depth[:, cols] > 50).all(axis=1)
    lost = {k: int(((clu == k) & ~covered).sum()) for k in R.CLUSTERS}
    total = {k: int((clu == k).sum()) for k in R.CLUSTERS}
    task = Task(
        id="cc03_rnr1_heteroplasmy",
        title="How common is mt-Rnr1 heteroplasmy in winner, intermediate and loser epiblast cells",
        paper_anchor="Fig. 6h and legend; Methods, analysis of heteroplasmy in single-cell RNA-seq",
        prompt=f"""{DESIGN}

The epiblast cells from CI-treated embryos fall into three clusters along a trajectory: winner cells, which are also common in DMSO-treated embryos, loser cells, found only when apoptosis is blocked, and intermediate cells between them. The lab asks whether mitochondrial DNA variants are more common toward the loser end.

`cc_mt_allele_counts.csv` gives, for each of the 274 epiblast cells from CI-treated embryos, the number of reads supporting A, C, G and T at every mtDNA position (columns `<position>_<base>`, positions 1 to 16299 of the mm10 mtDNA; reads uniquely mapped to the mtDNA, base quality above 33). `cc_ci_epiblast_cells.csv` gives each cell's cluster (`winner`, `intermediate` or `loser`), batch and losing score.

Heteroplasmy at a position is 1 minus the frequency of the most common base at that position in that cell, and it can only be estimated where the position is covered by more than 50 reads in that cell. Take the six mt-Rnr1 positions 300, 303, 304, 305, 326 and 327. For each cluster, report the fraction of cells whose mean heteroplasmy over these six positions is above 0.01, and report the total number of cells that your three fractions are based on.""",
        fields=[
            Field("frac_heteroplasmic_winner", "float",
                  f"share of winner cells with mean heteroplasmy over the six positions above 0.01, {FRACTION}",
                  tolerance={"value": 0.05}, digits=3),
            Field("frac_heteroplasmic_intermediate", "float",
                  f"the same share among intermediate cells, {FRACTION}", tolerance={"value": 0.06}, digits=3),
            Field("frac_heteroplasmic_loser", "float",
                  f"the same share among loser cells, {FRACTION}", tolerance={"value": 0.08}, digits=3),
            Field("n_cells_used", "int", "number of cells, all three clusters together, behind the three fractions",
                  tolerance={"value": 20}),
        ],
        data=["cc_mt_allele_counts.csv", "cc_ci_epiblast_cells.csv"],
        failure_modes=[
            "marker_reference_label_scope.dropout_panel_limit_or_absence_overclaim",
            "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        ],
        category="variant_calling",
        kit="scrna",
        intent="Test whether the agent keeps cells without enough reads out of a per-cell variant summary "
               "instead of counting them as cells without the variant.",
        trap=(
            f"Heteroplasmy is masked where coverage is too low, as the prompt says, but the per-cell mean is "
            f"then missing for cells without full coverage, a missing value compares as False against 0.01 "
            f"without any error, and those cells stay in their cluster's denominator: no data is read as no "
            f"variant. Coverage is incomplete for {lost['winner']} of {total['winner']} winner, "
            f"{lost['intermediate']} of {total['intermediate']} intermediate and {lost['loser']} of "
            f"{total['loser']} loser cells, so the dilution falls where heteroplasmy is common: the loser "
            f"fraction drops from {ref['frac_heteroplasmic_loser']:.2f} to "
            f"{trap['frac_heteroplasmic_loser']:.2f}, the intermediate fraction from "
            f"{ref['frac_heteroplasmic_intermediate']:.2f} to {trap['frac_heteroplasmic_intermediate']:.2f}, "
            f"and the cell count from {ref['n_cells_used']} to {trap['n_cells_used']}."
        ),
        tolerance_rationale=(
            f"The windows accept the paper's per-cell filter on covered mtDNA positions or none, a 50-read cut "
            f"taken as at least 50, and averaging over the covered subset of the six positions. They do not "
            f"accept keeping cells without full coverage in the denominator, whether their mean is missing or "
            f"their uncovered positions are set to 0 (intermediate {zero['frac_heteroplasmic_intermediate']:.2f}, "
            f"loser {zero['frac_heteroplasmic_loser']:.2f}, {zero['n_cells_used']} cells)."
        ),
        not_said="The prompt does not say what to do with cells that lack coverage at some of the six "
                 "positions, or that coverage differs between clusters.",
        provenance=(
            f"The paper reports about 5%, 55% and 87% (Fig. 6h, 215 cells). The allele counts are "
            f"raw_counts_allele_mt.csv from E-MTAB-8640, the table behind Fig. 6. Its rows are pileup files "
            f"named by older sequencing-file names, joined to cells through the authors' own name tables "
            f"(name_full.Rda and provo.Rda), as their figure code does. The join was checked by script: the "
            f"cells carrying inherited non-reference alleles at mtDNA positions 12217 and 13536 are the same "
            f"in this table and in the independent re-count shipped with the authors' MitoHEAR package "
            f"({chk['agree']} of {chk['cells']} cell-position calls agree, {chk['carriers']} carriers). "
            f"Cluster labels and losing scores are the authors' (MitoHEAR after_qc.Rda)."
        ),
        key_review=NUMERIC_REVIEW,
    )
    return Built(task, ref, controls, [files["cc_mt_allele_counts.csv"], files["cc_ci_epiblast_cells.csv"]])


# --- shared context ----------------------------------------------------------------


def context(d: Data) -> dict:
    """Trajectory and trend tests used by cc02, cc04 and the anchors."""
    from .prepare import ci_epiblast_cells, epiblast_cells
    from .trajectory import spearman

    epi = epiblast_cells(d)
    traj = R.epiblast_trajectory(d.tpm, d.hvg, epi)
    ci = ci_epiblast_cells(d)
    if list(ci.index) != traj["ci_cells"]:
        raise ValueError("CI epiblast cell order differs between tables")
    traj["rho_with_authors"] = spearman(traj["score"], ci["losing_score"].to_numpy())
    tpm_ci = d.tpm[list(ci.index)]
    tg = R.trend_genes(tpm_ci, ci["losing_score"].to_numpy(), R.expressed(tpm_ci, "gt15_in_gt10"))
    trend = {"n_tested": len(tg["ids"]), "n_de": len(tg["down"]) + len(tg["up"]),
             "n_down": len(tg["down"]), "n_up": len(tg["up"]), "genes": tg}
    return {"traj": traj, "trend": trend, "tpm_ci": tpm_ci, "ci": ci}


def _symbol_ids(d: Data) -> dict:
    """Gene symbol -> Ensembl id, for ids present in the TPM table."""
    present = set(d.tpm.index)
    out = {}
    for gid, sym in zip(d.genes["gene_id"], d.genes["gene_name"]):
        if gid in present and sym not in out:
            out[sym] = gid
    return out


# --- cc02 -------------------------------------------------------------------------

CC02_GENES = (("rho_samm50", "Samm50"), ("rho_mt_atp6", "mt-Atp6"), ("rho_atf4", "Atf4"), ("rho_sesn2", "Sesn2"))


def make_cc02(d: Data, files: dict, ctx: dict) -> Built:
    traj, ci = ctx["traj"], ctx["ci"]
    sym = _symbol_ids(d)
    genes = {field: sym[name] for field, name in CC02_GENES}
    lx = R.log_tpm(d.tpm[traj["ci_cells"]])

    def ans(score) -> dict:
        return R.gene_correlations(np.asarray(score, float), lx, genes)

    def end_nbrs(i_ci: int) -> int:
        return traj["dmso_nbrs"](traj["ci_idx"][i_ci])

    ref = ans(traj["score"])
    trap = ans(traj["reversed"])
    hvg_ci = traj["lx_hvg"][traj["ci_cells"]].to_numpy().T
    tpm_ci = d.tpm[traj["ci_cells"]]
    detected = tpm_ci[(tpm_ci > 0).sum(axis=1) >= 10]
    controls = [
        Control("destiny_style_dpt_dmso_oriented", "reference",
                what="Rebuild the authors' diffusion map of CI epiblast cells (Spearman distance over the lab's "
                     "variable genes) and root the pseudotime at the tip whose nearest neighbours include the "
                     "most cells from DMSO-treated embryos.",
                answer=ref),
        Control("authors_published_losing_score", "variant",
                what="Use the authors' own losing score for these cells.",
                answer=ans(ci["losing_score"].to_numpy())),
        Control("map_on_all_epiblast_cells", "variant",
                what="Build the diffusion map on all 512 epiblast cells, root it the same way, keep the CI cells.",
                answer=ans(traj["score_all_epiblast_map"])),
        Control("scanpy_pca_knn_dpt", "variant",
                what="scanpy: PCA on the lab's variable genes, 15-nearest-neighbour graph, diffmap and dpt, "
                     "rooted at the end of the first diffusion component next to DMSO cells.",
                answer=ans(R.scanpy_dpt(hvg_ci, end_nbrs))),
        Control("scanpy_own_variable_genes", "variant",
                what="The same scanpy pipeline on its own 2,000 variable genes from ln(TPM + 1) of the CI cells.",
                answer=ans(R.scanpy_own_hvg_dpt(d.tpm[traj["ci_cells"]], end_nbrs))),
        Control("scanpy_no_gene_selection", "variant",
                what="The same scanpy pipeline on ln(TPM + 1) of every gene detected in at least 10 CI cells, "
                     "with no variable-gene selection.",
                answer=ans(R.scanpy_dpt(np.log1p(detected.to_numpy().T), end_nbrs))),
        Control("root_at_ci_only_end", "trap",
                what="Root the pseudotime at the other end of the trajectory, the tip with no DMSO cell among "
                     "its nearest neighbours.",
                answer=trap),
    ]
    nbrs = traj["tip_dmso_neighbours"]
    tol = 0.3
    shift = max(abs(c.answer[k] - ref[k]) for c in controls if c.kind == "variant" for k in ref)
    outside = min(abs(trap[k] - ref[k]) - tol for k in ref)
    task = Task(
        id="cc02_losing_score_direction",
        title="Which way the winner-to-loser trajectory runs",
        paper_anchor="Fig. 2a, 2d-e, Fig. 3c and Extended Data Fig. 4d; Methods, identification of a "
                     "single-cell trajectory in the epiblast",
        prompt=f"""{DESIGN}

`cc_epiblast_tpm.csv` holds TPM values (one row per gene, Ensembl gene id in `gene_id`, one column per cell) for the 512 epiblast cells of the experiment, and `cc_epiblast_cells.csv` gives each cell's treatment (`CI` or `DMSO`) and batch. `cc_genes.csv` gives gene symbols; its `authors_hvg` column marks the 1,921 highly variable genes the lab selected on all cells.

In CI-treated embryos the epiblast cells form one continuum, from a state that is also common in DMSO-treated embryos to states seen only when apoptosis is blocked. Using the 274 CI-treated epiblast cells, order the cells along this continuum by diffusion pseudotime and express it as a losing score: 0 at the end where cells survive without the inhibitor, increasing toward the cells that would otherwise have been eliminated. Report the Spearman correlation, across the 274 CI-treated epiblast cells, between your losing score and log10(TPM + 1) of each of Samm50, mt-Atp6, Atf4 and Sesn2.""",
        fields=[
            Field("rho_samm50", "float", "Spearman correlation of the losing score with Samm50",
                  tolerance={"value": tol}, digits=3),
            Field("rho_mt_atp6", "float", "Spearman correlation of the losing score with mt-Atp6",
                  tolerance={"value": tol}, digits=3),
            Field("rho_atf4", "float", "Spearman correlation of the losing score with Atf4",
                  tolerance={"value": tol}, digits=3),
            Field("rho_sesn2", "float", "Spearman correlation of the losing score with Sesn2",
                  tolerance={"value": tol}, digits=3),
        ],
        data=["cc_epiblast_tpm.csv", "cc_epiblast_cells.csv", "cc_genes.csv"],
        failure_modes=["biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error"],
        category="trajectory",
        kit="scrna",
        intent="Test whether the agent decides which end of a pseudotime is its start from the experimental "
               "design (where cells from untreated embryos fall), instead of leaving it to the arbitrary sign "
               "of a diffusion component.",
        trap=(
            f"Diffusion pseudotime is rooted wherever the first diffusion component happens to put its "
            f"extreme. The tips of the CI map have {nbrs[0]}, {nbrs[1]} and {nbrs[2]} DMSO cells among their "
            f"{R.N_DMSO_NEIGHBOURS} nearest neighbours; rooting at the tip with none reverses the score "
            f"(Spearman {R.spearman(traj['reversed'], ci['losing_score'].to_numpy()):.3f} with the authors' "
            f"score) and every sign: Samm50 {trap['rho_samm50']:+.2f} instead of {ref['rho_samm50']:+.2f}, "
            f"Atf4 {trap['rho_atf4']:+.2f} instead of {ref['rho_atf4']:+.2f}. Read that way, mitochondrial "
            f"genes would rise and stress genes would fall as cells approach elimination."
        ),
        tolerance_rationale=(
            f"Windows of +/-{tol} accept the authors' own score, a map built on all epiblast cells, and scanpy "
            f"diffusion pseudotime on the lab's variable genes, on its own, or with no gene selection (largest "
            f"shift {shift:.2f}): the point is the direction, not the strength of each correlation. A reversed "
            f"score lands at least {outside:.2f} outside every window."
        ),
        not_said="The prompt does not say how to tell the two ends of the trajectory apart.",
        provenance=(
            f"The paper roots destiny's DPT in the winner cluster of CI epiblast cells and calls it the losing "
            f"score. The rebuild here (destiny-style kernel on Spearman distances over the lab's 1,921 variable "
            f"genes) gives Spearman {traj['rho_with_authors']:.4f} with the authors' published score (MitoHEAR "
            f"after_qc.Rda). Samm50 and mt-Atp6 fall along the trajectory in Fig. 3c; Atf4 and Sesn2 rise in "
            f"Extended Data Fig. 4d."
        ),
        key_review=NUMERIC_REVIEW,
    )
    return Built(task, ref, controls,
                 [files["cc_epiblast_tpm.csv"], files["cc_epiblast_cells.csv"], files["cc_genes.csv"]])


# --- cc04 -------------------------------------------------------------------------


def make_cc04(d: Data, files: dict, ctx: dict) -> Built:
    import re

    ci, tpm_ci = ctx["ci"], ctx["tpm_ci"]
    score = ci["losing_score"].to_numpy(float)
    pat = re.compile(R.RESP_CHAIN_PATTERN)
    id2sym = dict(zip(d.genes["gene_id"], d.genes["gene_name"]))

    def in_set(ids) -> np.ndarray:
        return np.array([bool(pat.match(str(id2sym.get(g, "")))) for g in ids])

    def run(rule: str, test: str = "linear") -> dict:
        tg = R.trend_genes(tpm_ci, score, R.expressed(tpm_ci, rule), test=test)
        return R.enrichment_answer(tg["ids"], tg["down"], in_set(tg["ids"]))

    tg = ctx["trend"]["genes"]
    ref = R.enrichment_answer(tg["ids"], tg["down"], in_set(tg["ids"]))
    genome = d.genes["gene_id"].to_numpy()
    trap = R.enrichment_answer(genome, tg["down"], in_set(genome))
    controls = [
        Control("tested_genes_background", "reference",
                what="Test genes with more than 15 TPM in more than 10 cells (the paper's filter) for a linear "
                     "trend on the losing score; use those tested genes as the background.",
                answer=ref),
        Control("spearman_trend_test", "variant",
                what="Same genes and background, Spearman correlation test instead of the linear trend.",
                answer=run("gt15_in_gt10", "spearman")),
        Control("lighter_expression_filter", "variant",
                what="Test genes above 1 TPM in at least 10% of cells; background = those genes.",
                answer=run("gt1_in_10pct")),
        Control("stricter_expression_filter", "variant",
                what="Test genes with mean TPM of at least 10; background = those genes.",
                answer=run("mean_ge10")),
        Control("whole_genome_background", "trap",
                what="Same decreasing genes, but the background is every gene in the annotation, as web "
                     "enrichment tools assume by default.",
                answer=trap),
    ]
    loose = run("detected_in_3")
    n_rc_tested = int(in_set(tg["ids"]).sum())
    n_rc_genome = int(in_set(genome).sum())
    task = Task(
        id="cc04_resp_chain_enrichment",
        title="Are respiratory chain genes over-represented among genes that fall along the losing score",
        paper_anchor="Fig. 3a-b; Methods, differential gene expression analysis along the trajectory",
        prompt=f"""{DESIGN}

`cc_epiblast_tpm.csv` holds TPM values (one row per gene, Ensembl gene id in `gene_id`, one column per cell) for the epiblast cells. `cc_ci_epiblast_cells.csv` lists the 274 epiblast cells from CI-treated embryos with the lab's losing score, a diffusion pseudotime that is 0 in cells common to both treatments and highest in cells found only when apoptosis is blocked. `cc_genes.csv` is the Ensembl gene annotation (gene id, symbol, chromosome, start, end).

The lab wants to know whether cells switch off the respiratory chain as they approach elimination. Using the 274 CI-treated epiblast cells, call a gene decreasing if a test of trend of log10(TPM + 1) against the losing score gives a negative slope (or correlation) at a Benjamini-Hochberg FDR below 0.01. Then test, with a one-sided Fisher exact test, whether respiratory chain and ATP synthase genes are over-represented among the decreasing genes. Define that gene set as the gene symbols matching the regular expression `{R.RESP_CHAIN_PATTERN}`.

Report: the share of the genes in your test's background that you call decreasing; the share of the respiratory chain genes in that background that you call decreasing; the fold enrichment (the second share divided by the first); and the one-sided Fisher p-value as -log10(p).""",
        fields=[
            Field("frac_genes_decreasing", "float",
                  f"share of all background genes called decreasing, {FRACTION}", tolerance={"value": 0.15},
                  digits=3),
            Field("frac_resp_chain_decreasing", "float",
                  f"share of respiratory chain genes in the background called decreasing, {FRACTION}",
                  tolerance={"value": 0.08}, digits=3),
            Field("fold_enrichment", "float", "the second share divided by the first",
                  tolerance={"value": 1.2}, digits=2),
            Field("neg_log10_p", "float", "-log10 of the one-sided Fisher exact test p-value",
                  tolerance={"value": 15}, digits=1),
        ],
        data=["cc_epiblast_tpm.csv", "cc_ci_epiblast_cells.csv", "cc_genes.csv"],
        failure_modes=["metric_statistic_mismatch.null_background_or_control_calibration_error"],
        category="pathway_enrichment",
        kit="scrna",
        intent="Test whether the agent takes the background of an over-representation test from the genes "
               "that could have been called, not from the whole genome.",
        trap=(
            f"The decreasing genes are compared with every gene in the annotation ({len(genome):,} genes) "
            f"instead of the {len(tg['ids']):,} genes that were expressed enough to be tested. Most of the "
            f"genome could never have been called, so the share of background genes called decreasing falls "
            f"from {ref['frac_genes_decreasing']:.3f} to {trap['frac_genes_decreasing']:.3f}, the fold "
            f"enrichment rises from {ref['fold_enrichment']:.2f} to {trap['fold_enrichment']:.2f}, and -log10 p "
            f"from {ref['neg_log10_p']:.0f} to {trap['neg_log10_p']:.0f}. The respiratory chain result holds "
            f"either way; its strength is what the background inflates."
        ),
        tolerance_rationale=(
            f"The windows accept a Spearman trend test and expression filters from about 1 TPM in 10% of "
            f"cells to a mean of 10 TPM, with the tested genes as the background. They do not accept a "
            f"background of the whole annotation, nor testing every gene detected in at least 3 cells and "
            f"using those as the background (share decreasing {loose['frac_genes_decreasing']:.3f}, fold "
            f"enrichment {loose['fold_enrichment']:.2f}, -log10 p {loose['neg_log10_p']:.0f})."
        ),
        not_said="The prompt does not say which genes to test or which genes form the background.",
        provenance=(
            f"The paper ran IPA, a commercial tool, on the genes changing along the trajectory with all tested "
            f"genes (more than 15 TPM in more than 10 cells) as the background, and reported Mitochondrial "
            f"Dysfunction (-log10 p = 21.1) and Oxidative Phosphorylation (-log10 p = 18.6) as the top "
            f"pathways. This task uses an open, explicitly defined gene set ({n_rc_tested} genes among the "
            f"tested genes, {n_rc_genome} in the annotation) and Fisher's exact test; it does not reproduce "
            f"IPA. With the paper's filter, {ctx['trend']['n_tested']:,} genes are tested and "
            f"{ctx['trend']['n_de']:,} change at FDR < 0.01 (the paper reports 5,311). The losing score is "
            f"the authors' (MitoHEAR after_qc.Rda)."
        ),
        key_review=NUMERIC_REVIEW,
    )
    return Built(task, ref, controls,
                 [files["cc_epiblast_tpm.csv"], files["cc_ci_epiblast_cells.csv"], files["cc_genes.csv"]])
