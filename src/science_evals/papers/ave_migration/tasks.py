"""Task specifications: prompts, answer fields, tolerances and notes.

Figures quoted in the notes come from the `facts` dicts computed in reference.py, so no
number in a task is typed by hand. The panel of genes in ave02 and the keyed option of ave04
are also chosen by code from the data.
"""
from __future__ import annotations

from ...spec import Field, Task
from . import inputs as io
from .choices import AVE04_CHOICES, ENRICHED, NULL_BAND

DECISION_REVIEW = "Decision key checked against the data by script. Not yet reviewed by a person."
NUMERIC_REVIEW = "Numeric key, computed from the data by script. No decision field. Not yet reviewed by a person."
FRACTION = "a fraction between 0 and 1"
AVE02_RHO_TOL = 0.25  # window on the Cer1 correlation in ave02


def _pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def ave01(f: dict) -> Task:
    return Task(
        id="ave01_plate_qc",
        title="Which wells of a Smart-seq2 plate hold one good cell",
        paper_anchor="STAR Methods, Transcript quantification, quality control and normalization "
                     "(255 good-quality cells at 5.5 dpc)",
        prompt=f"""Cells from 40 mouse embryos at 5.5 days post coitum were enriched for visceral endoderm (outer cells labelled with a membrane dye; live, labelled cells sorted by FACS into 384 wells) and processed with Smart-seq2 (paired-end 125 bp reads; ERCC spike-in mix added to every well at 1:10^7). You have:

- `{io.COUNTS_E55}`: Salmon gene-level counts (estimated, so not always whole numbers), one column per well and one row per Ensembl gene, plus 92 rows for the ERCC spike-ins (row names starting with `ERCC-`).
- `{io.WELLS_E55}`: one row per well with `well_id` (the column name in the counts file), `well_type` (the well's annotation in the submitters' sample sheet: `single cell`, `no cell` or `multiple cells`) and `read_pairs` (read pairs sequenced from the well).
- `{io.GENES}`: Ensembl gene ID to gene symbol (GRCm38; mitochondrial genes have symbols starting with `mt-`).

Before normalization and clustering, the lab needs the set of wells that each hold one good-quality cell. Decide which wells to keep, and report how many you keep and which they are.""",
        fields=[
            Field("n_wells_kept", "int", "number of wells kept for downstream analysis",
                  tolerance={"value": 20}, digits=0),
            Field("wells_kept", "labels", f"the `well_id` of every well kept, as in the columns of `{io.COUNTS_E55}`",
                  pass_threshold=0.85),
        ],
        data=[io.COUNTS_E55, io.WELLS_E55, io.GENES],
        failure_modes=[
            "qc_control_artifact.context_inappropriate_qc_threshold",
            "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        ],
        category="qc",
        kit="scrna",
        intent="Test whether the agent fits the QC to plate-based full-length data with spike-ins: in a "
               "FACS-sorted Smart-seq2 plate, the ERCC share of the reads, the mapping rate and the number "
               "of genes detected separate wells that hold a cell's RNA from wells that hold little more "
               "than the spike-in mix, and wells annotated as holding no cell or several cells are not single "
               "cells.",
        trap=f"Droplet-style defaults are applied (keep wells with at least 200 detected genes and a "
             f"mitochondrial share under 20%). Of the {f['n_single']} wells annotated as single cells, "
             f"{f['n_removed_single']} hold too little cellular RNA: their median ERCC share of read pairs is "
             f"{_pct(f['removed_median_ercc'])} (against {_pct(f['kept_median_ercc'])} in the kept wells), their "
             f"median mapped fraction {f['removed_median_mapped']:.2f} (against {f['kept_median_mapped']:.2f}) "
             f"and their median number of detected genes {f['removed_median_genes']:,.0f} (against "
             f"{f['kept_median_genes']:,.0f}). {f['removed_passing_droplet']} of them pass the default filters, "
             f"so {f['n_droplet']} wells are kept instead of {f['n_ref']}, {f['droplet_annotated']} wells annotated as no cell or multiple cells among them. No mitochondrial "
             f"filter can catch them: the highest mitochondrial share of read pairs in any well is "
             f"{_pct(f['max_mito_all'])}.",
        tolerance_rationale=f"The windows accept the lab's five thresholds with or without setting aside the "
                            f"wells annotated as no cell or multiple cells ({f['controls_passing']} of them pass "
                            f"the thresholds: {f['controls_passing_types'].get('no cell', 0)} no cell and "
                            f"{f['controls_passing_types'].get('multiple cells', 0)} multiple cells, with "
                            f"profiles like good single cells), "
                            f"the same thresholds on shares of counted reads when read pairs are not used, a gene "
                            f"floor of 2,000, and 3-MAD outlier calls among single-cell wells. They do not accept "
                            f"QC that ignores the spike-ins and the mapping rate: droplet defaults, or a floor of "
                            f"1,000 genes with a mitochondrial cap. The windows are a compromise: a strict floor "
                            f"of 3,000 genes above 10 rpm with depth and mitochondrial cuts but no spike-in or "
                            f"mapping cut keeps {f['strict_floor_no_spikein']} wells and fails only narrowly, "
                            f"and an ERCC cut taken as a share of read pairs without the mapping cut lets most of "
                            f"the poorly mapped wells through (their unmapped reads dilute the spike-in share).",
        not_said="The prompt does not mention the spike-in share, the mapping rate or any threshold, and does "
                 "not say that many wells annotated as single cells are nearly empty or that some wells "
                 "annotated as no cell look like good cells.",
        provenance=f"The paper's criteria (STAR Methods): more than 3,000 genes above 10 reads per million, "
                   f"log10 total reads above 4, mapped fraction above 0.5, mitochondrial fraction below 0.1 and "
                   f"ERCC fraction below 0.3, fractions taken of read pairs. On the public counts with read pairs "
                   f"from the ENA run report (ERP124375, read count / 2) they keep {f['n_ref'] + f['controls_passing']} "
                   f"wells; setting aside the {f['controls_passing']} annotated as no cell or multiple cells leaves the "
                   f"{f['n_ref']} "
                   f"wells the authors marked as passing in the ArrayExpress SDRF, well for well. Of the "
                   f"{f['n_removed_single']} single-cell wells removed, {f['fail_genes']} fail the gene floor, "
                   f"{f['fail_ercc']} the ERCC share and {f['fail_mapped']} the mapping rate. Well annotations are the "
                   f"SDRF's single cell well quality column ({f['n_no_cell']} no cell, {f['n_multi']} multiple "
                   f"cells). Gene symbols: the authors' normalized file lists symbols in Ensembl ID order, checked "
                   f"against the 726 genes it names by ID.",
        key_review=NUMERIC_REVIEW,
    )


def ave02(f: dict, panel: list[str]) -> Task:
    panel_text = ", ".join(panel)
    return Task(
        id="ave02_pseudotime_root",
        title="Which end of the embryonic VE trajectory is the AVE",
        paper_anchor="Fig. 1D and STAR Methods, Diffusion pseudotime analysis of AVE and emVE cells "
                     "(952 genes at 5.5 dpc; high-in-AVE and low-in-AVE groups)",
        prompt=f"""`{io.NORM_E55}` holds Smart-seq2 profiles of the 255 cells from 5.5 dpc mouse embryos that passed quality control (`X`: scran-normalized expression, not log-transformed; variable names are gene symbols). `{io.EMB_VE_E55}` lists the 173 of them that form the embryonic visceral endoderm, the VE over the epiblast. At this stage that population contains the anterior visceral endoderm (AVE), at the distal tip with high Cer1, Hhex and Lefty1, and the rest of the embryonic VE (emVE).

Work on these 173 cells only, with ln(x + 1) of the normalized values. The lab orders them along a diffusion pseudotime that starts at the AVE end of the trajectory, with this pipeline: the 3,000 most highly variable genes among the 173 cells; PCA; a k-nearest-neighbour graph (k = 15) on the first 10 principal components, with the distance sqrt((1 - r) / 2), r being the Spearman correlation between two cells; a diffusion map; diffusion pseudotime (DPT), which needs a root cell.

Then test every gene detected (value above 0) in at least 10 of the 173 cells for a change of its ln expression along pseudotime: regress it on a natural cubic spline of pseudotime with 4 degrees of freedom, F-test against an intercept-only model, Benjamini-Hochberg correction across genes. A gene with FDR < 0.01 is high-in-AVE if its Spearman correlation with pseudotime is negative (it is highest at the start) and low-in-AVE if it is positive.

Report the Spearman correlation between pseudotime and the ln expression of Cer1 across the 173 cells, the fraction of the FDR < 0.01 genes that are high-in-AVE, and which of these genes are high-in-AVE: {panel_text}.""",
        fields=[
            Field("rho_pseudotime_cer1", "float",
                  "Spearman correlation between pseudotime and ln(x + 1) of Cer1 across the 173 cells",
                  tolerance={"value": AVE02_RHO_TOL}, digits=4),
            Field("frac_high_in_ave", "float", f"share of the FDR < 0.01 genes that are high-in-AVE, {FRACTION}",
                  tolerance={"value": 0.09}, digits=4),
            Field("high_in_ave_panel", "labels", "the genes of the listed panel that are high-in-AVE",
                  pass_threshold=0.6),
        ],
        data=[io.NORM_E55, io.EMB_VE_E55],
        failure_modes=["biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error"],
        category="trajectory",
        kit="scrna",
        intent="Test whether the agent anchors the pseudotime in the AVE and checks which population its root "
               "cell belongs to, instead of trusting the extreme of a diffusion component or the first cell "
               "of the file.",
        trap=f"The root is taken at an extreme of the first diffusion component, or at the first cell of the "
             f"file, without checking where it sits. The sign of a diffusion component is arbitrary; here the "
             f"extreme taken and the first cell ({f['first_cell']}) are emVE cells, so the pseudotime runs from "
             f"the emVE to the AVE: the Cer1 correlation is {f['rev_rho_cer1']:+.2f} instead of "
             f"{f['rho_cer1']:+.2f}, {_pct(f['rev_frac_high'])} of the significant genes are called high-in-AVE "
             f"instead of {_pct(f['frac_high'])}, and every panel gene lands in the wrong group, with no error.",
        tolerance_rationale=f"The task is about where the pseudotime starts, not which AVE cell is the root, so "
                            f"the windows accept any root at an AVE extreme of the diffusion map: the AVE end of "
                            f"the first component, the end of the second component in the AVE (the paper's "
                            f"recipe), the cell the authors used ({f['authors_root']}), or the cell with the "
                            f"highest mean z-score of Cer1, Hhex and Lefty1; also the Euclidean distance, genes "
                            f"scaled before the PCA, and a linear, 3-df spline or Spearman test of each gene. They "
                            f"do not accept a pseudotime that starts in the emVE, nor one rooted at an AVE cell "
                            f"that sits among the emVE cells in the graph: used as roots one at a time, "
                            f"{f['n_ave_roots_outside']} of the {f['n_ave']} AVE cells give a Cer1 correlation "
                            f"outside the window (median over all AVE roots {f['median_ave_root_rho']:+.2f}).",
        not_said="The prompt says the pseudotime starts at the AVE end but not how to find that end. It does not "
                 "say that a diffusion component has no fixed sign or that the first cell of the file is an emVE "
                 "cell.",
        provenance=f"The paper roots the pseudotime at the cell with the minimum (5.5 dpc) or maximum (6.25 dpc) "
                   f"of the second diffusion component, which in its run lay in the AVE; the authors' coordinates "
                   f"(code repository, DATA/DC_55.csv) show that their first component also separates AVE from "
                   f"emVE, while a rerun returns the components in another order and sign, so the task asks for "
                   f"the AVE end and not for a component. Cells: the authors' AVE ({f['n_ave']}) and emVE "
                   f"({f['n_cells'] - f['n_ave']}) clusters from the ArrayExpress SDRF, merged and not labelled. "
                   f"With the reference root (the AVE end of the first diffusion component, chosen by its "
                   f"correlation with a Cer1, Hhex and Lefty1 score, rho {f['dc1_marker_rho']:.2f}), the AVE "
                   f"cells have mean pseudotime {f['mean_pt_ave']:.2f} and the emVE cells {f['mean_pt_emve']:.2f}; "
                   f"{f['n_sig']:,} of {f['n_tested']:,} genes have FDR < 0.01 ({f['n_high']} high-in-AVE), against "
                   f"952 in the paper with R gam's loess smoother. Panel: the 5 most negatively and the 5 most "
                   f"positively correlated FDR < 0.01 genes of the reference, Cer1, Hhex and Lefty1 left out.",
        key_review=NUMERIC_REVIEW,
    )


def ave03(f: dict) -> Task:
    return Task(
        id="ave03_ephrin_bootstrap_signs",
        title="Which ephrin and Eph genes carry the pathway's main axes, with PCA signs bootstrapped correctly",
        paper_anchor="Fig. S5A-D and STAR Methods, Analysis of genes in the Ephrin/Eph-signalling pathway "
                     "(14 and 13 genes kept; 6 and 8 significant)",
        prompt=f"""`{io.EPHRIN}` holds the scran-normalized expression (not log-transformed) of the 22 ephrin and Eph receptor genes (Efna1-5, Efnb1-3, Epha1-8, Epha10, Ephb1-4, Ephb6) in two Smart-seq2 datasets of mouse embryos enriched for visceral endoderm: {f['n_cells_e55']} cells at 5.5 dpc and {f['n_cells_e625']} cells at 6.25 dpc (column `stage`), each with its cluster (`cell_type`: AVE, emVE, exVE, Epi or ExE).

The lab wants to know which of these genes drive the main axes of variation of the Ephrin/Eph pathway at each stage. Analyse each stage separately, on ln(x + 1) values:

1. Keep the genes whose mean ln(x + 1) over that stage's cells is above 0.5.
2. Run a PCA on the kept genes (centred, not scaled) and take the loadings of the first 3 principal components.
3. Draw 1,000 bootstrap samples of the cells (with replacement) and rerun the PCA on each. Give each loading of the first 3 components a p-value: the fraction of bootstrap samples in which that loading has the opposite sign to the same loading from the full data.
4. Call a gene significant if at least one of its 3 loadings has p < 0.01.

Report, for each stage, how many genes you kept, how many are significant, and which.""",
        fields=[
            Field("n_genes_kept_e55", "int", "genes kept at 5.5 dpc (mean ln(x + 1) above 0.5)",
                  tolerance={"value": 1}, digits=0),
            Field("n_significant_e55", "int", "significant genes at 5.5 dpc", tolerance={"value": 2}, digits=0),
            Field("significant_genes_e55", "labels", "symbols of the significant genes at 5.5 dpc",
                  pass_threshold=0.7),
            Field("n_genes_kept_e625", "int", "genes kept at 6.25 dpc", tolerance={"value": 1}, digits=0),
            Field("n_significant_e625", "int", "significant genes at 6.25 dpc", tolerance={"value": 2}, digits=0),
            Field("significant_genes_e625", "labels", "symbols of the significant genes at 6.25 dpc",
                  pass_threshold=0.7),
        ],
        data=[io.EPHRIN],
        failure_modes=[
            "metric_statistic_mismatch.null_background_or_control_calibration_error",
            "biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error",
        ],
        category="dimensionality_reduction",
        kit="scrna",
        intent="Test whether the agent knows that a principal axis has no intrinsic sign, so each bootstrap "
               "axis must be aligned with the full-data axis (reflected when it points the other way) before "
               "the signs of its loadings are compared.",
        trap=f"Signs are compared as each PCA returns them. A bootstrap axis comes back pointing either way, so "
             f"every loading flips in a large share of the samples (at 5.5 dpc the median p-value is "
             f"{f['median_p_raw_e55']:.2f} and the smallest {f['min_p_raw_e55']:.2f}) and no gene is significant at "
             f"5.5 dpc, against {f['n_sig_ref_e55']} once each axis is reflected to match the full data. The "
             f"sign test then measures the solver's sign convention, not the stability of the loadings.",
        tolerance_rationale="The windows accept reflection with or without matching bootstrap axes to the "
                            "full-data axes by their loadings (axis reordering, which adds Ephb1 at 5.5 dpc), "
                            "another seed, 2,000 samples, and scikit-learn's PCA with reflection. They do not "
                            "accept unaligned signs, which leave no significant gene at 5.5 dpc with numpy's or "
                            "scikit-learn's PCA (scikit-learn's own sign convention keeps 6 of the 8 genes at "
                            "6.25 dpc, so the 5.5 dpc fields are the ones that catch it).",
        not_said="The prompt defines the p-value as a comparison of signs and does not say that principal axes "
                 "have arbitrary sign or mention reflection.",
        provenance=f"The paper used scanpy's PCA and 1,000 bootstrap samples with 'the fraction of sign "
                   f"inversions ... accounting for possible axis reflections, as described by Peres-Neto et al.', "
                   f"and reports 14 and 13 genes above the expression cut with 6 and 8 significant genes; the "
                   f"reference (numpy SVD, seed 0) reproduces all four numbers. The first 3 components explain "
                   f"{_pct(sum(f['var_top3_e55']))} of the variance at 5.5 dpc and {_pct(sum(f['var_top3_e625']))} "
                   f"at 6.25 dpc. Cells: the authors' labelled cells (the 3 unassigned cells of each stage left "
                   f"out), clusters from the ArrayExpress SDRF; expression from the authors' normalized files.",
        key_review=NUMERIC_REVIEW,
    )


def ave04(f: dict, key: int) -> Task:
    lengths = [len(c) for c in AVE04_CHOICES]
    if lengths[key - 1] in (max(lengths), min(lengths)):
        raise ValueError("ave04: the keyed option must be neither the longest nor the shortest")
    return Task(
        id="ave04_dkk1_timing",
        title="Is Dkk1 an AVE signal before or after migration",
        paper_anchor="Fig. S3C and Results (Dkk1 elevated in the anterior only at 6.25 dpc, arguing against "
                     "DKK1 guiding AVE migration)",
        prompt=f"""Two Smart-seq2 datasets of mouse embryos enriched for visceral endoderm (VE), both after quality control: `{io.NORM_E55}` (5.5 dpc, when the anterior visceral endoderm (AVE) is at the distal tip or migrating towards the future anterior) and `{io.NORM_E625}` (6.25 dpc, when the AVE has reached the boundary between epiblast and extra-embryonic ectoderm and its proximal migration is over). In both, `X` is scran-normalized expression, not log-transformed, with gene symbols as variable names. `{io.LABELS}` gives the stage and cluster of each cell (AVE, emVE, exVE, Epi or ExE); cells that are not listed there were left unassigned and are not used.

Earlier work proposed that DKK1 secreted by the AVE guides its migration, and a collaborator asks whether these data support that. Compare AVE cells with emVE cells, the embryonic VE through which the AVE migrates. For Dkk1, and for Cer1 as a reference AVE gene, report at each stage the AUC of ln(x + 1) expression for AVE against emVE cells: the probability that a randomly chosen AVE cell has higher expression than a randomly chosen emVE cell, ties counted as one half. Then choose the conclusion the data best support.""",
        fields=[
            Field("dkk1_auc_e55", "float", "AUC of Dkk1, AVE against emVE, at 5.5 dpc", tolerance={"value": 0.06}),
            Field("dkk1_auc_e625", "float", "AUC of Dkk1, AVE against emVE, at 6.25 dpc", tolerance={"value": 0.06}),
            Field("cer1_auc_e55", "float", "AUC of Cer1, AVE against emVE, at 5.5 dpc", tolerance={"value": 0.06}),
            Field("cer1_auc_e625", "float", "AUC of Cer1, AVE against emVE, at 6.25 dpc", tolerance={"value": 0.06}),
            Field("conclusion", "choice"),
        ],
        data=[io.NORM_E55, io.NORM_E625, io.LABELS],
        failure_modes=[
            "biological_interpretation.prior_driven_or_overclaimed_biological_interpretation",
            "subset_scope_compartment.wrong_comparator_contrast_or_baseline",
        ],
        category="differential_expression",
        kit="scrna",
        choice_question="Which conclusion do these data best support?",
        choices=list(AVE04_CHOICES),
        correct_choice=key,
        intent="Test whether the agent reads the timing of the AVE enrichment against the timing of the "
               "migration, rather than taking DKK1's published place among AVE signals, or an enrichment "
               "pooled over stages, as support.",
        trap=f"The prior-driven reading: Dkk1 is AVE-enriched at 6.25 dpc (AUC {f['auc_dkk1_e625']:.2f}) and, "
             f"pooled over both stages, still looks AVE-enriched (AUC {f['pooled_dkk1_auc']:.2f}, Wilcoxon "
             f"p = {f['pooled_dkk1_p']:.1e}); DKK1 is a published guidance candidate, so this is read as support. "
             f"At 5.5 dpc, while the AVE migrates, Dkk1 is no higher in the AVE than in the emVE: mean ln "
             f"expression {f['Dkk1_e55_mean_ave']:.2f} against {f['Dkk1_e55_mean_emve']:.2f} (AUC "
             f"{f['auc_dkk1_e55']:.2f}, p = {f['Dkk1_e55_p']:.2f}).",
        tolerance_rationale="AUCs are rank statistics, so the windows (+/- 0.06) accept any monotone transform "
                            "of the expression and any standard AUC routine. They do not accept AUCs from cells "
                            f"pooled over the two stages, which put Dkk1 at {f['pooled_dkk1_auc']:.2f} for both.",
        not_said="The prompt asks for the numbers by stage but does not say which stage bears on the "
                 "hypothesis or what an AUC near 0.5 implies for it.",
        provenance=f"The paper reaches the same timing from its pseudotime fits and HCR (Fig. S3) and "
                   f"concludes that it argues against a role for DKK1 in symmetry breaking and in guiding AVE "
                   f"migration. Key chosen by script from the measured pattern: Dkk1 AUC within 0.5 +/- "
                   f"{NULL_BAND} at 5.5 dpc ({f['auc_dkk1_e55']:.3f}) and at least {ENRICHED} at 6.25 dpc "
                   f"({f['auc_dkk1_e625']:.3f}, p = {f['Dkk1_e625_p']:.1e}); Dkk1 is detected in "
                   f"{_pct(f['Dkk1_e55_detected'])} of the 5.5 dpc AVE and emVE cells (so option 6 does not hold); "
                   f"Cer1 AUC {f['auc_cer1_e55']:.2f} and {f['auc_cer1_e625']:.2f}, while Dkk1 is not enriched at "
                   f"5.5 dpc (so option 3 does not hold) and neither Dkk1 AUC is at or below {1 - ENRICHED:.2f} (so "
                   f"option 5 does not hold). Cells: {f['n_ave_e55']} AVE "
                   f"and {f['n_emve_e55']} emVE at 5.5 dpc, {f['n_ave_e625']} and {f['n_emve_e625']} at 6.25 dpc, "
                   f"clusters from the ArrayExpress SDRF; expression from the authors' normalized files.",
        key_review=DECISION_REVIEW,
    )
