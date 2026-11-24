"""The four tasks: prompt, answer fields, tolerances, notes, and their scripted controls."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...controls import Control as _Control
from ...spec import Field, Task
from .. import Built
from . import inputs as I
from . import reference as R

RNA_FILES = [I.RNA_COUNTS, I.LENGTHS, I.SAMPLES, I.GENES]
RNA_CONTEXT = f"""Human iPSC-derived microglia from 4 isogenic lines were transplanted into the brains of AppNL-G-F mice: lines carrying a single APOE2, APOE3 or APOE4 allele, and an APOE knockout line (APOE-KO). After 12 months the human microglia were sorted from each mouse and profiled by bulk RNA-seq. Reads were quantified against GRCh38 transcripts with Salmon (with its GC-bias correction) and summarised to genes with tximport. You have:

- `{I.RNA_COUNTS}`: gene-level estimated counts as deposited in GEO (rows: Entrez gene IDs; columns: samples)
- `{I.LENGTHS}`: the matching tximport average transcript lengths (same rows and columns)
- `{I.SAMPLES}`: the genotype group, iPSC line and GEO accession of each sample
- `{I.GENES}`: gene symbols for the Entrez IDs that have one"""

TEST_RULES = ("fit one negative binomial GLM (DESeq2 or an equivalent) with genotype as the only factor, "
              "{test}, and adjust the p-values with Benjamini-Hochberg over every gene that receives a "
              "p-value (no independent filtering). A gene is differentially expressed at FDR < 0.05.")


def Control(name: str, kind: str, what: str, answer: dict) -> _Control:
    return _Control(name=name, kind=kind, answer=answer, what=what)


def _f(x, d=3):
    return f"{x:.{d}f}"


def _pct(x):
    return f"{100 * x:.0f}%"


# --- apoe01: failed knockouts ---------------------------------------------------


def t1_failed_knockouts(rna: R.Rna, files: dict, cache: Path) -> Built:
    ko_all = rna.columns(["KO"])
    failed = R.failed_knockouts(rna)
    every = list(rna.counts.columns)
    kept = [c for c in every if c not in failed]
    up = "n_higher_in_ko"

    def ans(r):
        return R.de_summary(r, up)

    ref_r = R.de(rna, kept, "KO", "E3", "offsets")
    ref = ans(ref_r)
    trap_r = R.de(rna, every, "KO", "E3", "offsets")
    controls = [
        Control("reference", "reference",
                "Drop the knockout grafts whose APOE level is at least a quarter of the APOE3 median, then "
                "DESeq2 with the tximport lengths as offsets, BH over all genes.", ref),
        Control("plain_counts", "variant", "Same samples; counts without the transcript lengths.",
                ans(R.de(rna, kept, "KO", "E3", "counts"))),
        Control("length_scaled_counts", "variant", "Same samples; lengths divided out of the counts first.",
                ans(R.de(rna, kept, "KO", "E3", "length_scaled"))),
        Control("independent_filtering", "variant",
                "Same samples and offsets; DESeq2's padj with independent filtering instead of BH over all genes.",
                ans(R.de(rna, kept, "KO", "E3", "offsets", adjust="independent_filtering"))),
        Control("pydeseq2", "variant", "Same samples; the pydeseq2 package on the counts.",
                ans(R.pydeseq2_de(rna, kept, "KO", "E3", "counts", cache))),
        Control("all_knockout_grafts", "trap",
                "Keep all five APOE-KO samples, two of which express APOE; otherwise the reference analysis.",
                ans(trap_r)),
    ]
    cpm = R.apoe_cpm(rna)
    e3 = cpm[rna.columns(["E3"])].median()
    rel = {c: cpm[c] / e3 for c in ko_all}
    good = [c for c in ko_all if c not in failed]
    sig = ref_r["fdr"] < 0.05
    task = Task(
        id="apoe01_failed_knockouts",
        title="APOE knockout against APOE3: confirm the knockout before comparing",
        paper_anchor="Fig. 1b and Supplementary Figs. 2-3 (APOE-KO against APOE3); Methods, RNA-seq QC and processing",
        prompt=RNA_CONTEXT + "\n\nThe lab wants the genes that change when microglia lose APOE. Compare APOE-KO "
        "microglia with APOE3 microglia, with all four genotype groups in one model: "
        + TEST_RULES.format(test="test APOE-KO against APOE3 with a Wald test")
        + " Report the number of differentially expressed genes, how many of them are higher in APOE-KO, and the "
        "log2 fold change of APOE (APOE-KO over APOE3).",
        fields=[
            Field("n_deg", "int", "number of genes with FDR < 0.05 for APOE-KO against APOE3",
                  tolerance={"lower": 45, "upper": 75}, digits=0),
            Field(up, "int", "how many of those genes have a positive log2 fold change (higher in APOE-KO)",
                  tolerance={"lower": 25, "upper": 35}, digits=0),
            Field("apoe_log2fc", "float", "log2 fold change of APOE (Entrez 348), APOE-KO over APOE3",
                  tolerance={"value": 0.4}, digits=3),
        ],
        data=list(RNA_FILES),
        failure_modes=[
            "assay_intervention_feature_linkage.barcode_feature_observation_or_metadata_linkage_error",
            "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        ],
        category="differential_expression",
        kit="bulk_rnaseq",
        intent="Test whether the agent checks that the knockout grafts lack APOE before it uses them as the "
               "knockout group.",
        trap=(f"All five APOE-KO samples are used. {' and '.join(failed)} carry APOE at "
              f"{' and '.join(_pct(rel[c]) for c in failed)} of the APOE3 median (counts per million), while "
              f"{', '.join(good)} sit at {_pct(min(rel[c] for c in good))} to {_pct(max(rel[c] for c in good))}; "
              f"they are not knockouts in this experiment. With them the APOE log2 fold change shrinks from "
              f"{_f(ref['apoe_log2fc'], 2)} to {_f(trap_r.loc[R.APOE, 'log2FoldChange'], 2)}, the number of "
              f"differentially expressed genes rises from {ref['n_deg']} to {int((trap_r['fdr'] < 0.05).sum())}, and "
              f"the number higher in APOE-KO from {ref[up]} to {int(((trap_r['fdr'] < 0.05) & (trap_r['log2FoldChange'] > 0)).sum())}."),
        tolerance_rationale="The windows accept the tximport lengths as offsets, plain counts, length-scaled "
                            "counts, pydeseq2 in place of DESeq2, and DESeq2's independent filtering in place of "
                            "BH over all genes. They do not accept keeping the two knockout grafts that express "
                            "APOE.",
        not_said="The prompt does not say that any graft failed or that the knockout should be confirmed from "
                 "the data. The sample sheet lists all five knockout samples under the same line.",
        provenance=(f"The paper kept the 3 knockout samples in which APOE was significantly lower (Fig. 1b, "
                    f"Supplementary Fig. 2) and reports 126 genes up and 206 down in APOE-KO against APOE3. The "
                    f"reference here, with the tximport lengths as DESeq2 offsets, finds {ref['n_deg']} "
                    f"({ref[up]} up, {int((sig & (ref_r['log2FoldChange'] < 0)).sum())} down) over every gene "
                    f"with a p-value; the authors' Supplementary Data 2 lists 335 (127 up, 208 down) over the "
                    f"genes that have a symbol. Counts are from GEO GSE271385; the lengths are the avgTxLength "
                    f"assay of the DESeqDataSet in the authors' repository."),
    )
    return Built(task, ref, controls, [files[n] for n in RNA_FILES])


# --- apoe02: tximport transcript lengths ------------------------------------------


def t2_length_offsets(rna: R.Rna, files: dict, cache: Path) -> Built:
    cols = rna.columns(["E2", "E3", "E4"])
    up = "n_higher_in_apoe4"

    def ans(r):
        return R.de_summary(r, up)

    ref_r = R.de(rna, cols, "E4", "E3", "offsets")
    trap_r = R.de(rna, cols, "E4", "E3", "counts")
    ref = ans(ref_r)
    paper_cols = [c for c in rna.counts.columns if c not in R.failed_knockouts(rna)]
    paper_r = R.de(rna, paper_cols, "E4", "E3", "offsets")
    controls = [
        Control("reference", "reference",
                "APOE2, APOE3 and APOE4 samples; DESeq2 with the tximport lengths as offsets; BH over all genes.",
                ref),
        Control("length_scaled_counts", "variant", "Lengths divided out of the counts, then DESeq2 with size factors.",
                ans(R.de(rna, cols, "E4", "E3", "length_scaled"))),
        Control("pydeseq2_length_scaled", "variant", "Length-scaled counts with the pydeseq2 package.",
                ans(R.pydeseq2_de(rna, cols, "E4", "E3", "length_scaled", cache))),
        Control("independent_filtering", "variant", "Offsets; DESeq2's padj with independent filtering.",
                ans(R.de(rna, cols, "E4", "E3", "offsets", adjust="independent_filtering"))),
        Control("paper_model", "variant",
                "The paper's model: offsets, all samples except the two knockout grafts that express APOE.",
                ans(paper_r)),
        Control("counts_without_lengths", "trap",
                "The gene-level counts tested as plain read counts, size factors only.", ans(trap_r)),
    ]
    both = (ref_r["fdr"] < 0.05) & (trap_r["fdr"] < 0.05)
    only_trap = (trap_r["fdr"] < 0.05) & ~(ref_r["fdr"] < 0.05)
    from scipy.stats import ttest_ind

    log_len = np.log2(rna.lengths)
    e3 = rna.columns(["E3"])
    ratio_e4 = log_len[rna.columns(["E4"])].mean(1) - log_len[e3].mean(1)
    expressed = rna.counts[cols].mean(1) >= 10
    shifted = {g: int(np.nansum(ttest_ind(log_len.loc[expressed, rna.columns([g])], log_len.loc[expressed, e3],
                                          axis=1).pvalue < 0.001)) for g in ("E2", "E4")}
    agree = np.mean(np.sign(trap_r.loc[only_trap, "log2FoldChange"]) == np.sign(ratio_e4[only_trap]))
    task = Task(
        id="apoe02_transcript_length_offsets",
        title="APOE4 against APOE3 from tximport counts",
        paper_anchor="Fig. 1d and Fig. 2b (APOE4 against APOE3); Methods, RNA-seq QC and processing",
        prompt=RNA_CONTEXT + "\n\nA colleague wants the list of genes that differ between APOE4 and APOE3 "
        "microglia. Use the APOE2, APOE3 and APOE4 samples only: " + TEST_RULES.format(
            test="test APOE4 against APOE3 with a Wald test") + " Report the number of differentially expressed "
        "genes, how many of them are higher in APOE4, and the log2 fold change of APOE (APOE4 over APOE3).",
        fields=[
            Field("n_deg", "int", "number of genes with FDR < 0.05 for APOE4 against APOE3",
                  tolerance={"lower": 50, "upper": 90}, digits=0),
            Field(up, "int", "how many of those genes have a positive log2 fold change (higher in APOE4)",
                  tolerance={"lower": 30, "upper": 50}, digits=0),
            Field("apoe_log2fc", "float", "log2 fold change of APOE (Entrez 348), APOE4 over APOE3",
                  tolerance={"value": 0.12}, digits=3),
        ],
        data=list(RNA_FILES),
        failure_modes=["qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology"],
        category="differential_expression",
        kit="bulk_rnaseq",
        intent="Test whether the agent uses tximport's per-sample transcript lengths when it tests gene-level "
               "counts summarised from Salmon.",
        trap=(f"The gene-level counts are tested as plain read counts, with size factors only. The APOE4 "
              f"libraries differ from the APOE3 ones in average transcript length far more consistently than the "
              f"APOE2 ones do: {shifted['E4']} of {int(expressed.sum())} genes with a mean of 10 or more counts "
              f"differ between APOE4 and APOE3 at t-test p < 0.001, against {shifted['E2']} between APOE2 and "
              f"APOE3. Without the lengths {int((trap_r['fdr'] < 0.05).sum())} genes come out as differentially "
              f"expressed instead of {ref['n_deg']}; of the {int(only_trap.sum())} genes called only without the "
              f"lengths, {_pct(agree)} change in the same direction as their APOE4/APOE3 length ratio. The APOE "
              f"log2 fold change moves from {_f(ref['apoe_log2fc'])} to {_f(trap_r.loc[R.APOE, 'log2FoldChange'])}."),
        tolerance_rationale="The windows accept the lengths as DESeq2 offsets, length-scaled counts with DESeq2 or "
                            "pydeseq2, DESeq2's independent filtering, and the paper's own model with the three "
                            "knockout grafts it kept. They do not accept the counts without the lengths.",
        not_said="The prompt describes the length matrix but does not say what to do with it.",
        provenance=(f"The authors built their DESeq2 object with DESeqDataSetFromTximport (manuscript_figures.qmd "
                    f"in their repository), so the lengths entered as offsets. GEO holds only the counts; the "
                    f"lengths here are the avgTxLength assay of the authors' saved DESeqDataSet. The paper reports "
                    f"332 genes (170 up, 158 down) for APOE4 against APOE3 from a model that also held the knockout "
                    f"grafts; that model here gives {paper_r['fdr'].lt(0.05).sum()}. {int(both.sum())} genes are "
                    f"called with and without the lengths."),
    )
    return Built(task, ref, controls, [files[n] for n in RNA_FILES])


# --- apoe03: VDR targets, direction of the stored comparison ----------------------


def t3_vdr_direction(degs: dict, vdr, files: dict) -> Built:
    sheet = degs["E4vsE2"]
    targets = R.vdr_induced(vdr)
    ref = R.vdr_test(sheet, targets)
    authors = R.vdr_test(sheet, targets, unit="authors")
    rows = dict(ref)
    rows["n_higher_in_apoe2"] = authors["n_higher_in_apoe2"]
    trap = R.vdr_test(sheet, targets, higher_in_apoe2_sign=1)
    controls = [
        Control("reference", "reference",
                "Genes with fdr < 0.05 and a negative log2 fold change in the APOE4 against APOE2 sheet; unique "
                "symbols; hypergeometric test inside the sheet's symbols.", ref),
        Control("rows_not_symbols", "variant", "Count table rows instead of unique symbols.", rows),
        Control("authors_universe", "variant",
                "The authors' code: all rows as the universe and every VDR target counted, inside the sheet or not.",
                authors),
        Control("reversed_direction", "trap",
                "Take a positive log2 fold change in the APOE4 against APOE2 sheet as higher in APOE2.", trap),
    ]
    shared = set(trap["vdr_targets_higher_in_apoe2"]) & set(ref["vdr_targets_higher_in_apoe2"])
    task = Task(
        id="apoe03_vdr_targets_direction",
        title="Are vitamin D receptor targets higher in APOE2 or in APOE4 microglia",
        paper_anchor="Fig. 7b (VDR target genes among genes up in APOE2); Results, vitamin D receptor binding",
        prompt=f"""Regions that open up in APOE2 microglia carry the vitamin D receptor (VDR) motif, so a lab wants to know whether VDR target genes are expressed more highly in APOE2 than in APOE4 microglia. You have differential expression results for human iPSC-derived microglia grafted into AppNL-G-F mice (`{I.DEGS}`: one sheet per comparison of DESeq2 results, with columns entrezgene_id, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, fdr, comparison) and the response of human monocytes to vitamin D after 24 h (`{I.VDR}`: GeneSymbol, log2FoldChange, pvalue, padj for vitamin D against vehicle).

Define VDR targets as the genes induced by vitamin D in the monocytes (log2FoldChange > 0 and padj < 0.05). Take the genes expressed more highly in APOE2 than in APOE4 microglia at FDR < 0.05 (the fdr column), and test whether VDR targets are over-represented among them with a one-sided hypergeometric test. Work with unique gene symbols: the universe is the set of symbols in the comparison's sheet, and only VDR targets inside that universe count. Report the number of genes higher in APOE2, how many of them are VDR targets, their symbols, and -log10 of the p-value.""",
        fields=[
            Field("n_higher_in_apoe2", "int",
                  "number of unique gene symbols expressed more highly in APOE2 than in APOE4 at FDR < 0.05",
                  tolerance={"value": 15}, digits=0),
            Field("n_vdr_targets", "int", "how many of those genes are VDR targets", tolerance={"value": 3},
                  digits=0),
            Field("vdr_targets_higher_in_apoe2", "labels", "the gene symbols of those VDR targets",
                  pass_threshold=0.9),
            Field("neglog10_p", "float", "-log10 of the one-sided hypergeometric p-value",
                  tolerance={"lower": 2.5, "upper": 1.0}, digits=2),
        ],
        data=[I.DEGS, I.VDR],
        failure_modes=["biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error"],
        category="enrichment",
        kit="bulk_rnaseq",
        intent="Test whether the agent reads which way a stored comparison runs before it selects the genes "
               "higher in one group.",
        trap=(f"The sheet for this pair (E4vE2_DEGs, comparison E4vsE2) holds APOE4 over APOE2, and a positive "
              f"log2 fold change is read as higher in APOE2. The genes taken are then those higher in APOE4: "
              f"{trap['n_higher_in_apoe2']} genes, {trap['n_vdr_targets']} of them VDR targets "
              f"(-log10 p {_f(trap['neglog10_p'], 2)}), and {len(shared)} of those targets are among the "
              f"{ref['n_vdr_targets']} that are higher in APOE2 (-log10 p {_f(ref['neglog10_p'], 2)})."),
        tolerance_rationale=(f"The windows accept counting table rows instead of unique symbols and the authors' "
                             f"universe (all rows, every VDR target counted), which moves -log10 p from "
                             f"{_f(ref['neglog10_p'], 2)} to {_f(authors['neglog10_p'], 2)}. They do not accept the "
                             f"reversed direction."),
        not_said="The prompt does not say which way the stored comparison runs; the sheet name and the comparison "
                 "column do.",
        provenance="The paper reports VDR targets over-represented among genes up in APOE2 against APOE3 and "
                   "against APOE4 (Fig. 7b); its code takes the genes down in the APOE4 against APOE2 table for the "
                   "second test. The DESeq2 table is the authors' Supplementary Data 2 (Zenodo); the monocyte table "
                   "is the one in their repository.",
    )
    return Built(task, ref, controls, [files[I.DEGS], files[I.VDR]])


# --- apoe04: promoter accessibility against expression -----------------------------


def t4_matched_concordance(degs: dict, peaks, files: dict) -> Built:
    ref = R.concordance(degs, peaks)
    trap = R.concordance(degs, peaks, pooled=True)
    controls = [
        Control("reference", "reference",
                "For each comparison, median promoter log2 fold change per gene from the same comparison; "
                "Pearson r over DEGs with a promoter peak; genes matched by Entrez ID.", ref),
        Control("mean_per_gene", "variant", "Mean instead of median of a gene's promoter peaks.",
                R.concordance(degs, peaks, agg="mean")),
        Control("match_by_symbol", "variant", "Genes matched by symbol instead of Entrez ID.",
                R.concordance(degs, peaks, key="symbol")),
        Control("pooled_comparisons", "trap",
                "Each gene's promoter summary taken over the rows of all four comparisons at once.", trap),
    ]
    tags = (("e2_vs_e3", "E2vsE3"), ("e4_vs_e3", "E4vsE3"), ("e4_vs_e2", "E4vsE2"), ("ko_vs_e3", "KOvsE3"))
    fields = []
    for tag, comp in tags:
        fields.append(Field(f"r_{tag}", "float",
                            f"Pearson r between expression and median promoter accessibility log2 fold changes, {comp}",
                            tolerance={"value": 0.035}, digits=4))
    for tag, comp in tags:
        fields.append(Field(f"n_{tag}", "int", f"number of genes in that correlation, {comp}",
                            tolerance={"value": 15 if tag == "e4_vs_e2" else 10}, digits=0))
    gaps = {comp: ref[f"r_{tag}"] - trap[f"r_{tag}"] for tag, comp in tags}
    matched = ", ".join(comp + " " + _f(ref["r_" + tag], 2) for tag, comp in tags)
    task = Task(
        id="apoe04_matched_contrast_concordance",
        title="Promoter accessibility against expression, one comparison at a time",
        paper_anchor="Supplementary Fig. 6 (expression and promoter accessibility log2 fold changes); Results, "
                     "concordance of ATAC-seq and RNA-seq",
        prompt=f"""Bulk ATAC-seq and RNA-seq were run on human iPSC-derived microglia (APOE2, APOE3, APOE4 and APOE-KO lines) grafted into AppNL-G-F mice. You have:

- `{I.DEGS}`: DESeq2 differential expression, one sheet per comparison (E2vsE3_DEGs, E4vsE3_DEGs, E4vE2_DEGs, KOvsE3_DEGs; columns entrezgene_id, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, fdr, comparison)
- `{I.PROMOTER_PEAKS}`: DESeq2 differential accessibility for every ATAC peak within 3 kb of an annotated transcription start site, one row per peak and comparison (comparison, Peak, Chr, Start, End, log2FoldChange, pvalue, fdr, annotation, distanceToTSS, geneId as Entrez ID, SYMBOL)

In both files the log2 fold change is the first-named genotype over the second (E2vsE3 is APOE2 over APOE3).

The lab wants to know whether promoter accessibility moves with expression in each comparison. For each of the 4 comparisons (E2vsE3, E4vsE3, E4vsE2, KOvsE3): take the genes differentially expressed at FDR < 0.05 (the fdr column), summarise each gene's promoter accessibility change as the median log2 fold change of the peaks annotated to it (match genes by Entrez ID), and compute the Pearson correlation between the expression log2 fold change and that median over the genes with at least one annotated peak. Report the correlation and the number of genes for each comparison.""",
        fields=fields,
        data=[I.DEGS, I.PROMOTER_PEAKS],
        failure_modes=["subset_scope_compartment.wrong_comparator_contrast_or_baseline"],
        category="multiomic_integration",
        kit="atac",
        intent="Test whether the agent takes the accessibility change from the same comparison as the expression "
               "change when both sit in long tables.",
        trap=(f"Each gene's accessibility change is summarised over the rows of all four comparisons at once, "
              f"instead of the comparison that matches the expression change. The median then mixes in the changes of the other three comparisons, and r drops by "
              f"{_f(gaps['E2vsE3'], 2)} for E2vsE3, {_f(gaps['KOvsE3'], 2)} for KOvsE3, "
              f"{_f(gaps['E4vsE2'], 2)} for E4vsE2 and {_f(gaps['E4vsE3'], 2)} for E4vsE3; the gene counts do not "
              f"change."),
        tolerance_rationale="The windows accept the mean instead of the median and matching genes by symbol "
                            "instead of Entrez ID (which loses up to 10 genes). They do not accept accessibility "
                            "summarised over all comparisons at once; that answer fails on E2vsE3, E4vsE2 and KOvsE3.",
        not_said="The prompt does not point out that the promoter table stacks the four comparisons in one file.",
        provenance=(f"The paper reports a strong correlation across all APOE comparisons between expression and "
                    f"accessibility changes at promoter peaks (Supplementary Fig. 6), and its text describes matched "
                    f"comparisons. Matched, the correlations are "
                    f"{matched}, all higher than pooled. Both tables are the authors' own: Supplementary Data 2 (Zenodo) "
                    f"and the promoter-peak table in their repository, reduced to the columns above."),
    )
    return Built(task, ref, controls, [files[I.DEGS], files[I.PROMOTER_PEAKS]])
