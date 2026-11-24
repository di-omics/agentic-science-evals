# Tasks: APOE microglia, Nat Commun 2025

Murphy KB, Hu D, Wolfs L, Rohde SK, Leguia Fauro G, Geric I, Mancuso R, De Strooper B, Marzi SJ. The APOE isoforms differentially shape the transcriptomic and epigenomic landscapes of human microglia xenografted into a mouse model of Alzheimer's disease. Nat Commun 16, 4883 (2025). doi:10.1038/s41467-025-60099-4

## `apoe01_failed_knockouts`: APOE knockout against APOE3: confirm the knockout before comparing

From the paper: Fig. 1b and Supplementary Figs. 2-3 (APOE-KO against APOE3); Methods, RNA-seq QC and processing

```
Human iPSC-derived microglia from 4 isogenic lines were transplanted into the brains of AppNL-G-F mice: lines carrying a single APOE2, APOE3 or APOE4 allele, and an APOE knockout line (APOE-KO). After 12 months the human microglia were sorted from each mouse and profiled by bulk RNA-seq. Reads were quantified against GRCh38 transcripts with Salmon (with its GC-bias correction) and summarised to genes with tximport. You have:

- `GSE271385_APOE_microglia_RNAseq_counts.csv.gz`: gene-level estimated counts as deposited in GEO (rows: Entrez gene IDs; columns: samples)
- `rnaseq_avg_tx_length.csv`: the matching tximport average transcript lengths (same rows and columns)
- `rnaseq_samples.csv`: the genotype group, iPSC line and GEO accession of each sample
- `rnaseq_genes.csv`: gene symbols for the Entrez IDs that have one

The lab wants the genes that change when microglia lose APOE. Compare APOE-KO microglia with APOE3 microglia, with all four genotype groups in one model: fit one negative binomial GLM (DESeq2 or an equivalent) with genotype as the only factor, test APOE-KO against APOE3 with a Wald test, and adjust the p-values with Benjamini-Hochberg over every gene that receives a p-value (no independent filtering). A gene is differentially expressed at FDR < 0.05. Report the number of differentially expressed genes, how many of them are higher in APOE-KO, and the log2 fold change of APOE (APOE-KO over APOE3).

Return EXACTLY:

<EVAL_ANSWER>
{"n_deg": <integer>, "n_higher_in_ko": <integer>, "apoe_log2fc": <float>}
</EVAL_ANSWER>

`n_deg`: number of genes with FDR < 0.05 for APOE-KO against APOE3
`n_higher_in_ko`: how many of those genes have a positive log2 fold change (higher in APOE-KO)
`apoe_log2fc`: log2 fold change of APOE (Entrez 348), APOE-KO over APOE3
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "n_deg": 338,
   "n_higher_in_ko": 133,
   "apoe_log2fc": -4.077
  },
  "tolerances": {
   "n_deg": {
    "type": "absolute",
    "lower": 45,
    "upper": 75
   },
   "n_higher_in_ko": {
    "type": "absolute",
    "lower": 25,
    "upper": 35
   },
   "apoe_log2fc": {
    "type": "absolute",
    "value": 0.4
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent checks that the knockout grafts lack APOE before it uses them as the knockout group.

2) Primary trap: All five APOE-KO samples are used. KO-1 and KO-2 carry APOE at 55% and 48% of the APOE3 median (counts per million), while KO-3, KO-4, KO-5 sit at 5% to 7%; they are not knockouts in this experiment. With them the APOE log2 fold change shrinks from -4.08 to -2.13, the number of differentially expressed genes rises from 338 to 458, and the number higher in APOE-KO from 133 to 231.

3) Scripted controls, run on the task data:
   reference reference: Drop the knockout grafts whose APOE level is at least a quarter of the APOE3 median, then DESeq2 with the tximport lengths as offsets, BH over all genes. -> n_deg=338, n_higher_in_ko=133, apoe_log2fc=-4.077 (passes)
   variant plain_counts: Same samples; counts without the transcript lengths. -> n_deg=333, n_higher_in_ko=136, apoe_log2fc=-4.21 (passes)
   variant length_scaled_counts: Same samples; lengths divided out of the counts first. -> n_deg=340, n_higher_in_ko=133, apoe_log2fc=-4.078 (passes)
   variant independent_filtering: Same samples and offsets; DESeq2's padj with independent filtering instead of BH over all genes. -> n_deg=401, n_higher_in_ko=164, apoe_log2fc=-4.077 (passes)
   variant pydeseq2: Same samples; the pydeseq2 package on the counts. -> n_deg=319, n_higher_in_ko=131, apoe_log2fc=-4.21 (passes)
   trap all_knockout_grafts: Keep all five APOE-KO samples, two of which express APOE; otherwise the reference analysis. -> n_deg=458, n_higher_in_ko=231, apoe_log2fc=-2.134 (fails)

4) Tolerance: The windows accept the tximport lengths as offsets, plain counts, length-scaled counts, pydeseq2 in place of DESeq2, and DESeq2's independent filtering in place of BH over all genes. They do not accept keeping the two knockout grafts that express APOE.

5) Not said in the prompt: The prompt does not say that any graft failed or that the knockout should be confirmed from the data. The sample sheet lists all five knockout samples under the same line.

6) Provenance: The paper kept the 3 knockout samples in which APOE was significantly lower (Fig. 1b, Supplementary Fig. 2) and reports 126 genes up and 206 down in APOE-KO against APOE3. The reference here, with the tximport lengths as DESeq2 offsets, finds 338 (133 up, 205 down) over every gene with a p-value; the authors' Supplementary Data 2 lists 335 (127 up, 208 down) over the genes that have a symbol. Counts are from GEO GSE271385; the lengths are the avgTxLength assay of the DESeqDataSet in the authors' repository.
```

## `apoe02_transcript_length_offsets`: APOE4 against APOE3 from tximport counts

From the paper: Fig. 1d and Fig. 2b (APOE4 against APOE3); Methods, RNA-seq QC and processing

```
Human iPSC-derived microglia from 4 isogenic lines were transplanted into the brains of AppNL-G-F mice: lines carrying a single APOE2, APOE3 or APOE4 allele, and an APOE knockout line (APOE-KO). After 12 months the human microglia were sorted from each mouse and profiled by bulk RNA-seq. Reads were quantified against GRCh38 transcripts with Salmon (with its GC-bias correction) and summarised to genes with tximport. You have:

- `GSE271385_APOE_microglia_RNAseq_counts.csv.gz`: gene-level estimated counts as deposited in GEO (rows: Entrez gene IDs; columns: samples)
- `rnaseq_avg_tx_length.csv`: the matching tximport average transcript lengths (same rows and columns)
- `rnaseq_samples.csv`: the genotype group, iPSC line and GEO accession of each sample
- `rnaseq_genes.csv`: gene symbols for the Entrez IDs that have one

A colleague wants the list of genes that differ between APOE4 and APOE3 microglia. Use the APOE2, APOE3 and APOE4 samples only: fit one negative binomial GLM (DESeq2 or an equivalent) with genotype as the only factor, test APOE4 against APOE3 with a Wald test, and adjust the p-values with Benjamini-Hochberg over every gene that receives a p-value (no independent filtering). A gene is differentially expressed at FDR < 0.05. Report the number of differentially expressed genes, how many of them are higher in APOE4, and the log2 fold change of APOE (APOE4 over APOE3).

Return EXACTLY:

<EVAL_ANSWER>
{"n_deg": <integer>, "n_higher_in_apoe4": <integer>, "apoe_log2fc": <float>}
</EVAL_ANSWER>

`n_deg`: number of genes with FDR < 0.05 for APOE4 against APOE3
`n_higher_in_apoe4`: how many of those genes have a positive log2 fold change (higher in APOE4)
`apoe_log2fc`: log2 fold change of APOE (Entrez 348), APOE4 over APOE3
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "n_deg": 270,
   "n_higher_in_apoe4": 153,
   "apoe_log2fc": -0.463
  },
  "tolerances": {
   "n_deg": {
    "type": "absolute",
    "lower": 50,
    "upper": 90
   },
   "n_higher_in_apoe4": {
    "type": "absolute",
    "lower": 30,
    "upper": 50
   },
   "apoe_log2fc": {
    "type": "absolute",
    "value": 0.12
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent uses tximport's per-sample transcript lengths when it tests gene-level counts summarised from Salmon.

2) Primary trap: The gene-level counts are tested as plain read counts, with size factors only. The APOE4 libraries differ from the APOE3 ones in average transcript length far more consistently than the APOE2 ones do: 327 of 13173 genes with a mean of 10 or more counts differ between APOE4 and APOE3 at t-test p < 0.001, against 23 between APOE2 and APOE3. Without the lengths 1220 genes come out as differentially expressed instead of 270; of the 993 genes called only without the lengths, 87% change in the same direction as their APOE4/APOE3 length ratio. The APOE log2 fold change moves from -0.463 to -0.726.

3) Scripted controls, run on the task data:
   reference reference: APOE2, APOE3 and APOE4 samples; DESeq2 with the tximport lengths as offsets; BH over all genes. -> n_deg=270, n_higher_in_apoe4=153, apoe_log2fc=-0.4634 (passes)
   variant length_scaled_counts: Lengths divided out of the counts, then DESeq2 with size factors. -> n_deg=273, n_higher_in_apoe4=154, apoe_log2fc=-0.4636 (passes)
   variant pydeseq2_length_scaled: Length-scaled counts with the pydeseq2 package. -> n_deg=254, n_higher_in_apoe4=142, apoe_log2fc=-0.4636 (passes)
   variant independent_filtering: Offsets; DESeq2's padj with independent filtering. -> n_deg=320, n_higher_in_apoe4=173, apoe_log2fc=-0.4634 (passes)
   variant paper_model: The paper's model: offsets, all samples except the two knockout grafts that express APOE. -> n_deg=333, n_higher_in_apoe4=179, apoe_log2fc=-0.4633 (passes)
   trap counts_without_lengths: The gene-level counts tested as plain read counts, size factors only. -> n_deg=1220, n_higher_in_apoe4=628, apoe_log2fc=-0.7255 (fails)

4) Tolerance: The windows accept the lengths as DESeq2 offsets, length-scaled counts with DESeq2 or pydeseq2, DESeq2's independent filtering, and the paper's own model with the three knockout grafts it kept. They do not accept the counts without the lengths.

5) Not said in the prompt: The prompt describes the length matrix but does not say what to do with it.

6) Provenance: The authors built their DESeq2 object with DESeqDataSetFromTximport (manuscript_figures.qmd in their repository), so the lengths entered as offsets. GEO holds only the counts; the lengths here are the avgTxLength assay of the authors' saved DESeqDataSet. The paper reports 332 genes (170 up, 158 down) for APOE4 against APOE3 from a model that also held the knockout grafts; that model here gives 333. 227 genes are called with and without the lengths.
```

## `apoe03_vdr_targets_direction`: Are vitamin D receptor targets higher in APOE2 or in APOE4 microglia

From the paper: Fig. 7b (VDR target genes among genes up in APOE2); Results, vitamin D receptor binding

```
Regions that open up in APOE2 microglia carry the vitamin D receptor (VDR) motif, so a lab wants to know whether VDR target genes are expressed more highly in APOE2 than in APOE4 microglia. You have differential expression results for human iPSC-derived microglia grafted into AppNL-G-F mice (`Supplementary_Table2_DEGs.xlsx`: one sheet per comparison of DESeq2 results, with columns entrezgene_id, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, fdr, comparison) and the response of human monocytes to vitamin D after 24 h (`VDR_target_genes_24h_monocytes_Warwick2021.csv`: GeneSymbol, log2FoldChange, pvalue, padj for vitamin D against vehicle).

Define VDR targets as the genes induced by vitamin D in the monocytes (log2FoldChange > 0 and padj < 0.05). Take the genes expressed more highly in APOE2 than in APOE4 microglia at FDR < 0.05 (the fdr column), and test whether VDR targets are over-represented among them with a one-sided hypergeometric test. Work with unique gene symbols: the universe is the set of symbols in the comparison's sheet, and only VDR targets inside that universe count. Report the number of genes higher in APOE2, how many of them are VDR targets, their symbols, and -log10 of the p-value.

Return EXACTLY:

<EVAL_ANSWER>
{"n_higher_in_apoe2": <integer>, "n_vdr_targets": <integer>, "vdr_targets_higher_in_apoe2": <list of strings>, "neglog10_p": <float>}
</EVAL_ANSWER>

`n_higher_in_apoe2`: number of unique gene symbols expressed more highly in APOE2 than in APOE4 at FDR < 0.05
`n_vdr_targets`: how many of those genes are VDR targets
`vdr_targets_higher_in_apoe2`: the gene symbols of those VDR targets
`neglog10_p`: -log10 of the one-sided hypergeometric p-value
```

Grader:

```json
{
 "type": "all_of",
 "config": {
  "children": [
   {
    "type": "numeric_tolerance",
    "config": {
     "ground_truth": {
      "n_higher_in_apoe2": 890,
      "n_vdr_targets": 82,
      "neglog10_p": 10.44
     },
     "tolerances": {
      "n_higher_in_apoe2": {
       "type": "absolute",
       "value": 15
      },
      "n_vdr_targets": {
       "type": "absolute",
       "value": 3
      },
      "neglog10_p": {
       "type": "absolute",
       "lower": 2.5,
       "upper": 1.0
      }
     }
    }
   },
   {
    "type": "label_set_jaccard",
    "config": {
     "ground_truth_labels": [
      "ABCC1",
      "ACO2",
      "ADAMTS10",
      "ADGRE5",
      "ARHGEF10L",
      "BCL9L",
      "CHCHD10",
      "CMTM7",
      "COX4I1",
      "CSF3R",
      "CTSB",
      "CTSD",
      "DNAJC4",
      "DNASE2",
      "EFNB1",
      "ELL",
      "EPB41L1",
      "FAM20C",
      "FAM53B",
      "FGFR1",
      "FTL",
      "FUCA1",
      "GIPC3",
      "GLIPR1",
      "GM2A",
      "GNPDA1",
      "GPI",
      "HBEGF",
      "HIGD2A",
      "HPCAL1",
      "IFNGR2",
      "ITGB2",
      "KLF11",
      "LAMC1",
      "LAMP1",
      "LGALS9",
      "LGI2",
      "LILRB4",
      "LRCH4",
      "MARCO",
      "MDH2",
      "MFSD12",
      "MGLL",
      "MRPL4",
      "MYO1F",
      "MYO9B",
      "NCF2",
      "NCF4",
      "NCOR2",
      "NDUFS7",
      "NFAM1",
      "NFKBIA",
      "NHSL2",
      "NINJ1",
      "PABPC1",
      "PIEZO1",
      "PILRA",
      "PLTP",
      "PRKACA",
      "RBM47",
      "RIN2",
      "RIN3",
      "SASH3",
      "SCARB1",
      "SEMA4B",
      "SHE",
      "SIPA1",
      "SLC27A1",
      "SLC37A2",
      "SMAP2",
      "SMPDL3A",
      "TEX264",
      "THBD",
      "THEMIS2",
      "TMC6",
      "TMEM127",
      "TMEM205",
      "UAP1L1",
      "VRK3",
      "XXYLT1",
      "XYLT1",
      "ZMIZ1"
     ],
     "scoring": {
      "pass_threshold": 0.9
     },
     "answer_field": "vdr_targets_higher_in_apoe2"
    }
   }
  ]
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent reads which way a stored comparison runs before it selects the genes higher in one group.

2) Primary trap: The sheet for this pair (E4vE2_DEGs, comparison E4vsE2) holds APOE4 over APOE2, and a positive log2 fold change is read as higher in APOE2. The genes taken are then those higher in APOE4: 735 genes, 30 of them VDR targets (-log10 p 0.20), and 0 of those targets are among the 82 that are higher in APOE2 (-log10 p 10.44).

3) Scripted controls, run on the task data:
   reference reference: Genes with fdr < 0.05 and a negative log2 fold change in the APOE4 against APOE2 sheet; unique symbols; hypergeometric test inside the sheet's symbols. -> n_higher_in_apoe2=890, n_vdr_targets=82, vdr_targets_higher_in_apoe2=82 labels, neglog10_p=10.44 (passes)
   variant rows_not_symbols: Count table rows instead of unique symbols. -> n_higher_in_apoe2=891, n_vdr_targets=82, vdr_targets_higher_in_apoe2=82 labels, neglog10_p=10.44 (passes)
   variant authors_universe: The authors' code: all rows as the universe and every VDR target counted, inside the sheet or not. -> n_higher_in_apoe2=891, n_vdr_targets=82, vdr_targets_higher_in_apoe2=82 labels, neglog10_p=8.294 (passes)
   trap reversed_direction: Take a positive log2 fold change in the APOE4 against APOE2 sheet as higher in APOE2. -> n_higher_in_apoe2=735, n_vdr_targets=30, vdr_targets_higher_in_apoe2=30 labels, neglog10_p=0.1966 (fails)

4) Tolerance: The windows accept counting table rows instead of unique symbols and the authors' universe (all rows, every VDR target counted), which moves -log10 p from 10.44 to 8.29. They do not accept the reversed direction.

5) Not said in the prompt: The prompt does not say which way the stored comparison runs; the sheet name and the comparison column do.

6) Provenance: The paper reports VDR targets over-represented among genes up in APOE2 against APOE3 and against APOE4 (Fig. 7b); its code takes the genes down in the APOE4 against APOE2 table for the second test. The DESeq2 table is the authors' Supplementary Data 2 (Zenodo); the monocyte table is the one in their repository.
```

## `apoe04_matched_contrast_concordance`: Promoter accessibility against expression, one comparison at a time

From the paper: Supplementary Fig. 6 (expression and promoter accessibility log2 fold changes); Results, concordance of ATAC-seq and RNA-seq

```
Bulk ATAC-seq and RNA-seq were run on human iPSC-derived microglia (APOE2, APOE3, APOE4 and APOE-KO lines) grafted into AppNL-G-F mice. You have:

- `Supplementary_Table2_DEGs.xlsx`: DESeq2 differential expression, one sheet per comparison (E2vsE3_DEGs, E4vsE3_DEGs, E4vE2_DEGs, KOvsE3_DEGs; columns entrezgene_id, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, fdr, comparison)
- `atac_promoter_peaks.csv`: DESeq2 differential accessibility for every ATAC peak within 3 kb of an annotated transcription start site, one row per peak and comparison (comparison, Peak, Chr, Start, End, log2FoldChange, pvalue, fdr, annotation, distanceToTSS, geneId as Entrez ID, SYMBOL)

In both files the log2 fold change is the first-named genotype over the second (E2vsE3 is APOE2 over APOE3).

The lab wants to know whether promoter accessibility moves with expression in each comparison. For each of the 4 comparisons (E2vsE3, E4vsE3, E4vsE2, KOvsE3): take the genes differentially expressed at FDR < 0.05 (the fdr column), summarise each gene's promoter accessibility change as the median log2 fold change of the peaks annotated to it (match genes by Entrez ID), and compute the Pearson correlation between the expression log2 fold change and that median over the genes with at least one annotated peak. Report the correlation and the number of genes for each comparison.

Return EXACTLY:

<EVAL_ANSWER>
{"r_e2_vs_e3": <float>, "r_e4_vs_e3": <float>, "r_e4_vs_e2": <float>, "r_ko_vs_e3": <float>, "n_e2_vs_e3": <integer>, "n_e4_vs_e3": <integer>, "n_e4_vs_e2": <integer>, "n_ko_vs_e3": <integer>}
</EVAL_ANSWER>

`r_e2_vs_e3`: Pearson r between expression and median promoter accessibility log2 fold changes, E2vsE3
`r_e4_vs_e3`: Pearson r between expression and median promoter accessibility log2 fold changes, E4vsE3
`r_e4_vs_e2`: Pearson r between expression and median promoter accessibility log2 fold changes, E4vsE2
`r_ko_vs_e3`: Pearson r between expression and median promoter accessibility log2 fold changes, KOvsE3
`n_e2_vs_e3`: number of genes in that correlation, E2vsE3
`n_e4_vs_e3`: number of genes in that correlation, E4vsE3
`n_e4_vs_e2`: number of genes in that correlation, E4vsE2
`n_ko_vs_e3`: number of genes in that correlation, KOvsE3
```

Grader:

```json
{
 "type": "numeric_tolerance",
 "config": {
  "ground_truth": {
   "r_e2_vs_e3": 0.4351,
   "r_e4_vs_e3": 0.6398,
   "r_e4_vs_e2": 0.4458,
   "r_ko_vs_e3": 0.6273,
   "n_e2_vs_e3": 444,
   "n_e4_vs_e3": 308,
   "n_e4_vs_e2": 1548,
   "n_ko_vs_e3": 316
  },
  "tolerances": {
   "r_e2_vs_e3": {
    "type": "absolute",
    "value": 0.035
   },
   "r_e4_vs_e3": {
    "type": "absolute",
    "value": 0.035
   },
   "r_e4_vs_e2": {
    "type": "absolute",
    "value": 0.035
   },
   "r_ko_vs_e3": {
    "type": "absolute",
    "value": 0.035
   },
   "n_e2_vs_e3": {
    "type": "absolute",
    "value": 10
   },
   "n_e4_vs_e3": {
    "type": "absolute",
    "value": 10
   },
   "n_e4_vs_e2": {
    "type": "absolute",
    "value": 15
   },
   "n_ko_vs_e3": {
    "type": "absolute",
    "value": 10
   }
  }
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent takes the accessibility change from the same comparison as the expression change when both sit in long tables.

2) Primary trap: Each gene's accessibility change is summarised over the rows of all four comparisons at once, instead of the comparison that matches the expression change. The median then mixes in the changes of the other three comparisons, and r drops by 0.29 for E2vsE3, 0.12 for KOvsE3, 0.05 for E4vsE2 and 0.03 for E4vsE3; the gene counts do not change.

3) Scripted controls, run on the task data:
   reference reference: For each comparison, median promoter log2 fold change per gene from the same comparison; Pearson r over DEGs with a promoter peak; genes matched by Entrez ID. -> r_e2_vs_e3=0.4351, r_e4_vs_e3=0.6398, r_e4_vs_e2=0.4458, r_ko_vs_e3=0.6273 (passes)
   variant mean_per_gene: Mean instead of median of a gene's promoter peaks. -> r_e2_vs_e3=0.4628, r_e4_vs_e3=0.6344, r_e4_vs_e2=0.4369, r_ko_vs_e3=0.6194 (passes)
   variant match_by_symbol: Genes matched by symbol instead of Entrez ID. -> r_e2_vs_e3=0.4344, r_e4_vs_e3=0.6402, r_e4_vs_e2=0.4387, r_ko_vs_e3=0.6304 (passes)
   trap pooled_comparisons: Each gene's promoter summary taken over the rows of all four comparisons at once. -> r_e2_vs_e3=0.1426, r_e4_vs_e3=0.6061, r_e4_vs_e2=0.3916, r_ko_vs_e3=0.5055 (fails)

4) Tolerance: The windows accept the mean instead of the median and matching genes by symbol instead of Entrez ID (which loses up to 10 genes). They do not accept accessibility summarised over all comparisons at once; that answer fails on E2vsE3, E4vsE2 and KOvsE3.

5) Not said in the prompt: The prompt does not point out that the promoter table stacks the four comparisons in one file.

6) Provenance: The paper reports a strong correlation across all APOE comparisons between expression and accessibility changes at promoter peaks (Supplementary Fig. 6), and its text describes matched comparisons. Matched, the correlations are E2vsE3 0.44, E4vsE3 0.64, E4vsE2 0.45, KOvsE3 0.63, all higher than pooled. Both tables are the authors' own: Supplementary Data 2 (Zenodo) and the promoter-peak table in their repository, reduced to the columns above.
```
