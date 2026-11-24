# Tasks: AVE migration, Dev Cell 2024

Thowfeequ S, Fiorentino J, Hu D, Solovey M, Ruane S, Whitehead M, Zhou F, Godwin J, Mateo-Otero Y, Vanhaesebroeck B, Scialdone A, Srinivas S. An integrated approach identifies the molecular underpinnings of murine anterior visceral endoderm migration. Developmental Cell 59, 2347-2363.e9 (2024). doi:10.1016/j.devcel.2024.05.014

## `ave01_plate_qc`: Which wells of a Smart-seq2 plate hold one good cell

From the paper: STAR Methods, Transcript quantification, quality control and normalization (255 good-quality cells at 5.5 dpc)

```
Cells from 40 mouse embryos at 5.5 days post coitum were enriched for visceral endoderm (outer cells labelled with a membrane dye; live, labelled cells sorted by FACS into 384 wells) and processed with Smart-seq2 (paired-end 125 bp reads; ERCC spike-in mix added to every well at 1:10^7). You have:

- `counts_e55.txt`: Salmon gene-level counts (estimated, so not always whole numbers), one column per well and one row per Ensembl gene, plus 92 rows for the ERCC spike-ins (row names starting with `ERCC-`).
- `wells_e55.tsv`: one row per well with `well_id` (the column name in the counts file), `well_type` (the well's annotation in the submitters' sample sheet: `single cell`, `no cell` or `multiple cells`) and `read_pairs` (read pairs sequenced from the well).
- `genes_grcm38.tsv`: Ensembl gene ID to gene symbol (GRCm38; mitochondrial genes have symbols starting with `mt-`).

Before normalization and clustering, the lab needs the set of wells that each hold one good-quality cell. Decide which wells to keep, and report how many you keep and which they are.

Return EXACTLY:

<EVAL_ANSWER>
{"n_wells_kept": <integer>, "wells_kept": <list of strings>}
</EVAL_ANSWER>

`n_wells_kept`: number of wells kept for downstream analysis
`wells_kept`: the `well_id` of every well kept, as in the columns of `counts_e55.txt`
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
      "n_wells_kept": 255
     },
     "tolerances": {
      "n_wells_kept": {
       "type": "absolute",
       "value": 20
      }
     }
    }
   },
   {
    "type": "label_set_jaccard",
    "config": {
     "ground_truth_labels": [
      "SS_sc7784898",
      "SS_sc7784899",
      "SS_sc7784900",
      "SS_sc7784902",
      "SS_sc7784903",
      "SS_sc7784904",
      "SS_sc7784905",
      "SS_sc7784907",
      "SS_sc7784909",
      "SS_sc7784911",
      "SS_sc7784912",
      "SS_sc7784913",
      "SS_sc7784914",
      "SS_sc7784915",
      "SS_sc7784917",
      "SS_sc7784920",
      "SS_sc7784921",
      "SS_sc7784922",
      "SS_sc7784923",
      "SS_sc7784924",
      "SS_sc7784926",
      "SS_sc7784927",
      "SS_sc7784928",
      "SS_sc7784930",
      "SS_sc7784932",
      "SS_sc7784933",
      "SS_sc7784934",
      "SS_sc7784935",
      "SS_sc7784936",
      "SS_sc7784937",
      "SS_sc7784938",
      "SS_sc7784939",
      "SS_sc7784940",
      "SS_sc7784942",
      "SS_sc7784945",
      "SS_sc7784947",
      "SS_sc7784948",
      "SS_sc7784949",
      "SS_sc7784951",
      "SS_sc7784953",
      "SS_sc7784954",
      "SS_sc7784955",
      "SS_sc7784957",
      "SS_sc7784959",
      "SS_sc7784961",
      "SS_sc7784965",
      "SS_sc7784966",
      "SS_sc7784967",
      "SS_sc7784968",
      "SS_sc7784970",
      "SS_sc7784971",
      "SS_sc7784972",
      "SS_sc7784973",
      "SS_sc7784974",
      "SS_sc7784976",
      "SS_sc7784977",
      "SS_sc7784978",
      "SS_sc7784979",
      "SS_sc7784980",
      "SS_sc7784981",
      "SS_sc7784982",
      "SS_sc7784983",
      "SS_sc7784984",
      "SS_sc7784985",
      "SS_sc7784986",
      "SS_sc7784987",
      "SS_sc7784988",
      "SS_sc7784989",
      "SS_sc7784990",
      "SS_sc7784991",
      "SS_sc7784993",
      "SS_sc7784994",
      "SS_sc7784995",
      "SS_sc7785001",
      "SS_sc7785002",
      "SS_sc7785003",
      "SS_sc7785005",
      "SS_sc7785006",
      "SS_sc7785007",
      "SS_sc7785011",
      "SS_sc7785012",
      "SS_sc7785013",
      "SS_sc7785014",
      "SS_sc7785015",
      "SS_sc7785016",
      "SS_sc7785018",
      "SS_sc7785019",
      "SS_sc7785020",
      "SS_sc7785021",
      "SS_sc7785022",
      "SS_sc7785024",
      "SS_sc7785026",
      "SS_sc7785027",
      "SS_sc7785028",
      "SS_sc7785030",
      "SS_sc7785031",
      "SS_sc7785032",
      "SS_sc7785034",
      "SS_sc7785035",
      "SS_sc7785036",
      "SS_sc7785037",
      "SS_sc7785038",
      "SS_sc7785039",
      "SS_sc7785043",
      "SS_sc7785044",
      "SS_sc7785045",
      "SS_sc7785047",
      "SS_sc7785048",
      "SS_sc7785049",
      "SS_sc7785050",
      "SS_sc7785051",
      "SS_sc7785052",
      "SS_sc7785053",
      "SS_sc7785055",
      "SS_sc7785056",
      "SS_sc7785057",
      "SS_sc7785059",
      "SS_sc7785060",
      "SS_sc7785061",
      "SS_sc7785062",
      "SS_sc7785063",
      "SS_sc7785065",
      "SS_sc7785066",
      "SS_sc7785068",
      "SS_sc7785069",
      "SS_sc7785070",
      "SS_sc7785071",
      "SS_sc7785072",
      "SS_sc7785074",
      "SS_sc7785075",
      "SS_sc7785076",
      "SS_sc7785077",
      "SS_sc7785088",
      "SS_sc7785089",
      "SS_sc7785090",
      "SS_sc7785091",
      "SS_sc7785092",
      "SS_sc7785093",
      "SS_sc7785094",
      "SS_sc7785097",
      "SS_sc7785099",
      "SS_sc7785101",
      "SS_sc7785103",
      "SS_sc7785105",
      "SS_sc7785106",
      "SS_sc7785107",
      "SS_sc7785110",
      "SS_sc7785111",
      "SS_sc7785112",
      "SS_sc7785113",
      "SS_sc7785114",
      "SS_sc7785115",
      "SS_sc7785116",
      "SS_sc7785117",
      "SS_sc7785118",
      "SS_sc7785119",
      "SS_sc7785120",
      "SS_sc7785121",
      "SS_sc7785123",
      "SS_sc7785124",
      "SS_sc7785126",
      "SS_sc7785128",
      "SS_sc7785129",
      "SS_sc7785130",
      "SS_sc7785131",
      "SS_sc7785132",
      "SS_sc7785135",
      "SS_sc7785136",
      "SS_sc7785137",
      "SS_sc7785139",
      "SS_sc7785140",
      "SS_sc7785141",
      "SS_sc7785142",
      "SS_sc7785143",
      "SS_sc7785145",
      "SS_sc7785147",
      "SS_sc7785148",
      "SS_sc7785149",
      "SS_sc7785153",
      "SS_sc7785154",
      "SS_sc7785155",
      "SS_sc7785156",
      "SS_sc7785157",
      "SS_sc7785158",
      "SS_sc7785160",
      "SS_sc7785162",
      "SS_sc7785164",
      "SS_sc7785165",
      "SS_sc7785166",
      "SS_sc7785168",
      "SS_sc7785170",
      "SS_sc7785172",
      "SS_sc7785173",
      "SS_sc7785174",
      "SS_sc7785175",
      "SS_sc7785176",
      "SS_sc7785178",
      "SS_sc7785179",
      "SS_sc7785180",
      "SS_sc7785181",
      "SS_sc7785183",
      "SS_sc7785184",
      "SS_sc7785185",
      "SS_sc7785187",
      "SS_sc7785188",
      "SS_sc7785190",
      "SS_sc7785192",
      "SS_sc7785193",
      "SS_sc7785195",
      "SS_sc7785197",
      "SS_sc7785199",
      "SS_sc7785200",
      "SS_sc7785203",
      "SS_sc7785204",
      "SS_sc7785205",
      "SS_sc7785206",
      "SS_sc7785208",
      "SS_sc7785209",
      "SS_sc7785210",
      "SS_sc7785211",
      "SS_sc7785212",
      "SS_sc7785213",
      "SS_sc7785214",
      "SS_sc7785215",
      "SS_sc7785216",
      "SS_sc7785217",
      "SS_sc7785218",
      "SS_sc7785220",
      "SS_sc7785221",
      "SS_sc7785222",
      "SS_sc7785223",
      "SS_sc7785224",
      "SS_sc7785226",
      "SS_sc7785227",
      "SS_sc7785228",
      "SS_sc7785231",
      "SS_sc7785233",
      "SS_sc7785234",
      "SS_sc7785235",
      "SS_sc7785237",
      "SS_sc7785239",
      "SS_sc7785241",
      "SS_sc7785242",
      "SS_sc7785243",
      "SS_sc7785245",
      "SS_sc7785255",
      "SS_sc7785257",
      "SS_sc7785259",
      "SS_sc7785260",
      "SS_sc7785262",
      "SS_sc7785265",
      "SS_sc7785266",
      "SS_sc7785267",
      "SS_sc7785268",
      "SS_sc7785269"
     ],
     "scoring": {
      "pass_threshold": 0.85
     },
     "answer_field": "wells_kept"
    }
   }
  ]
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent fits the QC to plate-based full-length data with spike-ins: in a FACS-sorted Smart-seq2 plate, the ERCC share of the reads, the mapping rate and the number of genes detected separate wells that hold a cell's RNA from wells that hold little more than the spike-in mix, and wells annotated as holding no cell or several cells are not single cells.

2) Primary trap: Droplet-style defaults are applied (keep wells with at least 200 detected genes and a mitochondrial share under 20%). Of the 364 wells annotated as single cells, 109 hold too little cellular RNA: their median ERCC share of read pairs is 40.9% (against 1.2% in the kept wells), their median mapped fraction 0.56 (against 0.76) and their median number of detected genes 793 (against 8,512). 105 of them pass the default filters, so 378 wells are kept instead of 255, 18 wells annotated as no cell or multiple cells among them. No mitochondrial filter can catch them: the highest mitochondrial share of read pairs in any well is 3.5%.

3) Scripted controls, run on the task data:
   reference reference: The paper's five thresholds (genes above 10 rpm > 3,000, log10 read pairs > 4, mapped > 0.5, mitochondrial < 0.1, ERCC < 0.3, as fractions of read pairs) on the wells annotated as single cells. -> n_wells_kept=255, wells_kept=255 labels (passes)
   variant thresholds_all_wells: The same thresholds, wells annotated as no cell or multiple cells not set aside. -> n_wells_kept=263, wells_kept=263 labels (passes)
   variant counted_reads_only: The thresholds on shares of counted reads, without read pairs or the mapping filter. -> n_wells_kept=259, wells_kept=259 labels (passes)
   variant gene_floor_2000: The thresholds with a gene floor of 2,000 instead of 3,000. -> n_wells_kept=261, wells_kept=261 labels (passes)
   variant mad_outliers: 3-MAD outliers on log read pairs, log genes detected, ERCC share and mitochondrial share, among single-cell wells. -> n_wells_kept=248, wells_kept=248 labels (passes)
   trap droplet_defaults: Droplet defaults: at least 200 genes detected and a mitochondrial share of counts under 20%. -> n_wells_kept=378, wells_kept=378 labels (fails)
   trap genes_and_mito_only: More than 1,000 genes detected and a mitochondrial share under 10%, on single-cell wells; spike-in share and mapping rate not used. -> n_wells_kept=306, wells_kept=306 labels (fails)

4) Tolerance: The windows accept the lab's five thresholds with or without setting aside the wells annotated as no cell or multiple cells (8 of them pass the thresholds: 7 no cell and 1 multiple cells, with profiles like good single cells), the same thresholds on shares of counted reads when read pairs are not used, a gene floor of 2,000, and 3-MAD outlier calls among single-cell wells. They do not accept QC that ignores the spike-ins and the mapping rate: droplet defaults, or a floor of 1,000 genes with a mitochondrial cap. The windows are a compromise: a strict floor of 3,000 genes above 10 rpm with depth and mitochondrial cuts but no spike-in or mapping cut keeps 282 wells and fails only narrowly, and an ERCC cut taken as a share of read pairs without the mapping cut lets most of the poorly mapped wells through (their unmapped reads dilute the spike-in share).

5) Not said in the prompt: The prompt does not mention the spike-in share, the mapping rate or any threshold, and does not say that many wells annotated as single cells are nearly empty or that some wells annotated as no cell look like good cells.

6) Provenance: The paper's criteria (STAR Methods): more than 3,000 genes above 10 reads per million, log10 total reads above 4, mapped fraction above 0.5, mitochondrial fraction below 0.1 and ERCC fraction below 0.3, fractions taken of read pairs. On the public counts with read pairs from the ENA run report (ERP124375, read count / 2) they keep 263 wells; setting aside the 8 annotated as no cell or multiple cells leaves the 255 wells the authors marked as passing in the ArrayExpress SDRF, well for well. Of the 109 single-cell wells removed, 82 fail the gene floor, 61 the ERCC share and 40 the mapping rate. Well annotations are the SDRF's single cell well quality column (16 no cell, 4 multiple cells). Gene symbols: the authors' normalized file lists symbols in Ensembl ID order, checked against the 726 genes it names by ID.
```

## `ave02_pseudotime_root`: Which end of the embryonic VE trajectory is the AVE

From the paper: Fig. 1D and STAR Methods, Diffusion pseudotime analysis of AVE and emVE cells (952 genes at 5.5 dpc; high-in-AVE and low-in-AVE groups)

```
`norm_data_55.h5ad` holds Smart-seq2 profiles of the 255 cells from 5.5 dpc mouse embryos that passed quality control (`X`: scran-normalized expression, not log-transformed; variable names are gene symbols). `embryonic_ve_e55.txt` lists the 173 of them that form the embryonic visceral endoderm, the VE over the epiblast. At this stage that population contains the anterior visceral endoderm (AVE), at the distal tip with high Cer1, Hhex and Lefty1, and the rest of the embryonic VE (emVE).

Work on these 173 cells only, with ln(x + 1) of the normalized values. The lab orders them along a diffusion pseudotime that starts at the AVE end of the trajectory, with this pipeline: the 3,000 most highly variable genes among the 173 cells; PCA; a k-nearest-neighbour graph (k = 15) on the first 10 principal components, with the distance sqrt((1 - r) / 2), r being the Spearman correlation between two cells; a diffusion map; diffusion pseudotime (DPT), which needs a root cell.

Then test every gene detected (value above 0) in at least 10 of the 173 cells for a change of its ln expression along pseudotime: regress it on a natural cubic spline of pseudotime with 4 degrees of freedom, F-test against an intercept-only model, Benjamini-Hochberg correction across genes. A gene with FDR < 0.01 is high-in-AVE if its Spearman correlation with pseudotime is negative (it is highest at the start) and low-in-AVE if it is positive.

Report the Spearman correlation between pseudotime and the ln expression of Cer1 across the 173 cells, the fraction of the FDR < 0.01 genes that are high-in-AVE, and which of these genes are high-in-AVE: Cd24a, Ctsb, Ctsz, Eps8, Gpc4, Grb10, Morc4, Pcyt2, Sfrp1, Tspan7.

Return EXACTLY:

<EVAL_ANSWER>
{"rho_pseudotime_cer1": <float>, "frac_high_in_ave": <float>, "high_in_ave_panel": <list of strings>}
</EVAL_ANSWER>

`rho_pseudotime_cer1`: Spearman correlation between pseudotime and ln(x + 1) of Cer1 across the 173 cells
`frac_high_in_ave`: share of the FDR < 0.01 genes that are high-in-AVE, a fraction between 0 and 1
`high_in_ave_panel`: the genes of the listed panel that are high-in-AVE
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
      "rho_pseudotime_cer1": -0.7887,
      "frac_high_in_ave": 0.3769
     },
     "tolerances": {
      "rho_pseudotime_cer1": {
       "type": "absolute",
       "value": 0.25
      },
      "frac_high_in_ave": {
       "type": "absolute",
       "value": 0.09
      }
     }
    }
   },
   {
    "type": "label_set_jaccard",
    "config": {
     "ground_truth_labels": [
      "Cd24a",
      "Gpc4",
      "Grb10",
      "Sfrp1",
      "Tspan7"
     ],
     "scoring": {
      "pass_threshold": 0.6
     },
     "answer_field": "high_in_ave_panel"
    }
   }
  ]
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent anchors the pseudotime in the AVE and checks which population its root cell belongs to, instead of trusting the extreme of a diffusion component or the first cell of the file.

2) Primary trap: The root is taken at an extreme of the first diffusion component, or at the first cell of the file, without checking where it sits. The sign of a diffusion component is arbitrary; here the extreme taken and the first cell (SS.sc7784898) are emVE cells, so the pseudotime runs from the emVE to the AVE: the Cer1 correlation is +0.81 instead of -0.79, 61.8% of the significant genes are called high-in-AVE instead of 37.7%, and every panel gene lands in the wrong group, with no error.

3) Scripted controls, run on the task data:
   reference reference: Root at the end of the first diffusion component where the Cer1, Hhex and Lefty1 score is high; natural spline F-test with 4 df. -> rho_pseudotime_cer1=-0.7887, frac_high_in_ave=0.3769, high_in_ave_panel=5 labels (passes)
   variant root_max_marker_score: Root at the cell with the highest Cer1, Hhex and Lefty1 score. -> rho_pseudotime_cer1=-0.8064, frac_high_in_ave=0.3834, high_in_ave_panel=5 labels (passes)
   variant euclidean_knn: Euclidean distance on the 10 principal components instead of the Spearman distance. -> rho_pseudotime_cer1=-0.8348, frac_high_in_ave=0.3863, high_in_ave_panel=5 labels (passes)
   variant scaled_pca: Genes scaled to unit variance before the PCA. -> rho_pseudotime_cer1=-0.7556, frac_high_in_ave=0.3898, high_in_ave_panel=5 labels (passes)
   variant root_dc2_extreme_in_ave: The paper's recipe, an extreme of the second diffusion component, taken at the end with the higher AVE marker score. -> rho_pseudotime_cer1=-0.5764, frac_high_in_ave=0.3411, high_in_ave_panel=5 labels (passes)
   variant authors_root_cell: The cell the authors rooted at (minimum of their second component, from DC_55.csv), placed in this graph. -> rho_pseudotime_cer1=-0.5528, frac_high_in_ave=0.4476, high_in_ave_panel=5 labels (passes)
   variant linear_trend_test: Linear regression on pseudotime instead of the spline. -> rho_pseudotime_cer1=-0.7887, frac_high_in_ave=0.379, high_in_ave_panel=5 labels (passes)
   variant spearman_test: Spearman correlation test instead of the spline. -> rho_pseudotime_cer1=-0.7887, frac_high_in_ave=0.4047, high_in_ave_panel=5 labels (passes)
   variant spline_3df: Natural spline with 3 degrees of freedom. -> rho_pseudotime_cer1=-0.7887, frac_high_in_ave=0.3848, high_in_ave_panel=5 labels (passes)
   trap root_at_emve_end: Root at the extreme of the first diffusion component that lies in the emVE, as an unchecked argmin or argmax gives when the sign falls that way. -> rho_pseudotime_cer1=0.8125, frac_high_in_ave=0.6178, high_in_ave_panel=5 labels (fails)
   trap root_first_cell: Root at the first cell of the file (iroot = 0), an emVE cell. -> rho_pseudotime_cer1=0.7197, frac_high_in_ave=0.6017, high_in_ave_panel=5 labels (fails)

4) Tolerance: The task is about where the pseudotime starts, not which AVE cell is the root, so the windows accept any root at an AVE extreme of the diffusion map: the AVE end of the first component, the end of the second component in the AVE (the paper's recipe), the cell the authors used (SS.sc7784978), or the cell with the highest mean z-score of Cer1, Hhex and Lefty1; also the Euclidean distance, genes scaled before the PCA, and a linear, 3-df spline or Spearman test of each gene. They do not accept a pseudotime that starts in the emVE, nor one rooted at an AVE cell that sits among the emVE cells in the graph: used as roots one at a time, 32 of the 95 AVE cells give a Cer1 correlation outside the window (median over all AVE roots -0.61).

5) Not said in the prompt: The prompt says the pseudotime starts at the AVE end but not how to find that end. It does not say that a diffusion component has no fixed sign or that the first cell of the file is an emVE cell.

6) Provenance: The paper roots the pseudotime at the cell with the minimum (5.5 dpc) or maximum (6.25 dpc) of the second diffusion component, which in its run lay in the AVE; the authors' coordinates (code repository, DATA/DC_55.csv) show that their first component also separates AVE from emVE, while a rerun returns the components in another order and sign, so the task asks for the AVE end and not for a component. Cells: the authors' AVE (95) and emVE (78) clusters from the ArrayExpress SDRF, merged and not labelled. With the reference root (the AVE end of the first diffusion component, chosen by its correlation with a Cer1, Hhex and Lefty1 score, rho 0.78), the AVE cells have mean pseudotime 0.47 and the emVE cells 0.76; 1,064 of 14,642 genes have FDR < 0.01 (401 high-in-AVE), against 952 in the paper with R gam's loess smoother. Panel: the 5 most negatively and the 5 most positively correlated FDR < 0.01 genes of the reference, Cer1, Hhex and Lefty1 left out.
```

## `ave03_ephrin_bootstrap_signs`: Which ephrin and Eph genes carry the pathway's main axes, with PCA signs bootstrapped correctly

From the paper: Fig. S5A-D and STAR Methods, Analysis of genes in the Ephrin/Eph-signalling pathway (14 and 13 genes kept; 6 and 8 significant)

```
`ephrin_family_expression.csv` holds the scran-normalized expression (not log-transformed) of the 22 ephrin and Eph receptor genes (Efna1-5, Efnb1-3, Epha1-8, Epha10, Ephb1-4, Ephb6) in two Smart-seq2 datasets of mouse embryos enriched for visceral endoderm: 252 cells at 5.5 dpc and 235 cells at 6.25 dpc (column `stage`), each with its cluster (`cell_type`: AVE, emVE, exVE, Epi or ExE).

The lab wants to know which of these genes drive the main axes of variation of the Ephrin/Eph pathway at each stage. Analyse each stage separately, on ln(x + 1) values:

1. Keep the genes whose mean ln(x + 1) over that stage's cells is above 0.5.
2. Run a PCA on the kept genes (centred, not scaled) and take the loadings of the first 3 principal components.
3. Draw 1,000 bootstrap samples of the cells (with replacement) and rerun the PCA on each. Give each loading of the first 3 components a p-value: the fraction of bootstrap samples in which that loading has the opposite sign to the same loading from the full data.
4. Call a gene significant if at least one of its 3 loadings has p < 0.01.

Report, for each stage, how many genes you kept, how many are significant, and which.

Return EXACTLY:

<EVAL_ANSWER>
{"n_genes_kept_e55": <integer>, "n_significant_e55": <integer>, "significant_genes_e55": <list of strings>, "n_genes_kept_e625": <integer>, "n_significant_e625": <integer>, "significant_genes_e625": <list of strings>}
</EVAL_ANSWER>

`n_genes_kept_e55`: genes kept at 5.5 dpc (mean ln(x + 1) above 0.5)
`n_significant_e55`: significant genes at 5.5 dpc
`significant_genes_e55`: symbols of the significant genes at 5.5 dpc
`n_genes_kept_e625`: genes kept at 6.25 dpc
`n_significant_e625`: significant genes at 6.25 dpc
`significant_genes_e625`: symbols of the significant genes at 6.25 dpc
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
      "n_genes_kept_e55": 14,
      "n_significant_e55": 6,
      "n_genes_kept_e625": 13,
      "n_significant_e625": 8
     },
     "tolerances": {
      "n_genes_kept_e55": {
       "type": "absolute",
       "value": 1
      },
      "n_significant_e55": {
       "type": "absolute",
       "value": 2
      },
      "n_genes_kept_e625": {
       "type": "absolute",
       "value": 1
      },
      "n_significant_e625": {
       "type": "absolute",
       "value": 2
      }
     }
    }
   },
   {
    "type": "label_set_jaccard",
    "config": {
     "ground_truth_labels": [
      "Efna3",
      "Efna5",
      "Efnb1",
      "Efnb2",
      "Epha4",
      "Ephb3"
     ],
     "scoring": {
      "pass_threshold": 0.7
     },
     "answer_field": "significant_genes_e55"
    }
   },
   {
    "type": "label_set_jaccard",
    "config": {
     "ground_truth_labels": [
      "Efna2",
      "Efna3",
      "Efna5",
      "Efnb1",
      "Efnb2",
      "Epha1",
      "Epha4",
      "Ephb3"
     ],
     "scoring": {
      "pass_threshold": 0.7
     },
     "answer_field": "significant_genes_e625"
    }
   }
  ]
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent knows that a principal axis has no intrinsic sign, so each bootstrap axis must be aligned with the full-data axis (reflected when it points the other way) before the signs of its loadings are compared.

2) Primary trap: Signs are compared as each PCA returns them. A bootstrap axis comes back pointing either way, so every loading flips in a large share of the samples (at 5.5 dpc the median p-value is 0.50 and the smallest 0.37) and no gene is significant at 5.5 dpc, against 6 once each axis is reflected to match the full data. The sign test then measures the solver's sign convention, not the stability of the loadings.

3) Scripted controls, run on the task data:
   reference reference: Each bootstrap axis reflected to match the full-data axis before signs are compared (numpy SVD, 1,000 samples, seed 0). -> n_genes_kept_e55=14, n_significant_e55=6, significant_genes_e55=6 labels, n_genes_kept_e625=13 (passes)
   variant axis_reordering_too: Bootstrap axes also matched to the full-data axes by their loadings (Peres-Neto reordering) before reflection. -> n_genes_kept_e55=14, n_significant_e55=7, significant_genes_e55=7 labels, n_genes_kept_e625=13 (passes)
   variant other_seed: Another random seed. -> n_genes_kept_e55=14, n_significant_e55=6, significant_genes_e55=6 labels, n_genes_kept_e625=13 (passes)
   variant sklearn_pca: scikit-learn PCA (full SVD) with reflection. -> n_genes_kept_e55=14, n_significant_e55=6, significant_genes_e55=6 labels, n_genes_kept_e625=13 (passes)
   variant two_thousand_samples: 2,000 bootstrap samples. -> n_genes_kept_e55=14, n_significant_e55=6, significant_genes_e55=6 labels, n_genes_kept_e625=13 (passes)
   trap signs_not_aligned: Loading signs compared as the SVD returns them, axes not aligned. -> n_genes_kept_e55=14, n_significant_e55=0, significant_genes_e55=0 labels, n_genes_kept_e625=13 (fails)
   trap sklearn_signs_not_aligned: Loading signs compared under scikit-learn's own sign convention, axes not aligned. -> n_genes_kept_e55=14, n_significant_e55=0, significant_genes_e55=0 labels, n_genes_kept_e625=13 (fails)

4) Tolerance: The windows accept reflection with or without matching bootstrap axes to the full-data axes by their loadings (axis reordering, which adds Ephb1 at 5.5 dpc), another seed, 2,000 samples, and scikit-learn's PCA with reflection. They do not accept unaligned signs, which leave no significant gene at 5.5 dpc with numpy's or scikit-learn's PCA (scikit-learn's own sign convention keeps 6 of the 8 genes at 6.25 dpc, so the 5.5 dpc fields are the ones that catch it).

5) Not said in the prompt: The prompt defines the p-value as a comparison of signs and does not say that principal axes have arbitrary sign or mention reflection.

6) Provenance: The paper used scanpy's PCA and 1,000 bootstrap samples with 'the fraction of sign inversions ... accounting for possible axis reflections, as described by Peres-Neto et al.', and reports 14 and 13 genes above the expression cut with 6 and 8 significant genes; the reference (numpy SVD, seed 0) reproduces all four numbers. The first 3 components explain 38.2% of the variance at 5.5 dpc and 46.1% at 6.25 dpc. Cells: the authors' labelled cells (the 3 unassigned cells of each stage left out), clusters from the ArrayExpress SDRF; expression from the authors' normalized files.
```

## `ave04_dkk1_timing`: Is Dkk1 an AVE signal before or after migration

From the paper: Fig. S3C and Results (Dkk1 elevated in the anterior only at 6.25 dpc, arguing against DKK1 guiding AVE migration)

```
Two Smart-seq2 datasets of mouse embryos enriched for visceral endoderm (VE), both after quality control: `norm_data_55.h5ad` (5.5 dpc, when the anterior visceral endoderm (AVE) is at the distal tip or migrating towards the future anterior) and `norm_data_625.h5ad` (6.25 dpc, when the AVE has reached the boundary between epiblast and extra-embryonic ectoderm and its proximal migration is over). In both, `X` is scran-normalized expression, not log-transformed, with gene symbols as variable names. `cell_labels.csv` gives the stage and cluster of each cell (AVE, emVE, exVE, Epi or ExE); cells that are not listed there were left unassigned and are not used.

Earlier work proposed that DKK1 secreted by the AVE guides its migration, and a collaborator asks whether these data support that. Compare AVE cells with emVE cells, the embryonic VE through which the AVE migrates. For Dkk1, and for Cer1 as a reference AVE gene, report at each stage the AUC of ln(x + 1) expression for AVE against emVE cells: the probability that a randomly chosen AVE cell has higher expression than a randomly chosen emVE cell, ties counted as one half. Then choose the conclusion the data best support.

Return EXACTLY:

<EVAL_ANSWER>
{"dkk1_auc_e55": <float>, "dkk1_auc_e625": <float>, "cer1_auc_e55": <float>, "cer1_auc_e625": <float>, "conclusion": <integer>}
</EVAL_ANSWER>

`dkk1_auc_e55`: AUC of Dkk1, AVE against emVE, at 5.5 dpc
`dkk1_auc_e625`: AUC of Dkk1, AVE against emVE, at 6.25 dpc
`cer1_auc_e55`: AUC of Cer1, AVE against emVE, at 5.5 dpc
`cer1_auc_e625`: AUC of Cer1, AVE against emVE, at 6.25 dpc

`conclusion`: Which conclusion do these data best support? Answer with the option number.
1. The data support DKK1 as a guidance cue: pooled over both stages, Dkk1 is clearly higher in the AVE than in the emVE, so the AVE releases DKK1 while it migrates.
2. The data argue against AVE-derived DKK1 starting the migration: Dkk1 is no higher in the AVE than in the emVE at 5.5 dpc and becomes AVE-enriched only at 6.25 dpc.
3. Dkk1 is enriched in the AVE at both stages about as strongly as Cer1, so it marks the AVE from its specification onwards and could guide the migration from the start.
4. Dkk1 is AVE-enriched at 6.25 dpc and DKK1 is an established AVE guidance cue, so this enrichment confirms that it directs the migration that comes before this stage.
5. Dkk1 is lower in the AVE than in the emVE at both stages, so any guidance by DKK1 would come from the emVE, pushing AVE cells away rather than being sent by them.
6. Dkk1 is barely detected in the visceral endoderm at 5.5 dpc, so the data say nothing about its role at that stage, and the question has to wait for other experiments.
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
      "dkk1_auc_e55": 0.5176,
      "dkk1_auc_e625": 0.7778,
      "cer1_auc_e55": 0.8278,
      "cer1_auc_e625": 0.8991
     },
     "tolerances": {
      "dkk1_auc_e55": {
       "type": "absolute",
       "value": 0.06
      },
      "dkk1_auc_e625": {
       "type": "absolute",
       "value": 0.06
      },
      "cer1_auc_e55": {
       "type": "absolute",
       "value": 0.06
      },
      "cer1_auc_e625": {
       "type": "absolute",
       "value": 0.06
      }
     }
    }
   },
   {
    "type": "numeric_tolerance",
    "config": {
     "ground_truth": {
      "conclusion": 2
     },
     "tolerances": {
      "conclusion": {
       "type": "absolute",
       "value": 0
      }
     }
    }
   }
  ]
 }
}
```

Notes:

```
1) Scientific intent: Test whether the agent reads the timing of the AVE enrichment against the timing of the migration, rather than taking DKK1's published place among AVE signals, or an enrichment pooled over stages, as support.

2) Primary trap: The prior-driven reading: Dkk1 is AVE-enriched at 6.25 dpc (AUC 0.78) and, pooled over both stages, still looks AVE-enriched (AUC 0.62, Wilcoxon p = 4.0e-04); DKK1 is a published guidance candidate, so this is read as support. At 5.5 dpc, while the AVE migrates, Dkk1 is no higher in the AVE than in the emVE: mean ln expression 2.70 against 2.70 (AUC 0.52, p = 0.69).

3) Scripted controls, run on the task data:
   reference reference: AUC of ln(x + 1), AVE against emVE, at each stage; conclusion chosen from the measured pattern. -> dkk1_auc_e55=0.5176, dkk1_auc_e625=0.7778, cer1_auc_e55=0.8278, cer1_auc_e625=0.8991 (passes)
   variant log2_scale: log2(x + 1) instead of ln(x + 1). -> dkk1_auc_e55=0.5176, dkk1_auc_e625=0.7778, cer1_auc_e55=0.8278, cer1_auc_e625=0.8991 (passes)
   variant roc_auc_score: AUC from scikit-learn's roc_auc_score. -> dkk1_auc_e55=0.5176, dkk1_auc_e625=0.7778, cer1_auc_e55=0.8278, cer1_auc_e625=0.8991 (passes)
   trap stages_pooled: AVE and emVE cells of both stages pooled, and the pooled enrichment read as support for guidance. -> dkk1_auc_e55=0.6246, dkk1_auc_e625=0.6246, cer1_auc_e55=0.8536, cer1_auc_e625=0.8536 (fails)
   trap prior_driven_reading: Correct AUCs, but the enrichment at 6.25 dpc and DKK1's published role taken as support. -> dkk1_auc_e55=0.5176, dkk1_auc_e625=0.7778, cer1_auc_e55=0.8278, cer1_auc_e625=0.8991 (fails)

4) Tolerance: AUCs are rank statistics, so the windows (+/- 0.06) accept any monotone transform of the expression and any standard AUC routine. They do not accept AUCs from cells pooled over the two stages, which put Dkk1 at 0.62 for both.

5) Not said in the prompt: The prompt asks for the numbers by stage but does not say which stage bears on the hypothesis or what an AUC near 0.5 implies for it.

6) Provenance: The paper reaches the same timing from its pseudotime fits and HCR (Fig. S3) and concludes that it argues against a role for DKK1 in symmetry breaking and in guiding AVE migration. Key chosen by script from the measured pattern: Dkk1 AUC within 0.5 +/- 0.1 at 5.5 dpc (0.518) and at least 0.65 at 6.25 dpc (0.778, p = 2.6e-06); Dkk1 is detected in 81.5% of the 5.5 dpc AVE and emVE cells (so option 6 does not hold); Cer1 AUC 0.83 and 0.90, while Dkk1 is not enriched at 5.5 dpc (so option 3 does not hold) and neither Dkk1 AUC is at or below 0.35 (so option 5 does not hold). Cells: 95 AVE and 78 emVE at 5.5 dpc, 63 and 39 at 6.25 dpc, clusters from the ArrayExpress SDRF; expression from the authors' normalized files.
```
