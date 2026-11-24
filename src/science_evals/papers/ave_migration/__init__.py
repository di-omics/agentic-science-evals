"""Tasks from Thowfeequ, Fiorentino et al., Dev Cell 2024: murine anterior visceral endoderm migration.

PAPER        citation, the repository owner's role, data, code, licence
DATA_FILES   every input, public (pinned URL) or derived here (pinned sha256)
build()      fetches and derives the inputs, runs every control, returns one Built per task
anchors()    statements from the paper recomputed from the public files

See README.md in this folder for the tasks and how to rebuild them.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...controls import Control
from .. import Built
from . import inputs as io
from . import tasks as specs

PAPER = {
    "short": "AVE migration, Dev Cell 2024",
    "citation": (
        "Thowfeequ S, Fiorentino J, Hu D, Solovey M, Ruane S, Whitehead M, Zhou F, Godwin J, "
        "Mateo-Otero Y, Vanhaesebroeck B, Scialdone A, Srinivas S. An integrated approach identifies the "
        "molecular underpinnings of murine anterior visceral endoderm migration. Developmental Cell 59, "
        "2347-2363.e9 (2024). doi:10.1016/j.devcel.2024.05.014"
    ),
    "doi": "10.1016/j.devcel.2024.05.014",
    "author_role": (
        "Di Hu (D.H.), from the paper's author contributions: \"D.H. and S.R. did sample collections for "
        "the transcriptomics and proteomics, respectively\". D.H. did the sample collections for the "
        "transcriptomics (the single-cell RNA-seq)."
    ),
    "data": (
        "ArrayExpress E-MTAB-9645 (processed counts and SDRF); ENA ERP124375 run report (read counts per "
        "run, no reads); the authors' normalized expression in github.com/ScialdoneLab/scAVE"
    ),
    "code": (
        f"https://github.com/ScialdoneLab/scAVE at commit {io.SCAVE_COMMIT} (tags 1.0.0 and 1.0.1; "
        "archived as Zenodo record 10829359, doi:10.5281/zenodo.10829359)"
    ),
    "license": (
        "Paper: CC BY 4.0. ArrayExpress and ENA data: public under EMBL-EBI terms of use. Code repository "
        "archive (Zenodo 10829359): CC BY 4.0."
    ),
}

DATA_FILES = [dict(e) for e in io.PUBLIC] + io.derived_entries()


def _controls(answers: dict[str, dict], kinds: dict[str, tuple[str, str]]) -> list[Control]:
    return [Control(name, kinds[name][0], answers[name], kinds[name][1]) for name in kinds]


def _files(data_dir: Path, names: list[str]) -> list[dict]:
    return [io.entry(data_dir, n) for n in names]


def build(data_dir: Path) -> list[Built]:
    from . import reference as ref

    data_dir = Path(data_dir)
    io.prepare(data_dir)
    built = []

    # ave01: plate QC
    f1 = ref.ave01_facts(data_dir)
    a1 = ref.ave01(data_dir)
    t1 = specs.ave01(f1)
    built.append(Built(t1, a1["reference"], _controls(a1, {
        "reference": ("reference", "The paper's five thresholds (genes above 10 rpm > 3,000, log10 read pairs > 4, "
                                   "mapped > 0.5, mitochondrial < 0.1, ERCC < 0.3, as fractions of read pairs) on "
                                   "the wells annotated as single cells."),
        "thresholds_all_wells": ("variant", "The same thresholds, wells annotated as no cell or multiple cells "
                                            "not set aside."),
        "counted_reads_only": ("variant", "The thresholds on shares of counted reads, without read pairs or "
                                          "the mapping filter."),
        "gene_floor_2000": ("variant", "The thresholds with a gene floor of 2,000 instead of 3,000."),
        "mad_outliers": ("variant", "3-MAD outliers on log read pairs, log genes detected, ERCC share and "
                                    "mitochondrial share, among single-cell wells."),
        "droplet_defaults": ("trap", "Droplet defaults: at least 200 genes detected and a mitochondrial share "
                                     "of counts under 20%."),
        "genes_and_mito_only": ("trap", "More than 1,000 genes detected and a mitochondrial share under 10%, on "
                                        "single-cell wells; spike-in share and mapping rate not used."),
    }), _files(data_dir, t1.data)))

    # ave02: pseudotime root
    f2 = ref.ave02_facts(data_dir)
    panel, _ = ref.ave02_panel(str(data_dir))
    a2 = ref.ave02(data_dir)
    rhos = f2.pop("ave_root_rhos")
    centre = round(a2["reference"]["rho_pseudotime_cer1"], 4)
    f2["n_ave_roots_outside"] = sum(abs(r - centre) > specs.AVE02_RHO_TOL for r in rhos)
    f2["median_ave_root_rho"] = float(np.median(rhos))
    t2 = specs.ave02(f2, panel)
    built.append(Built(t2, a2["reference"], _controls(a2, {
        "reference": ("reference", "Root at the end of the first diffusion component where the Cer1, Hhex and "
                                   "Lefty1 score is high; natural spline F-test with 4 df."),
        "root_max_marker_score": ("variant", "Root at the cell with the highest Cer1, Hhex and Lefty1 score."),
        "euclidean_knn": ("variant", "Euclidean distance on the 10 principal components instead of the "
                                     "Spearman distance."),
        "scaled_pca": ("variant", "Genes scaled to unit variance before the PCA."),
        "root_dc2_extreme_in_ave": ("variant", "The paper's recipe, an extreme of the second diffusion "
                                               "component, taken at the end with the higher AVE marker score."),
        "authors_root_cell": ("variant", "The cell the authors rooted at (minimum of their second component, "
                                         "from DC_55.csv), placed in this graph."),
        "linear_trend_test": ("variant", "Linear regression on pseudotime instead of the spline."),
        "spearman_test": ("variant", "Spearman correlation test instead of the spline."),
        "spline_3df": ("variant", "Natural spline with 3 degrees of freedom."),
        "root_at_emve_end": ("trap", "Root at the extreme of the first diffusion component that lies in the "
                                     "emVE, as an unchecked argmin or argmax gives when the sign falls that way."),
        "root_first_cell": ("trap", "Root at the first cell of the file (iroot = 0), an emVE cell."),
    }), _files(data_dir, t2.data)))

    # ave03: ephrin bootstrap
    f3 = ref.ave03_facts(data_dir)
    a3 = ref.ave03(data_dir)
    f3["n_sig_ref_e55"] = a3["reference"]["n_significant_e55"]
    t3 = specs.ave03(f3)
    built.append(Built(t3, a3["reference"], _controls(a3, {
        "reference": ("reference", "Each bootstrap axis reflected to match the full-data axis before signs are "
                                   "compared (numpy SVD, 1,000 samples, seed 0)."),
        "axis_reordering_too": ("variant", "Bootstrap axes also matched to the full-data axes by their loadings "
                                           "(Peres-Neto reordering) before reflection."),
        "other_seed": ("variant", "Another random seed."),
        "sklearn_pca": ("variant", "scikit-learn PCA (full SVD) with reflection."),
        "two_thousand_samples": ("variant", "2,000 bootstrap samples."),
        "signs_not_aligned": ("trap", "Loading signs compared as the SVD returns them, axes not aligned."),
        "sklearn_signs_not_aligned": ("trap", "Loading signs compared under scikit-learn's own sign convention, "
                                              "axes not aligned."),
    }), _files(data_dir, t3.data)))

    # ave04: Dkk1 timing
    f4 = ref.ave04_facts(data_dir)
    a4, key = ref.ave04(data_dir)
    r4 = a4["reference"]
    f4.update({"auc_dkk1_e55": r4["dkk1_auc_e55"], "auc_dkk1_e625": r4["dkk1_auc_e625"],
               "auc_cer1_e55": r4["cer1_auc_e55"], "auc_cer1_e625": r4["cer1_auc_e625"]})
    t4 = specs.ave04(f4, key)
    built.append(Built(t4, r4, _controls(a4, {
        "reference": ("reference", "AUC of ln(x + 1), AVE against emVE, at each stage; conclusion chosen from "
                                   "the measured pattern."),
        "log2_scale": ("variant", "log2(x + 1) instead of ln(x + 1)."),
        "roc_auc_score": ("variant", "AUC from scikit-learn's roc_auc_score."),
        "stages_pooled": ("trap", "AVE and emVE cells of both stages pooled, and the pooled enrichment read as "
                                  "support for guidance."),
        "prior_driven_reading": ("trap", "Correct AUCs, but the enrichment at 6.25 dpc and DKK1's published "
                                         "role taken as support."),
    }), _files(data_dir, t4.data)))
    return built


def anchors(data_dir: Path) -> list[dict]:
    from . import reference as ref

    data_dir = Path(data_dir)
    io.prepare(data_dir)
    out = []

    qc = ref.qc_both_stages(data_dir)
    e55, e625 = qc["e55"], qc["e625"]
    out.append({
        "paper_says": "255 good-quality cells at 5.5 dpc and 238 at 6.25 dpc after the five QC thresholds; "
                      "252 and 235 after 3 unassigned cells per stage are removed",
        "where": "STAR Methods, Transcript quantification, quality control and normalization; Cell clustering",
        "recomputed": f"The five thresholds on the public counts, with read pairs from the ENA run report, keep "
                      f"{e55['thresholds']} of {e55['wells']} wells at 5.5 dpc ({e55['single_cell_wells']} once "
                      f"the plate's control wells are set aside; {e55['agree_with_authors']} of {e55['wells']} "
                      f"wells then match the authors' pass or fail in the SDRF) and {e625['thresholds']} at "
                      f"6.25 dpc ({e625['agree_with_authors']} of {e625['wells']} match). The SDRF labels "
                      f"{e55['labelled']} and {e625['labelled']} cells.",
        "agrees": e55["single_cell_wells"] == 255 and e625["single_cell_wells"] == 238
                  and e55["labelled"] == 252 and e625["labelled"] == 235,
    })

    n_ave, n_emve, agree = ref.pam_split(data_dir)
    out.append({
        "paper_says": "At 5.5 dpc the embryonic VE cluster is split into AVE and emVE with pam (k = 2) on the "
                      "Spearman correlation distance",
        "where": "STAR Methods, Cell clustering and assignment of cluster identities",
        "recomputed": f"pam (k = 2) on sqrt((1 - r) / 2) over the 3,000 most variable genes of the 173 embryonic "
                      f"VE cells gives {n_ave} AVE and {n_emve} emVE cells; {100 * agree:.1f}% of cells match the "
                      f"authors' labels in the SDRF (95 AVE, 78 emVE).",
        "agrees": agree == 1.0,
    })

    g = ref.dpt_gene_counts(data_dir)
    close = all(abs(g[s]["n_sig"] - p) / p <= 0.15 for s, p in (("e55", 952), ("e625", 915)))
    out.append({
        "paper_says": "952 genes at 5.5 dpc and 915 at 6.25 dpc change along the diffusion pseudotime "
                      "(GAM, FDR < 0.01)",
        "where": "STAR Methods, Diffusion pseudotime analysis of AVE and emVE cells",
        "recomputed": f"{g['e55']['n_sig']} of {g['e55']['n_tested']} and {g['e625']['n_sig']} of "
                      f"{g['e625']['n_tested']} genes with a natural cubic spline F-test (4 df) in place of R "
                      f"gam's loess smoother, root at the AVE end of the first diffusion component; "
                      f"{g['e55']['n_high']} and {g['e625']['n_high']} of them high-in-AVE. Within 15% of the "
                      f"paper: {'yes' if close else 'no'}; the smoother and the test differ, so exact counts are "
                      f"not expected.",
        "agrees": close,
    })

    eph = ref.ave03(data_dir)["reference"]
    ok = (eph["n_genes_kept_e55"], eph["n_genes_kept_e625"], eph["n_significant_e55"],
          eph["n_significant_e625"]) == (14, 13, 6, 8)
    out.append({
        "paper_says": "14 ephrin/Eph genes at 5.5 dpc and 13 at 6.25 dpc have mean log-normalized expression "
                      "above 0.5; 6 and 8 of them have a significant PCA loading (bootstrap, p < 0.01)",
        "where": "STAR Methods, Analysis of genes in the Ephrin/Eph-signalling pathway; Fig. S5B",
        "recomputed": f"{eph['n_genes_kept_e55']} and {eph['n_genes_kept_e625']} genes; "
                      f"{eph['n_significant_e55']} significant at 5.5 dpc ({', '.join(eph['significant_genes_e55'])}) "
                      f"and {eph['n_significant_e625']} at 6.25 dpc ({', '.join(eph['significant_genes_e625'])}), "
                      f"with each bootstrap axis reflected to the full-data axis.",
        "agrees": ok,
    })

    f4 = ref.ave04_facts(data_dir)
    aucs = ref.stage_aucs(str(data_dir))
    out.append({
        "paper_says": "Unlike Cer1, Dkk1 became elevated in the anterior only at 6.25 dpc",
        "where": "Results, AP asymmetries; Fig. S3C",
        "recomputed": f"Dkk1 AUC, AVE against emVE: {aucs['dkk1_auc_e55']:.2f} at 5.5 dpc (Wilcoxon p = "
                      f"{f4['Dkk1_e55_p']:.2f}) and {aucs['dkk1_auc_e625']:.2f} at 6.25 dpc (p = "
                      f"{f4['Dkk1_e625_p']:.1e}); Cer1 {aucs['cer1_auc_e55']:.2f} and {aucs['cer1_auc_e625']:.2f}.",
        "agrees": abs(aucs["dkk1_auc_e55"] - 0.5) <= ref.NULL_BAND and aucs["dkk1_auc_e625"] >= ref.ENRICHED,
    })

    sp = ref.sema6d_plxna1(data_dir)
    s, p = sp["Sema6d"], sp["Plxna1"]
    out.append({
        "paper_says": "Sema6d is expressed at higher levels in the AVE, while its receptor Plxna1 is mainly "
                      "expressed in the surrounding emVE and Epi (5.5 dpc)",
        "where": "Results, Semaphorin signaling; Fig. 6C",
        "recomputed": "Mean ln(x + 1) at 5.5 dpc: Sema6d " + ", ".join(f"{k} {v:.2f}" for k, v in s.items())
                      + "; Plxna1 " + ", ".join(f"{k} {v:.2f}" for k, v in p.items()) + ".",
        "agrees": s["AVE"] > max(s["emVE"], s["Epi"]) and p["AVE"] < min(p["emVE"], p["Epi"]),
    })
    return out
