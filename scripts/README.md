# Scripts Map

This folder contains standalone workflow scripts. Most scripts are intentionally runnable from the repo root and write reports to existing `docs/*.md` paths for compatibility.

## Data And Label Construction

- `gdc_query_tcga_brca.py` - query GDC TCGA-BRCA slides and STAR-count RNA files.
- `build_tcga_brca_clinical_her2_labels.py` - build clinical HER2-positive/low/zero labels from TCGA clinical fields.
- `build_bcnb_her2_labels.py` - build BCNB HER2-positive/low/zero labels from the gated clinical workbook.
- `audit_bcnb_image_inputs.py` - inspect BCNB WSI/patch files and patient-ID mapping before model runs.
- `probe_xenium_breast_rna_validation.py` - feasibility probe for the 10x Xenium Human Breast RNA-validation cohort; downloads only the minimal artifacts (gene panel, H&E alignment, transcripts) and checks channel-gene coverage, alignment invertibility, per-channel transcript counts, and H&E/transcript extent overlap. Writes `docs/xenium_breast_rna_validation_probe.md`.
- `download_bcnb_wsis.py` - download approved BCNB full-WSI files into ignored `data/bcnb/WSIs`; supports SharePoint/OneDrive folder manifests, resumable file-by-file downloads, numeric patient-JPG filtering, and Google Drive fallback via `gdown`.
- `build_bcnb_wsi_slide_table.py` - build a patient/clinical-label slide table from locally downloaded BCNB full WSIs.
- `build_bcnb_patch_manifest.py` - build patient-linked manifests for BCNB precomputed patches, including deterministic capped patch sampling.
- `select_clinical_her2_cohort.py` - select balanced HER2 cohorts.
- `download_clinical_her2_cohort_slides.py` - download selected diagnostic slides.
- `download_selected_star_counts.py` - download selected STAR-count RNA files.
- `build_tcga_her2_trustworthy_slide_list.py` - label, metadata, file, and OpenSlide trust checks.

## Image Feature Extraction

- `run_gigatime_tcga_brca.py` - primary GigaTIME virtual mIF feature extraction; also supports BCNB flat `.jpg` WSIs via `--slide-backend pil` and copies selected slide-table clinical metadata into `slide_scores.csv`.
- `run_hoptimus_tcga_brca.py` - H0-mini/H-Optimus H&E embedding extraction.
- `run_virchow2_tcga_brca.py` - Virchow2 H&E embedding extraction (second embedding control).
- `run_bcnb_patch_embeddings.py` - extract patient-level BCNB patch embeddings from capped patch manifests using H-Optimus/H0-mini/Virchow2.

## Current Analysis And Sensitivity Checks

- `cleanup_gigatime_tile_features.py`
- `train_her2_cleaned_classifier_comparison.py`
- `analyze_high_trust_her2_sensitivity.py`
- `analyze_high_trust_case_drivers.py`
- `analyze_tissue_composition_sensitivity.py`
- `analyze_tumor_proxy_sensitivity.py`
- `analyze_classifier_permutation_sanity.py`
- `analyze_nested_classifier_model_selection.py`
- `analyze_clinical_covariate_sensitivity.py`
- `analyze_matched_low_zero_sensitivity.py`
- `analyze_source_site_generalization.py`
- `analyze_within_source_site_low_zero.py`
- `analyze_local_erbb2_expression_validation.py`
- `audit_her2_isoform_validation_feasibility.py`
- `analyze_hoptimus_embedding_control.py`
- `analyze_bcnb_patch_embedding_control.py` - patient-level BCNB HER2-low versus HER2-zero patch-embedding analysis with clinical covariate controls and shuffled-label null.
- `analyze_bcnb_patch_model_comparison.py` - paired BCNB H-Optimus-0 versus Virchow2 comparison, dual-model ensemble check, and cross-model score agreement.
- `analyze_bcnb_patch_stratified_performance.py` - BCNB clinical-slice robustness check for patient-mean out-of-fold image and clinical model scores across grade, ER/PR, subtype, nodal status, and Ki67 strata.
- `analyze_bcnb_patch_score_covariate_drivers.py` - explains BCNB image-model patient scores using HER2 label, clinical covariates, and patch QC, then tests clinical/patch-QC residual low-vs-zero signal.
- `validate_gigatime_xenium_rna.py` - within-slide validation of GigaTIME virtual channels against Xenium breast spatial RNA: tiles the post-Xenium H&E on GigaTIME's grid, bins transcripts onto the same tiles via the inverse H&E alignment affine, and reports per-channel virtual-vs-RNA Spearman with a spatial block-bootstrap CI, a channel-by-gene specificity matrix, and cellularity-controlled partial correlations. Includes a model-free `--alignment-check-only` tissue-vs-transcript sanity mode. Writes `docs/xenium_breast_rna_validation_results.md`.
- `make_gigatime_vs_rna_specificity_figure.py` - draft paper figure joining GigaTIME's published per-channel held-out-mIF agreement (Valanarasu et al., Cell 2026, Fig S5) with our Xenium RNA specificity (raw + cellularity-controlled partial r); shows that mIF agreement does not imply RNA channel specificity. Writes `docs/gigatime_vs_rna_specificity_comparison.md` assets.

## Historical Analysis

- `summarize_her2_gigatime.py`
- `summarize_clinical_her2_gigatime.py`
- `validate_gigatime_with_rna_signatures.py`
- `validate_gigatime_with_rna_programs.py`
- `train_her2_classifier_baseline.py`
- `compare_gigatime_run_agreement.py`

## Visuals And Reports

- `render_he_slide_images.py`
- `render_virtual_mif_channel_images.py`
- `render_virtual_mif_composites.py`
- `render_clinical_her2_visual_qc.py`
- `render_case_driver_visual_qc.py`
- `render_bcnb_patch_score_visual_qc.py` - render BCNB hash-capped patch montages for image-score extreme cases.
- `build_clinical_her2_findings_report.py`
- `build_her2_classifier_results_report.py`

## Model Smoke Tests

- `run_histoprism_one_vector_smoke.py`
- `run_deepspot_one_vector_smoke.py`
- `run_virchow2_one_slide_smoke.py`

## Cleanup Note

The next structural cleanup could move scripts into subfolders such as `scripts/data/`, `scripts/models/`, `scripts/analysis/`, and `scripts/reports/`. That should be done together with README command updates and any downstream automation changes.
