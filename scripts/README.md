# Scripts Map

This folder contains standalone workflow scripts. Most scripts are intentionally runnable from the repo root and write reports to existing `docs/*.md` paths for compatibility.

## Data And Label Construction

- `gdc_query_tcga_brca.py` - query GDC TCGA-BRCA slides and STAR-count RNA files.
- `build_tcga_brca_clinical_her2_labels.py` - build clinical HER2-positive/low/zero labels from TCGA clinical fields.
- `select_clinical_her2_cohort.py` - select balanced HER2 cohorts.
- `download_clinical_her2_cohort_slides.py` - download selected diagnostic slides.
- `download_selected_star_counts.py` - download selected STAR-count RNA files.
- `build_tcga_her2_trustworthy_slide_list.py` - label, metadata, file, and OpenSlide trust checks.

## Image Feature Extraction

- `run_gigatime_tcga_brca.py` - primary GigaTIME virtual mIF feature extraction.
- `run_hoptimus_tcga_brca.py` - H0-mini/H-Optimus H&E embedding extraction.
- `run_virchow2_tcga_brca.py` - Virchow2 H&E embedding extraction (second embedding control).

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
- `build_clinical_her2_findings_report.py`
- `build_her2_classifier_results_report.py`

## Model Smoke Tests

- `run_histoprism_one_vector_smoke.py`
- `run_deepspot_one_vector_smoke.py`
- `run_virchow2_one_slide_smoke.py`

## Cleanup Note

The next structural cleanup could move scripts into subfolders such as `scripts/data/`, `scripts/models/`, `scripts/analysis/`, and `scripts/reports/`. That should be done together with README command updates and any downstream automation changes.
