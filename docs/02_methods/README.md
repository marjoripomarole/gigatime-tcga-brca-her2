# Methods Map

This folder is a navigation layer for methods. The generated or hand-written reports still live at the top level of `docs/` for script compatibility.

## Data And Labels

- `tcga_her2_label_quality_assessment.md`
- `clinical_her2_trustworthy_slide_list.md`
- `clinical_her2_laptop_balanced61_selection.md`
- `clinical_her2_cohort_expanded20_selection.md`
- `clinical_her2_cohort_selection.md`

Core scripts:

- `scripts/gdc_query_tcga_brca.py`
- `scripts/build_tcga_brca_clinical_her2_labels.py`
- `scripts/select_clinical_her2_cohort.py`
- `scripts/build_tcga_her2_trustworthy_slide_list.py`
- `scripts/download_clinical_her2_cohort_slides.py`
- `scripts/download_selected_star_counts.py`

## Image Inference

- `scripts/run_gigatime_tcga_brca.py`
- `scripts/run_hoptimus_tcga_brca.py`

## Analysis And Validation

- `scripts/cleanup_gigatime_tile_features.py`
- `scripts/train_her2_cleaned_classifier_comparison.py`
- `scripts/analyze_tissue_composition_sensitivity.py`
- `scripts/analyze_tumor_proxy_sensitivity.py`
- `scripts/analyze_clinical_covariate_sensitivity.py`
- `scripts/analyze_matched_low_zero_sensitivity.py`
- `scripts/analyze_source_site_generalization.py`
- `scripts/analyze_within_source_site_low_zero.py`
- `scripts/analyze_local_erbb2_expression_validation.py`
- `scripts/audit_her2_isoform_validation_feasibility.py`

## Visual QC And Reports

- `scripts/render_case_driver_visual_qc.py`
- `scripts/render_clinical_her2_visual_qc.py`
- `scripts/render_virtual_mif_channel_images.py`
- `scripts/render_virtual_mif_composites.py`
- `scripts/build_clinical_her2_findings_report.py`
- `scripts/build_her2_classifier_results_report.py`
