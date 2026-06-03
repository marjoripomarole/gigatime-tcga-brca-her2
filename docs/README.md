# Documentation Index

Status: docs navigation map.

Start with `00_start_here.md`. This file is the compact index for the rest of the documentation.

## Start Here

- `00_start_here.md` - current project map and next scientific steps.
- `RUN_REGISTRY.md` - run-by-run evidence trail.
- `advisor_brief.md` - concise advisor-facing summary.
- `plain_language_methodology.md` - accessible methodology explanation.
- `paper_proposal_process_log.md` - process history and rationale.

## Current Results

Use `03_current_results/README.md` for the current-result reading order.

Core current reports:

- `clinical_her2_high_trust_tile128_results.md`
- `clinical_her2_high_trust_tile128_case_driver_analysis.md`
- `clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
- `clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`
- `clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`
- `clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`
- `clinical_her2_high_trust_tile128_source_site_generalization.md`
- `clinical_her2_high_trust_tile128_local_erbb2_validation.md`
- `her2_isoform_validation_feasibility.md`
- `her2_isoform_state_hypothesis.md`

## Methods

Use `02_methods/README.md` for the methods and script map.

Important method and cohort files:

- `tcga_her2_label_quality_assessment.md`
- `clinical_her2_trustworthy_slide_list.md`
- `clinical_her2_laptop_balanced61_selection.md`
- `virtual_mif_channel_outputs.md`

## Model Experiments

Use `04_model_experiments/README.md` for model-family tests beyond the primary GigaTIME result.

- `hoptimus_embedding_baseline.md`
- `histoprism_one_vector_smoke.md`
- `deepspot_one_vector_smoke.md`

## Archive

Use `90_archive/README.md` for historical 30-slide and 60-slide reports. These files remain useful for provenance but should not be cited as the latest result.

## Why Top-Level Reports Still Exist

Many scripts currently write directly to top-level `docs/*.md` paths. To avoid breaking reruns, this cleanup adds navigation folders while preserving existing generated-report paths.

A later refactor can physically move reports into subfolders after script output defaults and internal links are updated together.
