# Documentation Index

Status: docs navigation map.

Start with `00_start_here.md`. This file is the compact index for the rest of the documentation.

## Start Here

- `manuscript_draft.md` - **working draft of the cautionary-methods manuscript** (synthesizes the confound result + RNA-specificity audit + two-model field-level comparison; for author revision).
- `00_start_here.md` - current project map and next scientific steps.
- `RUN_REGISTRY.md` - run-by-run evidence trail.
- `advisor_brief.md` - concise advisor-facing summary.
- `plain_language_methodology.md` - accessible methodology explanation.
- `paper_proposal_process_log.md` - process history and rationale.

## Strategic / Future Directions

- `her2_new_directions_exploration.md` - 2026-06-06 brainstorm: pivot off the confounded H&E HER2-status direction toward the lab's HER2 isoform/ADC-resistance edge (Guardia et al. 2025), and a proposed cheap HER2 "targetability" marker (ECD-IV vs ICD epitope index). Hypothesis-generating; direction still open.

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
- `hest_rna_validation_summary.md` - cross-sample (9 sections, 2 platforms: HEST-1k Xenium + Visium) RNA-specificity generalization of the GigaTIME virtual channels; per-sample reports are `hest_rna_validation_<id>.md`.
- `gigatime_vs_rosie_field_level.md` - two-model field-level comparison (GigaTIME vs ROSIE) on the same 9 sections: the virtual-channel specificity ceiling generalizes across models, which disagree on which channels are reliable (concordance r=0.12).
- `gigatime_orion_finetune_results.md` - controlled in-domain experiment (Orion CRC): the specificity ceiling is **domain shift, not intrinsic** — in-domain fine-tuning rescues every channel (CD68 −0.16→+0.53; 9/9 specific). Out-of-domain baseline: `gigatime_orion_baseline_CRC01.md`.

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
