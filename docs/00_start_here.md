# Start Here

Status: current navigation spine for the project.

## One-Sentence Project

This project tests whether TCGA-BRCA H&E image features associate with clinically defined HER2-low versus HER2-zero breast cancer states, while explicitly stress-testing tissue-composition, source-site, slide-size, RNA, and model-family caveats.

## Current Claim

The current primary result is a hypothesis-generating tissue-context association:

- In the strict high-trust 171-slide TCGA-BRCA cohort, GigaTIME virtual immune/myeloid/checkpoint and CK-associated channels differ between HER2-low and HER2-zero.
- The signal survives several internal checks and shuffled-label sanity tests.
- It is not yet safe as independent HER2 biology because slide-size, TCGA source-site, and tissue-composition confounding remain strong.
- As of 2026-06-04, two independent generic foundation-model embeddings (H-Optimus-0 and Virchow2) reproduce the low-versus-zero separation and the same source-site collapse, so the GigaTIME virtual-immune framing is not required to explain it. The most parsimonious reading is generic morphology tracking TCGA acquisition structure. TCGA-internal evidence is now considered exhausted.
- BCNB is confirmed as the immediate external validation cohort: full clinical data preserve HER2 IHC 0/1+/2+/3+, yielding 127 HER2-zero and 654 HER2-low cases from one scanner, with grade/ER/PR/Ki67 available for covariate control. The precomputed patch archive is local and audited. H-Optimus-0 and Virchow2 hash-capped patient-level patch pilots are complete and agree: the signal is statistically non-null but modest (H-Optimus-0 PCA20 BA 0.585-0.597 / AUC 0.634-0.640 across two hash seeds; Virchow2 BA 0.600 / AUC 0.643; dual-model BA 0.609 / AUC 0.651), comparable to clinical covariates rather than a strong standalone classifier. Multi-seed patch-sampling, clinical-stratified, score-driver, and visual-QC checks change the effect size but not the conclusion: the signal is weak, partly clinical/patch-QC-associated, and enriched for aggressive morphology/covariate profiles rather than clean HER2 status.
- H0-mini has been probed only as an access smoke and is currently blocked by Hugging Face gated-model access/token visibility, not by a code or data error.

## Read First

1. `clinical_her2_high_trust_tile128_results.md` - current primary results.
2. `advisor_brief.md` - concise advisor-facing narrative.
3. `RUN_REGISTRY.md` - run-by-run evidence trail.
4. `plain_language_methodology.md` - accessible methodology explanation.
5. `her2_isoform_state_hypothesis.md` - careful biological framing and language guardrails.
6. `clinical_her2_high_trust_tile128_hoptimus_embedding_control.md` and `clinical_her2_high_trust_tile128_virchow2_embedding_control.md` - generic-embedding confound controls (the decisive 2026-06-04 result).
7. `external_validation_candidates.md` and `bcnb_exploration.md` - external cohort path off TCGA; BCNB is now the confirmed priority cohort.

## Navigation Folders

- `02_methods/README.md` - methods and rerun entry points.
- `03_current_results/README.md` - current evidence and caveats.
- `04_model_experiments/README.md` - H-Optimus, HistoPrism, DeepSpot, and related model tests.
- `90_archive/README.md` - historical 30-slide and 60-slide reports.

Existing report files currently remain at the top level of `docs/` because many scripts write those paths directly. The folder README files are curated maps over the existing report set.

## Next Best Scientific Steps

TCGA-internal evidence is exhausted (see the 2026-06-04 entry in `paper_proposal_process_log.md`). The generic-embedding controls are done and confirm the confound. The remaining productive moves are:

1. External / single-scanner validation in BCNB: full WSIs are now local and a balanced GigaTIME smoke has passed, so the next paper-grade step is the full 781-slide low/zero WSI run with clinical and tissue/tile controls. See `external_validation_candidates.md`, `bcnb_exploration.md`, `bcnb_gigatime_full_wsi_smoke.md`, `bcnb_paper_patches_audit.md`, `bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero.md`, `bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero.md`, `bcnb_patch_embedding_control_virchow2_hash_capped10_low_zero.md`, `bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero.md`, `bcnb_patch_sampling_sensitivity_hoptimus0.md`, `bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero.md`, `bcnb_patch_score_covariate_drivers_hoptimus0_virchow2_hash_capped10_low_zero.md`, and `bcnb_patch_score_visual_qc_hoptimus0_virchow2_hash_capped10_low_zero.md`.
2. Write a cautionary-methods paper; the two-model embedding control is the centerpiece figure and Valieris et al. 2024 is the key comparison.
3. Validate GigaTIME virtual immune channels against real multiplex IHC (e.g. IMPRESS) to close the RNA-validation gap.
4. Get a pathologist or tumor-region review loop around the case-driver tiles.
5. Do not pull more TCGA-BRCA slides: HER2-zero is capped at 61 cases (already fully used), so more data only tightens estimates around the confound.
6. Keep DeepSpot/HistoPrism as interpretive gene-expression-style follow-ups, not as primary validation; obtain transcript-level/junction-level RNA before any HER2 isoform claim.
7. Revisit H0-mini only after `bioptimus/H0-mini` access is approved and the token is visible to `conda run -n gigatime-tcga`.
