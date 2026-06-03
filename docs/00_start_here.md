# Start Here

Status: current navigation spine for the project.

## One-Sentence Project

This project tests whether TCGA-BRCA H&E image features associate with clinically defined HER2-low versus HER2-zero breast cancer states, while explicitly stress-testing tissue-composition, source-site, slide-size, RNA, and model-family caveats.

## Current Claim

The current primary result is a hypothesis-generating tissue-context association:

- In the strict high-trust 171-slide TCGA-BRCA cohort, GigaTIME virtual immune/myeloid/checkpoint and CK-associated channels differ between HER2-low and HER2-zero.
- The signal survives several internal checks and shuffled-label sanity tests.
- It is not yet safe as independent HER2 biology because slide-size, TCGA source-site, and tissue-composition confounding remain strong.

## Read First

1. `clinical_her2_high_trust_tile128_results.md` - current primary results.
2. `advisor_brief.md` - concise advisor-facing narrative.
3. `RUN_REGISTRY.md` - run-by-run evidence trail.
4. `plain_language_methodology.md` - accessible methodology explanation.
5. `her2_isoform_state_hypothesis.md` - careful biological framing and language guardrails.

## Navigation Folders

- `02_methods/README.md` - methods and rerun entry points.
- `03_current_results/README.md` - current evidence and caveats.
- `04_model_experiments/README.md` - H-Optimus, HistoPrism, DeepSpot, and related model tests.
- `90_archive/README.md` - historical 30-slide and 60-slide reports.

Existing report files currently remain at the top level of `docs/` because many scripts write those paths directly. The folder README files are curated maps over the existing report set.

## Next Best Scientific Steps

1. Get a pathologist or tumor-region review loop around the case-driver tiles.
2. Run H0-mini/H-Optimus embeddings as a generic H&E baseline when access/download is ready.
3. Compare model families on the same high-trust HER2-low/HER2-zero split.
4. Keep DeepSpot/HistoPrism as interpretive gene-expression-style follow-ups, not as primary validation.
5. Obtain transcript-level or junction-level RNA evidence before making any HER2 isoform claim.
