# Documentation Guide

Status: Start here. This is the navigation map for the docs folder.

This folder contains both current summaries and historical analysis reports. Start here if you are new to the project.

## Read First

1. `clinical_her2_expanded20_results.md`
   - Best current presentation summary.
   - Covers the expanded 60-slide run: 20 HER2-positive, 20 HER2-low, 20 HER2-zero.
   - Use this for the latest results.

2. `advisor_brief.md`
   - Short advisor-facing summary.
   - Good for meetings and high-level discussion.

3. `plain_language_methodology.md`
   - Explains the project for someone without genetics, pathology, or AI background.
   - Best teaching document.

4. `paper_proposal_process_log.md`
   - The living history/process log.
   - This is the document that records what we did over time and why.

5. `her2_isoform_state_hypothesis.md`
   - The paper-proposal framing around HER2 state, isoform hypotheses, targetability, and careful language.

## Current Latest Results

The latest top-line result is the expanded 20/20/20 clinical HER2 run:

- 60 TCGA-BRCA slides total.
- 20 HER2-positive, 20 HER2-low, 20 HER2-zero.
- Up to 256 tissue tiles per slide.
- 15,225 total GigaTIME tile predictions.
- STAR-count RNA-seq expression available for all 60 cases.

The strongest current finding:

- GigaTIME/H&E features continue to separate HER2-low from HER2-zero.
- Several HER2-low versus HER2-zero virtual immune/myeloid/checkpoint channel differences pass within-view BH correction in all-tissue or QC-cellular views.
- The HER2-low versus HER2-zero classifier remains around balanced accuracy 0.800 and macro AUC 0.820.
- RNA validation remains weak, so this is still hypothesis-generating and not clinical diagnosis.

## Current Detailed Expanded Reports

- `clinical_her2_expanded20_results.md`: current best findings summary.
- `clinical_her2_cohort_expanded20_selection.md`: current 20/20/20 cohort selection.
- `clinical_her2_expanded20_gigatime_data_cleanup.md`: current expanded cleanup/tile-filtering report.
- `clinical_her2_expanded20_cleaned_classifier_comparison.md`: current expanded cleaned classifier comparison.

## Historical 30-Slide Reports

These are still useful because they show how the project developed, but they should not be cited as the latest result:

- `clinical_her2_cohort_selection.md`: original 10/10/10 cohort.
- `clinical_her2_gigatime_run.md`: original 30-slide GigaTIME run.
- `clinical_her2_tile_sampling_robustness.md`: 30-slide 256-tile robustness check.
- `clinical_her2_rna_validation.md`: 30-slide marker-level RNA validation.
- `clinical_her2_rna_program_validation.md`: 30-slide broader RNA program validation.
- `clinical_her2_classifier_baseline.md`: 30-slide first classifier baseline.
- `clinical_her2_gigatime_data_cleanup.md`: 30-slide cleanup report.
- `clinical_her2_cleaned_classifier_comparison.md`: 30-slide cleaned classifier report.
- `clinical_her2_visual_qc.md`: initial 30-slide visual QC.

## Visual Explanation Files

- `virtual_mif_channel_outputs.md`: explains virtual mIF-style images and channel visualizations.
- `assets/`: tracked figures used by the markdown reports.

## How To Use This Folder

For a new researcher:

1. Read `clinical_her2_expanded20_results.md`.
2. Read `plain_language_methodology.md` if the biology or workflow is unfamiliar.
3. Read `paper_proposal_process_log.md` to understand the history.
4. Use historical reports only when you need details from a specific earlier analysis.

For a presentation:

1. Lead with `clinical_her2_expanded20_results.md`.
2. Use `advisor_brief.md` for a concise narrative.
3. Use `her2_isoform_state_hypothesis.md` for careful biological framing.
