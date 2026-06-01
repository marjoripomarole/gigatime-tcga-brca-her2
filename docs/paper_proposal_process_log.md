# Paper Proposal Process Log

Last updated: 2026-06-01

This document keeps a running record of the research process for a future paper or grant proposal. It is written to preserve both the scientific reasoning and the concrete computational steps used in this TCGA-BRCA GigaTIME HER2 project.

## Working Project Title

Evaluating virtual tumor immune microenvironment predictions from H&E histopathology across the HER2 axis in TCGA breast cancer.

## Central Research Question

Can an existing histopathology foundation model, GigaTIME, infer biologically interpretable virtual multiplex immunofluorescence features from TCGA-BRCA H&E slides, and do those predicted features differ across HER2-related breast cancer states?

The current project began with an ERBB2 RNA-expression pilot. The proposed next version should move toward clinically meaningful HER2 groups:

- HER2-positive
- HER2-low
- HER2-zero

This distinction matters because HER2-low and HER2-zero are defined by clinical protein-level testing, not simply by the amount of ERBB2 RNA measured in sequencing data.

## Plain-Language Framing

Breast cancer slides stained with H&E show tissue structure, but they do not directly show all immune-cell markers or molecular markers. GigaTIME is a released AI model that takes H&E image tiles and predicts virtual multiplex immunofluorescence-like marker maps. In this project, those predictions are used as research features that may describe the tumor immune microenvironment.

The study asks whether these predicted immune and tumor marker patterns vary along the HER2 axis. HER2 is clinically important in breast cancer because it affects biology and treatment options. The project is not trying to diagnose HER2 from H&E. Instead, it asks whether a previously released AI model produces tissue-microenvironment signals that correlate with known HER2 biology.

## Data Sources

### TCGA-BRCA

The data source is TCGA-BRCA, the breast invasive carcinoma project from The Cancer Genome Atlas. The files are accessed through the Genomic Data Commons (GDC).

The current workflow uses:

- Diagnostic H&E whole-slide images in `.svs` format.
- RNA-seq STAR-count files used to extract ERBB2 expression.
- Public clinical supplement files used to investigate HER2 IHC/ISH status.

Official references:

- GDC API overview: https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/
- GDC clinical supplement description: https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Supplement/

### GigaTIME

The project uses the official released GigaTIME implementation and model weights. The model predicts 23 virtual mIF channels from H&E pathology image tiles:

`DAPI`, `TRITC`, `Cy5`, `PD-1`, `CD14`, `CD4`, `T-bet`, `CD34`, `CD68`, `CD16`, `CD11c`, `CD138`, `CD20`, `CD3`, `CD8`, `PD-L1`, `CK`, `Ki67`, `Tryptase`, `Actin-D`, `Caspase3-D`, `PHH3-B`, and `Transgelin`.

These outputs are model predictions, not laboratory-measured mIF.

## Methods Completed So Far

### 1. Queried TCGA-BRCA Files From GDC

The script `scripts/gdc_query_tcga_brca.py` was created to query GDC for:

- TCGA-BRCA diagnostic H&E slide images.
- TCGA-BRCA STAR-count RNA-seq expression files.

Main local outputs:

- `data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv`
- `data/tcga_brca/tcga_brca_diagnostic_slides_files.csv`
- `data/tcga_brca/tcga_brca_star_counts_manifest.tsv`
- `data/tcga_brca/tcga_brca_star_counts_files.csv`
- `data/tcga_brca/file_metadata_slides.json`
- `data/tcga_brca/file_metadata_star_counts.json`

### 2. Extracted ERBB2 Expression

The workflow downloaded selected TCGA-BRCA STAR-count files and extracted ERBB2 expression from each file.

The ERBB2 gene was identified as:

- Gene symbol: `ERBB2`
- Ensembl gene ID: `ENSG00000141736`

Main local output:

- `data/tcga_brca/erbb2_expression.csv`

Current status:

- ERBB2 expression was extracted for 80 TCGA-BRCA cases.
- These expression values were used as the first HER2-biology proxy.

### 3. Selected ERBB2-Extreme Cases

The script `scripts/select_her2_extremes.py` selected the top and bottom ERBB2 TPM cases from the current expression pilot.

Main local output:

- `data/tcga_brca/her2_extreme_cases.csv`

Current selected cohort:

- 20 ERBB2-high cases.
- 20 ERBB2-low cases.

Important caveat: these labels are expression-based research labels. They do not mean confirmed clinical HER2-positive, HER2-low, or HER2-zero status.

### 4. Downloaded and Processed a Slide Subset

The workflow attempted to download H&E whole-slide images for the selected ERBB2-extreme cases. Slide downloads were slow and occasionally unstable, so the current processed subset is smaller than the target cohort.

Current processed pilot:

- 12 TCGA-BRCA diagnostic slides.
- 7 ERBB2-high slides.
- 5 ERBB2-low slides.
- 64 random tissue tiles per slide.
- 768 total tile predictions.
- CPU inference.

Main local outputs:

- `results/gigatime_tcga_brca_extremes/slide_scores.csv`
- `results/gigatime_tcga_brca_extremes/tile_scores.csv`
- `results/gigatime_tcga_brca_extremes/heatmaps/`

### 5. Ran GigaTIME on H&E Tiles

The script `scripts/run_gigatime_tcga_brca.py` tiles each slide, filters for tissue-containing tiles, runs GigaTIME, and aggregates virtual mIF channel predictions.

For each processed slide, the workflow stores:

- Tile-level channel activations.
- Slide-level mean activations.
- Slide-level fraction-positive or thresholded summaries for each channel.
- Heatmap-style visual outputs.

### 6. Summarized Virtual mIF Features by ERBB2 Group

The script `scripts/summarize_her2_gigatime.py` joins slide-level GigaTIME predictions to ERBB2 expression and compares ERBB2-high versus ERBB2-low groups.

Main local outputs:

- `results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/advisor_summary.md`
- Summary figures in `results/gigatime_tcga_brca_extremes/advisor_summary/`

This analysis is exploratory because the current processed subset contains only 12 slides.

### 7. Generated Virtual mIF Channel Figures

The script `scripts/render_virtual_mif_channel_images.py` renders documentation-facing all-channel GigaTIME figures.

Main local outputs:

- `docs/assets/virtual_mif_channels/virtual_mif_all_channel_group_means.png`
- `docs/assets/virtual_mif_channels/virtual_mif_slide_channel_matrix.png`
- `docs/assets/virtual_mif_channels/her2_high_reference_all_virtual_mif_channels.png`
- `docs/assets/virtual_mif_channels/her2_low_reference_all_virtual_mif_channels.png`

These figures show predicted channels across groups, slides, and spatial tile positions.

### 8. Generated Fluorescence-Style Virtual mIF Composites

The script `scripts/render_virtual_mif_composites.py` reruns GigaTIME on selected tiles and combines predicted channel maps into black-background fluorescence-style composites.

Main local outputs:

- `docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_virtual_mif_montage.png`
- `docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_virtual_mif_montage.png`
- `docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_he_vs_virtual_mif.png`
- `docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_he_vs_virtual_mif.png`
- Additional tumor/proliferation and myeloid/B-cell virtual panels.

These images are closer in appearance to real mIF images than the dot-grid figures, but they remain virtual predictions from H&E.

### 9. Built a Reproducible Clinical HER2 Label Table

The script `scripts/build_tcga_brca_clinical_her2_labels.py` now queries the GDC TCGA-BRCA clinical supplement, downloads the patient-level BCR Biotab, extracts HER2 IHC/ISH fields, and writes a traceable clinical HER2 label table.

Command:

```bash
conda run -n gigatime-tcga python scripts/build_tcga_brca_clinical_her2_labels.py
```

Main local outputs:

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt`

The generated label table includes one row per TCGA-BRCA clinical case, the raw HER2 IHC/ISH values, the assigned clinical HER2 group, and the exact rule used to assign that group.

### 10. Selected a Balanced Clinical HER2 Pilot Cohort

The script `scripts/select_clinical_her2_cohort.py` now joins the clinical HER2 label table with ERBB2 expression and slide metadata, then selects a deterministic balanced cohort for the next GigaTIME run.

Command:

```bash
conda run -n gigatime-tcga python scripts/select_clinical_her2_cohort.py
```

Main local outputs:

- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_summary.json`
- `docs/clinical_her2_cohort_selection.md`

Selection priority:

- Clinical HER2 group must be one of HER2-positive, HER2-low, or HER2-zero.
- Direct clinical HER2 labels are preferred over inferred labels.
- Already-downloaded slides are preferred.
- Smaller slide files are preferred to make the next pilot more practical.
- Case IDs are used for deterministic tie-breaking.

Selected cohort:

| Cohort group | Selected cases | Already downloaded slides |
|---|---:|---:|
| HER2-positive | 10 | 4 |
| HER2-low | 10 | 3 |
| HER2-zero | 10 | 1 |

This gives the first clean 30-case clinical HER2 pilot cohort for running GigaTIME across HER2-positive, HER2-low, and HER2-zero disease.

## Initial Biological Findings From the ERBB2-Extreme Pilot

The current processed dataset is too small for strong claims. The main result so far is that the workflow is feasible and produces interpretable tables and figures.

Current working interpretation:

- GigaTIME can be run on TCGA-BRCA H&E slides.
- Its 23 predicted virtual mIF channels can be aggregated by slide and compared with ERBB2 expression.
- The current figures are useful for advisor discussion and proposal planning.
- The current analysis should not yet be framed as a definitive HER2 biological result.

The correct proposal language is "exploratory pilot" rather than "validated biomarker analysis."

## Investigation of Clinical HER2-Zero, HER2-Low, and HER2-Positive Labels

The advisor-facing research direction requires true clinical HER2 groups, not only ERBB2 RNA extremes.

HER2-low is generally defined as:

- IHC `1+`, or
- IHC `2+` with ISH negative.

The FDA describes HER2-low breast cancer using this IHC/ISH definition in the context of trastuzumab deruxtecan eligibility. Reference: https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-fam-trastuzumab-deruxtecan-nxki-unresectable-or-metastatic-hr-positive-her2-low-or-her2

### Clinical Fields Found in TCGA-BRCA

The TCGA-BRCA GDC clinical patient Biotab contains HER2-related fields, including:

- `lab_proc_her2_neu_immunohistochemistry_receptor_status`
- `her2_erbb_pos_finding_cell_percent_category`
- `her2_immunohistochemistry_level_result`
- `pos_finding_her2_erbb2_other_measurement_scale_text`
- `lab_procedure_her2_neu_in_situ_hybrid_outcome_type`
- `her2_neu_chromosone_17_signal_ratio_value`

This means TCGA-BRCA can support a better clinical HER2 grouping than the original ERBB2-expression-only pilot.

### Implemented Clinical HER2 Mapping

The implemented mapping is:

- `HER2-positive`: IHC `3+`, ISH positive, or receptor status positive when detailed fields are missing.
- `HER2-low`: IHC `1+`, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH.
- `HER2-unknown`: missing, not evaluated, contradictory, or incomplete HER2 data.

This mapping should be documented carefully in the methods section because TCGA clinical supplement fields are not always complete.

### Counts From the Clinical HER2 Label Table

Using the implemented mapping above:

| Dataset | HER2-positive | HER2-low | HER2-zero | HER2-unknown |
|---|---:|---:|---:|---:|
| TCGA-BRCA clinical rows | 174 | 407 | 61 | 455 |
| Cases with current slide metadata and ERBB2 expression | 13 | 35 | 10 | 22 |
| Current 40 ERBB2-extreme selected cases | 13 | 10 | 8 | 9 |
| Current 12 GigaTIME-processed slides | 6 | 3 | 1 | 2 |

Interpretation:

- TCGA-BRCA appears useful for clinical HER2 grouping.
- The current processed GigaTIME subset is too small for a three-group HER2-zero versus HER2-low versus HER2-positive comparison.
- The next scientific step is to select and process a balanced clinical HER2 cohort rather than continuing only with ERBB2 RNA extremes.

## Multiplex Immunofluorescence Comparison Question

We investigated whether TCGA has real multiplex immunofluorescence results matched to the H&E slides used here. The current working conclusion is that this project does not yet have matched real mIF for those TCGA slides.

Therefore:

- We cannot currently validate GigaTIME virtual mIF predictions by direct TCGA matched mIF comparison.
- The virtual mIF images should be presented as model-generated predictions, not as experimental ground truth.
- Trustworthiness needs to be assessed using indirect and external validation strategies.

Potential validation strategies:

- Compare GigaTIME-predicted immune channels with bulk expression signatures from RNA-seq.
- Compare predicted epithelial/tumor channels with pathology/tumor-purity annotations where available.
- Compare HER2 group trends with known breast cancer biology.
- Use an external dataset with paired H&E and real mIF, if available.
- Perform manual pathology review of selected H&E and virtual mIF panels.
- Check whether tile-level high-signal regions correspond to plausible tissue structures.

## Proposed Next Analyses

### Analysis 1: Re-select Cases by Clinical HER2 Group

This is now completed for a first balanced 10/10/10 clinical HER2 pilot. Instead of selecting only ERBB2-high and ERBB2-low cases, the workflow now creates a clinical HER2 cohort:

- HER2-positive cases.
- HER2-low cases.
- HER2-zero cases.

The selection balances slide availability and chooses one primary-tumor slide per case whenever possible.

### Analysis 2: Run GigaTIME on a Larger Clinical HER2 Cohort

The current 12-slide pilot is not enough. The proposal should target a larger, balanced run.

Possible first target:

- 30 HER2-positive slides.
- 30 HER2-low slides.
- 30 HER2-zero slides.

If slide downloads are unstable, start with a smaller but balanced pilot:

- 10 HER2-positive slides.
- 10 HER2-low slides.
- 10 HER2-zero slides.

### Analysis 3: Compare Virtual mIF Features Across Clinical HER2 Groups

For each GigaTIME channel, compare slide-level mean activations across HER2 groups.

Candidate statistical tests:

- Kruskal-Wallis test for three-group comparison.
- Pairwise Wilcoxon rank-sum tests for HER2-positive versus HER2-low, HER2-low versus HER2-zero, and HER2-positive versus HER2-zero.
- Benjamini-Hochberg FDR correction across channels.
- Effect sizes with confidence intervals where possible.

Primary endpoints should be pre-specified before scaling:

- Immune checkpoint channels: `PD-1`, `PD-L1`.
- T-cell channels: `CD3`, `CD8`, `CD4`.
- Myeloid/macrophage channels: `CD68`, `CD11c`, `CD14`.
- Tumor/proliferation channels: `CK`, `Ki67`.

### Analysis 4: Compare Virtual mIF With RNA-Seq Immune Signatures

Because matched real mIF is not currently available, use RNA-seq as an indirect validation layer.

Examples:

- Compare virtual `CD3` or `CD8` channels with T-cell gene signatures.
- Compare virtual `PD-L1` with `CD274` expression.
- Compare virtual macrophage-like channels with macrophage-related genes or signatures.
- Compare virtual `Ki67` with proliferation-related genes such as `MKI67`.

This does not prove the virtual mIF is correct, but it helps determine whether predicted tissue signals are directionally consistent with molecular data.

### Analysis 5: Visual QC and Trustworthiness Review

For selected cases from each HER2 group:

- Show source H&E tiles.
- Show GigaTIME virtual mIF composites.
- Show all-channel spatial maps.
- Flag artifacts, blank tissue, necrosis, folds, staining variation, and suspicious tile predictions.

This step is important because model outputs can look polished while still being wrong. The proposal should explicitly state that image-level QC is part of the methodology.

## Paper Proposal Structure

### Background

HER2 is a clinically important axis in breast cancer. The emergence of HER2-low as a therapeutically relevant category creates a need to understand whether tumor morphology and microenvironmental patterns differ across HER2-positive, HER2-low, and HER2-zero disease.

### Gap

Traditional H&E slides are widely available, but they do not directly provide multiplex immune-marker information. Real mIF is informative but expensive and not routinely available for large public cohorts. Virtual mIF models may provide a scalable way to generate hypotheses about immune microenvironment differences from existing pathology slides.

### Objective

Evaluate whether GigaTIME-derived virtual mIF features from TCGA-BRCA H&E slides differ across clinical HER2 groups and whether those predicted features show consistency with molecular and clinical annotations.

### Methods Overview

1. Retrieve TCGA-BRCA H&E slides, RNA-seq data, and clinical HER2 supplement fields from GDC.
2. Assign clinical HER2-positive, HER2-low, HER2-zero, or unknown labels using IHC/ISH fields.
3. Run GigaTIME on tissue-containing H&E tiles from selected diagnostic slides.
4. Aggregate virtual mIF channel predictions at tile and slide levels.
5. Compare GigaTIME channels across HER2 clinical groups.
6. Perform indirect validation against RNA-seq immune and proliferation signatures.
7. Generate visual QC panels and virtual mIF composites for interpretability.

### Expected Contribution

This study would not claim that GigaTIME diagnoses HER2 status. Instead, it would evaluate whether a released virtual mIF model can produce biologically interpretable immune microenvironment features from public breast cancer H&E slides and whether those features vary across clinically meaningful HER2 categories.

## Current Limitations to State Clearly

- The current processed cohort has only 12 slides.
- The current completed comparison uses ERBB2 RNA-expression extremes, not clinical HER2 groups.
- Clinical HER2 fields in TCGA are incomplete for many cases.
- TCGA clinical supplement files may contain missing, not evaluated, or inconsistent fields.
- No matched real mIF validation data is currently present in this project.
- GigaTIME predictions are research features, not clinical measurements.
- Whole-slide sampling, tile quality, tumor purity, and batch/stain variation need stronger QC.

## Reproducibility Checklist

Current workflow scripts:

- `scripts/gdc_query_tcga_brca.py`
- `scripts/build_tcga_brca_clinical_her2_labels.py`
- `scripts/select_clinical_her2_cohort.py`
- `scripts/select_her2_extremes.py`
- `scripts/run_gigatime_tcga_brca.py`
- `scripts/summarize_her2_gigatime.py`
- `scripts/render_he_slide_images.py`
- `scripts/render_virtual_mif_channel_images.py`
- `scripts/render_virtual_mif_composites.py`

Current documentation:

- `README.md`
- `docs/current_pilot_run.md`
- `docs/advisor_brief.md`
- `docs/plain_language_methodology.md`
- `docs/virtual_mif_channel_outputs.md`
- `docs/paper_proposal_process_log.md`
- `docs/clinical_her2_cohort_selection.md`

Current key result files:

- `data/tcga_brca/erbb2_expression.csv`
- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_summary.json`
- `data/tcga_brca/her2_extreme_cases.csv`
- `results/gigatime_tcga_brca_extremes/slide_scores.csv`
- `results/gigatime_tcga_brca_extremes/tile_scores.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv`

## Next Immediate Step

Use the generated clinical HER2 cohort manifest to download the missing selected slides and run GigaTIME on the 30-case cohort:

- 10 HER2-positive cases,
- 10 HER2-low cases,
- 10 HER2-zero cases.

The current selected cohort already has 8 slides downloaded locally and 22 slides still needing download. Once those missing slides are downloaded, run GigaTIME on this clinical HER2 cohort and summarize the virtual mIF channels by clinical group.
