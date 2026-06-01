# TCGA-BRCA GigaTIME HER2 Workflow

This project wraps the official GigaTIME implementation for a TCGA-BRCA pilot focused on HER2/ERBB2 expression.

The goal is to generate virtual multiplex immunofluorescence (mIF) features from TCGA-BRCA diagnostic H&E whole-slide images, join those features to clinical HER2 labels and ERBB2 RNA expression, and produce advisor-ready summary tables and figures.

## What This Contains

- `external/GigaTIME/`: cloned official GigaTIME code from `prov-gigatime/GigaTIME`.
- `scripts/gdc_query_tcga_brca.py`: queries GDC for TCGA-BRCA diagnostic slides and STAR-count RNA-seq files, writes GDC manifests, and can extract ERBB2 expression.
- `scripts/build_tcga_brca_clinical_her2_labels.py`: queries the GDC TCGA-BRCA clinical supplement and builds reproducible clinical HER2-positive/HER2-low/HER2-zero labels.
- `scripts/select_clinical_her2_cohort.py`: selects a balanced clinical HER2-positive/HER2-low/HER2-zero cohort and writes a slide manifest for the next GigaTIME run.
- `scripts/download_clinical_her2_cohort_slides.py`: downloads the selected clinical HER2 cohort slide files from GDC by file ID and writes a resumable status JSON.
- `scripts/run_gigatime_tcga_brca.py`: tiles TCGA-BRCA `.svs` slides, runs the official GigaTIME model, and aggregates virtual mIF channels per slide.
- `scripts/summarize_her2_gigatime.py`: joins GigaTIME slide scores with ERBB2 expression and makes HER2-high/HER2-low summary figures.
- `scripts/summarize_clinical_her2_gigatime.py`: compares GigaTIME virtual mIF outputs across clinical HER2-positive/HER2-low/HER2-zero groups.
- `scripts/validate_gigatime_with_rna_signatures.py`: compares GigaTIME virtual channels with matched RNA-seq marker signatures as an indirect validation check.
- `scripts/validate_gigatime_with_rna_programs.py`: compares GigaTIME virtual composite programs with broader RNA immune and tissue programs.
- `scripts/cleanup_gigatime_tile_features.py`: builds pre-classifier cleaned GigaTIME feature views from cellular and CK-enriched tile subsets.
- `scripts/train_her2_classifier_baseline.py`: trains first slide-level HER2 classifier baselines from GigaTIME features with leave-one-out cross-validation.
- `scripts/train_her2_cleaned_classifier_comparison.py`: reruns HER2 classifiers across all-tissue, cellular-tissue, and CK-enriched GigaTIME feature views.
- `scripts/render_virtual_mif_channel_images.py`: renders all-channel virtual mIF figures from GigaTIME tile and slide predictions.
- `scripts/render_virtual_mif_composites.py`: reruns GigaTIME on selected tiles and renders fluorescence-style virtual mIF composites from the full predicted channel maps.
- `scripts/render_clinical_her2_visual_qc.py`: renders clinical HER2 visual QC panels for cases driving high virtual `CD68`/`PD-L1`/`CD11c` signal.
- `scripts/build_clinical_her2_findings_report.py`: builds a simple display notebook and HTML report for the current clinical HER2 findings.
- `docs/virtual_mif_channel_outputs.md`: explains the generated virtual mIF channel images and how to interpret them.
- `docs/plain_language_methodology.md`: detailed non-specialist explanation of the study background, methodology, outputs, and current limitations.
- `docs/paper_proposal_process_log.md`: living process log for turning the pilot into a paper or grant proposal.
- `docs/clinical_her2_cohort_selection.md`: selected 30-case clinical HER2 pilot cohort and selection counts.
- `docs/clinical_her2_gigatime_run.md`: selected-cohort GigaTIME run status and full 30-slide clinical HER2 summary.
- `docs/clinical_her2_rna_validation.md`: first RNA-seq validation check for the clinical HER2 GigaTIME pilot.
- `docs/clinical_her2_visual_qc.md`: first visual/spatial QC pass for the clinical HER2 virtual immune-channel signal.
- `docs/clinical_her2_tile_sampling_robustness.md`: 256-tile robustness check showing whether the 64-tile HER2-zero versus HER2-low signal persists with denser sampling.
- `docs/clinical_her2_rna_program_validation.md`: broader RNA immune/tissue program validation after the 256-tile robustness run.
- `docs/clinical_her2_gigatime_data_cleanup.md`: pre-classifier tile cleanup using cellular tissue and virtual CK-enriched GigaTIME views.
- `docs/clinical_her2_classifier_baseline.md`: first diagnostic-model style classifier baseline for HER2-positive/negative, HER2-low/zero, and three-class HER2 prediction.
- `docs/clinical_her2_cleaned_classifier_comparison.md`: classifier comparison after GigaTIME tile cleanup and CK-enriched feature selection.
- `docs/advisor_brief.md`: concise project framing and discussion points.
- `docs/current_pilot_run.md`: current two-case run status and advisor-facing caveats.
- `configs/tcga_brca_her2.yaml`: default paths and pilot settings.
- `notebooks/clinical_her2_findings_simple.ipynb` and `notebooks/clinical_her2_findings_simple.html`: simple presentation-ready summary of the findings so far.

## Requirements

GigaTIME requires access to the gated Hugging Face model. Accept the terms on the model card, then set a read-only token:

```bash
export HF_TOKEN=<huggingface_read_token>
```

Create the working environment:

```bash
conda env create -f envs/gigatime-tcga.yml
conda activate gigatime-tcga
```

The workflow is research-only and not for clinical decision-making, matching the GigaTIME model license notice.

## 1. Query TCGA-BRCA and Extract ERBB2

Start with a pilot subset before running all BRCA slides:

```bash
python scripts/gdc_query_tcga_brca.py \
  --out-dir data/tcga_brca \
  --case-limit 25 \
  --download-expression
```

This writes:

- `data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv`
- `data/tcga_brca/tcga_brca_star_counts_manifest.tsv`
- `data/tcga_brca/erbb2_expression.csv`
- `data/tcga_brca/file_metadata_*.json`

To download slide files, either use the manifest with the GDC Data Transfer Tool:

```bash
gdc-client download \
  -m data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv \
  -d data/tcga_brca/slides
```

or download a very small pilot directly:

```bash
python scripts/gdc_query_tcga_brca.py \
  --out-dir data/tcga_brca \
  --case-limit 5 \
  --download-slides \
  --max-slide-downloads 5 \
  --slide-download-order smallest
```

To pull one specific case, add `--slide-case-id TCGA-3C-AALI`.

## 2. Build Clinical HER2 Labels

For analyses that compare HER2-positive, HER2-low, and HER2-zero disease, build labels from the TCGA-BRCA clinical supplement:

```bash
conda run -n gigatime-tcga python scripts/build_tcga_brca_clinical_her2_labels.py
```

This writes:

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt`

The label rules are:

- `HER2-positive`: IHC `3+`, ISH positive, or positive IHC receptor status when detailed score/ISH are missing.
- `HER2-low`: IHC `1+` with no positive ISH, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH.
- `HER2-unknown`: missing, not evaluated, equivocal without definitive ISH, or otherwise incomplete HER2 fields.

## 3. Select a Balanced Clinical HER2 Cohort

After clinical labels are available, select a balanced 10/10/10 pilot cohort:

```bash
conda run -n gigatime-tcga python scripts/select_clinical_her2_cohort.py
```

This writes:

- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_summary.json`

The default selector chooses 10 cases per clinical group, prioritizing direct clinical HER2 labels, already-downloaded slides, smaller slide files, and deterministic case IDs.

To download the selected slide files with the GDC Data Transfer Tool:

```bash
gdc-client download \
  -m data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv \
  -d data/tcga_brca/slides
```

If `gdc-client` is not installed, use the project downloader:

```bash
conda run -n gigatime-tcga python scripts/download_clinical_her2_cohort_slides.py \
  --only-missing
```

This downloads each selected slide by GDC file ID into `data/tcga_brca/slides/<case>/` and writes:

- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`

## 4. Run GigaTIME on TCGA-BRCA Slides

```bash
python scripts/run_gigatime_tcga_brca.py \
  --slides-dir data/tcga_brca/slides \
  --out-dir results/gigatime_tcga_brca \
  --tile-limit 512 \
  --tile-order random \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

Key output:

- `results/gigatime_tcga_brca/slide_scores.csv`
- `results/gigatime_tcga_brca/tile_scores.csv`
- `results/gigatime_tcga_brca/heatmaps/*.png`

For the first advisor meeting, `--tile-limit 512` is enough to demonstrate the pipeline. Increase or remove it for the full run.

To run GigaTIME only on the selected clinical HER2 cohort and skip slides that have not been downloaded yet:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/tcga_brca/clinical_her2_cohort_slides_files.csv \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2 \
  --tile-limit 64 \
  --tile-order random \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

The first robustness rerun used the same selected slides with denser sampling:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/tcga_brca/clinical_her2_cohort_slides_files.csv \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2_tile256 \
  --tile-limit 256 \
  --tile-order random \
  --random-seed 42 \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

## 5. Summarize Results

```bash
python scripts/summarize_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca/slide_scores.csv \
  --expression data/tcga_brca/erbb2_expression.csv \
  --out-dir results/gigatime_tcga_brca/advisor_summary
```

This writes joined data, channel-level HER2-high versus HER2-low summaries, figures, and `advisor_summary.md`.

For the clinical HER2-positive/HER2-low/HER2-zero cohort:

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

To run the first indirect RNA-seq validation layer:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_signatures.py
```

This compares GigaTIME channels such as `CD68`, `PD-L1`, `CD11c`, and `Ki67` with simple matched RNA marker signatures from the available STAR-count files.

To run broader RNA program validation after the 256-tile clinical HER2 rerun:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_programs.py
```

This compares virtual composite programs such as myeloid/checkpoint and T-cell/checkpoint with broader RNA programs such as cytotoxic T-cell, checkpoint/IFNG, myeloid/macrophage, B-cell, stromal, endothelial, epithelial, and proliferation signatures.

To run the first slide-level HER2 classifier baseline:

```bash
conda run -n gigatime-tcga python scripts/train_her2_classifier_baseline.py
```

This trains regularized logistic and nearest-centroid baselines with leave-one-out cross-validation. It reports HER2-positive versus negative, HER2-low versus zero, and full three-class HER2 prediction performance.

To build cleaned pre-classifier GigaTIME feature views from the 256-tile output:

```bash
conda run -n gigatime-tcga python scripts/cleanup_gigatime_tile_features.py
```

This creates cellular-tissue and virtual CK-enriched slide feature tables under `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/` and tracked figures under `docs/assets/clinical_her2_gigatime_cleanup/`.

To rerun HER2 classifiers across those cleaned feature views:

```bash
conda run -n gigatime-tcga python scripts/train_her2_cleaned_classifier_comparison.py
```

This writes cleaned-view classifier metrics under `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/` and tracked figures under `docs/assets/clinical_her2_cleaned_classifier/`.

## 6. Render All Virtual mIF Channel Images

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_channel_images.py
```

This writes documentation-facing figures to `docs/assets/virtual_mif_channels/`, including all-channel group means, a slide-by-channel activation matrix, and HER2-high/HER2-low reference grids for the 23 GigaTIME virtual mIF channels. See `docs/virtual_mif_channel_outputs.md` for interpretation.

To create fluorescence-style virtual mIF images that look closer to real multiplex immunofluorescence panels:

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_composites.py
```

This writes H&E-versus-virtual-mIF panels and marker-composite montages to `docs/assets/virtual_mif_composites/`. These are still GigaTIME predictions, not experimental mIF data.

To render the clinical HER2 visual QC panels for cases driving high virtual `CD68`, `PD-L1`, and `CD11c`:

```bash
conda run -n gigatime-tcga python scripts/render_clinical_her2_visual_qc.py
```

This writes tracked QC panels and selected-case tables to `docs/assets/clinical_her2_visual_qc/`.

To rebuild the simple presentation notebook and HTML report:

```bash
conda run -n gigatime-tcga python scripts/build_clinical_her2_findings_report.py
```

This writes:

- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## Notes for the Advisor Discussion

- HER2 is represented here by `ERBB2` RNA expression from TCGA-BRCA STAR-count files.
- GigaTIME outputs virtual mIF maps for 23 channels, including immune markers such as `CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `PD-1`, and `PD-L1`.
- The first deliverable is a replication/adaptation pilot, not a new model: run the released model on TCGA-BRCA H&E slides and ask whether virtual TIME signatures differ across clinical HER2 groups.
- The initial run should be treated as exploratory until tissue QC, slide-level aggregation, and HER2 clinical annotations are reviewed.
- The current clinical HER2 pilot has processed 30 selected slides: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero. The strongest pilot signal is higher GigaTIME-predicted CD68, PD-L1, and CD11c in HER2-zero versus HER2-low, but these are hypothesis-generating and not FDR-significant after pairwise correction.
- The 256-tile robustness rerun reproduced the same HER2-zero greater than HER2-low direction for CD68, PD-L1, and CD11c. The leading pairwise q values improved to about 0.113 but remained above 0.05.
- The first RNA-seq validation check did not strongly confirm the virtual immune-channel signal; correlations between matched RNA marker signatures and GigaTIME channels were weak and not FDR-significant.
- Broader RNA program validation also did not positively confirm the virtual immune/checkpoint signal. The strongest FDR-significant associations were negative correlations between virtual immune/checkpoint programs and endothelial RNA signal.
- The first classifier baseline suggests possible GigaTIME signal for HER2-low versus HER2-zero, but not reliable HER2-positive/negative or three-class diagnosis. This is not clinically usable.
- The pre-classifier cleanup shows that the HER2-zero greater than HER2-low CD68/PD-L1/CD11c signal persists after cellular-tissue filtering, but weakens under strict CK-enriched tile selection. This suggests the original signal may depend partly on broader tissue context, not only tumor-rich tiles.
- The cleaned-view classifier comparison preserves HER2-low versus HER2-zero balanced accuracy at 0.800 after cellular-tissue filtering, but drops to 0.650 in CK-enriched views. This supports a microenvironment/tissue-context interpretation more than a purely tumor-epithelial HER2 classifier.
- The first visual QC pass found that high virtual CD68/PD-L1/CD11c tiles were tissue-containing and cellular rather than obvious blank background, but this still does not validate the virtual marker biology.
