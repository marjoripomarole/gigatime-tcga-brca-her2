# TCGA-BRCA GigaTIME HER2 Workflow

This project wraps the official GigaTIME implementation for a TCGA-BRCA pilot focused on HER2/ERBB2 expression.

The goal is to generate virtual multiplex immunofluorescence (mIF) features from TCGA-BRCA diagnostic H&E whole-slide images, join those features to ERBB2 RNA expression, and produce advisor-ready summary tables and figures.

## What This Contains

- `external/GigaTIME/`: cloned official GigaTIME code from `prov-gigatime/GigaTIME`.
- `scripts/gdc_query_tcga_brca.py`: queries GDC for TCGA-BRCA diagnostic slides and STAR-count RNA-seq files, writes GDC manifests, and can extract ERBB2 expression.
- `scripts/build_tcga_brca_clinical_her2_labels.py`: queries the GDC TCGA-BRCA clinical supplement and builds reproducible clinical HER2-positive/HER2-low/HER2-zero labels.
- `scripts/run_gigatime_tcga_brca.py`: tiles TCGA-BRCA `.svs` slides, runs the official GigaTIME model, and aggregates virtual mIF channels per slide.
- `scripts/summarize_her2_gigatime.py`: joins GigaTIME slide scores with ERBB2 expression and makes HER2-high/HER2-low summary figures.
- `scripts/render_virtual_mif_channel_images.py`: renders all-channel virtual mIF figures from GigaTIME tile and slide predictions.
- `scripts/render_virtual_mif_composites.py`: reruns GigaTIME on selected tiles and renders fluorescence-style virtual mIF composites from the full predicted channel maps.
- `docs/virtual_mif_channel_outputs.md`: explains the generated virtual mIF channel images and how to interpret them.
- `docs/plain_language_methodology.md`: detailed non-specialist explanation of the study background, methodology, outputs, and current limitations.
- `docs/paper_proposal_process_log.md`: living process log for turning the pilot into a paper or grant proposal.
- `docs/advisor_brief.md`: concise project framing and discussion points.
- `docs/current_pilot_run.md`: current two-case run status and advisor-facing caveats.
- `configs/tcga_brca_her2.yaml`: default paths and pilot settings.

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

## 3. Run GigaTIME on TCGA-BRCA Slides

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

## 4. Summarize by HER2/ERBB2 Expression

```bash
python scripts/summarize_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca/slide_scores.csv \
  --expression data/tcga_brca/erbb2_expression.csv \
  --out-dir results/gigatime_tcga_brca/advisor_summary
```

This writes joined data, channel-level HER2-high versus HER2-low summaries, figures, and `advisor_summary.md`.

## 5. Render All Virtual mIF Channel Images

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_channel_images.py
```

This writes documentation-facing figures to `docs/assets/virtual_mif_channels/`, including all-channel group means, a slide-by-channel activation matrix, and HER2-high/HER2-low reference grids for the 23 GigaTIME virtual mIF channels. See `docs/virtual_mif_channel_outputs.md` for interpretation.

To create fluorescence-style virtual mIF images that look closer to real multiplex immunofluorescence panels:

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_composites.py
```

This writes H&E-versus-virtual-mIF panels and marker-composite montages to `docs/assets/virtual_mif_composites/`. These are still GigaTIME predictions, not experimental mIF data.

## Notes for the Advisor Discussion

- HER2 is represented here by `ERBB2` RNA expression from TCGA-BRCA STAR-count files.
- GigaTIME outputs virtual mIF maps for 23 channels, including immune markers such as `CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `PD-1`, and `PD-L1`.
- The first deliverable is a replication/adaptation pilot, not a new model: run the released model on TCGA-BRCA H&E slides and ask whether virtual TIME signatures differ across ERBB2 expression strata.
- The initial run should be treated as exploratory until tissue QC, slide-level aggregation, and HER2 clinical annotations are reviewed.
