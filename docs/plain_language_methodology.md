# Plain-Language Methodology: TCGA-BRCA GigaTIME HER2 Pilot

This document explains what has been done in this project so far for a reader who does not have a genetics, pathology, or computational biology background.

The short version is: this project takes public breast cancer microscope images, runs an existing artificial intelligence model called GigaTIME on small pieces of those images, and compares the model's immune-marker predictions with the activity level of a breast cancer gene called `ERBB2`, which is the gene connected to HER2 biology.

This is an early research pilot. It is not a clinical test, not a diagnostic tool, and not a new AI model.

## 1. The Research Question

The current working question is:

Can the released GigaTIME model generate interpretable immune-environment features from public TCGA breast cancer H&E pathology slides, and do those features look different between tumors with high versus low `ERBB2` expression?

In plainer terms:

- We have breast cancer tissue images.
- We have a separate measurement of how active the `ERBB2` gene is in those tumors.
- We use an existing AI model to estimate immune and tumor marker patterns from the tissue images.
- We ask whether the image-derived marker patterns differ between tumors with high `ERBB2` activity and tumors with low `ERBB2` activity.

This is a replication/adaptation pilot. The project does not train a new model. It applies a previously released model to a new breast cancer use case.

## 2. Biological Background Without Assuming Genetics Knowledge

### What is DNA, a gene, and gene expression?

DNA is the instruction library inside cells. A gene is one instruction unit inside that library.

Cells do not use every gene equally. Some genes are very active in a cell, while others are quiet. When a gene is active, the cell makes RNA copies from that gene. Measuring RNA is one common way to estimate how active a gene is. This is called gene expression.

In this project, the relevant gene is `ERBB2`.

### What is HER2, and how is it related to ERBB2?

HER2 is a protein found on the surface of some breast cancer cells. In some tumors, there is too much HER2 activity. This matters clinically because HER2 can influence tumor behavior and treatment choices.

`ERBB2` is the gene that contains the instructions for making the HER2 protein.

Important distinction:

- `ERBB2` is the gene.
- HER2 is the protein and clinical biology associated with that gene.
- This project currently uses `ERBB2` RNA expression as a proxy for HER2-related biology.
- This project does not yet use official clinical HER2 labels from immunohistochemistry or FISH testing.

So when this project says "HER2-high" and "HER2-low", it currently means high or low `ERBB2` RNA expression in the available TCGA data. It does not mean a confirmed clinical HER2-positive or HER2-negative diagnosis.

### What is RNA-seq and TPM?

RNA-seq is a laboratory method that measures RNA molecules in a tissue sample. It gives an estimate of which genes are active and how active they are.

TPM means "transcripts per million." It is a normalized number used to compare gene expression levels. A higher TPM for `ERBB2` means the tumor sample had more RNA from the `ERBB2` gene.

In this project, the file `data/tcga_brca/erbb2_expression.csv` stores the extracted `ERBB2` TPM values.

## 3. Pathology Background Without Assuming Medical Training

### What is an H&E slide?

An H&E slide is a standard pathology slide. The tissue is stained with two dyes:

- Hematoxylin stains cell nuclei, usually blue or purple.
- Eosin stains other tissue structures, usually pink.

Pathologists use H&E slides every day to look at tissue structure under a microscope.

In this project, the H&E slides are digital whole-slide images. A whole-slide image is a very large scanned image of the entire tissue section.

### Why split a slide into tiles?

Whole-slide images are huge. An AI model usually cannot process the entire image at once.

So the workflow cuts each slide into many small square image patches called tiles. In this project, the GigaTIME inference script uses 256 by 256 pixel tiles by default.

Not every tile is useful. Some areas may be blank background or have little tissue. The script estimates the tissue fraction of each tile and keeps tiles that pass a minimum tissue threshold.

For the current ERBB2-extreme pilot run, the workflow processed 64 random tissue tiles per slide.

## 4. What is GigaTIME?

GigaTIME is an existing released model. This project uses the official GigaTIME implementation and model weights.

The key idea is that GigaTIME takes ordinary H&E pathology image tiles as input and predicts virtual multiplex immunofluorescence-like marker maps as output.

That phrase has several parts:

- Multiplex immunofluorescence, often shortened to mIF, is a lab technique that can stain tissue for many biological markers at the same time.
- A marker is a biological signal associated with a cell type or process. For example, `CD3` and `CD8` are often used as immune-cell markers.
- "Virtual" means the model predicts marker-like signal from the H&E image computationally. The tissue was not actually stained for those markers in this project.

GigaTIME predicts 23 channels:

`DAPI`, `TRITC`, `Cy5`, `PD-1`, `CD14`, `CD4`, `T-bet`, `CD34`, `CD68`, `CD16`, `CD11c`, `CD138`, `CD20`, `CD3`, `CD8`, `PD-L1`, `CK`, `Ki67`, `Tryptase`, `Actin-D`, `Caspase3-D`, `PHH3-B`, and `Transgelin`.

For the summary analysis, the current project focuses on a smaller set of interpretable channels:

`CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `CD11c`, `PD-1`, `PD-L1`, `CK`, and `Ki67`.

Very roughly:

- `CD3`, `CD4`, and `CD8` relate to T cells, a type of immune cell.
- `CD20` relates to B cells, another immune cell type.
- `CD68` relates to macrophage-like immune cells.
- `PD-1` and `PD-L1` relate to immune checkpoint biology.
- `CK` relates to epithelial/tumor-cell structure.
- `Ki67` relates to cell proliferation.

These are model-predicted research features. They should not be interpreted as direct laboratory measurements without validation.

## 5. Data Source

The data source is TCGA-BRCA.

TCGA means The Cancer Genome Atlas. It is a large public cancer research dataset.

BRCA is TCGA's breast cancer project. In this context, BRCA means breast invasive carcinoma; it is not the same thing as the `BRCA1` or `BRCA2` genes.

This project uses two kinds of TCGA-BRCA data from the Genomic Data Commons, or GDC:

- Diagnostic H&E whole-slide images, stored as `.svs` slide files.
- RNA-seq gene expression files, specifically STAR-count files, used to extract `ERBB2` expression.

The main query/download script is:

```bash
scripts/gdc_query_tcga_brca.py
```

## 6. What Has Been Done So Far

### Step 1: Query TCGA-BRCA files from GDC

The workflow queried GDC for two open-access file types:

- TCGA-BRCA slide images where `data_type` is `Slide Image` and `data_format` is `SVS`.
- TCGA-BRCA RNA-seq expression files where the workflow type is `STAR - Counts`.

The script writes metadata and manifests under:

```text
data/tcga_brca/
```

Important files include:

```text
data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv
data/tcga_brca/tcga_brca_diagnostic_slides_files.csv
data/tcga_brca/tcga_brca_star_counts_manifest.tsv
data/tcga_brca/tcga_brca_star_counts_files.csv
data/tcga_brca/file_metadata_slides.json
data/tcga_brca/file_metadata_star_counts.json
```

The manifests are download instructions. The CSV and JSON files are metadata tables describing which files were found.

### Step 2: Download RNA-seq files and extract ERBB2 expression

The same GDC script can download STAR-count RNA-seq files and extract the row corresponding to `ERBB2`.

The gene is identified by:

```text
Gene symbol: ERBB2
Ensembl gene ID: ENSG00000141736
```

The extracted expression table is:

```text
data/tcga_brca/erbb2_expression.csv
```

Current status:

- `ERBB2` expression was extracted for 80 TCGA-BRCA cases.
- In those 80 cases, the observed `ERBB2` TPM range is approximately 6.7 to 3236.8.

### Step 3: Select ERBB2-high and ERBB2-low extreme cases

Instead of starting with all available cases, the project selected the extremes:

- The 20 cases with the highest `ERBB2` TPM values.
- The 20 cases with the lowest `ERBB2` TPM values.

This creates a clearer first comparison than using all cases immediately.

The script is:

```bash
scripts/select_her2_extremes.py
```

The output is:

```text
data/tcga_brca/her2_extreme_cases.csv
```

Current selected groups:

- 20 `HER2-high` cases by `ERBB2` expression, with `ERBB2` TPM from about 219.2 to 3236.8.
- 20 `HER2-low` cases by `ERBB2` expression, with `ERBB2` TPM from about 6.7 to 78.2.

Again, these labels are expression-based research labels, not clinical HER2 labels.

### Step 4: Download slide images for selected cases

The workflow then tried to download diagnostic H&E whole-slide images for the selected cases.

Slide downloads are large and can be slow or unstable. The current workspace notes that GDC slide downloads repeatedly dropped connections, so only a subset of the selected cases has been processed so far.

Current pilot status:

- Selected target cohort: 40 ERBB2-extreme cases, 20 high and 20 low.
- Successfully processed so far: 12 slides.
- Processed groups so far: 7 HER2-high slides and 5 HER2-low slides.

### Step 5: Run GigaTIME on the H&E slides

The GigaTIME inference script is:

```bash
scripts/run_gigatime_tcga_brca.py
```

For each slide, the script:

1. Opens the digital pathology slide.
2. Divides the slide into 256 by 256 pixel tiles.
3. Estimates whether each tile contains enough tissue.
4. Keeps tissue-containing tiles.
5. Randomly samples tiles when a tile limit is used.
6. Normalizes the tile image in the same general style used for ImageNet-based models.
7. Runs the tile through the GigaTIME model.
8. Gets predicted marker maps for the 23 GigaTIME channels.
9. Summarizes each tile into marker scores.
10. Aggregates tile scores into one row per slide.

For the current ERBB2-extreme pilot:

- Slides processed: 12.
- Tiles per slide: 64 random tissue tiles.
- Total tile predictions: 768.
- Device used: CPU.

The main outputs are:

```text
results/gigatime_tcga_brca_extremes/slide_scores.csv
results/gigatime_tcga_brca_extremes/tile_scores.csv
results/gigatime_tcga_brca_extremes/heatmaps/
```

### Step 6: Summarize the GigaTIME predictions by ERBB2 group

The summary script is:

```bash
scripts/summarize_her2_gigatime.py
```

This script combines:

- The slide-level GigaTIME output.
- The `ERBB2` RNA expression table.
- The optional explicit HER2-high/HER2-low grouping file.

It writes a joined table:

```text
results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv
```

It also writes a channel comparison table:

```text
results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv
```

For each marker channel, the summary compares the average virtual activation in the HER2-high group versus the HER2-low group.

The summary includes:

- Mean marker activation in HER2-high slides.
- Mean marker activation in HER2-low slides.
- Difference between the two means.
- Spearman correlation between `ERBB2` TPM and marker activation.
- Exploratory statistical tests.
- Effect-size estimates.

The current strongest mean differences in the 12-slide pilot are small and exploratory. The largest absolute mean difference among the summarized channels is currently for `CK`, followed by `PD-1`, `CD3`, `CD4`, and `CD20`. These are not definitive findings because the current sample is small.

The generated figures include:

```text
results/gigatime_tcga_brca_extremes/advisor_summary/erbb2_tpm_distribution.png
results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_deltas.png
results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_boxplots.png
results/gigatime_tcga_brca_extremes/advisor_summary/erbb2_vs_virtual_mif_scatter.png
docs/assets/virtual_mif_channels/virtual_mif_all_channel_group_means.png
docs/assets/virtual_mif_channels/virtual_mif_slide_channel_matrix.png
docs/assets/virtual_mif_channels/her2_high_reference_all_virtual_mif_channels.png
docs/assets/virtual_mif_channels/her2_low_reference_all_virtual_mif_channels.png
docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_virtual_mif_montage.png
docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_virtual_mif_montage.png
```

The files under `docs/assets/virtual_mif_channels/` are documentation-facing images of the GigaTIME virtual mIF outputs. They show all 23 predicted channels, including representative spatial tile maps for one ERBB2-high slide and one ERBB2-low slide. See `docs/virtual_mif_channel_outputs.md` for a figure-by-figure explanation.

The files under `docs/assets/virtual_mif_composites/` look closer to real multiplex immunofluorescence images. They are made by rerunning GigaTIME on selected H&E tiles, keeping the full predicted channel maps, and compositing marker colors on a black background. They are still virtual predictions, not real mIF measurements.

### Step 7: Render visual examples from the H&E slides

The image-rendering script is:

```bash
scripts/render_he_slide_images.py
```

This script makes visual material for discussion and presentation:

- A low-resolution thumbnail of the whole H&E slide.
- An overlay showing where sampled tiles came from.
- Colored overlays showing virtual marker activation for selected markers.
- Montages of the top-scoring H&E tiles for selected virtual markers.

The current output directory is:

```text
results/gigatime_tcga_brca_extremes/he_examples/
```

These images are useful because they connect the abstract model outputs back to tissue regions that a human can inspect.

## 7. What the Output Tables Mean

### `slide_scores.csv`

This table has one row per processed slide.

Each row contains:

- The slide file path.
- The slide ID.
- The TCGA case ID.
- The number of tiles analyzed.
- The average tissue fraction.
- For each GigaTIME channel, the average predicted activation across tiles.
- For each GigaTIME channel, the fraction of tile area above an activation threshold.

This is the main slide-level model output.

### `tile_scores.csv`

This table has one row per tile.

Each row contains:

- The tile coordinate on the slide.
- The estimated tissue fraction.
- The slide ID.
- The TCGA case ID.
- The GigaTIME marker scores for that tile.

This table is more detailed than `slide_scores.csv` and is useful for heatmaps or spatial inspection.

### `joined_slide_her2_gigatime.csv`

This table joins model output to gene expression.

Each row represents a processed slide with:

- GigaTIME slide-level marker scores.
- The matching TCGA case ID.
- `ERBB2` TPM expression.
- The HER2-high or HER2-low expression group.

This is the table used for the HER2/ERBB2 comparison.

### `her2_group_channel_summary.csv`

This table summarizes group differences by marker channel.

For each marker, it asks:

- Is the average virtual marker activation higher in the ERBB2-high group or the ERBB2-low group?
- How large is the difference?
- Is there a simple monotonic relationship between ERBB2 expression and the marker score?
- Are the exploratory statistics suggestive enough to prioritize follow-up?

This table should be treated as a guide for discussion, not as final evidence.

## 8. What This Study Has Not Done Yet

The current pilot has not yet:

- Processed all 40 selected ERBB2-extreme cases.
- Performed full tissue quality control.
- Confirmed clinical HER2 status using IHC or FISH annotations.
- Compared `ERBB2` expression groups against official pathology HER2 labels.
- Validated GigaTIME predictions against real multiplex immunofluorescence staining in these TCGA slides.
- Trained or fine-tuned a new model.
- Produced a clinically deployable classifier.

These limitations are important. The current goal is proof of workflow and early biological exploration, not a final scientific claim.

## 9. Why This Is Still Useful

This pilot is useful because it proves several practical pieces:

- TCGA-BRCA H&E slides can be queried and downloaded through GDC.
- Matching RNA-seq expression files can be queried and used to extract `ERBB2`.
- ERBB2-high and ERBB2-low expression groups can be selected reproducibly.
- The released GigaTIME model can be run on TCGA-BRCA pathology tiles.
- The GigaTIME outputs can be aggregated into slide-level marker features.
- Those marker features can be joined back to `ERBB2` expression.
- The results can be summarized in tables, plots, and visual examples.

In other words, the technical pipeline works end to end on a small pilot subset.

## 10. Current Scientific Interpretation

The safest interpretation is:

This is an exploratory feasibility run showing that an existing H&E-to-virtual-mIF model can be applied to TCGA-BRCA breast cancer slides and connected to `ERBB2` expression data.

The current 12-slide result is too small for strong biological conclusions. It is useful for checking whether the workflow is plausible and for deciding whether to scale up to the full selected cohort.

The next scientific step is to complete the 20 high / 20 low ERBB2-extreme cohort, rerun the same pipeline with the same settings, add tissue quality review, and compare the expression-based HER2 grouping with clinical HER2 annotations if available.

## 11. Reproducible Workflow Summary

The workflow is organized around these scripts:

```text
scripts/gdc_query_tcga_brca.py
scripts/select_her2_extremes.py
scripts/run_gigatime_tcga_brca.py
scripts/summarize_her2_gigatime.py
scripts/render_he_slide_images.py
```

The typical flow is:

```text
Query GDC files
  -> download RNA-seq files
  -> extract ERBB2 expression
  -> select ERBB2-high and ERBB2-low cases
  -> download H&E slides
  -> run GigaTIME on slide tiles
  -> aggregate tile predictions to slide scores
  -> join slide scores to ERBB2 expression
  -> summarize and visualize group differences
  -> render H&E examples for human inspection
```

## 12. Practical Notes for a New Reader

If you are reading the project for the first time, start with:

1. `README.md` for commands and file locations.
2. `docs/current_pilot_run.md` for the current run status.
3. This document for the conceptual explanation.
4. `docs/advisor_brief.md` for the short advisor-facing summary.
5. The executed notebook in `notebooks/` for a presentation-style view of the current pilot.

The most important caution is that "HER2-high" and "HER2-low" currently mean high and low `ERBB2` RNA expression in this pilot dataset. They do not yet mean clinically validated HER2 status.
