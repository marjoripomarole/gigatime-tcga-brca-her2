# Plain-Language Methodology: TCGA-BRCA GigaTIME HER2 Pilot

This document explains what has been done in this project so far for a reader who does not have a genetics, pathology, or computational biology background.

The short version is: this project takes public breast cancer microscope images, runs an existing artificial intelligence model called GigaTIME on small pieces of those images, and compares the model's predicted immune-marker patterns across clinically defined HER2 groups.

This is an early research pilot. It is not a clinical test, not a diagnostic tool, and not a new AI model.

## 1. The Research Question

The current working question is:

Can the released GigaTIME model generate interpretable immune-environment features from public TCGA breast cancer H&E pathology slides, and do those predicted features look different between HER2-positive, HER2-low, and HER2-zero tumors?

In plainer terms:

- We have breast cancer tissue images.
- We have clinical HER2 information from TCGA clinical records.
- We use an existing AI model to estimate immune and tumor marker patterns from the tissue images.
- We ask whether those image-derived marker patterns differ across HER2-positive, HER2-low, and HER2-zero breast cancers.

This is a replication/adaptation pilot. The project does not train a new model. It applies a previously released model to a breast cancer research question.

## 2. Biological Background Without Assuming Genetics Knowledge

### What is HER2?

HER2 is a protein that can be present on the surface of breast cancer cells. In some tumors, HER2 is strongly present or amplified. That can affect tumor biology and treatment options.

HER2 status is usually assessed in clinical pathology using protein or gene-copy tests:

- IHC, or immunohistochemistry, estimates how much HER2 protein is visible in the tumor tissue.
- ISH/FISH, or in situ hybridization, checks whether the HER2 gene region is amplified.

### What are HER2-positive, HER2-low, and HER2-zero?

For this project, the groups are defined from TCGA clinical HER2 fields:

- `HER2-positive`: IHC `3+`, ISH positive, or a positive HER2 receptor status when more detailed fields are missing.
- `HER2-low`: IHC `1+`, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH evidence.
- `HER2-unknown`: missing, not evaluated, equivocal, contradictory, or incomplete fields.

This distinction matters because HER2-low and HER2-zero can be biologically and clinically different, even though both are not classic HER2-positive disease.

### What is ERBB2?

`ERBB2` is the gene that contains the instructions for making the HER2 protein.

At the start of this project, `ERBB2` RNA expression was used as a first HER2-biology proxy. That was useful for building the pipeline, but RNA expression is not the same as clinical HER2 status. The current better analysis uses clinical IHC/ISH-derived HER2 groups.

## 3. Pathology Background Without Assuming Medical Training

### What is an H&E slide?

An H&E slide is a standard pathology slide. The tissue is stained with two dyes:

- Hematoxylin stains cell nuclei, usually blue or purple.
- Eosin stains other tissue structures, usually pink.

Pathologists use H&E slides every day to look at tissue structure under a microscope.

In this project, the H&E slides are digital whole-slide images. A whole-slide image is a very large scanned image of the entire tissue section.

### Why split a slide into tiles?

Whole-slide images are huge. An AI model usually cannot process the entire image at once.

So the workflow cuts each slide into many small image patches called tiles. The GigaTIME inference script uses 256 by 256 pixel tiles by default.

Not every tile is useful. Some areas are blank background or have little tissue. The script estimates the tissue fraction of each tile and keeps tiles that pass a minimum tissue threshold.

For the current clinical HER2 pilot, the workflow processed 64 random tissue tiles per slide.

## 4. What is GigaTIME?

GigaTIME is an existing released model. This project uses the official GigaTIME implementation and model weights.

The key idea is that GigaTIME takes ordinary H&E pathology image tiles as input and predicts virtual multiplex immunofluorescence-like marker maps as output.

That phrase has several parts:

- Multiplex immunofluorescence, often shortened to mIF, is a lab technique that can stain tissue for many biological markers at the same time.
- A marker is a biological signal associated with a cell type or process. For example, `CD3` and `CD8` are often used as immune-cell markers.
- "Virtual" means the model predicts marker-like signal from the H&E image computationally. The tissue was not actually stained for those markers in this project.

GigaTIME predicts 23 channels:

`DAPI`, `TRITC`, `Cy5`, `PD-1`, `CD14`, `CD4`, `T-bet`, `CD34`, `CD68`, `CD16`, `CD11c`, `CD138`, `CD20`, `CD3`, `CD8`, `PD-L1`, `CK`, `Ki67`, `Tryptase`, `Actin-D`, `Caspase3-D`, `PHH3-B`, and `Transgelin`.

For the clinical HER2 summary analysis, the project focuses on a smaller set of interpretable channels:

`CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `CD11c`, `PD-1`, `PD-L1`, `CK`, and `Ki67`.

Very roughly:

- `CD3`, `CD4`, and `CD8` relate to T cells, a type of immune cell.
- `CD20` relates to B cells.
- `CD68` and `CD11c` relate to macrophage/myeloid immune biology.
- `PD-1` and `PD-L1` relate to immune checkpoint biology.
- `CK` relates to epithelial/tumor-cell structure.
- `Ki67` relates to cell proliferation.

These are model-predicted research features. They should not be interpreted as direct laboratory measurements without validation.

## 5. Data Source

The data source is TCGA-BRCA.

TCGA means The Cancer Genome Atlas. It is a large public cancer research dataset.

BRCA is TCGA's breast cancer project. In this context, BRCA means breast invasive carcinoma; it is not the same thing as the `BRCA1` or `BRCA2` genes.

This project uses three kinds of TCGA-BRCA data from the Genomic Data Commons, or GDC:

- Diagnostic H&E whole-slide images, stored as `.svs` slide files.
- RNA-seq gene expression files, used to extract `ERBB2` expression.
- Clinical supplement fields, used to assign HER2-positive, HER2-low, HER2-zero, or HER2-unknown labels.

## 6. What Has Been Done So Far

### Step 1: Query TCGA-BRCA files from GDC

The workflow queried GDC for:

- TCGA-BRCA slide images where `data_type` is `Slide Image` and `data_format` is `SVS`.
- TCGA-BRCA RNA-seq expression files where the workflow type is `STAR - Counts`.
- TCGA-BRCA clinical supplement files containing HER2-related fields.

Important local metadata files include:

```text
data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv
data/tcga_brca/tcga_brca_diagnostic_slides_files.csv
data/tcga_brca/tcga_brca_star_counts_manifest.tsv
data/tcga_brca/tcga_brca_star_counts_files.csv
data/tcga_brca/file_metadata_slides.json
data/tcga_brca/file_metadata_star_counts.json
```

### Step 2: Extract ERBB2 expression

The GDC script downloaded selected STAR-count RNA-seq files and extracted the row corresponding to `ERBB2`.

The extracted expression table is:

```text
data/tcga_brca/erbb2_expression.csv
```

This was used first to make an ERBB2-high versus ERBB2-low pilot. That pilot proved the workflow could run, but it is not the main clinical HER2 comparison.

### Step 3: Build clinical HER2 labels

The clinical HER2 labeling script is:

```bash
scripts/build_tcga_brca_clinical_her2_labels.py
```

It downloads the TCGA-BRCA patient-level clinical supplement and extracts HER2 IHC/ISH fields.

Main outputs:

```text
data/tcga_brca/clinical_her2_labels.csv
data/tcga_brca/clinical_her2_labels_metadata.json
data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt
```

The resulting label table found:

| Clinical HER2 group | TCGA-BRCA clinical rows |
|---|---:|
| HER2-positive | 174 |
| HER2-low | 407 |
| HER2-zero | 61 |
| HER2-unknown | 455 |

### Step 4: Select a balanced clinical HER2 pilot cohort

The cohort selection script is:

```bash
scripts/select_clinical_her2_cohort.py
```

It joins clinical HER2 labels, ERBB2 expression, and slide metadata, then selects a balanced pilot:

| Clinical HER2 group | Selected cases |
|---|---:|
| HER2-positive | 10 |
| HER2-low | 10 |
| HER2-zero | 10 |

The selection prefers direct clinical labels, local slide availability, smaller slide files, and deterministic case IDs.

Main outputs:

```text
data/tcga_brca/clinical_her2_cohort_cases.csv
data/tcga_brca/clinical_her2_cohort_slides_files.csv
data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv
data/tcga_brca/clinical_her2_cohort_summary.json
docs/clinical_her2_cohort_selection.md
```

### Step 5: Download the selected slides

The selected clinical HER2 cohort needed 30 diagnostic H&E slides. Eight were already present locally, and 22 missing slides were downloaded with:

```bash
scripts/download_clinical_her2_cohort_slides.py
```

The downloader uses GDC file IDs and writes a status file:

```text
data/tcga_brca/clinical_her2_cohort_slide_download_status.json
```

Current selected-slide status:

| Clinical HER2 group | Selected slides | Downloaded slides |
|---|---:|---:|
| HER2-positive | 10 | 10 |
| HER2-low | 10 | 10 |
| HER2-zero | 10 | 10 |

### Step 6: Run GigaTIME on the clinical HER2 cohort

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
6. Normalizes the tile image.
7. Runs the tile through the GigaTIME model.
8. Gets predicted marker maps for the GigaTIME channels.
9. Summarizes each tile into marker scores.
10. Aggregates tile scores into one row per slide.

Current clinical HER2 run:

- Slides processed: 30.
- Cases processed: 30.
- Groups processed: 10 HER2-positive, 10 HER2-low, 10 HER2-zero.
- Tiles per slide: 64 random tissue tiles.
- Total tile predictions: about 1,920.
- Device used: Apple MPS in the current local run.

Main outputs:

```text
results/gigatime_tcga_brca_clinical_her2/slide_scores.csv
results/gigatime_tcga_brca_clinical_her2/tile_scores.csv
results/gigatime_tcga_brca_clinical_her2/heatmaps/
```

### Step 7: Summarize the GigaTIME predictions by clinical HER2 group

The clinical summary script is:

```bash
scripts/summarize_clinical_her2_gigatime.py
```

It combines:

- The slide-level GigaTIME output.
- The selected clinical HER2 cohort table.
- The clinical HER2 group labels.

It writes:

```text
results/gigatime_tcga_brca_clinical_her2/clinical_summary/joined_slide_clinical_her2_gigatime.csv
results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv
results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv
results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md
results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_boxplots.png
results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_group_mean_heatmap.png
results/gigatime_tcga_brca_clinical_her2/clinical_summary/erbb2_tpm_by_clinical_her2_group.png
```

## 7. What the Current Full-Pilot Results Show

The full clinical HER2 pilot includes 30 joined slides:

- 10 HER2-positive.
- 10 HER2-low.
- 10 HER2-zero.

The strongest three-group differences were:

| Channel | Kruskal p | Highest mean group | Lowest mean group |
|---|---:|---|---|
| CD68 | 0.0242 | HER2-zero | HER2-low |
| PD-L1 | 0.0423 | HER2-zero | HER2-low |
| CD11c | 0.0494 | HER2-zero | HER2-low |
| CD4 | 0.0794 | HER2-zero | HER2-low |
| Ki67 | 0.0920 | HER2-zero | HER2-low |

The pattern is that HER2-zero had higher predicted mean signal than HER2-low for several immune and checkpoint-related virtual channels. HER2-positive was usually between those two groups rather than clearly separated from HER2-low.

The strongest pairwise comparisons were HER2-low versus HER2-zero:

| Channel | Direction | Mann-Whitney p | BH q |
|---|---|---:|---:|
| CD68 | HER2-zero higher than HER2-low | 0.0091 | 0.2113 |
| CD11c | HER2-zero higher than HER2-low | 0.0173 | 0.2113 |
| PD-L1 | HER2-zero higher than HER2-low | 0.0211 | 0.2113 |
| CD4 | HER2-zero higher than HER2-low | 0.0312 | 0.2258 |
| Ki67 | HER2-zero higher than HER2-low | 0.0376 | 0.2258 |

No pairwise comparison remained statistically significant after multiple-testing correction. This means the result is a signal worth following, not proof of a biological conclusion.

Plain-language interpretation:

> In this first 30-slide pilot, GigaTIME predicted more immune/checkpoint-like signal in HER2-zero tumors than in HER2-low tumors, especially for CD68, PD-L1, and CD11c. The sample is small, so this should be treated as a hypothesis for validation.

## 8. What the Output Tables Mean

### `slide_scores.csv`

This table has one row per processed slide.

Each row contains:

- The slide file path.
- The slide ID.
- The TCGA case ID.
- The number of tiles analyzed.
- The average tissue fraction.
- For each GigaTIME channel, the average predicted activation across tiles.
- For each GigaTIME channel, thresholded summaries.

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

### `joined_slide_clinical_her2_gigatime.csv`

This table joins model output to clinical HER2 labels.

Each row represents a processed slide with:

- GigaTIME slide-level marker scores.
- The matching TCGA case ID.
- The clinical HER2 group.
- ER, PR, HER2 IHC/ISH, and ERBB2 expression context where available.

This is the table used for the clinical HER2 comparison.

### `clinical_her2_channel_summary.csv`

This table summarizes three-group differences by marker channel.

For each marker, it asks:

- Which HER2 group has the highest average virtual marker activation?
- Which HER2 group has the lowest average virtual marker activation?
- Is the three-group difference suggestive by Kruskal-Wallis testing?
- How large is the difference between the highest and lowest group means?

### `clinical_her2_pairwise_tests.csv`

This table compares pairs of HER2 groups:

- HER2-positive versus HER2-low.
- HER2-positive versus HER2-zero.
- HER2-low versus HER2-zero.

It includes Mann-Whitney p values and Benjamini-Hochberg corrected q values.

## 9. Visual Outputs

The project also generated documentation-facing virtual mIF images:

```text
docs/assets/virtual_mif_channels/
docs/assets/virtual_mif_composites/
```

The files under `docs/assets/virtual_mif_channels/` show all 23 predicted channels, including group-level channel means and slide-by-channel matrices.

The files under `docs/assets/virtual_mif_composites/` look closer to real multiplex immunofluorescence images. They are made by rerunning GigaTIME on selected H&E tiles, keeping the full predicted channel maps, and compositing marker colors on a black background.

These images are still virtual predictions, not real mIF measurements.

## 10. What This Study Has Not Done Yet

The current pilot has not yet:

- Processed a large clinical HER2 cohort beyond the 30 selected cases.
- Used more exhaustive whole-slide sampling.
- Performed formal tissue quality control on every sampled tile.
- Validated GigaTIME predictions against real multiplex immunofluorescence staining in these TCGA slides.
- Compared GigaTIME immune channels with RNA-seq immune signatures.
- Trained or fine-tuned a new model.
- Produced a clinically deployable classifier.

These limitations are important. The current goal is proof of workflow and early biological exploration, not a final scientific claim.

## 11. Why This Is Still Useful

This pilot is useful because it proves several practical pieces:

- TCGA-BRCA H&E slides can be queried and downloaded through GDC.
- TCGA-BRCA clinical HER2 IHC/ISH fields can be used to create HER2-positive, HER2-low, and HER2-zero groups.
- A balanced 10/10/10 clinical HER2 pilot cohort can be selected reproducibly.
- The released GigaTIME model can be run on TCGA-BRCA pathology tiles.
- The GigaTIME outputs can be aggregated into slide-level marker features.
- Those marker features can be compared across clinical HER2 groups.
- The results can be summarized in tables, plots, and visual examples.

In other words, the technical pipeline now works end to end for the clinically meaningful HER2 grouping.

## 12. Current Scientific Interpretation

The safest interpretation is:

This is an exploratory feasibility run showing that an existing H&E-to-virtual-mIF model can be applied to TCGA-BRCA breast cancer slides and connected to clinical HER2 labels.

The current 30-slide result is too small for strong biological conclusions. It is useful because it identifies a specific hypothesis: HER2-zero tumors may show higher GigaTIME-predicted immune/checkpoint signal than HER2-low tumors in this selected TCGA-BRCA pilot.

The next scientific step is to test whether this pattern remains when:

- More tiles per slide are sampled.
- More cases are included if reliable HER2-zero cases are available.
- The virtual immune channels are compared with RNA-seq immune marker expression or immune signatures.
- A human reviews representative H&E tiles and virtual mIF composites for plausibility.

## 13. Reproducible Workflow Summary

The workflow is organized around these scripts:

```text
scripts/gdc_query_tcga_brca.py
scripts/build_tcga_brca_clinical_her2_labels.py
scripts/select_clinical_her2_cohort.py
scripts/download_clinical_her2_cohort_slides.py
scripts/run_gigatime_tcga_brca.py
scripts/summarize_clinical_her2_gigatime.py
scripts/render_he_slide_images.py
scripts/render_virtual_mif_channel_images.py
scripts/render_virtual_mif_composites.py
```

The current clinical HER2 flow is:

```text
Query GDC files
  -> extract ERBB2 expression
  -> build clinical HER2 labels from IHC/ISH fields
  -> select balanced HER2-positive / HER2-low / HER2-zero cases
  -> download selected H&E slides
  -> run GigaTIME on slide tiles
  -> aggregate tile predictions to slide scores
  -> join slide scores to clinical HER2 groups
  -> summarize and visualize group differences
  -> use visual QC and RNA validation as next checks
```

## 14. Practical Notes for a New Reader

If you are reading the project for the first time, start with:

1. `README.md` for commands and file locations.
2. `docs/clinical_her2_cohort_selection.md` for the selected 30-case cohort.
3. `docs/clinical_her2_gigatime_run.md` for the current full clinical HER2 result.
4. This document for the conceptual explanation.
5. `docs/advisor_brief.md` for the short advisor-facing summary.

The most important caution is that GigaTIME outputs are predicted virtual mIF research features. They are not real multiplex immunofluorescence measurements and should be validated before making biological or clinical claims.
