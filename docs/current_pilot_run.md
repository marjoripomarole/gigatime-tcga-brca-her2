# Current Pilot Run

## Run Status

The current workspace contains a completed TCGA-BRCA GigaTIME pilot on a balanced clinical HER2 cohort:

- HER2-positive cases processed: 10
- HER2-low cases processed: 10
- HER2-zero cases processed: 10
- Slides processed: 30
- GigaTIME tiles per slide: 64 random tissue tiles
- Total tile predictions: about 1,920
- Device used in the current run: Apple MPS

The HER2 group labels come from `data/tcga_brca/clinical_her2_labels.csv`, which was built from TCGA-BRCA clinical HER2 IHC/ISH fields. The selected 30-case cohort is recorded in `data/tcga_brca/clinical_her2_cohort_cases.csv` and summarized in `docs/clinical_her2_cohort_selection.md`.

## Main Outputs

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`
- `results/gigatime_tcga_brca_clinical_her2/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/heatmaps/`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/*.png`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/rna_validation_summary.md`
- `docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv`
- `docs/assets/clinical_her2_visual_qc/*_he_vs_virtual_mif_qc.png`

The earlier ERBB2-high versus ERBB2-low pilot outputs are still present under `results/gigatime_tcga_brca_extremes/`, and the documentation-facing virtual mIF images are still under:

- `docs/assets/virtual_mif_channels/`
- `docs/assets/virtual_mif_composites/`

## Current Finding

The full 30-slide clinical HER2 pilot suggests a possible HER2-zero versus HER2-low immune-microenvironment signal.

The strongest three-group virtual-channel differences were:

| Channel | Kruskal p | Highest mean group | Lowest mean group |
|---|---:|---|---|
| CD68 | 0.0242 | HER2-zero | HER2-low |
| PD-L1 | 0.0423 | HER2-zero | HER2-low |
| CD11c | 0.0494 | HER2-zero | HER2-low |

Pairwise HER2-low versus HER2-zero comparisons were strongest for CD68, CD11c, PD-L1, CD4, and Ki67, but none remained significant after Benjamini-Hochberg correction.

## RNA Validation Check

The first indirect validation layer compared GigaTIME virtual channels with matched RNA-seq marker signatures from the same 30 cases.

Result:

- No channel had an FDR-significant GigaTIME-versus-RNA signature correlation.
- `Ki67` had the strongest positive trend, Spearman rho 0.294.
- `CD68`, `PD-L1`, and `CD11c` did not show strong positive correlations with their matching RNA signatures.

Interpretation: the clinical HER2 virtual immune signal is interesting, but not yet validated. It needs visual QC, more tiles per slide, and stronger orthogonal validation before making biological claims.

## Visual QC Check

The first visual QC pass selected the top `CD68` + `PD-L1` + `CD11c` case from each HER2 group:

| Clinical HER2 group | Selected case | Combined signal |
|---|---|---:|
| HER2-positive | TCGA-A2-A0EQ | 0.115 |
| HER2-low | TCGA-A2-A04Q | 0.086 |
| HER2-zero | TCGA-A2-A0T2 | 0.126 |

The high-scoring tiles were tissue-containing and cellular rather than obvious blank regions. This supports continued investigation, but it does not validate the virtual marker biology. The selected HER2-positive case also had visually plausible high-signal tiles, so the current result should be framed as a slide-level pilot trend rather than a clean visual separation.

## Commands

The clinical HER2 run used:

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

The clinical HER2 summary used:

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

## Caveat

This is still a pilot, not a definitive biological result. It is stronger than the first ERBB2-expression proof-of-work because it uses clinical HER2 groups and a balanced 10/10/10 design. The next scientific step is to test robustness with more tiles per slide, stronger tissue QC, and indirect validation against RNA-seq immune signatures or another dataset with real mIF.
