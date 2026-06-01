# Clinical HER2 Tile-Sampling Robustness Check

Status: Historical 30-slide robustness report. For the latest expanded 60-slide result, use `docs/clinical_her2_expanded20_results.md`.

This document records the first robustness check after the 64-tile clinical HER2 pilot. The goal was to test whether the HER2-zero versus HER2-low virtual immune/checkpoint signal remained visible when each slide was sampled more densely.

## Why This Was Done

The 30-slide clinical HER2 pilot used 64 random tissue tiles per slide. That was enough for a fast feasibility run, but a whole-slide image is large and heterogeneous. A 64-tile sample can miss important tissue regions or overrepresent a small part of the slide.

The next question was:

> If we rerun the exact same 30 slides with more tissue tiles per slide, does the same HER2-zero greater than HER2-low signal remain?

This is a robustness check. It does not validate the virtual markers against real mIF, but it helps determine whether the signal was mainly a random tile-sampling artifact.

## Command

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

Then the same summary, RNA validation, and visual QC scripts were rerun against the `tile256` output directory.

## Outputs

Local result outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/rna_validation_summary.md`

Tracked figure outputs:

- `docs/assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png`
- `docs/assets/clinical_her2_tile256/clinical_her2_channel_boxplots.png`
- `docs/assets/clinical_her2_tile256/gigatime_rna_correlation_heatmap.png`
- `docs/assets/clinical_her2_visual_qc_tile256/*_he_vs_virtual_mif_qc.png`

## Main Result

The denser 256-tile run reproduced the same main pattern from the 64-tile pilot:

- `HER2-zero` had the highest mean virtual signal.
- `HER2-low` had the lowest mean virtual signal.
- The clearest channels were again `CD68`, `PD-L1`, and `CD11c`.

| Channel | 64-tile Kruskal p | 256-tile Kruskal p | 64-tile max-min | 256-tile max-min | 256-tile highest group | 256-tile lowest group |
|---|---:|---:|---:|---:|---|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero | HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero | HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero | HER2-low |
| CD3 | 0.1876 | 0.0686 | 0.02661 | 0.03289 | HER2-zero | HER2-low |
| CD4 | 0.0794 | 0.0819 | 0.02739 | 0.03277 | HER2-zero | HER2-low |
| Ki67 | 0.0920 | 0.1033 | 0.00416 | 0.00400 | HER2-zero | HER2-low |

![256-tile clinical HER2 group mean heatmap](assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png)

## Pairwise Results

The strongest 256-tile pairwise tests were again HER2-low versus HER2-zero:

| Channel | Comparison | Delta mean | Mann-Whitney p | BH q |
|---|---|---:|---:|---:|
| CD68 | HER2-low vs HER2-zero | -0.01044 | 0.00729 | 0.1133 |
| PD-L1 | HER2-low vs HER2-zero | -0.02061 | 0.00911 | 0.1133 |
| CD11c | HER2-low vs HER2-zero | -0.00504 | 0.01133 | 0.1133 |
| CD3 | HER2-low vs HER2-zero | -0.03289 | 0.03121 | 0.1560 |
| PD-1 | HER2-low vs HER2-zero | -0.05282 | 0.02575 | 0.1560 |
| CD4 | HER2-low vs HER2-zero | -0.03277 | 0.03121 | 0.1560 |

The corrected q values improved compared with the 64-tile run, but they still did not pass the usual 0.05 FDR threshold. This is stronger pilot evidence, not a definitive result.

## RNA Validation

The 256-tile rerun did not strongly validate the virtual immune-channel signal against matched bulk RNA-seq marker signatures.

| Channel | Spearman rho | p | BH q |
|---|---:|---:|---:|
| Ki67 | 0.252 | 0.179 | 0.5968 |
| PD-L1 | 0.103 | 0.5897 | 0.7371 |
| CD68 | -0.036 | 0.851 | 0.9284 |
| CD11c | -0.184 | 0.3292 | 0.6038 |
| CD20 | -0.337 | 0.0686 | 0.5968 |

Interpretation: the tile-sampling result is more robust, but the virtual marker biology is still not orthogonally validated by bulk RNA-seq.

![256-tile RNA correlation heatmap](assets/clinical_her2_tile256/gigatime_rna_correlation_heatmap.png)

## Broader RNA Program Follow-Up

After this robustness check, the RNA validation was expanded to broader immune and tissue programs. That analysis also did not positively confirm the virtual immune/checkpoint signal.

Key result:

- Virtual myeloid/checkpoint remained higher in HER2-zero than HER2-low, but was not FDR-significant: p 0.0176, BH q 0.0878.
- No broad RNA immune program showed an FDR-significant HER2-group difference.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal.

See `docs/clinical_her2_rna_program_validation.md`.

## Visual QC

The 256-tile visual QC selected the same representative cases as the 64-tile QC:

| Clinical HER2 group | Selected case | Combined CD68 + PD-L1 + CD11c | mean CD68 | mean PD-L1 | mean CD11c |
|---|---|---:|---:|---:|---:|
| HER2-positive | TCGA-A2-A0EQ | 0.110 | 0.028 | 0.069 | 0.013 |
| HER2-low | TCGA-A2-A04Q | 0.080 | 0.017 | 0.054 | 0.009 |
| HER2-zero | TCGA-A2-A0T2 | 0.140 | 0.041 | 0.077 | 0.022 |

![256-tile HER2-zero H&E versus virtual mIF QC](assets/clinical_her2_visual_qc_tile256/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png)

The selected HER2-zero case had a clearer combined signal in the 256-tile run. The top tiles remained tissue-containing and cellular, supporting continued investigation. This still does not prove that the predicted channels correspond to real `CD68`, `PD-L1`, or `CD11c` protein staining.

## Interpretation

What became stronger:

- The same HER2-zero greater than HER2-low signal persisted after denser sampling.
- The top three channels stayed stable: `CD68`, `PD-L1`, and `CD11c`.
- The group mean gaps for those channels increased slightly.
- The leading pairwise FDR values improved from about 0.21 to about 0.11.

What remains unresolved:

- The results still do not pass FDR correction.
- RNA-seq validation remains weak.
- TCGA does not provide matched real mIF for these exact slides in this project.
- A pathologist has not yet reviewed whether the high-signal H&E regions look biologically plausible.

## Proposal Language

A careful way to describe the robustness result:

> After increasing slide sampling from 64 to up to 256 tissue tiles per slide, the same HER2-zero greater than HER2-low GigaTIME virtual immune/checkpoint signal persisted for CD68, PD-L1, and CD11c. This supports robustness to tile sampling, but the signal remains hypothesis-generating because RNA-seq marker-signature validation was weak and no matched real mIF ground truth is available.

## Next Step

The next step after this robustness check was broader RNA-program validation and a first held-out classifier baseline. Both are now complete.

The next scientific step should be validation and better classifier input design rather than another simple rerun:

1. Ask an advisor/pathologist to review the high-signal H&E tiles and virtual mIF panels.
2. Restrict the next classifier to tumor-rich tiles rather than all sampled tissue tiles.
3. Add tile distribution features and, if available, GigaTIME/pathology embeddings.
4. Add tumor purity or immune deconvolution covariates if available.
5. Check whether endothelial/stromal/tissue-composition differences might explain part of the virtual signal.
6. Look for an external dataset with paired H&E and real mIF to directly evaluate GigaTIME predictions.

See `docs/clinical_her2_classifier_baseline.md` for the first classifier baseline.
