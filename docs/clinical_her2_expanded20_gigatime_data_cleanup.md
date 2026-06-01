# Clinical HER2 GigaTIME Data Cleanup

This cleanup step goes back before classifier training. It asks whether the GigaTIME input features should be aggregated from all sampled tissue tiles, or from cleaner tile subsets that are more cellular and more epithelial/tumor enriched.

## Cleanup Views

- `all_sampled_tissue`: all tissue tiles selected by the original GigaTIME run.
- `qc_cellular_tissue`: tiles with H&E tissue fraction >= 0.70 and virtual DAPI mean >= 0.05.
- `ck_enriched_top50`: the top 50% virtual CK tiles within each slide after cellular-tissue QC.
- `ck_enriched_top25`: the top 25% virtual CK tiles within each slide after cellular-tissue QC.

Important: CK and DAPI are still GigaTIME virtual predictions from H&E, not laboratory stains. These filters create tumor-enriched research feature views, not confirmed tumor masks.

## Tile Retention

| Cleanup view | Median retained tiles | Median retained fraction | Median DAPI | Median CK |
| --- | --- | --- | --- | --- |
| All sampled tissue | 256.0 | 1.000 | 0.327 | 0.209 |
| QC cellular tissue | 191.5 | 0.748 | 0.356 | 0.215 |
| CK-enriched top 50% | 96.5 | 0.377 | 0.444 | 0.326 |
| CK-enriched top 25% | 48.5 | 0.189 | 0.498 | 0.415 |

![Tile retention by cleanup view](assets/clinical_her2_expanded20_gigatime_cleanup/cleanup_retained_tiles_by_filter.png)

## Tile-Level Cleanup Map

This plot shows how the cleanup rules move from broad tissue tiles toward cellular, CK-enriched tiles.

![CK and DAPI tile distribution](assets/clinical_her2_expanded20_gigatime_cleanup/cleanup_ck_dapi_distribution.png)

## Group Means After Cleanup

![Cleaned GigaTIME group mean heatmap](assets/clinical_her2_expanded20_gigatime_cleanup/cleanup_key_channel_heatmap.png)

## Top Three-Group Signals Across Cleanup Views

| Cleanup view | Channel | Kruskal p | BH q within view | Highest group | Lowest group | Max-min mean |
| --- | --- | --- | --- | --- | --- | --- |
| All sampled tissue | CD4 | 0.0068 | 0.0383 | HER2-positive | HER2-low | 0.0506 |
| All sampled tissue | CD3 | 0.0085 | 0.0383 | HER2-positive | HER2-low | 0.0495 |
| All sampled tissue | CD11c | 0.0136 | 0.0407 | HER2-positive | HER2-low | 0.0137 |
| All sampled tissue | CD68 | 0.0208 | 0.0467 | HER2-positive | HER2-low | 0.0218 |
| All sampled tissue | CD20 | 0.0306 | 0.0550 | HER2-positive | HER2-low | 0.0355 |
| QC cellular tissue | CD4 | 0.0051 | 0.0277 | HER2-positive | HER2-low | 0.0676 |
| QC cellular tissue | CD3 | 0.0061 | 0.0277 | HER2-positive | HER2-low | 0.0668 |
| QC cellular tissue | CD11c | 0.0110 | 0.0331 | HER2-positive | HER2-low | 0.0192 |
| QC cellular tissue | CD20 | 0.0223 | 0.0427 | HER2-positive | HER2-low | 0.0454 |
| QC cellular tissue | CD68 | 0.0237 | 0.0427 | HER2-positive | HER2-low | 0.0297 |
| CK-enriched top 50% | CD3 | 0.0072 | 0.0391 | HER2-positive | HER2-low | 0.0634 |
| CK-enriched top 50% | CD4 | 0.0087 | 0.0391 | HER2-positive | HER2-low | 0.0644 |
| CK-enriched top 50% | CD68 | 0.0217 | 0.0652 | HER2-positive | HER2-low | 0.0293 |
| CK-enriched top 50% | CD11c | 0.0378 | 0.0850 | HER2-positive | HER2-low | 0.0189 |
| CK-enriched top 50% | CK | 0.0563 | 0.1013 | HER2-zero | HER2-positive | 0.0781 |
| CK-enriched top 25% | CD68 | 0.0559 | 0.2030 | HER2-positive | HER2-low | 0.0291 |
| CK-enriched top 25% | CD4 | 0.0627 | 0.2030 | HER2-positive | HER2-low | 0.0585 |
| CK-enriched top 25% | CK | 0.0678 | 0.2030 | HER2-zero | HER2-positive | 0.0808 |
| CK-enriched top 25% | CD3 | 0.0902 | 0.2030 | HER2-positive | HER2-low | 0.0573 |
| CK-enriched top 25% | CD11c | 0.2076 | 0.3737 | HER2-positive | HER2-low | 0.0179 |

## HER2-Low Versus HER2-Zero Focus

Negative delta means HER2-low is lower than HER2-zero.

| Cleanup view | Channel | HER2-low minus HER2-zero | Mann-Whitney p | BH q within view |
| --- | --- | --- | --- | --- |
| All sampled tissue | CD68 | -0.0070 | 0.0051 | 0.0326 |
| All sampled tissue | PD-L1 | -0.0145 | 0.0123 | 0.0556 |
| All sampled tissue | CD11c | -0.0051 | 0.0026 | 0.0252 |
| QC cellular tissue | CD68 | -0.0079 | 0.0066 | 0.0320 |
| QC cellular tissue | PD-L1 | -0.0180 | 0.0071 | 0.0320 |
| QC cellular tissue | CD11c | -0.0062 | 0.0023 | 0.0206 |
| CK-enriched top 50% | CD68 | -0.0076 | 0.0084 | 0.0752 |
| CK-enriched top 50% | PD-L1 | -0.0122 | 0.0315 | 0.0946 |
| CK-enriched top 50% | CD11c | -0.0040 | 0.0114 | 0.0772 |
| CK-enriched top 25% | CD68 | -0.0065 | 0.0179 | 0.2433 |
| CK-enriched top 25% | PD-L1 | -0.0059 | 0.1136 | 0.3066 |
| CK-enriched top 25% | CD11c | -0.0024 | 0.0720 | 0.2948 |

![Key channel boxplots after cleanup](assets/clinical_her2_expanded20_gigatime_cleanup/cleanup_key_channel_boxplots.png)

## Interpretation

This cleanup does not validate the virtual markers, but it makes the next classifier input more biologically defensible. The original baseline averaged all sampled tissue tiles. The cleaned tables let us rerun summaries or classifiers using cellular tissue tiles and CK-enriched tumor-context tiles.

If HER2-low versus HER2-zero signal remains after CK enrichment, it is less likely to be explained only by blank tissue or broad non-cellular sampling. If the signal disappears, then the original classifier may have been leaning on non-tumor tissue composition rather than tumor-region biology.

## Outputs

- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/tile_qc_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/cleaned_slide_features.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/filter_retention_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/cleanup_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/cleanup_pairwise_tests.csv`

## Next Step

Use `cleaned_slide_features.csv` to rerun the classifier separately for each cleanup view, especially `qc_cellular_tissue`, `ck_enriched_top50`, and `ck_enriched_top25`. The comparison should show whether tumor-enriched GigaTIME features improve HER2 prediction or expose the current model as tissue-composition driven.
