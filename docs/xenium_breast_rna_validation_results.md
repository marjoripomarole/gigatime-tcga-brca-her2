# Xenium Breast RNA-Validation Results

Status: within-slide validation of GigaTIME virtual channels against Xenium spatial RNA. Sample `Xenium_FFPE_Human_Breast_Cancer_Rep1`.

## Method

- H&E full resolution: 24241 x 30786 px; 6460 tissue tiles at 256 px (stride 256).
- Transcripts: 42,638,083 total; binned to the tile grid via the H&E alignment affine (direction `he_to_morph`, in-bounds fraction 1.000).
- Per channel: within-slide Spearman correlation of virtual-channel mean activation vs transcript density across tiles, with a spatial block-bootstrap 95% CI.

## Alignment Sanity (model-free)

Spearman(tile tissue fraction, total transcript density) = **0.181** (p=1.2e-48, 95% CI [0.120, 0.249]).
A strongly positive value confirms the transcript-to-H&E coordinate mapping is correct before interpreting channels.

## Channel Correlations (virtual channel vs RNA)

| Channel | Gene(s) | Spearman r | 95% CI | p | Transcripts on grid |
|---|---|---:|---|---:|---:|
| CD3 | CD3D, CD3E, CD3G | 0.428 | [0.365, 0.490] | 1.3e-286 | 204,300 |
| CD4 | CD4 | 0.428 | [0.366, 0.485] | 2.0e-286 | 177,629 |
| CD8 | CD8A, CD8B | 0.358 | [0.303, 0.409] | 1.7e-194 | 113,636 |
| Ki67 | MKI67 | 0.337 | [0.266, 0.402] | 2.6e-171 | 58,975 |
| Tryptase | TPSAB1 | 0.331 | [0.286, 0.378] | 2.5e-165 | 21,714 |
| CD11c | ITGAX | 0.330 | [0.273, 0.384] | 2.8e-164 | 52,613 |
| PD-L1 | CD274 | 0.223 | [0.171, 0.272] | 1.8e-73 | 9,099 |
| CD16 | FCGR3A | 0.221 | [0.164, 0.277] | 2.8e-72 | 51,735 |
| PD-1 | PDCD1 | 0.214 | [0.181, 0.246] | 6.0e-68 | 1,219 |
| CD20 | MS4A1 | 0.214 | [0.154, 0.272] | 1.5e-67 | 28,637 |
| CD14 | CD14 | 0.167 | [0.114, 0.224] | 2.1e-41 | 87,746 |
| CK | KRT8, KRT7, EPCAM | 0.146 | [0.078, 0.210] | 6.0e-32 | 2,706,696 |
| CD68 | CD68 | 0.132 | [0.072, 0.196] | 1.6e-26 | 180,150 |

### Scatter plots

![CD3 scatter](assets/gigatime_xenium_rna_validation/scatter_CD3.png)
![CD4 scatter](assets/gigatime_xenium_rna_validation/scatter_CD4.png)
![CD8 scatter](assets/gigatime_xenium_rna_validation/scatter_CD8.png)
![Ki67 scatter](assets/gigatime_xenium_rna_validation/scatter_Ki67.png)
![Tryptase scatter](assets/gigatime_xenium_rna_validation/scatter_Tryptase.png)
![CD11c scatter](assets/gigatime_xenium_rna_validation/scatter_CD11c.png)
![PD-L1 scatter](assets/gigatime_xenium_rna_validation/scatter_PD-L1.png)
![CD16 scatter](assets/gigatime_xenium_rna_validation/scatter_CD16.png)
![PD-1 scatter](assets/gigatime_xenium_rna_validation/scatter_PD-1.png)
![CD20 scatter](assets/gigatime_xenium_rna_validation/scatter_CD20.png)
![CD14 scatter](assets/gigatime_xenium_rna_validation/scatter_CD14.png)
![CK scatter](assets/gigatime_xenium_rna_validation/scatter_CK.png)
![CD68 scatter](assets/gigatime_xenium_rna_validation/scatter_CD68.png)

## Channel Specificity (is the signal channel-specific, not just cellularity?)

Two tests beyond the raw correlation. (1) Row-max: for each virtual channel, is its own gene the most correlated gene-set among all channels? Own-gene is the row maximum for **2/13** channels. (2) Partial correlation: does the virtual-vs-own-gene correlation survive partialling out total transcript density per tile (a per-tile cellularity control)? It stays positive (95% CI > 0) for **7/13** channels.

| Channel | Own-gene r | Partial r (control total tx) | Partial 95% CI | Own-gene row-max? | Closest other channel |
|---|---:|---:|---|:--:|---|
| CK | 0.146 | 0.308 | [0.258, 0.355] | yes | Ki67 (0.091) |
| CD3 | 0.428 | 0.256 | [0.192, 0.324] | no | CD11c (0.440) |
| CD8 | 0.358 | 0.237 | [0.177, 0.293] | no | CD3 (0.391) |
| CD4 | 0.428 | 0.211 | [0.143, 0.279] | no | CD11c (0.451) |
| CD11c | 0.330 | 0.147 | [0.098, 0.199] | yes | CD4 (0.310) |
| CD20 | 0.214 | 0.085 | [0.031, 0.136] | no | CD11c (0.363) |
| PD-1 | 0.214 | 0.053 | [0.029, 0.077] | no | CK (0.404) |
| Tryptase | 0.331 | 0.028 | [-0.018, 0.075] | no | CD11c (0.488) |
| Ki67 | 0.337 | -0.020 | [-0.055, 0.019] | no | CK (0.388) |
| CD16 | 0.221 | -0.036 | [-0.078, 0.009] | no | CD11c (0.258) |
| CD14 | 0.167 | -0.048 | [-0.097, 0.009] | no | CD11c (0.305) |
| PD-L1 | 0.223 | -0.072 | [-0.106, -0.034] | no | CK (0.368) |
| CD68 | 0.132 | -0.333 | [-0.366, -0.300] | no | CK (0.403) |

![Specificity matrix](assets/gigatime_xenium_rna_validation/specificity_matrix.png)

Read the heatmap diagonal: a channel-specific model has its brightest cell on the diagonal (virtual-X tracks gene-X more than other genes). Off-diagonal brightness is expected among co-localized cell types (e.g. T-cell markers travel together).

## Interpretation

- Raw within-slide correlations are positive and significant for all 13 channels (r about 0.13 to 0.43), so the virtual channels do carry real, spatially-localized signal that tracks RNA. This is the first RNA check of these channels. But raw correlation is not the same as channel specificity, and the specificity tests qualify it sharply.
- Specificity is limited. Own-gene is the most-correlated gene-set for only 2/13 channels (CK, CD11c); for the rest, some other channel's gene correlates as well or better, and the immune channels (CD3/CD4/CD8/CD20) collectively track lymphocyte-dense regions rather than their specific cell type.
- After partialling out total per-tile transcript density (a cellularity control), channel-specific signal survives (95% CI > 0) for 7/13 channels and is meaningful for only a few: CK 0.31 (epithelium, the most specific), then the T-cell channels CD3 0.26 / CD8 0.24 / CD4 0.21. Ki67, CD14, CD16 and PD-L1 collapse to about zero, and CD68 goes strongly negative (about -0.33) - i.e. virtual CD68 tracks cellularity/epithelium, not macrophages.
- Note the CK inversion: CK had the weakest raw correlation (0.15) but the strongest cellularity-controlled correlation (0.31), because epithelium-rich tiles are immune-poor, which suppresses the raw number until cellularity is removed. This is exactly why the specificity control matters.
- Takeaway: GigaTIME virtual channels mostly reflect a broad epithelial-versus-immune/cellularity contrast rather than faithful per-marker stains. Only the epithelial (CK) and aggregate T-cell channels are even modestly marker-specific. Use GigaTIME as interpretive context, not as a quantitative cell-type readout and not as load-bearing biological evidence.
- Caveats: single section (repeat across Xenium breast replicates and/or HEST-1k for generalization); sparse channels are exploratory (PD-1 n=1,219; PD-L1 n=9,099); GigaTIME predicts protein (IF) so RNA is a proxy with a concordance ceiling, a partial excuse for low coefficients but not for the failed specificity.

## Output Files

- `results/gigatime_xenium_rna_validation/xenium_rna_validation_report.json`
- `docs/assets/gigatime_xenium_rna_validation/`
