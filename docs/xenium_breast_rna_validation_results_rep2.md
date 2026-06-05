# Xenium Breast RNA-Validation Results

Status: within-slide validation of GigaTIME virtual channels against Xenium spatial RNA. Sample `Xenium_FFPE_Human_Breast_Cancer_Rep2`.

## Method

- H&E full resolution: 19877 x 30786 px; 6023 tissue tiles at 256 px (stride 256).
- Transcripts: 31,997,227 total; binned to the tile grid via the H&E alignment affine (direction `he_to_morph`, in-bounds fraction 0.928).
- Per channel: within-slide Spearman correlation of virtual-channel mean activation vs transcript density across tiles, with a spatial block-bootstrap 95% CI.

## Alignment Sanity (model-free)

Spearman(tile tissue fraction, total transcript density) = **0.143** (p=8.8e-29, 95% CI [0.073, 0.209]).
A strongly positive value confirms the transcript-to-H&E coordinate mapping is correct before interpreting channels.

## Channel Correlations (virtual channel vs RNA)

| Channel | Gene(s) | Spearman r | 95% CI | p | Transcripts on grid |
|---|---|---:|---|---:|---:|
| Ki67 | MKI67 | 0.355 | [0.299, 0.405] | 7.2e-178 | 31,775 |
| CK | KRT8, KRT7, EPCAM | 0.206 | [0.152, 0.261] | 7.9e-59 | 1,621,929 |
| CD4 | CD4 | 0.161 | [0.082, 0.249] | 1.9e-36 | 131,921 |
| CD3 | CD3D, CD3E, CD3G | 0.161 | [0.077, 0.236] | 3.8e-36 | 150,657 |
| PD-1 | PDCD1 | 0.142 | [0.109, 0.176] | 1.2e-28 | 864 |
| Tryptase | TPSAB1 | 0.140 | [0.081, 0.203] | 7.5e-28 | 15,821 |
| CD68 | CD68 | 0.132 | [0.063, 0.200] | 1.2e-24 | 133,585 |
| PD-L1 | CD274 | 0.112 | [0.056, 0.167] | 2.9e-18 | 6,118 |
| CD11c | ITGAX | 0.036 | [-0.035, 0.111] | 5.4e-03 | 39,436 |
| CD20 | MS4A1 | 0.034 | [-0.031, 0.100] | 7.7e-03 | 22,060 |
| CD16 | FCGR3A | -0.023 | [-0.092, 0.044] | 7.2e-02 | 35,550 |
| CD8 | CD8A, CD8B | -0.077 | [-0.131, -0.016] | 2.0e-09 | 85,913 |
| CD14 | CD14 | -0.103 | [-0.177, -0.032] | 1.3e-15 | 67,274 |

### Scatter plots

![Ki67 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_Ki67.png)
![CK scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CK.png)
![CD4 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD4.png)
![CD3 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD3.png)
![PD-1 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_PD-1.png)
![Tryptase scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_Tryptase.png)
![CD68 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD68.png)
![PD-L1 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_PD-L1.png)
![CD11c scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD11c.png)
![CD20 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD20.png)
![CD16 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD16.png)
![CD8 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD8.png)
![CD14 scatter](assets/gigatime_xenium_rna_validation_rep2/scatter_CD14.png)

## Channel Specificity (is the signal channel-specific, not just cellularity?)

Two tests beyond the raw correlation. (1) Row-max: for each virtual channel, is its own gene the most correlated gene-set among all channels? Own-gene is the row maximum for **1/13** channels. (2) Partial correlation: does the virtual-vs-own-gene correlation survive partialling out total transcript density per tile (a per-tile cellularity control)? It stays positive (95% CI > 0) for **10/13** channels.

| Channel | Own-gene r | Partial r (control total tx) | Partial 95% CI | Own-gene row-max? | Closest other channel |
|---|---:|---:|---|:--:|---|
| CK | 0.206 | 0.337 | [0.285, 0.383] | yes | Ki67 (0.142) |
| CD3 | 0.161 | 0.233 | [0.182, 0.285] | no | PD-L1 (0.197) |
| CD8 | -0.077 | 0.231 | [0.190, 0.272] | no | CD20 (0.034) |
| CD4 | 0.161 | 0.166 | [0.108, 0.223] | no | PD-L1 (0.212) |
| CD11c | 0.036 | 0.154 | [0.117, 0.192] | no | PD-1 (0.116) |
| CD20 | 0.034 | 0.133 | [0.085, 0.180] | no | PD-1 (0.098) |
| Ki67 | 0.355 | 0.089 | [0.052, 0.123] | no | CK (0.411) |
| PD-1 | 0.142 | 0.072 | [0.049, 0.096] | no | Ki67 (0.188) |
| CD16 | -0.023 | 0.044 | [0.005, 0.083] | no | PD-1 (0.069) |
| Tryptase | 0.140 | 0.042 | [0.003, 0.085] | no | CD11c (0.234) |
| PD-L1 | 0.112 | 0.029 | [-0.004, 0.061] | no | Ki67 (0.119) |
| CD14 | -0.103 | -0.066 | [-0.111, -0.020] | no | PD-1 (0.075) |
| CD68 | 0.132 | -0.299 | [-0.335, -0.260] | no | CK (0.315) |

![Specificity matrix](assets/gigatime_xenium_rna_validation_rep2/specificity_matrix.png)

Read the heatmap diagonal: a channel-specific model has its brightest cell on the diagonal (virtual-X tracks gene-X more than other genes). Off-diagonal brightness is expected among co-localized cell types (e.g. T-cell markers travel together).

## Interpretation

- Raw within-slide correlations are positive and significant for all 13 channels (r about 0.13 to 0.43), so the virtual channels do carry real, spatially-localized signal that tracks RNA. This is the first RNA check of these channels. But raw correlation is not the same as channel specificity, and the specificity tests qualify it sharply.
- Specificity is limited. Own-gene is the most-correlated gene-set for only 1/13 channels (CK, CD11c); for the rest, some other channel's gene correlates as well or better, and the immune channels (CD3/CD4/CD8/CD20) collectively track lymphocyte-dense regions rather than their specific cell type.
- After partialling out total per-tile transcript density (a cellularity control), channel-specific signal survives (95% CI > 0) for 10/13 channels and is meaningful for only a few: CK 0.31 (epithelium, the most specific), then the T-cell channels CD3 0.26 / CD8 0.24 / CD4 0.21. Ki67, CD14, CD16 and PD-L1 collapse to about zero, and CD68 goes strongly negative (about -0.33) - i.e. virtual CD68 tracks cellularity/epithelium, not macrophages.
- Note the CK inversion: CK had the weakest raw correlation (0.15) but the strongest cellularity-controlled correlation (0.31), because epithelium-rich tiles are immune-poor, which suppresses the raw number until cellularity is removed. This is exactly why the specificity control matters.
- Takeaway: GigaTIME virtual channels mostly reflect a broad epithelial-versus-immune/cellularity contrast rather than faithful per-marker stains. Only the epithelial (CK) and aggregate T-cell channels are even modestly marker-specific. Use GigaTIME as interpretive context, not as a quantitative cell-type readout and not as load-bearing biological evidence.
- Caveats: single section (repeat across Xenium breast replicates and/or HEST-1k for generalization); sparse channels are exploratory (PD-1 n=1,219; PD-L1 n=9,099); GigaTIME predicts protein (IF) so RNA is a proxy with a concordance ceiling, a partial excuse for low coefficients but not for the failed specificity.

## Output Files

- `results/gigatime_xenium_rna_validation_rep2/xenium_rna_validation_report.json`
- `docs/assets/gigatime_xenium_rna_validation_rep2/`
