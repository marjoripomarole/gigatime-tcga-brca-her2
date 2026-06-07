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
| CK | KRT8, KRT7, EPCAM | 0.328 | [0.268, 0.383] | 1.3e-161 | 2,706,696 |
| Ki67 | MKI67 | 0.195 | [0.146, 0.243] | 1.2e-56 | 58,975 |
| CD8 | CD8A, CD8B | 0.194 | [0.140, 0.250] | 5.9e-56 | 113,636 |
| CD68 | CD68 | 0.182 | [0.131, 0.236] | 3.0e-49 | 180,150 |
| PD-L1 | CD274 | 0.145 | [0.106, 0.184] | 1.1e-31 | 9,099 |
| CD14 | CD14 | 0.119 | [0.072, 0.168] | 5.8e-22 | 87,746 |
| CD4 | CD4 | 0.050 | [0.003, 0.099] | 5.9e-05 | 177,629 |
| CD3 | CD3D, CD3E, CD3G | 0.043 | [-0.003, 0.091] | 4.8e-04 | 204,300 |
| CD11c | ITGAX | 0.042 | [-0.013, 0.094] | 6.4e-04 | 52,613 |
| CD20 | MS4A1 | -0.006 | [-0.057, 0.044] | 6.4e-01 | 28,637 |
| PD-1 | PDCD1 | -0.048 | [-0.080, -0.016] | 1.2e-04 | 1,219 |

### Scatter plots

![CK scatter](assets/rosie_xenium_rna_validation/scatter_CK.png)
![Ki67 scatter](assets/rosie_xenium_rna_validation/scatter_Ki67.png)
![CD8 scatter](assets/rosie_xenium_rna_validation/scatter_CD8.png)
![CD68 scatter](assets/rosie_xenium_rna_validation/scatter_CD68.png)
![PD-L1 scatter](assets/rosie_xenium_rna_validation/scatter_PD-L1.png)
![CD14 scatter](assets/rosie_xenium_rna_validation/scatter_CD14.png)
![CD4 scatter](assets/rosie_xenium_rna_validation/scatter_CD4.png)
![CD3 scatter](assets/rosie_xenium_rna_validation/scatter_CD3.png)
![CD11c scatter](assets/rosie_xenium_rna_validation/scatter_CD11c.png)
![CD20 scatter](assets/rosie_xenium_rna_validation/scatter_CD20.png)
![PD-1 scatter](assets/rosie_xenium_rna_validation/scatter_PD-1.png)

## Channel Specificity (is the signal channel-specific, not just cellularity?)

Two tests beyond the raw correlation. (1) Row-max: for each virtual channel, is its own gene the most correlated gene-set among all channels? Own-gene is the row maximum for **2/11** channels. (2) Partial correlation: does the virtual-vs-own-gene correlation survive partialling out total transcript density per tile (a per-tile cellularity control)? It stays positive (95% CI > 0) for **7/11** channels.

| Channel | Own-gene r | Partial r (control total tx) | Partial 95% CI | Own-gene row-max? | Closest other channel |
|---|---:|---:|---|:--:|---|
| CD4 | 0.050 | 0.336 | [0.295, 0.374] | no | CD20 (0.123) |
| CD11c | 0.042 | 0.296 | [0.259, 0.330] | no | CD20 (0.190) |
| CD14 | 0.119 | 0.245 | [0.210, 0.279] | yes | CD20 (0.115) |
| PD-L1 | 0.145 | 0.167 | [0.137, 0.199] | no | CD3 (0.207) |
| CD3 | 0.043 | 0.099 | [0.060, 0.139] | no | CD14 (0.057) |
| CD8 | 0.194 | 0.090 | [0.039, 0.142] | no | CD3 (0.207) |
| CD68 | 0.182 | 0.047 | [-0.004, 0.100] | no | CD11c (0.247) |
| PD-1 | -0.048 | 0.033 | [0.004, 0.063] | no | CD20 (0.149) |
| CK | 0.328 | -0.030 | [-0.067, 0.009] | yes | Ki67 (0.291) |
| CD20 | -0.006 | -0.119 | [-0.153, -0.081] | no | CK (0.173) |
| Ki67 | 0.195 | -0.133 | [-0.162, -0.104] | no | CD11c (0.217) |

![Specificity matrix](assets/rosie_xenium_rna_validation/specificity_matrix.png)

Read the heatmap diagonal: a channel-specific model has its brightest cell on the diagonal (virtual-X tracks gene-X more than other genes). Off-diagonal brightness is expected among co-localized cell types (e.g. T-cell markers travel together).

## Interpretation

- Raw within-slide correlations are positive and significant for all 13 channels (r about 0.13 to 0.43), so the virtual channels do carry real, spatially-localized signal that tracks RNA. This is the first RNA check of these channels. But raw correlation is not the same as channel specificity, and the specificity tests qualify it sharply.
- Specificity is limited. Own-gene is the most-correlated gene-set for only 2/11 channels (CK, CD11c); for the rest, some other channel's gene correlates as well or better, and the immune channels (CD3/CD4/CD8/CD20) collectively track lymphocyte-dense regions rather than their specific cell type.
- After partialling out total per-tile transcript density (a cellularity control), channel-specific signal survives (95% CI > 0) for 7/11 channels and is meaningful for only a few: CK 0.31 (epithelium, the most specific), then the T-cell channels CD3 0.26 / CD8 0.24 / CD4 0.21. Ki67, CD14, CD16 and PD-L1 collapse to about zero, and CD68 goes strongly negative (about -0.33) - i.e. virtual CD68 tracks cellularity/epithelium, not macrophages.
- Note the CK inversion: CK had the weakest raw correlation (0.15) but the strongest cellularity-controlled correlation (0.31), because epithelium-rich tiles are immune-poor, which suppresses the raw number until cellularity is removed. This is exactly why the specificity control matters.
- Takeaway: GigaTIME virtual channels mostly reflect a broad epithelial-versus-immune/cellularity contrast rather than faithful per-marker stains. Only the epithelial (CK) and aggregate T-cell channels are even modestly marker-specific. Use GigaTIME as interpretive context, not as a quantitative cell-type readout and not as load-bearing biological evidence.
- Caveats: single section (repeat across Xenium breast replicates and/or HEST-1k for generalization); sparse channels are exploratory (PD-1 n=1,219; PD-L1 n=9,099); GigaTIME predicts protein (IF) so RNA is a proxy with a concordance ceiling, a partial excuse for low coefficients but not for the failed specificity.

## Output Files

- `results/rosie_xenium_rna_validation/xenium_rna_validation_report.json`
- `docs/assets/rosie_xenium_rna_validation/`
