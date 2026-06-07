# Orion CRC in-domain specificity baseline — CRC01 (released GigaTIME)

Status: in-domain (CRC), in-modality (protein-vs-protein) specificity audit of the RELEASED (lung-trained) GigaTIME on Orion CRC. Per 256px H&E tile, virtual channel activation is correlated against the density of marker-positive cells (Otsu-gated single-cell intensities), with total cells/tile as the cellularity control — the same audit used for the breast RNA work, reusing its stats core.

- H&E: 57360 x 78417 px (0.325 um/px); 1,620,375 segmented cells; 2900 tiles (>=1 cell) used.

## Per-channel specificity (cellularity-controlled)

| Channel (Orion marker) | raw Spearman | partial r \| total cells | partial 95% CI | own-marker row-max? |
|---|---:|---:|---|:--:|
| CK (Pan-CK) | 0.855 | 0.654 | [0.623, 0.683] | yes |
| PD-1 (PD-1) | 0.693 | 0.046 | [0.005, 0.084] | yes |
| Ki67 (Ki67) | 0.676 | 0.224 | [0.185, 0.263] | no |
| CD3 (CD3e) | 0.664 | 0.520 | [0.483, 0.553] | yes |
| CD4 (CD4) | 0.603 | 0.509 | [0.474, 0.542] | no |
| CD8 (CD8a) | 0.386 | -0.022 | [-0.059, 0.013] | no |
| PD-L1 (PD-L1) | 0.301 | -0.113 | [-0.153, -0.066] | no |
| CD20 (CD20) | 0.249 | 0.196 | [0.154, 0.238] | no |
| CD68 (CD68) | 0.196 | -0.159 | [-0.198, -0.122] | no |

Own-marker is the row-max for **3/9** channels; cellularity-controlled partial r stays >0 for **6/9**.

## Output Files

- `CRC01` baseline JSON in `results/gigatime_orion_baseline/`
