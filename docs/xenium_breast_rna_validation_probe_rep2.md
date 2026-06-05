# Xenium Breast RNA-Validation Feasibility Probe

Status: feasibility probe for validating GigaTIME virtual channels against Xenium breast RNA. Sample `Xenium_FFPE_Human_Breast_Cancer_Rep2`.

## Purpose

Confirm that the 10x Xenium Human Breast dataset can serve as a within-slide RNA ground truth for
GigaTIME virtual immune/tumor channels before wiring GigaTIME inference to it. Within-slide transcript
co-localization sidesteps the TCGA composition/batch confound, because the correlation is computed inside
one section rather than across patients.

## Verdict

- Feasible for RNA validation: **NO**
- Every GigaTIME channel maps to a panel gene: True
- Every channel gene has transcripts: True
- H&E<->Xenium alignment invertible: True
- Alignment extent cross-check: False

## Channel-Gene Coverage

Panel genes: 313. Channels covered: 6 / 6.

| GigaTIME channel | Panel gene(s) present | Transcript count |
|---|---|---:|
| CD3 | CD3D, CD3E, CD3G | 157,492 |
| CD8 | CD8A, CD8B | 91,353 |
| PD-L1 | CD274 | 6,891 |
| CK | KRT8, KRT7, EPCAM | 1,698,279 |
| Ki67 | MKI67 | 34,143 |
| myeloid | CD68, CD163, LYZ, ITGAX | 566,204 |

## Transcripts

- Total transcripts: 31,997,227 across 313 detected genes (146,838 control-probe transcripts).
- Feature column: `feature_name`.
- Coordinate extent (microns): x [0.0, 7522.0], y [4.2, 5472.4].

## Alignment

- Matrix shape: [3, 3], invertible: True, determinant: -2.9329007894366765.
- Matrix (H&E pixel -> Xenium morphology pixel, homogeneous):

```
 1.71255   0.00835507  -8115.57
 0.00835507  -1.71255   19334.6
 0   0   1
```

## H&E Image

- Full-resolution shape: [19877, 30786, 3], pyramid levels: 7, microns/pixel tag: 0.36378800370496256.
- Cross-check: H&E frame maps to transcript extent with overlap x=0.67, y=0.48 (applies_cleanly=False, xenium_mpp=0.2125).

## Next Steps

1. If feasible: tile the post-Xenium H&E at GigaTIME's expected microns/pixel and run virtual-channel inference.
2. Bin transcripts to the same tile grid via the alignment transform; sum each channel gene per tile.
3. Within-slide Spearman correlation of virtual-channel intensity vs transcript density per tile, per channel,
   with a block-bootstrap CI over tiles. Lead on CD8A/CD3D/MKI67/keratins; report CD274/PD-L1 with the known
   RNA-protein concordance caveat.
4. For multi-slide breadth, repeat across Xenium breast replicates and/or HEST-1k breast Visium samples.

## Sources

- Janesick et al., Nat Commun 2023 (Xenium FFPE Human Breast Cancer).
- 10x Xenium Human Breast Dataset Explorer: https://www.10xgenomics.com/products/xenium-in-situ/human-breast-dataset-explorer
- CC BY 4.0 broader-panel alternative: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer
