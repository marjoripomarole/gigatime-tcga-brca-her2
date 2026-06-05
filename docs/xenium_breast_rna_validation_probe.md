# Xenium Breast RNA-Validation Feasibility Probe

Status: feasibility probe for validating GigaTIME virtual channels against Xenium breast RNA. Sample `Xenium_FFPE_Human_Breast_Cancer_Rep1`.

## Purpose

Confirm that the 10x Xenium Human Breast dataset can serve as a within-slide RNA ground truth for
GigaTIME virtual immune/tumor channels before wiring GigaTIME inference to it. Within-slide transcript
co-localization sidesteps the TCGA composition/batch confound, because the correlation is computed inside
one section rather than across patients.

## Verdict

- Feasible for RNA validation: **YES**
- Every GigaTIME channel maps to a panel gene: True
- Every channel gene has transcripts: True
- H&E<->Xenium alignment invertible: True
- Alignment extent cross-check: None

## Channel-Gene Coverage

Panel genes: 313. Channels covered: 6 / 6.

| GigaTIME channel | Panel gene(s) present | Transcript count |
|---|---|---:|
| CD3 | CD3D, CD3E, CD3G | 209,815 |
| CD8 | CD8A, CD8B | 116,827 |
| PD-L1 | CD274 | 9,348 |
| CK | KRT8, KRT7, EPCAM | 2,729,111 |
| Ki67 | MKI67 | 59,855 |
| myeloid | CD68, CD163, LYZ, ITGAX | 662,371 |

## Transcripts

- Total transcripts: 42,638,083 across 313 detected genes (204,548 control-probe transcripts).
- Feature column: `feature_name`.
- Coordinate extent (microns): x [-1.9, 7522.7], y [4.4, 5473.5].

## Alignment

- Matrix shape: [3, 3], invertible: True, determinant: -2.9327295160409075.
- Matrix (H&E pixel -> Xenium morphology pixel, homogeneous):

```
 1.7125   0.00834508  -10397.2
 0.00834508  -1.7125   37211.5
 0   0   1
```

## H&E Image

- H&E not downloaded (use --include-he).

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
