# Clinical HER2 RNA Program Validation

This document records the next validation step after the 256-tile clinical HER2 robustness run.

## Why This Was Done

The marker-level RNA validation compared one GigaTIME virtual channel with a small matching RNA marker set. That is useful, but it can be too narrow. For example, a macrophage-like tissue signal may not be captured well by only a few genes, and a checkpoint-like pattern may reflect a broader inflamed program rather than a single marker.

So this step tested broader RNA programs:

- T cell / cytotoxic
- Checkpoint / IFNG
- Myeloid / macrophage
- Dendritic / antigen-presenting cell
- B cell
- Proliferation
- Epithelial / tumor
- Stromal / fibroblast
- Endothelial

It also created broader GigaTIME virtual programs:

- Virtual myeloid/checkpoint: `CD68`, `CD11c`, `PD-L1`
- Virtual T cell/checkpoint: `CD3`, `CD4`, `CD8`, `PD-1`
- Virtual all immune/checkpoint: `CD3`, `CD4`, `CD8`, `CD20`, `CD68`, `CD11c`, `PD-1`, `PD-L1`
- Virtual proliferation: `Ki67`
- Virtual epithelial: `CK`

The goal was to ask whether the GigaTIME virtual immune/checkpoint pattern aligns with broader RNA evidence.

## Command

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_programs.py
```

Inputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/joined_slide_clinical_her2_gigatime.csv`
- `data/tcga_brca/expression_files/`

Local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/case_rna_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/case_virtual_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/joined_virtual_rna_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_rna_program_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/rna_program_group_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_program_group_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/rna_program_validation_summary.md`

Tracked figures:

- `docs/assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png`
- `docs/assets/clinical_her2_rna_program_validation/rna_programs_by_her2_group.png`
- `docs/assets/clinical_her2_rna_program_validation/virtual_programs_by_her2_group.png`

## Main Result

The broader RNA program validation still did not positively validate the GigaTIME virtual immune/checkpoint signal.

The strongest virtual-vs-RNA associations were:

| Virtual program | RNA program | Spearman rho | p | BH q |
|---|---|---:|---:|---:|
| Virtual T cell/checkpoint | Endothelial | -0.585 | 0.0007 | 0.0309 |
| Virtual all immune/checkpoint | Endothelial | -0.556 | 0.0014 | 0.0320 |
| Virtual T cell/checkpoint | B cell | -0.399 | 0.0288 | 0.3769 |
| Virtual proliferation | Endothelial | -0.375 | 0.0410 | 0.3769 |
| Virtual all immune/checkpoint | B cell | -0.361 | 0.0503 | 0.3769 |
| Virtual epithelial | Proliferation | 0.349 | 0.0590 | 0.3769 |

Only two associations were FDR-significant, and both were negative correlations with the endothelial RNA program. This is not the validation pattern we would want if the virtual immune/checkpoint outputs were cleanly reflecting bulk immune RNA programs.

![Virtual vs RNA program correlation heatmap](assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png)

## RNA Programs Across HER2 Groups

No broad RNA immune program showed an FDR-significant difference across the three clinical HER2 groups.

The most suggestive RNA group pattern was proliferation:

| RNA program | Kruskal p | BH q | Highest group | Lowest group | Max-min mean |
|---|---:|---:|---|---|---:|
| Proliferation | 0.1063 | 0.9569 | HER2-zero | HER2-low | 0.922 |
| Myeloid / macrophage | 0.4116 | 0.9987 | HER2-low | HER2-positive | 0.393 |
| B cell | 0.4317 | 0.9987 | HER2-low | HER2-zero | 0.570 |
| Dendritic / APC | 0.5697 | 0.9987 | HER2-low | HER2-zero | 0.407 |
| Checkpoint / IFNG | 0.9031 | 0.9987 | HER2-zero | HER2-positive | 0.177 |
| T cell / cytotoxic | 0.9351 | 0.9987 | HER2-low | HER2-positive | 0.281 |

This means the RNA data itself does not currently reproduce a clear HER2-zero greater than HER2-low immune-program pattern.

![RNA programs by clinical HER2 group](assets/clinical_her2_rna_program_validation/rna_programs_by_her2_group.png)

## Virtual Programs Across HER2 Groups

The GigaTIME virtual composite programs still showed the same general HER2-zero greater than HER2-low direction.

| Virtual program | Kruskal p | BH q | Highest group | Lowest group | Max-min mean |
|---|---:|---:|---|---|---:|
| Virtual myeloid/checkpoint | 0.0176 | 0.0878 | HER2-zero | HER2-low | 0.0120 |
| Virtual all immune/checkpoint | 0.0568 | 0.1292 | HER2-zero | HER2-low | 0.0235 |
| Virtual T cell/checkpoint | 0.0992 | 0.1292 | HER2-zero | HER2-low | 0.0309 |
| Virtual proliferation | 0.1033 | 0.1292 | HER2-zero | HER2-low | 0.0040 |
| Virtual epithelial | 0.2883 | 0.2883 | HER2-zero | HER2-low | 0.0518 |

The virtual myeloid/checkpoint composite was the strongest virtual program, but it still did not pass the usual 0.05 FDR threshold.

![Virtual programs by clinical HER2 group](assets/clinical_her2_rna_program_validation/virtual_programs_by_her2_group.png)

## Interpretation

This result sharpens the project story.

What became stronger:

- The virtual signal is internally consistent: the same GigaTIME immune/checkpoint channels combine into a HER2-zero greater than HER2-low virtual program.
- The 256-tile sampling result is not just a single-channel accident.

What became more concerning:

- Broader RNA immune programs did not validate the virtual immune/checkpoint signal.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal.
- RNA immune programs did not show the same HER2-zero greater than HER2-low group pattern.

This does not prove GigaTIME is wrong. Bulk RNA-seq and selected H&E image tiles can disagree because they are measured from different tissue material and at different biological scales. But it does mean we should be very careful.

## Proposal Language

A careful way to describe this step:

> We extended validation from single marker genes to broader RNA immune and tissue programs. The GigaTIME virtual myeloid/checkpoint composite retained the HER2-zero greater than HER2-low direction, but broader RNA programs did not confirm this immune pattern. The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal, reinforcing that the GigaTIME output should be treated as hypothesis-generating and prioritized for pathology review and external validation.

## Next Step

The next step after this validation check was to train a first held-out classifier baseline. That classifier is now complete and is documented in `docs/clinical_her2_classifier_baseline.md`.

The next scientific step should combine trustworthiness review with better classifier inputs:

1. Ask an advisor/pathologist to review the high virtual immune/checkpoint H&E regions.
2. Restrict the next classifier to tumor-rich tiles rather than all sampled tissue tiles.
3. Add tile distribution features and, if available, GigaTIME/pathology embeddings.
4. Add tumor purity or immune deconvolution covariates if available.
5. Check whether endothelial/stromal/tissue-composition differences might explain part of the virtual signal.
6. Search for an external dataset with paired H&E and real mIF.
7. Only then expand to a larger TCGA cohort.
