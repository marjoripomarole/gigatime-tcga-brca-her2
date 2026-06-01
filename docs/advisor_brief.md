# Advisor Brief: TCGA-BRCA GigaTIME HER2 Pilot

## Working Question

Can the released GigaTIME model generate biologically interpretable virtual tumor immune microenvironment features from TCGA-BRCA H&E slides, and do those features vary across clinically defined HER2-positive, HER2-low, and HER2-zero breast cancers?

## Concrete Pilot

1. Query TCGA-BRCA diagnostic H&E slides, RNA-seq STAR-count files, and clinical HER2 supplement fields from GDC.
2. Assign clinical HER2-positive, HER2-low, HER2-zero, or unknown labels using IHC/ISH fields.
3. Select a balanced 30-case pilot: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero.
4. Run the official GigaTIME model on diagnostic slide tiles.
5. Aggregate virtual mIF activations per slide for the GigaTIME channels.
6. Compare virtual TIME markers across clinical HER2 groups and render virtual mIF figures for review.

## Current Full-Pilot Finding

The first balanced clinical HER2 pilot processed all 30 selected slides using 64 random tissue tiles per slide.

The strongest pilot signal was not HER2-positive versus HER2-low. Instead, HER2-zero showed higher mean GigaTIME-predicted immune/checkpoint signals than HER2-low, especially:

- `CD68`
- `PD-L1`
- `CD11c`
- `CD4`
- `Ki67`

The top unadjusted three-group tests were CD68, PD-L1, and CD11c. Pairwise HER2-low versus HER2-zero tests were the strongest, but none remained significant after Benjamini-Hochberg correction. This should be framed as a hypothesis-generating signal.

## First Validation Check

We compared the GigaTIME virtual channels with matched TCGA RNA-seq marker signatures for the same 30 cases.

This did not strongly confirm the virtual immune-channel signal:

- `Ki67` had the strongest positive correlation with its RNA signature, but it was weak and not FDR-significant.
- `CD68`, `PD-L1`, and `CD11c`, the main virtual channels driving the HER2-zero versus HER2-low signal, did not show strong positive RNA-signature correlations.
- This means the current result should remain hypothesis-generating until visual QC, more tile sampling, and additional validation are done.

## Visual QC Update

We rendered H&E-versus-virtual-mIF QC panels for the top `CD68` + `PD-L1` + `CD11c` case in each HER2 group.

The high-scoring tiles were real tissue-containing, cellular H&E regions rather than obvious blank background. That supports continuing the analysis. However, the high-signal tiles were not visually unique to HER2-zero; the selected HER2-positive case also had strong high-signal tiles. The result remains a slide-level pilot trend, not a clean single-case visual phenotype.

## 256-Tile Robustness Update

We reran the same 30 selected clinical HER2 slides with up to 256 random tissue tiles per slide. This was done to test whether the original 64-tile signal was mainly a sparse-sampling artifact.

The main result persisted:

| Channel | 64-tile p | 256-tile p | 64 max-min | 256 max-min | 256 direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |

Pairwise HER2-low versus HER2-zero q values improved for CD68, PD-L1, and CD11c to 0.1133, but they still did not meet the usual 0.05 FDR threshold.

The 256-tile RNA validation remained weak. No virtual channel had an FDR-significant correlation with its matched RNA marker signature. Therefore, the current interpretation is stronger sampling robustness, not biological validation.

## RNA Program Validation Update

We then tested broader RNA immune and tissue programs rather than only single marker-channel signatures. This compared GigaTIME virtual composite programs with RNA programs for T-cell/cytotoxic, checkpoint/IFNG, myeloid/macrophage, dendritic/APC, B-cell, proliferation, epithelial, stromal, and endothelial biology.

The broader validation still did not positively confirm the virtual immune/checkpoint signal.

Key findings:

- The virtual myeloid/checkpoint composite retained the HER2-zero greater than HER2-low direction, but was not FDR-significant: Kruskal p 0.0176, BH q 0.0878.
- No broad RNA immune program showed an FDR-significant HER2-group difference.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal:
  - Virtual T cell/checkpoint versus endothelial RNA: Spearman rho -0.585, BH q 0.0309.
  - Virtual all immune/checkpoint versus endothelial RNA: Spearman rho -0.556, BH q 0.0320.

This is a cautionary result. It suggests the virtual signal is reproducible inside GigaTIME, but not yet validated against orthogonal RNA evidence. It also raises the possibility that tissue composition, stromal/endothelial context, or slide sampling may be influencing the predictions.

## First Classifier Baseline

We then moved from group-average comparisons to a first diagnostic-model style classifier. The classifier used slide-level GigaTIME features from the 256-tile run and leave-one-out cross-validation.

Three tasks were tested:

- HER2-positive versus HER2-negative.
- HER2-low versus HER2-zero.
- Full three-class HER2-positive versus HER2-low versus HER2-zero.

Main result:

| Task | Best GigaTIME/H&E feature set | Balanced accuracy | Macro AUC |
|---|---|---:|---:|
| HER2-low vs HER2-zero | GigaTIME mean + fraction channels | 0.800 | 0.870 |
| HER2-positive vs HER2-negative | GigaTIME mean + fraction channels | 0.475 | 0.430 |
| Three-class HER2 group | GigaTIME mean + fraction channels | 0.333 | 0.555 |

Interpretation:

- The HER2-low versus HER2-zero result is promising but very small-sample and potentially unstable.
- GigaTIME/H&E features do not currently classify HER2-positive status reliably.
- Full three-class prediction is at chance.
- ERBB2 RNA, included as a non-H&E reference, classified HER2-positive versus HER2-negative better than GigaTIME/H&E features. This means the labels contain molecular signal, but the current image-derived features are not capturing the clinical HER2-positive signal reliably.

## Why This Is a Good First Step

- It is replication-first: the model is not retrained, only applied to public TCGA-BRCA data.
- TCGA-BRCA gives paired histology, transcriptomic context, and clinical HER2 IHC/ISH fields suitable for an exploratory HER2 axis.
- The output is easy to inspect with an advisor: slide-level CSVs, channel summaries, and figures.
- The all-channel virtual mIF figures in `docs/assets/virtual_mif_channels/` show group-level channel means, slide-by-channel relative activation, and representative HER2-high/HER2-low spatial channel grids.

## Caveats to Discuss

- Clinical HER2 labels are derived from TCGA clinical supplement IHC/ISH fields, which are incomplete and must be described carefully.
- The first balanced clinical HER2 run is still small: 10 cases per group.
- The 64-tile-per-slide run is a practical pilot, not a final whole-slide sampling strategy.
- The 256-tile rerun supports robustness to denser sampling, but it is still not exhaustive whole-slide analysis.
- Bulk RNA-seq is an indirect validation layer and did not strongly validate the current GigaTIME immune-channel pattern, even with broader RNA programs.
- The first classifier baseline is not clinically usable; it is a feasibility and failure-mode analysis.
- Visual QC supports that the signal is not just blank background, but it is not biological validation.
- TCGA slide quality, tissue sampling, and tumor purity need QC before strong biological claims.
- GigaTIME is research-only and not a clinical HER2 classifier.
- The virtual mIF channel images are GigaTIME predictions from H&E tiles, not real multiplex immunofluorescence measurements from TCGA.
