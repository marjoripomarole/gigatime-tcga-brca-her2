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

## Why This Is a Good First Step

- It is replication-first: the model is not retrained, only applied to public TCGA-BRCA data.
- TCGA-BRCA gives paired histology, transcriptomic context, and clinical HER2 IHC/ISH fields suitable for an exploratory HER2 axis.
- The output is easy to inspect with an advisor: slide-level CSVs, channel summaries, and figures.
- The all-channel virtual mIF figures in `docs/assets/virtual_mif_channels/` show group-level channel means, slide-by-channel relative activation, and representative HER2-high/HER2-low spatial channel grids.

## Caveats to Discuss

- Clinical HER2 labels are derived from TCGA clinical supplement IHC/ISH fields, which are incomplete and must be described carefully.
- The first balanced clinical HER2 run is still small: 10 cases per group.
- The 64-tile-per-slide run is a practical pilot, not a final whole-slide sampling strategy.
- TCGA slide quality, tissue sampling, and tumor purity need QC before strong biological claims.
- GigaTIME is research-only and not a clinical HER2 classifier.
- The virtual mIF channel images are GigaTIME predictions from H&E tiles, not real multiplex immunofluorescence measurements from TCGA.
