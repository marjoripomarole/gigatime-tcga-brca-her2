# Advisor Brief: TCGA-BRCA GigaTIME HER2 Pilot

## Working Question

Can the released GigaTIME model generate biologically interpretable virtual tumor immune microenvironment features from TCGA-BRCA H&E slides, and do those features vary with HER2/ERBB2 expression?

## Concrete Pilot

1. Query TCGA-BRCA diagnostic H&E slides and matched open RNA-seq STAR-count files from GDC.
2. Extract ERBB2 TPM from each TCGA-BRCA RNA-seq file.
3. Run the official GigaTIME model on diagnostic slide tiles.
4. Aggregate virtual mIF activations per slide for the 23 GigaTIME channels.
5. Compare virtual TIME markers across ERBB2-high and ERBB2-low strata.
6. Render all-channel virtual mIF figures for visual inspection and advisor discussion.

## Why This Is a Good First Step

- It is replication-first: the model is not retrained, only applied to public TCGA-BRCA data.
- TCGA-BRCA gives paired histology and transcriptomic context suitable for an exploratory HER2 axis.
- The output is easy to inspect with an advisor: slide-level CSVs, channel summaries, and figures.
- The all-channel virtual mIF figures in `docs/assets/virtual_mif_channels/` show group-level channel means, slide-by-channel relative activation, and representative HER2-high/HER2-low spatial channel grids.

## Caveats to Discuss

- HER2 status is proxied by ERBB2 RNA expression unless clinical IHC/FISH annotations are added.
- TCGA slide quality, tissue sampling, and tumor purity need QC before strong biological claims.
- GigaTIME is research-only and not a clinical HER2 classifier.
- The virtual mIF channel images are GigaTIME predictions from H&E tiles, not real multiplex immunofluorescence measurements from TCGA.
- The first run should be a small reproducible pilot before scaling to all TCGA-BRCA diagnostic slides.
