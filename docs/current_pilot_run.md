# Current Pilot Run

## Run Status

The current workspace contains a completed TCGA-BRCA GigaTIME pilot on ERBB2 extreme-expression cases:

- HER2-high cases processed: 7
- HER2-low cases processed: 5
- Slides processed: 12
- GigaTIME tiles per slide: 64 random tissue tiles
- Total tile predictions: 768
- Device used: CPU

The HER2 group labels come from `data/tcga_brca/her2_extreme_cases.csv`, which selects the top 20 and bottom 20 ERBB2 TPM cases available in the current 80-case TCGA-BRCA metadata/expression pilot. GDC slide downloads were slow and repeatedly dropped connections, so 12 of the 40 selected slides are currently processed.

## Main Outputs

- `data/tcga_brca/erbb2_expression.csv`
- `data/tcga_brca/her2_extreme_cases.csv`
- `results/gigatime_tcga_brca_extremes/slide_scores.csv`
- `results/gigatime_tcga_brca_extremes/tile_scores.csv`
- `results/gigatime_tcga_brca_extremes/heatmaps/`
- `results/gigatime_tcga_brca_extremes/advisor_summary/advisor_summary.md`
- `results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/*.png`
- `docs/assets/virtual_mif_channels/virtual_mif_all_channel_group_means.png`
- `docs/assets/virtual_mif_channels/virtual_mif_slide_channel_matrix.png`
- `docs/assets/virtual_mif_channels/her2_high_reference_all_virtual_mif_channels.png`
- `docs/assets/virtual_mif_channels/her2_low_reference_all_virtual_mif_channels.png`
- `docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_virtual_mif_montage.png`
- `docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_virtual_mif_montage.png`

## Virtual mIF Channel Figures

The virtual mIF figures are generated from the GigaTIME tile and slide prediction tables with:

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_channel_images.py
```

These figures show all 23 GigaTIME channels. They include a group-level HER2-high versus HER2-low channel comparison, a slide-by-channel matrix, and representative spatial channel grids for one ERBB2-high and one ERBB2-low slide. See `docs/virtual_mif_channel_outputs.md` for captions and interpretation guardrails.

The fluorescence-style virtual mIF composites can be generated with:

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_composites.py
```

Those composites look closer to real mIF panels because they combine full predicted channel maps into colored black-background images. They are still model predictions from H&E, not real experimental mIF.

## Caveat

This is still a pilot, not a definitive biological result. It is stronger than the first two-case proof-of-work because it uses explicit ERBB2 top/bottom groups and case-level statistics, but the next scientific step is to finish downloading the remaining selected slides, rerun the same commands on the full 20 high / 20 low cohort, and add stronger tissue QC.
