# Expanded 20/20/20 Clinical HER2 Results

Status: Historical expanded 60-slide presentation summary. For the current primary result, use `docs/clinical_her2_high_trust_tile128_results.md`; use this report as a 256-tile comparator/provenance record.

## What Was Run

We expanded the TCGA-BRCA clinical HER2 GigaTIME pilot from 30 slides to 60 slides:

| Group | Slides/cases |
|---|---:|
| HER2-positive | 20 |
| HER2-low | 20 |
| HER2-zero | 20 |

The expanded run used the same denser sampling strategy as the robustness run:

- Up to 256 tissue tiles per slide.
- 60 total slides.
- 15,225 total tile predictions.
- Median 256 tiles per slide, range 176-256.
- Matched STAR-count RNA-seq expression was downloaded for all 60 selected cases.

## Strongest Result To Present

The strongest presentation result is that the HER2-low versus HER2-zero separation became stronger after expansion.

In the original 30-slide pilot, the leading HER2-low versus HER2-zero pairwise q values improved but did not pass 0.05. In the expanded 60-slide run, several HER2-low versus HER2-zero virtual-channel differences pass BH correction within the all-tissue and QC-cellular cleanup views.

### HER2-Low Versus HER2-Zero Pairwise Signal

Negative delta means HER2-low has lower mean virtual activation than HER2-zero.

| Cleanup view | Channel | HER2-low minus HER2-zero | Mann-Whitney p | BH q within view |
|---|---:|---:|---:|---:|
| All sampled tissue | CD4 | -0.0372 | 0.0016 | 0.0252 |
| All sampled tissue | CD11c | -0.0051 | 0.0026 | 0.0252 |
| All sampled tissue | CD3 | -0.0371 | 0.0028 | 0.0252 |
| All sampled tissue | CD68 | -0.0070 | 0.0051 | 0.0326 |
| QC cellular tissue | CD4 | -0.0448 | 0.0012 | 0.0200 |
| QC cellular tissue | CD3 | -0.0449 | 0.0015 | 0.0200 |
| QC cellular tissue | CD11c | -0.0062 | 0.0023 | 0.0206 |
| QC cellular tissue | CD68 | -0.0079 | 0.0066 | 0.0320 |
| QC cellular tissue | PD-L1 | -0.0180 | 0.0071 | 0.0320 |

This is stronger than the earlier 30-slide result because the expanded run now has FDR-passing HER2-low versus HER2-zero differences for multiple virtual immune/myeloid/checkpoint-associated channels.

## Classifier Result

The HER2-low versus HER2-zero classifier signal also persisted in the expanded 40-case binary comparison.

| Input view | Best feature set | N | Balanced accuracy | Macro AUC | Sensitivity | Specificity |
|---|---|---:|---:|---:|---:|---:|
| All sampled tissue | Mean + fraction channels | 40 | 0.800 | 0.820 | 0.950 | 0.650 |
| QC cellular tissue | Mean + fraction channels | 40 | 0.775 | 0.820 | 0.900 | 0.650 |
| CK-enriched top 50% | Mean channels | 40 | 0.750 | 0.807 | 0.850 | 0.650 |
| CK-enriched top 25% | Mean + fraction channels | 40 | 0.800 | 0.820 | 0.850 | 0.750 |

This is useful because the result did not disappear when the cohort doubled from 10/10/10 to 20/20/20. It is still not a diagnostic model, but it is a stronger feasibility signal.

## Three-Group Pattern

The expanded run adds an important nuance.

In the 30-slide run, the clearest group-average story was HER2-zero greater than HER2-low for `CD68`, `PD-L1`, and `CD11c`.

In the 60-slide run, HER2-low still tends to be the lowest group, but HER2-positive often becomes the highest group for broader immune/checkpoint channels:

| Virtual program | Kruskal p | BH q | Highest group | Lowest group |
|---|---:|---:|---|---|
| Virtual all immune/checkpoint | 0.0256 | 0.0488 | HER2-positive | HER2-low |
| Virtual T cell/checkpoint | 0.0284 | 0.0488 | HER2-positive | HER2-low |
| Virtual myeloid/checkpoint | 0.0293 | 0.0488 | HER2-positive | HER2-low |
| Virtual epithelial | 0.0401 | 0.0501 | HER2-zero | HER2-low |

So the cleaner biological statement is:

> The expanded run supports an image-derived HER2-low versus HER2-zero difference, with HER2-low often showing lower virtual immune/checkpoint and myeloid-associated activation. The three-group pattern is not simply HER2-zero highest for every channel; HER2-positive becomes high for several broader virtual immune programs.

## What Did Not Validate Strongly

The RNA validation is still weak.

Marker-level RNA validation:

- 60 cases had paired GigaTIME and RNA-seq data.
- No virtual channel had an FDR-significant positive correlation with its matching bulk RNA marker signature.
- The strongest positive correlations were very weak, for example `CD11c` rho 0.058 and `CD68` rho 0.050.

Broader RNA program validation:

- No RNA immune program showed an FDR-significant HER2-group difference.
- The strongest virtual-vs-RNA program associations did not pass FDR correction.

This means the virtual signal is stronger statistically inside the GigaTIME/image feature space, but it is still not molecularly validated by bulk RNA-seq.

## HER2-Positive Classification

GigaTIME/H&E still does not reliably classify HER2-positive disease.

Best GigaTIME-only HER2-positive versus HER2-negative result:

- Balanced accuracy: 0.575 with all sampled tissue.
- Balanced accuracy: 0.613 after QC-cellular cleanup.
- Sensitivity remained low, around 0.25.

ERBB2 RNA did better as expected:

- HER2-positive versus HER2-negative balanced accuracy: 0.750.
- Specificity: 1.000.

This is important because it shows the clinical labels contain molecular signal, but the current image-derived features are still weak for identifying HER2-positive disease.

## Best Presentation Framing

The strongest presentation claim is:

> After expanding the pilot to 60 TCGA-BRCA slides, GigaTIME-derived H&E features continued to separate HER2-low from HER2-zero tumors. Several virtual immune/myeloid/checkpoint channels now show FDR-corrected HER2-low versus HER2-zero differences, and a cross-validated HER2-low versus HER2-zero classifier maintained balanced accuracy around 0.80. However, RNA validation remains weak and GigaTIME/H&E does not yet reliably classify HER2-positive disease, so these results should be presented as hypothesis-generating image-derived HER2-state associations, not as clinical HER2 diagnosis or direct isoform detection.

## Outputs

- Expanded cohort selection: `docs/clinical_her2_cohort_expanded20_selection.md`
- Expanded GigaTIME run: `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/`
- Clinical summary: `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/clinical_summary/clinical_her2_summary.md`
- RNA marker validation: `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/rna_validation/rna_validation_summary.md`
- RNA program validation: `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/rna_program_validation/rna_program_validation_summary.md`
- Cleanup summary: `docs/clinical_her2_expanded20_gigatime_data_cleanup.md`
- Cleaned classifier comparison: `docs/clinical_her2_expanded20_cleaned_classifier_comparison.md`
