# Cleaned GigaTIME HER2 Classifier Comparison

This analysis reruns the slide-level HER2 classifier after cleaning the GigaTIME tile inputs. It compares all sampled tissue against cellular and virtual CK-enriched feature views.

Every prediction is leave-one-out cross-validated. This remains a small 60-case pilot, not a clinical model.

## Feature Views

- All sampled tissue: the original 256-tile slide aggregation.
- QC cellular tissue: tissue fraction >= 0.70 and virtual DAPI mean >= 0.05.
- CK-enriched top 50%: top half of virtual CK tiles within each slide after QC.
- CK-enriched top 25%: top quarter of virtual CK tiles within each slide after QC.

Virtual CK and DAPI are GigaTIME predictions, not real stains or pathologist tumor annotations.

## Feature Sets

| Feature set | Number of features |
| --- | --- |
| Mean channels | 23 |
| Mean + fraction channels | 46 |
| Interpretable means | 10 |
| Interpretable distribution features | 50 |
| Virtual programs | 8 |
| ERBB2 RNA reference | 1 |

## Best GigaTIME/H&E Result Per View and Task

| Cleanup view | Task | Best feature set | N | Accuracy | Balanced accuracy | Macro AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| All sampled tissue | HER2-low vs HER2-zero | Mean + fraction channels | 40 | 0.800 | 0.800 | 0.820 | 0.950 | 0.650 |
| All sampled tissue | HER2-positive vs HER2-negative | Mean + fraction channels | 60 | 0.683 | 0.575 | 0.537 | 0.250 | 0.900 |
| All sampled tissue | HER2-positive vs HER2-low vs HER2-zero | Interpretable distribution features | 60 | 0.567 | 0.567 | 0.643 |  |  |
| QC cellular tissue | HER2-low vs HER2-zero | Mean + fraction channels | 40 | 0.775 | 0.775 | 0.820 | 0.900 | 0.650 |
| QC cellular tissue | HER2-positive vs HER2-negative | Interpretable means | 60 | 0.733 | 0.613 | 0.281 | 0.250 | 0.975 |
| QC cellular tissue | HER2-positive vs HER2-low vs HER2-zero | Interpretable means | 60 | 0.550 | 0.550 | 0.635 |  |  |
| CK-enriched top 50% | HER2-low vs HER2-zero | Mean channels | 40 | 0.750 | 0.750 | 0.807 | 0.850 | 0.650 |
| CK-enriched top 50% | HER2-positive vs HER2-negative | Interpretable distribution features | 60 | 0.667 | 0.587 | 0.631 | 0.350 | 0.825 |
| CK-enriched top 50% | HER2-positive vs HER2-low vs HER2-zero | Mean + fraction channels | 60 | 0.567 | 0.567 | 0.741 |  |  |
| CK-enriched top 25% | HER2-low vs HER2-zero | Mean + fraction channels | 40 | 0.800 | 0.800 | 0.820 | 0.850 | 0.750 |
| CK-enriched top 25% | HER2-positive vs HER2-negative | Mean + fraction channels | 60 | 0.700 | 0.600 | 0.659 | 0.300 | 0.900 |
| CK-enriched top 25% | HER2-positive vs HER2-low vs HER2-zero | Mean channels | 60 | 0.567 | 0.567 | 0.726 |  |  |

![Best classifier by cleanup view](assets/clinical_her2_expanded20_cleaned_classifier/cleaned_classifier_best_by_view.png)

## HER2-Low Versus HER2-Zero Focus

| Cleanup view | Best feature set | Accuracy | Balanced accuracy | Macro AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| All sampled tissue | Mean + fraction channels | 0.800 | 0.800 | 0.820 | 0.950 | 0.650 |
| QC cellular tissue | Mean + fraction channels | 0.775 | 0.775 | 0.820 | 0.900 | 0.650 |
| CK-enriched top 50% | Mean channels | 0.750 | 0.750 | 0.807 | 0.850 | 0.650 |
| CK-enriched top 25% | Mean + fraction channels | 0.800 | 0.800 | 0.820 | 0.850 | 0.750 |

![HER2-low versus HER2-zero feature-set comparison](assets/clinical_her2_expanded20_cleaned_classifier/cleaned_classifier_low_zero_feature_sets.png)

![HER2-low versus HER2-zero confusion matrices](assets/clinical_her2_expanded20_cleaned_classifier/cleaned_classifier_low_zero_confusions.png)

## Main Result

- All sampled tissue HER2-low versus HER2-zero balanced accuracy: 0.800, macro AUC: 0.820.
- QC cellular tissue preserved the HER2-low versus HER2-zero result: balanced accuracy 0.775, macro AUC 0.820.
- CK-enriched top 50% reduced HER2-low versus HER2-zero performance to balanced accuracy 0.750.
- CK-enriched top 25% HER2-low versus HER2-zero balanced accuracy was 0.800.
- CK-enriched top 25% modestly improved HER2-positive versus HER2-negative balanced accuracy to 0.600, but sensitivity remained low at 0.300.

## ERBB2 RNA Reference

ERBB2 RNA is repeated as a non-H&E reference. It does not depend on the cleanup view and should not be interpreted as image-derived performance.

| Cleanup view | Task | Accuracy | Balanced accuracy | Macro AUC |
| --- | --- | --- | --- | --- |
| All sampled tissue | HER2-low vs HER2-zero | 0.500 | 0.500 | 0.430 |
| All sampled tissue | HER2-positive vs HER2-negative | 0.833 | 0.750 | 0.806 |
| All sampled tissue | HER2-positive vs HER2-low vs HER2-zero | 0.200 | 0.200 | 0.571 |
| QC cellular tissue | HER2-low vs HER2-zero | 0.500 | 0.500 | 0.430 |
| QC cellular tissue | HER2-positive vs HER2-negative | 0.833 | 0.750 | 0.806 |
| QC cellular tissue | HER2-positive vs HER2-low vs HER2-zero | 0.200 | 0.200 | 0.571 |
| CK-enriched top 50% | HER2-low vs HER2-zero | 0.500 | 0.500 | 0.430 |
| CK-enriched top 50% | HER2-positive vs HER2-negative | 0.833 | 0.750 | 0.806 |
| CK-enriched top 50% | HER2-positive vs HER2-low vs HER2-zero | 0.200 | 0.200 | 0.571 |
| CK-enriched top 25% | HER2-low vs HER2-zero | 0.500 | 0.500 | 0.430 |
| CK-enriched top 25% | HER2-positive vs HER2-negative | 0.833 | 0.750 | 0.806 |
| CK-enriched top 25% | HER2-positive vs HER2-low vs HER2-zero | 0.200 | 0.200 | 0.571 |

## Interpretation

The cleaned-view comparison asks whether the classifier signal is stronger in tumor-enriched tile views or in broader tissue context.

In this run, cellular-tissue cleanup preserved the HER2-low versus HER2-zero classifier signal, which argues against the result being only blank/background artifact. However, the signal weakened when restricted to the most CK-enriched tiles. That suggests the current GigaTIME signal may be capturing broader tissue or microenvironment context more than a purely epithelial tumor-cell HER2 phenotype.

The practical next step is not to claim diagnosis. It is to inspect the cases that change classification between all-tissue/QC-cellular and CK-enriched views, because those flips can reveal whether GigaTIME is using tumor regions, stromal context, immune infiltrates, or tile-selection artifacts.

## Outputs

- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/cleaned_classifier_comparison/cleaned_classifier_predictions.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/cleaned_classifier_comparison/cleaned_classifier_metrics.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/cleaned_classifier_comparison/cleaned_classifier_confusion_matrices.csv`
- `docs/assets/clinical_her2_expanded20_cleaned_classifier/`
