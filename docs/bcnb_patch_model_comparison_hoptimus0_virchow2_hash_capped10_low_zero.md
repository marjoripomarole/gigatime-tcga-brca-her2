# BCNB Patch Model Comparison: H-Optimus-0 Versus Virchow2

Status: paired external-cohort comparison using the same BCNB patients and the same deterministic patch sample.

## Method

- Cohort: 781 BCNB patients (654 HER2-low, 127 HER2-zero).
- Inputs: patient-level means of the same 10 deterministic hash-sampled 256x256 tumor-region patches per patient.
- Models: H-Optimus-0 (1536-d) and Virchow2 (2560-d).
- Classifier: class-balanced regularized logistic regression with repeated stratified 5-fold CV (5 repeats).
- Dimensionality control: PCA is fit inside each training fold only; dual-model rows use 20 components per embedding model.
- Sanity: 200 shuffled-label permutations for the dual-model embedding classifier.

## Results

| Feature set | Features after PCA | Balanced accuracy | AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- |
| Clinical covariates | 13 | 0.643 | 0.627 | 0.532 | 0.753 |
| H-Optimus-0 embedding | 20 | 0.597 | 0.640 | 0.539 | 0.655 |
| Virchow2 embedding | 20 | 0.600 | 0.643 | 0.545 | 0.654 |
| H-Optimus-0 + Virchow2 | 40 | 0.609 | 0.651 | 0.534 | 0.684 |
| H-Optimus-0 + Virchow2 + clinical | 53 | 0.615 | 0.661 | 0.532 | 0.698 |
| Average probability ensemble | 2 | 0.610 | 0.650 | 0.550 | 0.669 |

![BCNB patch model comparison](assets/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_comparison_metrics.png)

## Cross-Model Agreement

| Agreement measure | Value |
| --- | --- |
| OOF probability Pearson r | 0.787 |
| OOF probability Spearman rho | 0.778 |
| Patient-mean probability Pearson r | 0.804 |
| Patient-mean probability Spearman rho | 0.793 |
| Patient threshold disagreement fraction | 0.186 |
| Patient mean absolute probability difference | 0.086 |

![BCNB patch model probability agreement](assets/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_probability_agreement.png)

## Dual-Model Shuffled-Label Sanity

| Metric | Observed | Null mean | Null 95% | Empirical p |
| --- | --- | --- | --- | --- |
| Balanced accuracy | 0.609 | 0.498 | 0.536 | 0.0050 |
| AUC | 0.651 | 0.496 | 0.548 | 0.0050 |

## Interpretation

- H-Optimus-0 and Virchow2 produce concordant but not identical patient scores on the same BCNB patches.
- The dual-model embedding does not create a large jump over either model alone, which argues against a hidden strong signal missed by one encoder.
- Clinical covariates remain at least as strong by balanced accuracy, so the current BCNB result should be framed as a weak, reproducible morphology/covariate-associated signal rather than a clinically deployable HER2-low versus zero classifier.
- The next escalation, if needed for a manuscript, is broader patch-sampling sensitivity, more patches per patient, or full-WSI processing to test whether broader tissue context changes the effect size.

## Output Files

- `docs/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero.md`
- `results/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_comparison_metrics.csv`
- `results/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_oof_predictions.csv`
- `results/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_patient_predictions.csv`
- `results/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/bcnb_patch_model_agreement.json`
- `docs/assets/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero/`
