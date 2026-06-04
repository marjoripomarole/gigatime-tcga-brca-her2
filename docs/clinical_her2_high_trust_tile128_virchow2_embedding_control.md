# Generic H&E Embedding Control (Virchow2)

Status: control experiment for the HER2-low versus HER2-zero result.

## Why This Control

The primary result uses GigaTIME virtual mIF channels and is presented as immune/myeloid/checkpoint biology. This control asks whether that interpretation is necessary. Virchow2 is a generic pathology foundation-model embedding with no immune-channel meaning. If a generic embedding separates HER2-low from HER2-zero about as well as GigaTIME, and if it also collapses under source-site holdout while slide-size covariates stay strong, then the low-versus-zero axis is better described as generic morphology/tissue composition tracking TCGA acquisition structure than as GigaTIME-specific virtual immune biology.

## Method

- Cohort: 118 strict high-trust slides with Virchow2 embeddings (57 HER2-low, 61 HER2-zero).
- Embedding: `paige-ai/Virchow2`, 2560-d, mean-pooled over 128 random tissue tiles per slide (same slide list and 128-tile sampling as the GigaTIME primary run).
- Classifier: regularized logistic regression on identical folds as the GigaTIME analyses; PCA (20 components) fit inside each training fold only.
- Comparators: slide-size covariates, source-site covariates, and GigaTIME mean channels (`qc_cellular_tissue` whole-slide view).
- Validation: repeated stratified 5-fold CV (3 repeats) and leave-one-source-site-out.
- Sanity: 200 shuffled-label permutations for the embedding under repeated CV.

## Head-to-Head Results

| Feature set | Validation | Features | Balanced accuracy | AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| Slide-size covariates | Repeated stratified CV | 3 | 0.888 | 0.921 | 0.869 | 0.906 |
| Slide-size covariates | Leave source site out | 3 | 0.882 | 0.915 | 0.869 | 0.895 |
| Source-site covariates | Repeated stratified CV | 20 | 0.873 | 0.924 | 0.956 | 0.789 |
| Source-site covariates | Leave source site out | 20 | 0.500 | 0.053 | 0.000 | 1.000 |
| GigaTIME mean channels | Repeated stratified CV | 23 | 0.710 | 0.742 | 0.760 | 0.661 |
| GigaTIME mean channels | Leave source site out | 23 | 0.617 | 0.637 | 0.672 | 0.561 |
| Virchow2 embedding (PCA) | Repeated stratified CV | 2560 | 0.693 | 0.730 | 0.738 | 0.649 |
| Virchow2 embedding (PCA) | Leave source site out | 2560 | 0.551 | 0.622 | 0.541 | 0.561 |
| Virchow2 + slide-size | Repeated stratified CV | 2563 | 0.699 | 0.734 | 0.732 | 0.667 |
| Virchow2 + slide-size | Leave source site out | 2563 | 0.568 | 0.631 | 0.574 | 0.561 |

![Embedding vs GigaTIME and confound baselines](assets/clinical_her2_high_trust_tile128_virchow2_embedding_control/embedding_control_balanced_accuracy.png)

## Embedding PCA-Component Robustness (Repeated CV)

| PCA components | Balanced accuracy | AUC |
| --- | --- | --- |
| 10 | 0.693 | 0.751 |
| 20 | 0.693 | 0.730 |
| 30 | 0.659 | 0.690 |

## Shuffled-Label Sanity (Embedding, Repeated CV)

| Observed bal acc | Null mean | Null 95% | Empirical p |
| --- | --- | --- | --- |
| 0.693 | 0.489 | 0.562 | 0.0050 |

## Interpretation

- The generic Virchow2 embedding reaches balanced accuracy 0.693 under repeated stratified CV (shuffled-label null mean 0.489, empirical p 0.0050), versus 0.710 for GigaTIME mean channels and 0.888 for slide-size covariates.
- Under leave-source-site-out validation the embedding moves to 0.551 (GigaTIME 0.617, slide-size 0.882).
- Read together, a generic morphology embedding with no immune interpretation separates HER2-low from HER2-zero about as well as GigaTIME and also loses ground under source-site holdout while slide-size covariates stay strong. This supports the confound reading: the low-versus-zero axis is better explained as generic tissue/morphology that tracks TCGA acquisition structure than as GigaTIME-specific virtual immune biology.

- This is an internal control on TCGA, not external validation. It cannot prove the GigaTIME signal is biologically meaningless; it tests whether the GigaTIME-specific virtual-immune framing is required to reproduce the low-versus-zero separation.

## Output Files

- `docs/clinical_her2_high_trust_tile128_virchow2_embedding_control.md`
- `results/virchow2_tcga_brca_high_trust_tile128/embedding_low_zero_control/embedding_control_metrics.csv`
- `results/virchow2_tcga_brca_high_trust_tile128/embedding_low_zero_control/embedding_control_pca_robustness.csv`
- `results/virchow2_tcga_brca_high_trust_tile128/embedding_low_zero_control/embedding_control_permutation.csv`
- `docs/assets/clinical_her2_high_trust_tile128_virchow2_embedding_control/`
