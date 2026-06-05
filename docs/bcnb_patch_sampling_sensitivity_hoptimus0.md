# BCNB H-Optimus-0 Patch Sampling Sensitivity

Status: multi-seed patch-selection sensitivity check for the BCNB external patch pilot.

## Method

- Cohort: 781 BCNB low/zero patients (654 HER2-low, 127 HER2-zero).
- Model: H-Optimus-0 patient-level mean embeddings from 10 precomputed 256x256 patches per patient.
- Comparison: older lexicographic capped sample versus two deterministic hash-capped samples (`--sampling-method hash`, seeds `20260604` and `20260605`).
- Classifier: same class-balanced regularized logistic regression, repeated stratified 5-fold CV with 5 repeats.
- Default embedding dimensionality: PCA20 fit inside each training fold only.
- Sanity: 200 shuffled-label permutations for each embedding run.
- The two hash seeds are a meaningful perturbation: across all 1,058 BCNB patients, the 20260605 sample shares a mean 3.09/10 patches per patient with the 20260604 sample; only 19 patients have identical 10-patch sets.

## Results

| Feature set | Patch sample | PCA | Balanced accuracy | AUC | Shuffled-label p |
|---|---|---:|---:|---:|---:|
| H-Optimus-0 embedding | Lexicographic capped10 | 20 | 0.575 | 0.633 | 0.005 |
| H-Optimus-0 embedding | Hash capped10 seed 20260604 | 20 | 0.597 | 0.640 | 0.005 |
| H-Optimus-0 embedding | Hash capped10 seed 20260605 | 20 | 0.585 | 0.634 | 0.005 |
| H-Optimus-0 + clinical covariates | Lexicographic capped10 | 20 | 0.580 | 0.636 |  |
| H-Optimus-0 + clinical covariates | Hash capped10 seed 20260604 | 20 | 0.595 | 0.641 |  |
| H-Optimus-0 + clinical covariates | Hash capped10 seed 20260605 | 20 | 0.582 | 0.635 |  |
| Clinical covariates | Same patients |  | 0.643 | 0.627 |  |

PCA-grid sensitivity for the image-only H-Optimus-0 embedding:

| Patch sample | PCA5 BA/AUC | PCA10 BA/AUC | PCA20 BA/AUC | PCA30 BA/AUC | PCA50 BA/AUC |
|---|---|---|---|---|---|
| Lexicographic capped10 | 0.531 / 0.556 | 0.588 / 0.618 | 0.575 / 0.633 | 0.643 / 0.686 | 0.630 / 0.665 |
| Hash capped10 seed 20260604 | 0.573 / 0.590 | 0.590 / 0.611 | 0.597 / 0.640 | 0.626 / 0.680 | 0.610 / 0.663 |
| Hash capped10 seed 20260605 | 0.542 / 0.562 | 0.573 / 0.598 | 0.585 / 0.634 | 0.622 / 0.666 | 0.621 / 0.661 |

## Interpretation

- The predeclared PCA20 effect is somewhat patch-sample-sensitive by balanced accuracy (0.575 to 0.597), but AUC is stable (0.633 to 0.640), and all three runs beat shuffled-label nulls.
- This supports the same conclusion as the paired H-Optimus-0/Virchow2 comparison: BCNB contains a weak external low/zero-associated morphology signal, not a strong standalone image classifier.
- The PCA grid shows that choosing components after looking at the data could inflate the apparent balanced accuracy. Treat the grid as sensitivity evidence only; the default PCA20 result remains the conservative headline.
- The second independent hash seed makes it less likely that the result is a one-sample patch artifact. A manuscript-grade sampling check would still be stronger with more seeds, more patches per patient, or full WSI processing.

## Commands

```bash
conda run -n gigatime-tcga python scripts/analyze_bcnb_patch_embedding_control.py \
  --embeddings results/bcnb_patch_embeddings_hoptimus0_capped10_low_zero/patient_embeddings.csv \
  --model-label H-Optimus-0 \
  --model-id bioptimus/H-optimus-0 \
  --out-dir results/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero \
  --asset-dir docs/assets/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero \
  --out-markdown docs/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero.md \
  --folds 5 --repeats 5 --permutations 200 \
  --pca-components 20 --pca-grid 5,10,20,30,50
```

Preferred hash-sampled comparator:

```bash
conda run -n gigatime-tcga python scripts/analyze_bcnb_patch_embedding_control.py \
  --embeddings results/bcnb_patch_embeddings_hoptimus0_hash_capped10_low_zero/patient_embeddings.csv \
  --model-label H-Optimus-0 \
  --model-id bioptimus/H-optimus-0 \
  --out-dir results/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero \
  --asset-dir docs/assets/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero \
  --out-markdown docs/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero.md \
  --folds 5 --repeats 5 --permutations 200 \
  --pca-components 20 --pca-grid 5,10,20,30,50
```

Second independent hash-sampled comparator:

```bash
conda run -n gigatime-tcga python scripts/build_bcnb_patch_manifest.py \
  --max-patches-per-patient 10 \
  --sampling-method hash \
  --sampling-seed 20260605 \
  --output data/bcnb/bcnb_patch_manifest_hash20260605_capped10.csv

conda run -n gigatime-tcga python scripts/run_bcnb_patch_embeddings.py \
  --model hoptimus0 \
  --patch-manifest data/bcnb/bcnb_patch_manifest_hash20260605_capped10.csv \
  --groups HER2-zero,HER2-low \
  --max-patches-per-patient 10 \
  --out-dir results/bcnb_patch_embeddings_hoptimus0_hash20260605_capped10_low_zero \
  --batch-size 16 \
  --resume

conda run -n gigatime-tcga python scripts/analyze_bcnb_patch_embedding_control.py \
  --embeddings results/bcnb_patch_embeddings_hoptimus0_hash20260605_capped10_low_zero/patient_embeddings.csv \
  --model-label H-Optimus-0 \
  --model-id bioptimus/H-optimus-0 \
  --out-dir results/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero \
  --asset-dir docs/assets/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero \
  --out-markdown docs/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero.md \
  --folds 5 --repeats 5 --permutations 200 \
  --pca-components 20 --pca-grid 5,10,20,30,50
```

## Output Files

- `docs/bcnb_patch_sampling_sensitivity_hoptimus0.md`
- `docs/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero.md`
- `results/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero/`
- `docs/assets/bcnb_patch_embedding_control_hoptimus0_lexicographic_capped10_low_zero/`
- `docs/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero.md`
- `results/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero/`
- `docs/assets/bcnb_patch_embedding_control_hoptimus0_hash20260605_capped10_low_zero/`
