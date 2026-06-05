# BCNB H0-mini Access Probe

Status: blocked pending Hugging Face gated-model access or token visibility.

## Probe

Attempted a minimal BCNB smoke with H0-mini: two patients maximum, one patch per patient, using the audited hash-capped BCNB patch manifest. This was an access/tensor-shape probe only, not a biological run.

```bash
conda run -n gigatime-tcga python scripts/run_bcnb_patch_embeddings.py \
  --model h0-mini \
  --patch-manifest data/bcnb/bcnb_patch_manifest_hash_capped10.csv \
  --groups HER2-zero,HER2-low \
  --max-patients-per-group 1 \
  --max-patches-per-patient 1 \
  --out-dir results/bcnb_patch_embeddings_h0mini_smoke2 \
  --batch-size 2 \
  --save-patch-csv \
  --resume
```

## Result

The run failed before any BCNB patch inference:

```text
HF_TOKEN/HUGGING_FACE_HUB_TOKEN is not set; gated H-Optimus download may fail.
Using device: mps
Could not load hf-hub:bioptimus/H0-mini. H-Optimus models are gated on Hugging Face; accept the model terms and set HF_TOKEN/HUGGING_FACE_HUB_TOKEN before running.
```

## Interpretation

- H0-mini is not currently a runnable BCNB model avenue in this workspace.
- This is an access/token gate, not evidence about model performance.
- Re-run only after Hugging Face access to `bioptimus/H0-mini` is approved and the token is visible inside the `gigatime-tcga` conda execution environment.
- H-Optimus-0 remains runnable in this workspace because the existing model access/cache path succeeded for the BCNB pilot.
