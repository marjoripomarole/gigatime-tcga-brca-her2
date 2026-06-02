# H-Optimus Embedding Baseline

Status: ready to run once Hugging Face access is available.

## Why We Are Trying This

The current GigaTIME result is interesting but confounded. H-Optimus gives a generic H&E foundation-model embedding baseline, which asks a cleaner next question:

> Is there a broader H&E morphology signal around HER2-low versus HER2-zero that is not tied to GigaTIME virtual mIF channels?

This should be treated as a research embedding baseline, not as clinical HER2 diagnosis.

## Model Choice

The extractor supports both:

- `h0-mini`: laptop-friendlier distilled H-Optimus-family model. This is the default first run.
- `hoptimus0`: full `bioptimus/H-optimus-0` model. This is larger and likely slower/heavier.

Both models are gated on Hugging Face. Accept the model terms in the browser and set a token before inference:

```bash
export HF_TOKEN=<huggingface_read_token>
```

## Smoke Test

First check that the environment can import the runtime and find the local high-trust slide list:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py --dry-run
```

Then inspect one slide without loading the model or running inference:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --dry-run \
  --inspect-slide \
  --max-slides 1
```

Then run a tiny two-slide H0-mini smoke test:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset h0-mini \
  --out-dir results/hoptimus_tcga_brca_high_trust_tile224_smoke \
  --max-slides 2 \
  --tile-limit 16 \
  --batch-size 4
```

If that works, scale to the strict high-trust list:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset h0-mini \
  --out-dir results/hoptimus_tcga_brca_high_trust_tile224 \
  --tile-limit 64 \
  --batch-size 8 \
  --resume
```

To try full H-Optimus-0 instead:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset hoptimus0 \
  --out-dir results/hoptimus0_tcga_brca_high_trust_tile224 \
  --tile-limit 64 \
  --batch-size 4 \
  --resume
```

## Output

The main output is:

- `slide_embeddings.csv`: one row per slide with mean tile embedding features.
- `tile_embeddings.csv`: optional, only written with `--save-tile-csv`; can become large.
- `hoptimus_embedding_summary.json`: model/run settings and embedding dimension.

The default tile extraction uses 224-pixel model inputs at target 0.5 microns per pixel, matching the H-Optimus convention as closely as TCGA slide metadata allows.

## First Analysis Plan

1. Generate slide-level H-Optimus/H0-mini embeddings for the 171-slide high-trust list.
2. Join embeddings to the existing HER2 labels and source-site/slide-size covariates.
3. Repeat the key low-versus-zero checks with strict caveats:
   - ordinary cross-validation,
   - source-site held-out validation,
   - within-source-site sensitivity,
   - source-site and slide-size covariate comparison.
4. Present the result as a generic H&E embedding baseline against GigaTIME, not as a diagnostic model.
