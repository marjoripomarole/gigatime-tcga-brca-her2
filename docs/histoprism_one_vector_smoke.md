# HistoPrism one-vector smoke test

This smoke test checks that the published HistoPrism checkpoint can be loaded
and can produce gene-level predictions from one HistoPrism-format UNI feature
vector.

It is not yet a TCGA raw-tile run. HistoPrism split 0 expects 1024-dimensional
UNI features, and `MahmoodLab/UNI` is gated for the current Hugging Face account.
Until UNI access is available, this script uses one official precomputed UNI
feature vector bundled with the upstream HistoPrism repository.

## Local setup

Clone the upstream repository locally:

```bash
git clone https://github.com/susuhu/HistoPrism.git external/HistoPrism
```

The local clone is ignored by this repository, matching the existing
`external/GigaTIME` pattern.

## Command

```bash
conda run -n gigatime-tcga python scripts/run_histoprism_one_vector_smoke.py \
  --out-dir results/histoprism_one_vector_smoke \
  --device auto \
  --top-n 25
```

The script downloads the open `HuSusu/HistoPrism` split-0 checkpoint if it is
not already cached locally.

## Smoke-test output

The initial run used:

- sample: `ZEN48`
- barcode: `AAACAGCTTTCAGAAG-1`
- oncotree conditioning: `READ`
- input kind: one precomputed 1024-dimensional UNI feature vector
- output size: `38,982` predicted gene values
- ERBB2 predicted value: `0.5776057243`

Ignored result files are written under:

- `results/histoprism_one_vector_smoke/all_gene_predictions.csv`
- `results/histoprism_one_vector_smoke/top_predicted_genes.csv`
- `results/histoprism_one_vector_smoke/bottom_predicted_genes.csv`
- `results/histoprism_one_vector_smoke/selected_marker_predictions.csv`
- `results/histoprism_one_vector_smoke/histoprism_one_vector_summary.json`

## Next step

After `MahmoodLab/UNI` access is approved, the same pattern can be extended to
extract a UNI embedding from one TCGA H&E tile and feed that vector into
HistoPrism for a true raw-tile smoke test.
