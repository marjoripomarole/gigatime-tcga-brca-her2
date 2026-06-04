# BCNB Paper Patches Audit

Status: completed 2026-06-04. This records the patch-first BCNB image-input triage. The raw archive is local under ignored `data/bcnb/` and is not redistributed.

## Source And Local Artifact

- Source: BALNMP GitHub README "processed WSI patches" Google Drive link (`paper_patches.zip`, file id `1wY5KIVixdwzZZq2m0IoqmBLp0YlwBAz6`).
- Local path: `data/bcnb/paper_patches.zip` (gitignored).
- Zip integrity: `zipfile.testzip()` returned `None`.
- Size on disk: ~1.8 GB.

## Archive Structure

The archive contains precomputed 256x256 RGB `.jpg` patches under patient folders:

```text
patches/938/938_9_256_256.jpg
patches/938/938_13_0_0.jpg
patches/938/938_1_0_256.jpg
...
```

Counts from the zip central directory:

| Metric | Count |
|---|---:|
| Total zip entries including directories | 77,637 |
| Image files | 76,578 |
| Metadata files inside zip | 0 |
| Patients mapped to BCNB label table | 1,058 / 1,058 |
| Ambiguous image filenames | 0 |
| Unmatched image filenames | 0 |

The mapping is clean because the first folder below `patches/` is the numeric `Patient ID` used by `data/bcnb/bcnb_her2_labels.csv`.

Two local manifests were generated from the zip central directory:

| Manifest | Rows | Purpose |
|---|---:|---|
| `data/bcnb/bcnb_patch_manifest.csv` | 76,578 | Full patch inventory |
| `data/bcnb/bcnb_patch_manifest_capped10.csv` | 10,580 | Deterministic smoke/pilot input with 10 patches per patient |

## HER2 Group Coverage

| Group | Patients | Patches | Median patches / patient | Min | Max | Mean |
|---|---:|---:|---:|---:|---:|---:|
| HER2-zero | 127 | 9,644 | 52 | 10 | 421 | 75.9 |
| HER2-low | 654 | 49,868 | 37 | 10 | 2,915 | 76.3 |
| HER2-positive | 277 | 17,066 | 34 | 10 | 1,241 | 61.6 |

Every patient has at least 10 patches. Patch counts are highly uneven, so analysis must aggregate at patient level or use a fixed/capped patch count per patient. A naive patch-level split would leak patient identity and overweight large-patch patients.

## Visual QC

One patch per HER2 group was extracted to `data/bcnb/patch_samples/` for local inspection. Each sample is 256x256 RGB and contains plausible H&E tissue. Some sampled patches include white masked/background regions, so the model pipeline should retain tissue/background QC even though the patches are precomputed.

## Interpretation

`paper_patches.zip` is usable for a fast BCNB smoke/pilot because it has complete patient coverage and clean patient-ID mapping. It is not a full substitute for WSIs in the strongest paper-grade analysis, because the patches are precomputed from annotated tumor regions and do not preserve whole-slide acquisition features such as slide size, tissue area, or full tissue-composition context. The appropriate first experiment is therefore:

1. Build a patch manifest from the zip central directory.
2. Run a one-patient or tiny balanced low/zero embedding smoke from extracted patches.
3. For a pilot, use `bcnb_patch_manifest_capped10.csv` or another capped manifest and aggregate embeddings to patient level before classification.
4. Treat full WSI download as still preferable if the patch pilot finds a signal worth testing with slide-level controls.

## Commands Run

```bash
conda run -n gigatime-tcga python scripts/audit_bcnb_image_inputs.py \
  --report-json data/bcnb/paper_patches_audit.json
conda run -n gigatime-tcga python scripts/build_bcnb_patch_manifest.py
conda run -n gigatime-tcga python scripts/build_bcnb_patch_manifest.py \
  --max-patches-per-patient 10 \
  --output data/bcnb/bcnb_patch_manifest_capped10.csv
```
