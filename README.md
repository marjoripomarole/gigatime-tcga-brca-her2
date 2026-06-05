# BRCA HER2 Pathology AI

Computational pathology workspace for studying HER2-related breast cancer states in TCGA-BRCA diagnostic H&E whole-slide images.

This repo is not a clinical HER2 classifier. It is an exploratory research workspace asking whether image-derived H&E features, starting with GigaTIME virtual mIF outputs, associate with clinically defined HER2-positive, HER2-low, and HER2-zero breast cancer states.

## Current Research Question

Can H&E-derived image features from TCGA-BRCA associate with the HER2-low versus HER2-zero boundary, and do those signals survive tissue-composition, source-site, slide-size, RNA, and model-family validation checks?

The current honest interpretation is:

- The strongest image signal is HER2-low versus HER2-zero, not HER2-positive detection.
- GigaTIME-derived virtual immune, myeloid, checkpoint, and CK-associated signals reproducibly differ between HER2-low and HER2-zero groups.
- The signal is hypothesis-generating and heavily caveated by tissue composition, slide size, and TCGA source-site imbalance; H-Optimus-0 and Virchow2 reproduce the same low-versus-zero separation and source-site collapse, so TCGA-internal evidence is now exhausted.
- Local STAR-count RNA supports broad HER2-positive status through ERBB2 expression, but does not strongly separate HER2-low from HER2-zero.
- Current local RNA files do not support direct HER2 isoform or junction-level validation.
- BCNB is now the confirmed external validation priority: full clinical data preserve HER2 IHC 0/1+/2+/3+, with 127 HER2-zero and 654 HER2-low single-scanner cases plus grade/ER/PR/Ki67 covariates. H-Optimus-0 and Virchow2 hash-capped patch pilots show a statistically non-null but modest signal (single-model BA ~0.60, AUC ~0.64; dual-model BA 0.609, AUC 0.651), comparable to clinical covariates rather than a strong standalone classifier. Patch-sampling, clinical-stratified, score-driver, and visual-QC checks change the effect size but not the conclusion: the signal is weak, partly clinical/patch-QC-associated, enriched for aggressive morphology/covariate profiles, and not HER2-low/zero classifier-grade.

## Start Here

Read these in order:

1. [docs/00_start_here.md](docs/00_start_here.md) for the current project map.
2. [docs/clinical_her2_high_trust_tile128_results.md](docs/clinical_her2_high_trust_tile128_results.md) for the current primary result.
3. [docs/advisor_brief.md](docs/advisor_brief.md) for the advisor-facing narrative.
4. [docs/RUN_REGISTRY.md](docs/RUN_REGISTRY.md) for the run-by-run evidence trail.
5. [docs/external_validation_candidates.md](docs/external_validation_candidates.md) and [docs/bcnb_exploration.md](docs/bcnb_exploration.md) for the confirmed BCNB external-validation path.
6. [scripts/README.md](scripts/README.md) for script groups and rerun entry points.

## Current Primary Result

The primary analysis is the strict high-trust 171-slide TCGA-BRCA GigaTIME run:

| Group | High-trust slides |
|---|---:|
| HER2-positive | 53 |
| HER2-low | 57 |
| HER2-zero | 61 |

Primary run:

- Model: GigaTIME virtual mIF from H&E tiles.
- Cohort: female TCGA-BRCA diagnostic H&E slides with HER2 label, file-integrity, and OpenSlide QC.
- Tile sampling: 128 random tissue tiles per slide.
- Main association: HER2-low has lower virtual immune/myeloid/checkpoint and CK-associated signal than HER2-zero.
- Main caveat: slide-size/source-site and tissue-composition confounding remain strong; generic H&E embeddings reproduce the separation, so the signal should be treated as morphology/acquisition entanglement until external validation.

Use [docs/clinical_her2_high_trust_tile128_results.md](docs/clinical_her2_high_trust_tile128_results.md) as the current presentation summary. Older 30-slide and 60-slide reports are historical.

## Repo Map

```text
configs/      Default path and run configuration.
data/         Local TCGA manifests, clinical labels, slides, and RNA files. Ignored except placeholders.
docs/         Human-readable reports, navigation, and tracked figure assets.
envs/         Conda environment specs.
external/     Ignored upstream model repositories such as GigaTIME, HistoPrism, and DeepSpot.
notebooks/    Presentation notebooks and rendered HTML summaries.
results/      Ignored machine-readable outputs from runs and smoke tests.
scripts/      Reproducible workflow scripts.
```

Important note: this cleanup keeps existing generated report paths stable because many scripts write directly to current top-level `docs/*.md` files. The new `docs/*/README.md` files are navigation overlays, not a breaking move of generated outputs.

## Setup

Create the working environment:

```bash
conda env create -f envs/gigatime-tcga.yml
conda activate gigatime-tcga
```

Gated model access may require a Hugging Face token:

```bash
export HF_TOKEN=<huggingface_read_token>
```

Gated or access-sensitive models used in this project include:

- GigaTIME
- `bioptimus/H0-mini`
- `bioptimus/H-optimus-0`
- `MahmoodLab/UNI`

## Main Workflows

Build HER2 clinical labels:

```bash
conda run -n gigatime-tcga python scripts/build_tcga_brca_clinical_her2_labels.py
conda run -n gigatime-tcga python scripts/build_bcnb_her2_labels.py
```

Audit BCNB image inputs before any model run:

```bash
conda run -n gigatime-tcga python scripts/audit_bcnb_image_inputs.py
conda run -n gigatime-tcga python scripts/build_bcnb_patch_manifest.py \
  --max-patches-per-patient 10 \
  --sampling-method hash \
  --sampling-seed 20260604 \
  --output data/bcnb/bcnb_patch_manifest_hash_capped10.csv
```

Run a patient-level BCNB patch embedding control:

```bash
conda run -n gigatime-tcga python scripts/run_bcnb_patch_embeddings.py \
  --model hoptimus0 \
  --patch-manifest data/bcnb/bcnb_patch_manifest_hash_capped10.csv \
  --groups HER2-zero,HER2-low \
  --max-patches-per-patient 10 \
  --out-dir results/bcnb_patch_embeddings_hoptimus0_hash_capped10_low_zero \
  --batch-size 16 \
  --resume

conda run -n gigatime-tcga python scripts/analyze_bcnb_patch_embedding_control.py \
  --embeddings results/bcnb_patch_embeddings_hoptimus0_hash_capped10_low_zero/patient_embeddings.csv \
  --model-label H-Optimus-0 \
  --model-id bioptimus/H-optimus-0 \
  --out-dir results/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero \
  --asset-dir docs/assets/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero \
  --out-markdown docs/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero.md

conda run -n gigatime-tcga python scripts/analyze_bcnb_patch_model_comparison.py \
  --hoptimus-embeddings results/bcnb_patch_embeddings_hoptimus0_hash_capped10_low_zero/patient_embeddings.csv \
  --virchow2-embeddings results/bcnb_patch_embeddings_virchow2_hash_capped10_low_zero/patient_embeddings.csv
```

Select/download a clinical HER2 cohort:

```bash
conda run -n gigatime-tcga python scripts/select_clinical_her2_cohort.py
conda run -n gigatime-tcga python scripts/download_clinical_her2_cohort_slides.py --only-missing
```

Run the current GigaTIME-style slide inference:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv \
  --slide-path-column slide_local_path \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2_high_trust_tile128 \
  --tile-limit 128 \
  --tile-order random \
  --batch-size 16 \
  --device auto \
  --save-tile-csv \
  --resume
```

Run key current sensitivity checks:

```bash
conda run -n gigatime-tcga python scripts/analyze_tissue_composition_sensitivity.py
conda run -n gigatime-tcga python scripts/analyze_tumor_proxy_sensitivity.py
conda run -n gigatime-tcga python scripts/analyze_clinical_covariate_sensitivity.py
conda run -n gigatime-tcga python scripts/analyze_source_site_generalization.py
conda run -n gigatime-tcga python scripts/analyze_local_erbb2_expression_validation.py
```

Run current model-family smoke tests:

```bash
conda run -n gigatime-tcga python scripts/run_histoprism_one_vector_smoke.py
conda run -n gigatime-tcga python scripts/run_deepspot_one_vector_smoke.py \
  --checkpoint /private/tmp/deepspot_weights_lung_visium/final_model.pkl \
  --genes-csv /private/tmp/deepspot_weights_lung_visium/info_highly_variable_genes.csv \
  --config /private/tmp/deepspot_weights_lung_visium/top_param_overall.yaml
```

See [scripts/README.md](scripts/README.md) for more commands.

## Model Strategy

The project is broader than one model:

- **GigaTIME** is the current primary model because it generates virtual mIF/TIME channels from H&E.
- **H-Optimus-0 and Virchow2** are the completed generic H&E embedding controls; both reproduce the TCGA low-versus-zero separation and source-site collapse.
- **BCNB H-Optimus-0 and Virchow2 patch embeddings** are the first external pilots: both show a modest non-null low-versus-zero signal. The paired model comparison shows high patient-score agreement and only a small dual-model gain. Sampling, clinical-stratified, score-driver, and visual-QC analyses do not overturn this, so neither model family is hiding a strong standalone classifier.
- **H0-mini** is not currently runnable for BCNB in this workspace because the Hugging Face gated model load fails before inference; see `docs/bcnb_h0mini_access_probe.md`.
- **H0-mini** remains a possible smaller gated Bioptimus follow-up if access is granted.
- **Phikon** is an open fallback for tile embeddings when gated model access blocks progress.
- **HistoPrism** and **DeepSpot** are interpretive follow-ups for tile/spot-level virtual gene-expression style outputs; they should not be treated as primary biological validation yet.

## Data Hygiene

`data/` and `results/` are ignored because they contain large local artifacts.

An unrelated GATK germline variant-calling scaffold (hg38 reference, BQSR known-sites VCFs, and an `NA12878` test BAM) was left over from the repo's initial template and removed on 2026-06-04, reclaiming ~19 GB locally. These were never git-tracked and not referenced by the pathology workflow:

- `data/reference/` (hg38 FASTA + BWA/GATK index) — removed
- `data/known_sites/` (dbSNP138, Mills, known-indels VCFs) — removed
- `data/raw/NA12878.bam` (Genome-in-a-Bottle germline sample) — removed

If germline variant calling is ever needed again, re-download the reference FASTA and known-sites from the GATK resource bundle.

## Guardrails

- Do not claim clinical diagnosis.
- Do not claim image AI detects HER2 isoforms.
- Do not present the TCGA HER2-low versus HER2-zero classifier as source-independent biology.
- Treat GigaTIME virtual mIF channels as model predictions, not measured multiplex immunofluorescence.
- Keep the strongest claim at the level of hypothesis-generating tissue-context association until pathologist/tumor-rich and external validation are added.
