# BCNB Exploration: Confirmed Clinical External Cohort

Status: BCNB clinical access has been resolved. The full clinical file confirms that BCNB can support the HER2-low versus HER2-zero question that TCGA cannot; the remaining task is WSI or patch acquisition and manifest construction.

## Why BCNB

- Single institution, single scanner (iScan Coreo, 200x), 1,058 patients, China, enrolled 2010-2020. A single-scanner cohort removes the slide-size / source-site acquisition confound that dominates the TCGA result by construction.
- Free with a registration form (non-commercial).
- Has clinical HER2, ER, PR, and (in the full release) reportedly HER2 expression, histological grading, Ki67, and molecular subtype.
- Tissue type is core-needle biopsy, not surgical resection (a difference from TCGA-BRCA diagnostic resections; relevant to tissue amount and tile sampling).

## Access Map: What Is Public vs Gated

There are two tiers of BCNB clinical data, and they differ critically:

1. Public (GitHub `bupt-ai-cz/BALNMP`, `code/dataset/clinical_data/`): five preprocessed `.xlsx` files for the lymph-node-metastasis task. Downloaded locally to `data/bcnb/clinical_public/` (2026-06-04). These contain only:

   | Column | Meaning |
   |---|---|
   | `p_id` | patient id (1-1058) |
   | `age` | standardized (z-scored) |
   | `tumor_size` | standardized (z-scored) |
   | `is_er_negative` | binary |
   | `is_pr_negative` | binary |
   | `is_her2_negative` | binary |

   Public-data distribution (1,058 patients): HER2-negative 781 / HER2-positive 277; ER-negative 227 / positive 831; PR-negative 268 / positive 790. This matches Table 1 of the BCNB paper (Front Oncol 2021, PMC8551965). HER2 here is binary only; there is no IHC score and no grade in the public files.

2. Gated (full dataset, via the registration form at https://bupt-ai-cz.github.io/BCNB/): the WSIs plus the full clinical `.xlsx`, which the dataset page lists as including `HER2`, a separate `HER2 expression` field, `histological grading`, `Ki67`, `molecular subtype`, tumor type, surgical info, and lymph-node fields.

## The Make-Or-Break Question, And Its Resolution

For the HER2-low versus HER2-zero question we need IHC score 0 to be separable from 1+ (and ideally grade). Before full access, the evidence looked like this:

- Public data: binary HER2 only. Cannot recover HER2-zero.
- BCNB paper (PMC8551965): Table 1 records HER2 only as Positive/Negative, and does not report histological grade. So the paper does not confirm score-level granularity.
- Dataset page: lists `HER2` and a separate `HER2 expression` field, plus `histological grading`. The naming strongly suggests `HER2 expression` is the IHC score level (0 / 1+ / 2+ / 3+), which would make HER2-zero recoverable, and that grade is available. But the actual values are not documented publicly.

Historical conclusion before full access: BCNB support for low-versus-zero was plausible but had to be settled by inspecting the `HER2 expression` and `histological grading` columns of the full clinical `.xlsx`. That gating fact is now resolved below.

## RESOLVED 2026-06-04: Full Clinical Data Obtained — BCNB Is Viable For Low-vs-Zero

Registration was approved and the full clinical file obtained (`patient-clinical-data.xlsx`, 76 KB, downloaded to `data/bcnb/`, gitignored and not redistributed per the non-commercial license). The gating question is answered: the `HER2 Expression` column IS the IHC score, and HER2-zero is recoverable.

Full clinical columns (1,058 patients): Patient ID, Age, Tumour Size, Tumour Type, ER, PR, HER2 (status), HER2 Expression (IHC score), Histological grading, Surgical, Ki67, Molecular subtype, Number of lymph node metastases, ALN status.

HER2 Expression (IHC score) distribution: 0 → 127, 1+ → 242, 2+ → 483, 3+ → 206. The binary `HER2` status already encodes ISH for 2+ cases (412 of the 2+ are HER2-negative/ISH−, 71 are HER2-positive/ISH+).

Derived clinical HER2 groups (written to `data/bcnb/bcnb_her2_labels.csv`):

| Group | Definition | N |
|---|---|---:|
| HER2-zero | IHC 0 | 127 |
| HER2-low | IHC 1+, or 2+/ISH-negative | 654 |
| HER2-positive | IHC 3+, or 2+/ISH-positive | 277 |

This is a single-scanner external HER2-low-versus-zero cohort of 654 vs 127. The HER2-zero group alone (127) is more than double all of TCGA-BRCA (61), and the whole low/zero set (781) is ~6.6x the TCGA high-trust low/zero set (118).

Grade is present and shows the literature-expected pattern: among graded low/zero cases, HER2-zero is proportionally more often grade 3 (42/90 = 47%) than HER2-low (181/577 = 31%), consistent with HER2-zero being more aggressive (see `external_validation_candidates.md` literature context). Because grade is available, the grade confounder that TCGA could not provide can now be tested and adjusted, exactly as slide-size/ER/PR were handled internally. ER also differs between groups (HER2-low 86% ER+, HER2-zero 65% ER+), and is likewise available as a covariate.

Image-input update: `paper_patches.zip` has now been downloaded and audited locally (see `bcnb_paper_patches_audit.md`). It maps cleanly to all 1,058 patients and can support a fast patch-based smoke/pilot. Full WSIs are now also downloaded locally from the approved BCNB SharePoint/OneDrive mirror and map cleanly to all 1,058 patients (see `bcnb_gigatime_full_wsi_smoke.md`).

## FIRST PATCH PILOTS 2026-06-04, SENSITIVITY UPDATED 2026-06-05: H-Optimus-0 And Virchow2 Find A Modest Non-Null Signal

The first BCNB external patch pilots are complete (`bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero.md`, `bcnb_patch_embedding_control_virchow2_hash_capped10_low_zero.md`, and `bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero.md`). They used deterministic hash-sampled capped patches (`10` patches per patient), patient-level mean foundation-model embeddings, class-balanced logistic regression, repeated stratified 5-fold CV with 5 repeats, and 200 shuffled-label permutations.

Key low-versus-zero results:

| Feature set | Balanced accuracy | AUC |
|---|---:|---:|
| H-Optimus-0 patch embedding | 0.597 | 0.640 |
| Virchow2 patch embedding | 0.600 | 0.643 |
| H-Optimus-0 + Virchow2 | 0.609 | 0.651 |
| H-Optimus-0 + Virchow2 + clinical covariates | 0.615 | 0.661 |
| Average probability ensemble | 0.610 | 0.650 |
| H-Optimus-0 + clinical covariates | 0.595 | 0.641 |
| Virchow2 + clinical covariates | 0.603 | 0.646 |
| Clinical covariates | 0.643 | 0.627 |
| Grade only | 0.595 | 0.604 |

Both single-model embedding results beat the shuffled-label null for balanced accuracy and AUC (empirical p=0.005 with 200 permutations), and the paired dual-model embedding also beats its shuffled-label null (BA 0.609 / AUC 0.651, empirical p=0.005). But the effect size is modest, the dual-model gain is small, and neither image-only model beats clinical covariates by balanced accuracy. H-Optimus-0 and Virchow2 patient-mean probabilities are highly concordant (Pearson r=0.804), arguing against a hidden strong classifier missed by one encoder.

The H-Optimus-0 patch-selection sensitivity check is also complete (`bcnb_patch_sampling_sensitivity_hoptimus0.md`): the older lexicographic capped10 sample gives BA 0.575 / AUC 0.633 at PCA20, the first hash-capped seed gives BA 0.597 / AUC 0.640, and a second independent hash seed gives BA 0.585 / AUC 0.634. This shifts the effect size slightly but preserves the same conclusion.

The clinical-stratified performance check is now complete (`bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero.md`). It scores patient-mean out-of-fold predictions inside grade, ER/PR, subtype, nodal, Ki67, and grade-by-ER slices. The pooled dual-model patient-mean result is still modest (BA 0.613 / AUC 0.660), and performance is uneven in clinically meaningful strata: ER-negative dual AUC 0.598, Luminal A dual AUC 0.557, Grade 2 dual AUC 0.617, Grade 3 dual AUC 0.668.

The score-driver analysis is now complete (`bcnb_patch_score_covariate_drivers_hoptimus0_virchow2_hash_capped10_low_zero.md`). It asks what explains the image-model score itself. HER2 label alone explains little dual-model score variance (R2 0.036), while clinical covariates + patch QC explain more (R2 0.161). After residualizing the dual image score against clinical covariates + patch QC, low/zero AUC drops from 0.660 to 0.592. That means the measured covariates explain part, but not all, of the weak image signal.

The visual score-extreme QC is now complete (`bcnb_patch_score_visual_qc_hoptimus0_virchow2_hash_capped10_low_zero.md`). It renders hash-sampled patch montages for true HER2-zero scored zero-like, HER2-low scored zero-like, and HER2-low scored low-like cases. The zero-like HER2-low false positives are not blank/low-tissue artifacts, and the zero-like score extremes are enriched for grade 3 / ER-negative / PR-negative / triple-negative profiles, matching the quantitative covariate-driver story.

This is the first real evidence that a low/zero-associated morphology signal exists outside TCGA, but it is not a strong standalone HER2-low-versus-zero classifier from the patch pilot. The result currently supports a careful interpretation: BCNB contains weak image-readable morphology/covariate signal, plausibly grade/receptor/tissue-context related, with some residual morphology not captured by the measured covariates. The next decision is whether the paper needs broader patch-sampling sensitivity, more patches per patient, or full WSI processing for stronger slide/tissue-area controls.

## Why BCNB Is Now The Priority External Cohort

- The clinical label gate is solved: HER2-zero, HER2-low, and HER2-positive can be derived from preserved IHC score plus binary HER2/ISH status.
- The acquisition confound that dominates TCGA is reduced by design: one institution and one scanner.
- Grade, ER, PR, Ki67, molecular subtype, and nodal status are available, so grade and receptor-status confounding can be modeled rather than hand-waved.
- The low/zero cohort is large enough for real robustness checks: 781 low/zero cases total, versus 118 high-trust low/zero TCGA slides.
- The tissue type is core-needle biopsy, which is not identical to TCGA diagnostic resections, but that difference is scientifically manageable and should be documented as an external-validation domain shift.

## Recommended Next Steps

1. Use the full-WSI path for the next paper-grade GigaTIME experiment now that the 1,058 patient JPG WSIs are local. Keep the patch pilots as fast foundation-model controls and sampling-sensitivity evidence.
2. Build or refresh the BCNB patch manifest with `scripts/build_bcnb_patch_manifest.py`, keeping restricted data under ignored `data/bcnb/`.
3. Extend patch-sampling/PCA sensitivity with more seeds or more patches per patient, or launch heavier WSI processing only if the paper needs stronger slide/tissue-area controls than the precomputed patch pilot can provide.
4. Run the full BCNB low/zero GigaTIME cohort from `data/bcnb/bcnb_wsi_slide_table_low_zero.csv` after confirming the balanced WSI smoke outputs.
5. Reuse the existing confound discipline: compare image embeddings against grade, ER/PR, Ki67, molecular subtype, nodal status, and tissue/slide-size features. In BCNB, slide-size/source-site should not classify low-vs-zero well; if it does, that is itself a warning sign.
6. Treat H-Optimus-0/Virchow2 as primary foundation-model controls; keep GigaTIME/DeepSpot/HistoPrism as interpretive follow-ups unless the BCNB signal survives clinical and acquisition controls.

## Local Workspace And Environment Notes

- Public clinical files: `data/bcnb/clinical_public/preprocessed-type-{0..4}.xlsx` (gitignored under `data/`).
- Full clinical file: `data/bcnb/patient-clinical-data.xlsx` (gitignored; non-commercial dataset file, not redistributed).
- Derived label table: `data/bcnb/bcnb_her2_labels.csv` (gitignored; local derivative used for analysis), reproducibly built by `scripts/build_bcnb_her2_labels.py`.
- Paper patch archive: `data/bcnb/paper_patches.zip` (gitignored; 76,578 256x256 RGB `.jpg` patches; clean patient-ID mapping; see `bcnb_paper_patches_audit.md`).
- Full WSI directory: `data/bcnb/WSIs/BCNB/WSIs/` is present locally and contains 1,058 numeric patient `.jpg` WSIs, plus 21 metadata `.json` files left from the earlier Google Drive partial download. The local WSI audit maps all 1,058 images to known BCNB patients.
- Full WSI slide tables: `data/bcnb/bcnb_wsi_slide_table_low_zero.csv` contains 781 HER2-low/zero rows with no missing patients; `data/bcnb/bcnb_wsi_slide_table_low_zero_balanced20.csv` is a deterministic 10 low / 10 zero smoke table.
- Patch manifests: `data/bcnb/bcnb_patch_manifest.csv`, `data/bcnb/bcnb_patch_manifest_capped10.csv`, preferred pilot `data/bcnb/bcnb_patch_manifest_hash_capped10.csv`, and second hash replicate `data/bcnb/bcnb_patch_manifest_hash20260605_capped10.csv` (gitignored; built by `scripts/build_bcnb_patch_manifest.py`).
- Image-input audit: `scripts/audit_bcnb_image_inputs.py` checks for `data/bcnb/WSIs/`, `data/bcnb/paper_patches.zip`, and `data/bcnb/paper_patches/` without extracting or running models.
- `openpyxl` was installed into the `gigatime-tcga` conda env on 2026-06-04 to read `.xlsx` (the full BCNB clinical file is also `.xlsx`).

## Caveats

- Core-needle biopsy tissue differs from TCGA diagnostic resections (less tissue per slide, different sampling); tile-count and tissue-fraction settings may need adjustment.
- HER2-low is a post-2019 clinical category; BCNB preserves the original IHC score, but the study should still phrase the low/zero grouping as a derived contemporary clinical grouping from historical IHC/ISH data.
- BCNB is non-commercial use only.

## Sources

- BCNB dataset page: https://bupt-ai-cz.github.io/BCNB/
- BALNMP repo (public preprocessed clinical data): https://github.com/bupt-ai-cz/BALNMP
- BCNB paper (Front Oncol 2021): https://pmc.ncbi.nlm.nih.gov/articles/PMC8551965/
