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

Image-input update: `paper_patches.zip` has now been downloaded and audited locally (see `bcnb_paper_patches_audit.md`). It maps cleanly to all 1,058 patients and can support a fast patch-based smoke/pilot. Full WSIs are still the stronger paper-grade input if a patch pilot finds a signal worth testing with slide-level controls.

## Why BCNB Is Now The Priority External Cohort

- The clinical label gate is solved: HER2-zero, HER2-low, and HER2-positive can be derived from preserved IHC score plus binary HER2/ISH status.
- The acquisition confound that dominates TCGA is reduced by design: one institution and one scanner.
- Grade, ER, PR, Ki67, molecular subtype, and nodal status are available, so grade and receptor-status confounding can be modeled rather than hand-waved.
- The low/zero cohort is large enough for real robustness checks: 781 low/zero cases total, versus 118 high-trust low/zero TCGA slides.
- The tissue type is core-needle biopsy, which is not identical to TCGA diagnostic resections, but that difference is scientifically manageable and should be documented as an external-validation domain shift.

## Recommended Next Steps

1. Decide the image input path:
   - Patch pilot: ready for the next smoke; use `paper_patches.zip` with capped patches per patient and patient-level aggregation.
   - Full WSIs: strongest and cleanest for a paper-grade analysis, because the same tile-sampling, tissue-fraction, and slide-size controls can be reused.
2. Build or refresh the BCNB patch manifest with `scripts/build_bcnb_patch_manifest.py`, keeping restricted data under ignored `data/bcnb/`.
3. Run a one-patient patch smoke first, then a small balanced low/zero pilot, before launching a full 781-patient embedding run.
4. Run `scripts/audit_bcnb_image_inputs.py` again after any WSI download to confirm which files are present and whether patient IDs map cleanly.
5. Reuse the existing confound discipline: compare image embeddings against grade, ER/PR, Ki67, molecular subtype, nodal status, and tissue/slide-size features. In BCNB, slide-size/source-site should not classify low-vs-zero well; if it does, that is itself a warning sign.
6. Treat H-Optimus-0/Virchow2 as primary foundation-model controls; keep GigaTIME/DeepSpot/HistoPrism as interpretive follow-ups unless the BCNB signal survives clinical and acquisition controls.

## Local Workspace And Environment Notes

- Public clinical files: `data/bcnb/clinical_public/preprocessed-type-{0..4}.xlsx` (gitignored under `data/`).
- Full clinical file: `data/bcnb/patient-clinical-data.xlsx` (gitignored; non-commercial dataset file, not redistributed).
- Derived label table: `data/bcnb/bcnb_her2_labels.csv` (gitignored; local derivative used for analysis), reproducibly built by `scripts/build_bcnb_her2_labels.py`.
- Paper patch archive: `data/bcnb/paper_patches.zip` (gitignored; 76,578 256x256 RGB `.jpg` patches; clean patient-ID mapping; see `bcnb_paper_patches_audit.md`).
- Patch manifests: `data/bcnb/bcnb_patch_manifest.csv` and `data/bcnb/bcnb_patch_manifest_capped10.csv` (gitignored; built by `scripts/build_bcnb_patch_manifest.py`).
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
