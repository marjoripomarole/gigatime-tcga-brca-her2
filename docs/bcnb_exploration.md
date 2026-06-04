# BCNB Exploration: Access Map And Clinical-Data Reconnaissance

Status: working reconnaissance of the BCNB (Early Breast Cancer Core-Needle Biopsy) dataset as an external validation candidate. Compiled 2026-06-04. The goal was to settle, before committing effort, whether BCNB can support the HER2-low versus HER2-zero question that TCGA cannot.

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

## The Make-Or-Break Question, And Its Honest Status

For the HER2-low versus HER2-zero question we need IHC score 0 to be separable from 1+ (and ideally grade). Status:

- Public data: binary HER2 only. Cannot recover HER2-zero.
- BCNB paper (PMC8551965): Table 1 records HER2 only as Positive/Negative, and does not report histological grade. So the paper does not confirm score-level granularity.
- Dataset page: lists `HER2` and a separate `HER2 expression` field, plus `histological grading`. The naming strongly suggests `HER2 expression` is the IHC score level (0 / 1+ / 2+ / 3+), which would make HER2-zero recoverable, and that grade is available. But the actual values are not documented publicly.

Conclusion: whether BCNB supports the low-versus-zero split is plausible but UNCONFIRMED, and can only be settled by registering and inspecting the `HER2 expression` and `histological grading` columns of the full clinical `.xlsx`. This is the single gating fact.

## Why It Is Worth Registering Anyway

- If `HER2 expression` encodes IHC 0/1+/2+/3+: BCNB becomes a single-scanner external HER2-low-versus-zero cohort. With 781 HER2-negative patients, a typical 0-versus-low split would yield on the order of a few hundred HER2-zero and a few hundred HER2-low cases, far beyond the 61 HER2-zero ceiling of all TCGA-BRCA. This is the ideal test of whether any low-versus-zero signal survives once acquisition is controlled.
- If `HER2 expression` is only binary after all: BCNB still provides a clean single-scanner external test of HER2-positive versus HER2-negative reproducibility (277 vs 781), and a check on whether acquisition-style confounds (which should be absent in a single-scanner cohort) reappear.
- Grade: if `histological grading` is present, it lets us test the literature-motivated confounder (grade differs low vs zero; see `external_validation_candidates.md`) that TCGA could not provide.

## Recommended Next Steps

1. Register at https://bupt-ai-cz.github.io/BCNB/ (name, email, institution, country) and download the full clinical `.xlsx`. This is a manual, human step.
2. Once the full clinical file is local, inspect the `HER2 expression` and `histological grading` column values to confirm (a) IHC-0-vs-1+ separability and (b) grade availability/encoding.
3. If HER2-zero is recoverable, build a single-scanner BCNB HER2-low/HER2-zero cohort and run the existing pipeline. The embedding runners (`scripts/run_hoptimus_tcga_brca.py`, `scripts/run_virchow2_tcga_brca.py`) and the control (`scripts/analyze_hoptimus_embedding_control.py`) all take an arbitrary slide table, so BCNB can be processed with the same machinery; only a small BCNB slide-list/label builder is needed.
4. Apply the same confound discipline: because BCNB is single-scanner, slide-size/source-site baselines should NOT classify low-vs-zero well; if they do, that is itself an important negative finding.

## Local Workspace And Environment Notes

- Public clinical files: `data/bcnb/clinical_public/preprocessed-type-{0..4}.xlsx` (gitignored under `data/`).
- `openpyxl` was installed into the `gigatime-tcga` conda env on 2026-06-04 to read `.xlsx` (the full BCNB clinical file is also `.xlsx`).

## Caveats

- Core-needle biopsy tissue differs from TCGA diagnostic resections (less tissue per slide, different sampling); tile-count and tissue-fraction settings may need adjustment.
- HER2-low is a post-2019 clinical category; a 2010-2020 cohort scored for it requires that the `HER2 expression` field preserves the original IHC score, not a re-derived binary.
- BCNB is non-commercial use only.

## Sources

- BCNB dataset page: https://bupt-ai-cz.github.io/BCNB/
- BALNMP repo (public preprocessed clinical data): https://github.com/bupt-ai-cz/BALNMP
- BCNB paper (Front Oncol 2021): https://pmc.ncbi.nlm.nih.gov/articles/PMC8551965/
