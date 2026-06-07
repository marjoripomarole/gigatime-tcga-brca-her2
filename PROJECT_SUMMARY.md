# Project Summary

Date: 2026-06-05. One-page capture of the whole arc and the defensible conclusions.
For navigation see `docs/00_start_here.md`; for the run-by-run trail see `docs/RUN_REGISTRY.md`.

## One sentence

Tests whether image features from breast-cancer H&E slides can distinguish clinically
defined **HER2-low vs HER2-zero** disease, using GigaTIME (an H&E to virtual
multiplex-immunofluorescence model) plus generic pathology foundation models, while
explicitly stress-testing every confound.

## The arc

1. **An apparent signal.** On TCGA-BRCA, GigaTIME's virtual immune/myeloid/checkpoint
   and CK channels differed between HER2-low and HER2-zero. (It never detected
   HER2-*positive* — that stayed at chance. The entire result is low-vs-zero.)

2. **The signal is mostly a confound.** It is most parsimoniously TCGA
   acquisition/composition batch, not biology:
   - Raw **slide file size** alone separated low/zero better than GigaTIME (~0.88 vs ~0.75 BA).
   - Holding out the **source site** collapsed GigaTIME toward chance.
   - Two **generic** foundation models (H-Optimus-0, Virchow2) reproduced the same
     separation and the same source-site collapse — so the "virtual immune" framing is
     not needed to explain it. TCGA-internal evidence is exhausted.

3. **A real but weak, non-HER2 signal underneath.** Literature confirms HER2-zero tumors
   are genuinely higher-grade and more immune-rich, and both are visible on H&E. So the
   image signal is **real morphology (grade + immune context) entangled with batch, but
   NOT HER2 protein.**

4. **External validation (BCNB, single-scanner — removes the batch confound by design).**
   127 HER2-zero vs 654 HER2-low, with grade/ER/PR/Ki67. The image signal survives but is
   **modest** (AUC ~0.60–0.66), comparable to clinical covariates and partly explained by
   them — a weak association, not a standalone HER2 classifier.

5. **RNA validation of GigaTIME itself (2026-06-05).** First check of the virtual channels
   against ground-truth RNA, within-slide on a 10x Xenium Human Breast section (spatial RNA
   + matched H&E). Raw per-channel correlations are positive (Spearman r 0.13–0.43), but a
   cellularity-controlled **specificity** analysis shows most of that is generic tissue
   density: own-gene is the top match for only **2/13** channels, and after partialling out
   per-tile transcript density only **CK 0.31** (epithelium) and the **T-cell channels
   CD3 0.26 / CD8 0.24 / CD4 0.21** retain channel-specific signal. Ki67/CD14/CD16/PD-L1
   collapse to ~0 and **CD68 ("macrophage") goes negative**. GigaTIME channels reflect a
   broad epithelial-vs-immune/cellularity contrast, not faithful per-marker stains.

6. **Cross-sample generalization (2026-06-06).** The audit was replicated across 9 within-slide
   sections on two platforms — 4 HEST-1k Xenium IDC patients + 3 HEST-1k Visium (IDC+ILC) + the two
   Janesick sections. Only the **T-cell channels (CD3/CD8/CD4)** are consistently marker-specific
   (CD3/CD8 positive in 8/9), **CK** is specific in 6/9 tumors, and **CD68/CD14/CD16/PD-L1/Ki67 are
   never specific** (CD68 0/9). Specificity is tissue-variable — one Xenium patient and the ILC case
   show essentially none — confirming and sharpening point 5. See `docs/hest_rna_validation_summary.md`.

7. **Two-model field-level confirmation (2026-06-06).** A second, independent H&E→virtual-mIF model — ROSIE
   (Wu et al. 2025, ConvNeXt, 50 markers) — was run through the identical RNA-specificity audit on the same 9
   sections. ROSIE also shows only weak, tissue-variable marker specificity, and the two models **disagree on which
   channels are trustworthy** (per-measurement concordance Pearson r=0.12; 44/83 channel-specific calls differ).
   Only the T-cell channels (CD8/CD4) are reliably shared; GigaTIME recovers CD3/CD11c/CK while ROSIE instead
   recovers CD14/CD68. So the limited, tissue-dependent specificity is a property of the **H&E→virtual-mIF approach**,
   not of GigaTIME — no single model's virtual channels can be trusted as quantitative readouts.
   See `docs/gigatime_vs_rosie_field_level.md`.

## Strongest findings (ranked)

1. **The flagship HER2 imaging signal is a confound, not biology** — demonstrated cleanly by
   two independent generic-embedding controls and source-site holdout. Most rigorous and
   publishable result.
2. **The low-vs-zero morphology signal is real but weak and non-specific** — survives only
   modestly (~0.60 AUC) in a clean single-scanner external cohort, consistent with grade /
   immune-context, not HER2 status.
3. **GigaTIME's virtual channels are only weakly marker-specific** — even as interpretive
   overlays they mostly track cellularity/epithelium; only CK and aggregate T-cell channels
   hold up against RNA. They cannot rescue a biological claim.

## What the evidence does NOT support

- HER2 status (low vs zero, or positive detection) as a reliable H&E imaging biomarker.
- Treating GigaTIME virtual channels as quantitative cell-type readouts or load-bearing
  biological evidence.
- Any HER2-isoform / transcript-junction claim (no transcript-level RNA in hand).

## Recommended path

Write the **cautionary-methods paper**: an apparent HER2 imaging biomarker is explained by
scanner/acquisition batch plus tumor grade and immune context — shown via generic-embedding
controls, external single-scanner validation (BCNB), and a within-slide RNA-specificity
audit of the interpretive model (Xenium). Centerpiece: the two-model embedding control;
key comparison: Valieris et al. 2024. GigaTIME appears as honest interpretive context with
its specificity caveats, not as primary evidence.

## Key artifacts

- Confound controls: `docs/clinical_her2_high_trust_tile128_hoptimus_embedding_control.md`,
  `..._virchow2_embedding_control.md`.
- External cohort: `docs/bcnb_exploration.md` and the BCNB patch/WSI reports.
- RNA validation: `docs/xenium_breast_rna_validation_probe.md`,
  `docs/xenium_breast_rna_validation_results.md`
  (scripts `scripts/probe_xenium_breast_rna_validation.py`,
  `scripts/validate_gigatime_xenium_rna.py`).
- Environment note: `conda run -n gigatime-tcga` is broken on this machine; use
  `~/miniconda3/envs/gigatime-tcga/bin/python` directly.
