# GigaTIME mIF Agreement vs Our RNA Specificity (draft figure + notes)

Status: draft paper material. The candidate centerpiece figure for the virtual-mIF
specificity-audit paper. Fully reproducible: the RNA specificity via
`scripts/validate_gigatime_xenium_rna.py`, the figure via
`scripts/make_gigatime_vs_rna_specificity_figure.py`, and GigaTIME's own mIF agreement
**recomputed** (not eyeballed) on its released test patches via
`scripts/recompute_gigatime_mif_pearson.py`.

## The one-line claim this figure supports

**Agreement with the training modality (measured mIF) does not imply channel
specificity against an independent modality (spatial RNA) once per-tile cellularity
is controlled.** GigaTIME validates each channel only as virtual-X vs measured-X mIF;
several channels that pass that test fail an independent RNA specificity test.

## The data

GigaTIME mIF agreement is GigaTIME's own 8x8-box activation-count Pearson (virtual vs
measured mIF), **recomputed on the released 50-patch sample test set with the released
model** (`prov-gigatime/GigaTIME`), sorted descending. RNA columns are our within-slide
Xenium audit (Rep1).

| Channel | GigaTIME vs measured mIF (Pearson, recomputed) | Our RNA raw Spearman | Our RNA partial r \| cellularity | Own-gene row-max? |
|---|---:|---:|---:|:--:|
| CD11c | 0.56 | 0.330 | **0.147** | yes |
| **PD-L1** | **0.55** | 0.223 | **-0.072** | no |
| **CK** | **0.53** | 0.146 | **0.308** | yes |
| CD16 | 0.52 | 0.221 | -0.036 | no |
| **CD68** | **0.52** | 0.132 | **-0.333** | no |
| CD4 | 0.47 | 0.428 | 0.211 | no |
| CD3 | 0.45 | 0.428 | 0.256 | no |
| Ki67 | 0.36 | 0.337 | -0.020 | no |
| CD14 | 0.30 | 0.167 | -0.048 | no |
| CD8 | 0.23 | 0.358 | 0.237 | no |
| PD-1 | 0.17 | 0.214 | 0.053 | no |
| Tryptase | 0.13 | 0.331 | 0.028 | no |
| CD20 | n/a (B cells too sparse) | 0.214 | 0.085 | no |

The released 50-patch sample is **not** the paper's full test set — the paper's Fig S5
reports higher full-set values (e.g. CK ~0.96, CD3 ~0.89) — but even GigaTIME's own
released-sample agreement (<=0.56 for markers) does not predict RNA specificity: at
~0.5 mIF agreement the cellularity-controlled RNA partial spans the full range from
**CD68 -0.33** to **CK +0.31**. CD20 had no measurable mIF agreement on this sample.

Context: GigaTIME's single best channel is **DAPI (recomputed Pearson 0.72)** — the
nuclear/cellularity marker. It cannot appear in a specificity test because it *is*
cellularity; that it tops the agreement ranking is itself the tell.

Figures: `assets/gigatime_vs_rna_specificity/mif_agreement_vs_rna_specificity_scatter.png`
(primary) and `assets/gigatime_vs_rna_specificity/per_channel_bars.png`.

## How to read it

- **Hold up (specific):** CK (0.31), CD3 (0.26), CD8 (0.24), CD4 (0.21), CD11c (0.15)
  keep positive cellularity-controlled signal. CK is the cleanest — weakest raw RNA
  correlation but strongest after cellularity control, because epithelium-rich tiles are
  immune-poor.
- **Fail (cellularity, not marker):** CD68 (-0.33), PD-L1 (-0.07), CD16 (-0.04),
  CD14 (-0.05), Ki67 (-0.02) collapse to or below zero. Several of these have *comparable*
  mIF agreement to the channels that hold (CD68 0.52, PD-L1 0.55, vs CK 0.53) — exactly the
  dissociation the figure shows.
- **Off-diagonal specificity failures:** own-gene is the row-max for only 2/13 channels
  (CK, CD11c). The immune channels (CD3/CD4/CD8/CD20) collectively track lymphocyte-dense
  regions; e.g. virtual-CD3 correlates with CD11c-RNA slightly more than CD3-RNA.

## Metric definitions (theirs vs ours) — for the methods contrast

GigaTIME (Valanarasu et al., Cell 2026, STAR Methods), all virtual-vs-measured-**mIF**:
- Dice (pixel): `Dice_c = 2*sum(y_c*yhat_c) / (sum(y_c)+sum(yhat_c))`, mIF binarized by Otsu.
- Pearson (cell): activation counts in 8x8 windows -> per-channel Pearson, virtual vs measured.
- Spearman (slide): activation density per 256x256 patch -> Spearman, virtual vs measured.
- External "validation": virtual-vs-virtual subtype-level agreement on TCGA (r=0.88); TCGA has no mIF.

Our recompute of GigaTIME's mIF Pearson (`scripts/recompute_gigatime_mif_pearson.py`) uses
GigaTIME's **own** released model, released test patches, dataset class and the identical 8x8-box
metric above — so the x-axis is GigaTIME's own claim measured on released data, not our reinterpretation.

Ours (`scripts/validate_gigatime_xenium_rna.py`), virtual-vs-independent-**RNA**, within-slide:
- Raw: per-tile Spearman of virtual-channel activation vs binned transcript density of the channel's gene(s).
- Specificity matrix: channel x gene-set Spearman; is own-gene the row maximum?
- Partial r: Spearman partial correlation controlling for total per-tile transcript density (cellularity).
- All with spatial block-bootstrap 95% CIs; alignment auto-validated (100% in-bounds, 36.5x tissue enrichment).

## Fairness caveats and how to defend (a reviewer WILL raise these)

GigaTIME predicts protein (mIF); our reference is RNA. A reviewer will say low RNA
specificity could be RNA-protein discordance, not model failure. Defenses, strongest first:

1. **Row-max is internal to RNA.** virtual-CD3 tracking CD11c-RNA above CD3-RNA is a
   specificity failure with no cross-modality excuse.
2. **The divergence is channel-specific, not uniform.** CK and T cells hold; CD68/PD-L1/Ki67
   collapse. A uniform RNA-protein ceiling cannot produce that pattern.
3. **Partial correlation isolates cellularity**, which is modality-agnostic.
4. **Lead with CD68 (-0.33), not PD-L1.** Macrophage RNA-protein concordance is fine, so its
   negative partial is hard to dismiss; PD-L1 has known RNA-protein discordance, so it is the
   weakest example to foreground.
5. **The dissociation is also model-level, not just GigaTIME.** A second virtual-mIF model
   (ROSIE) audited identically disagrees with GigaTIME on which channels are specific
   (see `gigatime_vs_rosie_field_level.md`) — so the unreliability is a property of the approach.

## Limitations of this draft

- GigaTIME mIF agreement is now **recomputed** on the released sample (no longer eyeballed). But
  the released 50-patch sample is not the paper's full test set, so the recomputed agreements are
  lower than the paper's reported Fig S5 values; they are used here only to show that GigaTIME's own
  measured agreement does not predict RNA specificity. (To use the paper's exact full-set numbers
  would require the article's source data, which is behind the Cell paywall as a figure.)
- Generalization confirmed across 9 sections / 2 platforms (4 HEST-1k Xenium patients + 3 Visium
  + Janesick Rep1/Rep2); see `hest_rna_validation_summary.md`. The single-section result holds and
  sharpens: specificity is tissue-variable, with the aggregate T-cell contrast the most reproducible.
- GigaTIME was trained on lung adenocarcinoma mIF; breast was only a generalization TMA check (mIF only).
  Our audit is the first breast + independent-modality + specificity test, but the train/test organ shift
  should be stated, not hidden.
- RNA is a proxy for protein with a concordance ceiling; report it as a lower bound on protein-level agreement.

## Draft caption (scatter)

"GigaTIME channel agreement with its own measured mIF (Pearson, GigaTIME's 8x8-box metric recomputed
on the released test patches with the released model) versus our independent within-slide RNA
specificity (Spearman partial correlation controlling for per-tile cellularity, Xenium breast).
GigaTIME's own measured-mIF agreement does not predict RNA specificity: PD-L1 and CD68 agree ~0.5
with measured mIF yet have non-positive cellularity-controlled RNA correlation, while CK retains the
strongest channel-specific RNA signal."

## Provenance

- GigaTIME paper: Valanarasu et al., "Multimodal AI generates virtual population for tumor
  microenvironment modeling," Cell 189:386-400 (Jan 22, 2026); doi:10.1016/j.cell.2025.11.016.
  Full-set per-channel mIF Pearson in Fig S5 / Fig 2C; metric definitions from STAR Methods.
- GigaTIME mIF agreement (x-axis): **recomputed** on the released sample test patches
  (`prov-gigatime/GigaTIME` model + Dropbox sample data linked in the GigaTIME README) with
  `scripts/recompute_gigatime_mif_pearson.py`; output `results/gigatime_mif_recompute/per_channel_pearson.csv`.
  Pickle masks were verified non-executing (pickletools) before loading.
- Our RNA audit: `docs/xenium_breast_rna_validation_results.md`,
  `results/gigatime_xenium_rna_validation/xenium_rna_validation_report.json`
  (10x Xenium Human Breast Rep1; Janesick et al. 2023).
- Regenerate: `~/miniconda3/envs/gigatime-tcga/bin/python scripts/make_gigatime_vs_rna_specificity_figure.py`

## TODO before this is paper-ready

1. DONE (2026-06-06): replicated across 7 more sections on a second platform — 4 HEST-1k Xenium
   IDC patients + 3 HEST-1k Visium whole-transcriptome (IDC+ILC); aggregated with Janesick Rep1/Rep2
   in `hest_rna_validation_summary.md` (9 sections, 2 platforms). T-cell channels (CD3/CD8/CD4)
   consistently specific (8/9), CK variable (6/9), CD68/CD14/CD16/PD-L1/Ki67 never specific (CD68 0/9).
2. DONE (2026-06-06): audited ROSIE (Wu et al. 2025, ConvNeXt, 50-marker; HF ericwu09/ROSIE) on the same 9 sections
   via the identical pipeline. Both models show weak, tissue-variable specificity AND disagree on which channels are
   trustworthy (per-measurement concordance Pearson r=0.12; 44/83 channel-specific calls differ): only the T-cell
   channels (CD8/CD4) are reliably shared; GigaTIME recovers CD3/CD11c/CK, ROSIE instead recovers CD14/CD68. Field-level
   claim established. See `gigatime_vs_rosie_field_level.md`.
3. DONE (2026-06-06): replaced the eyeballed Fig S5 transcription with GigaTIME's own mIF Pearson
   **recomputed** on the released test patches (released model + its 8x8-box metric;
   `scripts/recompute_gigatime_mif_pearson.py`). The figure x-axis is now reproducible. Caveat: the
   released sample is not the paper's full set (recomputed agreements are lower), but the dissociation
   between mIF agreement and RNA specificity holds (PD-L1/CD68 ~0.5 mIF yet fail RNA; CK holds).
