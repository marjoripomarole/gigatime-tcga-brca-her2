# GigaTIME mIF Agreement vs Our RNA Specificity (draft figure + notes)

Status: draft paper material. The candidate centerpiece figure for the virtual-mIF
specificity-audit paper. Numbers are reproducible via
`scripts/make_gigatime_vs_rna_specificity_figure.py`.

## The one-line claim this figure supports

**Agreement with the training modality (held-out mIF) does not imply channel
specificity against an independent modality (spatial RNA) once per-tile cellularity
is controlled.** GigaTIME validates each channel only as virtual-X vs measured-X mIF;
several channels that pass that test fail an independent RNA specificity test.

## The data

| Channel | GigaTIME vs held-out mIF (Pearson, Fig S5) | Our RNA raw Spearman | Our RNA partial r \| cellularity | Own-gene row-max? |
|---|---:|---:|---:|:--:|
| CK | 0.96 | 0.146 | **0.308** | yes |
| CD3 | 0.89 | 0.428 | 0.256 | no (tracks CD11c-RNA) |
| CD20 | 0.88 | 0.214 | 0.085 | no |
| CD4 | 0.86 | 0.428 | 0.211 | no |
| **PD-L1** | **0.74** | 0.223 | **-0.072** | no |
| CD16 | 0.56 | 0.221 | -0.036 | no |
| **CD68** | 0.54 | 0.132 | **-0.333** | no |
| CD14 | 0.51 | 0.167 | -0.048 | no |
| Tryptase | 0.49 | 0.331 | 0.028 | no |
| PD-1 | 0.47 | 0.214 | 0.053 | no |
| CD11c | 0.47 | 0.330 | 0.147 | yes |
| CD8 | 0.30 | 0.358 | 0.237 | no |
| Ki67 | 0.28 | 0.337 | -0.020 | no |

Context not in the join: GigaTIME's single best channel is **DAPI (Pearson 0.98)** —
the nuclear/cellularity marker. It cannot appear in a specificity test because it *is*
cellularity; that it tops their agreement ranking is itself the tell.

Figures: `assets/gigatime_vs_rna_specificity/mif_agreement_vs_rna_specificity_scatter.png`
(primary) and `assets/gigatime_vs_rna_specificity/per_channel_bars.png`.

## How to read it

- **Hold up (specific):** CK (0.31), CD3 (0.26), CD8 (0.24), CD4 (0.21), CD11c (0.15)
  keep positive cellularity-controlled signal. CK is the cleanest — weakest raw RNA
  correlation but strongest after cellularity control, because epithelium-rich tiles are
  immune-poor.
- **Fail (cellularity, not marker):** CD68 (-0.33), PD-L1 (-0.07), CD16 (-0.04),
  CD14 (-0.05), Ki67 (-0.02) collapse to or below zero. Several of these have *high*
  mIF agreement (CD68 0.54, PD-L1 0.74) — exactly the divergence the figure shows.
- **Off-diagonal specificity failures:** own-gene is the row-max for only 2/13 channels
  (CK, CD11c). The immune channels (CD3/CD4/CD8/CD20) collectively track lymphocyte-dense
  regions; e.g. virtual-CD3 correlates with CD11c-RNA slightly more than CD3-RNA.

## Metric definitions (theirs vs ours) — for the methods contrast

GigaTIME (Valanarasu et al., Cell 2026, STAR Methods), all virtual-vs-measured-**mIF**:
- Dice (pixel): `Dice_c = 2*sum(y_c*yhat_c) / (sum(y_c)+sum(yhat_c))`, mIF binarized by Otsu.
- Pearson (cell): activation counts in 8x8 windows -> per-channel Pearson, virtual vs measured.
- Spearman (slide): activation density per 256x256 patch -> Spearman, virtual vs measured.
- External "validation": virtual-vs-virtual subtype-level agreement on TCGA (r=0.88); TCGA has no mIF.

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

## Limitations of this draft

- GigaTIME mIF Pearson values are **transcribed from the published Fig S5 / 2C panels** (+/- ~0.01).
  Verify against the source figure (or request underlying values) before publication.
- Generalization now confirmed across 9 sections / 2 platforms (4 HEST-1k Xenium patients + 3 Visium
  + Janesick Rep1/Rep2); see `hest_rna_validation_summary.md`. The single-section result holds and
  sharpens: specificity is tissue-variable, with the aggregate T-cell contrast the most reproducible.
- GigaTIME was trained on lung adenocarcinoma mIF; breast was only a generalization TMA check (mIF only).
  Our audit is the first breast + independent-modality + specificity test, but the train/test organ shift
  should be stated, not hidden.
- RNA is a proxy for protein with a concordance ceiling; report it as a lower bound on protein-level agreement.

## Draft caption (scatter)

"GigaTIME channel agreement with its held-out training modality (measured mIF, Pearson;
Valanarasu et al. 2026, Fig S5) versus our independent within-slide RNA specificity (Spearman
partial correlation controlling for per-tile cellularity, Xenium breast). High mIF agreement
does not predict RNA specificity: CD68 and PD-L1 agree well with held-out mIF yet have
non-positive cellularity-controlled RNA correlation, while only CK and aggregate T-cell channels
retain channel-specific RNA signal."

## Provenance

- GigaTIME paper: Valanarasu et al., "Multimodal AI generates virtual population for tumor
  microenvironment modeling," Cell 189:386-400 (Jan 22, 2026); doi:10.1016/j.cell.2025.11.016.
  Per-channel mIF Pearson from Fig S5 / Fig 2C; metric definitions from STAR Methods.
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
3. Confirm or replace the transcribed GigaTIME Fig S5 values with source values.
