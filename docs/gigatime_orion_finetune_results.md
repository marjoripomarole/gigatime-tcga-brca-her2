# Does in-domain training rescue virtual-mIF specificity? GigaTIME × Orion CRC

Status: result of the CRC-Orion controlled experiment. Answers the open question left by the
cautionary paper — is GigaTIME's per-channel specificity ceiling **domain shift** (fixable by
in-domain training) or **intrinsic** to predicting protein-mIF from H&E?

## Design

Orion CRC (Lin et al. 2023; same-section H&E + 17-plex CyCIF) lets us test entirely **in-domain
and in-modality** (protein-vs-protein), with no RNA proxy and no organ-transfer confound. Targets
were built **without the 60–150 GB raw CyCIF image**: the release ships a single-cell table (per-cell
marker intensity + X/Y centroid), so each Otsu-gated marker-positive cell is painted as a disk at its
centroid → cell-level binary marker maps, coarsened /8 (GigaTIME-style). The same within-slide audit
as the breast work (per-256px-tile virtual activation vs marker-positive-cell density, partial
correlation controlling for total cells) is reused verbatim.

- **Train:** CRC02–06 (3,741 tissue tiles), warm-started from the released GigaTIME UNet++ with a
  fresh 9-channel head, BCE+Dice, 40 epochs, MPS, native 256 px (matching GigaTIME inference).
- **Held-out test:** CRC01 (2,900 tiles), never seen in training.
- **Baseline:** the released (lung-trained) GigaTIME on the same held-out CRC01 (out-of-domain).
- **Channels:** the 9 markers overlapping GigaTIME's panel — CD3, CD8, CD4, CD20, CD68, PD-L1, PD-1, Ki67, CK.

## Result — in-domain training rescues every failing channel

Cellularity-controlled partial correlation (virtual channel vs own-marker cell density) on held-out CRC01:

| Channel | Released (out-of-domain) | Fine-tuned (in-domain) | shift |
|---|---:|---:|---:|
| CD68 (macrophage) | **−0.16** | **+0.53** | **+0.69** |
| PD-L1 (checkpoint) | **−0.11** | **+0.47** | **+0.58** |
| CD8 (T cell) | −0.02 | **+0.43** | +0.45 |
| Ki67 (proliferation) | 0.22 | 0.46 | +0.24 |
| CK (epithelium) | 0.65 | 0.75 | +0.10 |
| CD3 (T cell) | 0.52 | 0.61 | +0.09 |
| CD4 (T cell) | 0.51 | 0.60 | +0.09 |
| CD20 (B cell) | 0.20 | 0.22 | +0.02 |
| PD-1 | 0.05 | 0.07 | +0.02 |
| **partial r > 0 (95% CI)** | **6/9** | **9/9** | |

The channels that were *negative* out-of-domain — CD68, PD-L1, CD8 — become *strongly positive* after
in-domain training, and **all 9 channels now have channel-specific signal surviving the cellularity
control (9/9 vs 6/9).** CK also reconverged (0.75).

## Verdict

**The specificity ceiling is domain shift, not intrinsic.** The same GigaTIME architecture, trained
in-domain, *does* learn marker-specific channels — including the macrophage/checkpoint channels that
failed on every out-of-domain test (breast-vs-RNA and lung-model-on-CRC). So the failures documented
in the cautionary paper are a property of applying the model **out of domain**, not a fundamental
limit of H&E→virtual-mIF.

**Implication for breast:** acquiring breast H&E + multiplex protein (mIF/CODEX) and fine-tuning would
very likely yield marker-specific breast virtual channels. This makes Option A (true breast
fine-tuning) scientifically worthwhile — gated only on **data**, which the dataset survey found is not
public for breast (so internal Sírio-Libanês mIF, or newly generated Orion/CODEX+H&E, would be needed).

## Honest caveats

- **CRC, not breast.** This isolates the domain-shift-vs-intrinsic question (organ-general), but it is
  a CRC demonstration; breast confirmation requires breast data.
- **Per-marker signal ≠ perfect channel isolation.** Own-marker is the single top correlate (row-max)
  for only 3/9 channels even fine-tuned (CK, CD3, PD-L1); co-localized immune markers (e.g. CD68 vs CD4,
  CD8 vs PD-1) still cross-correlate because those cell types share immune-rich regions. The rescue is
  in the (cellularity-controlled) *own-marker* signal, which is now positive for all 9.
- The fine-tune is modest (5 training specimens, 40 epochs, MPS); more data/compute would only
  strengthen it. The released model and fine-tuned model are evaluated by the identical pipeline at the
  identical native resolution (256 px) for a fair comparison.

## Reproduce

- `scripts/download_orion_specimens.py` (H&E + single-cell table per specimen; skips the raw image)
- `scripts/validate_gigatime_orion.py` (released-model in-domain baseline)
- `scripts/train_gigatime_orion.py` (warm-start fine-tune + held-out specificity eval)
- Reports: `results/gigatime_orion_baseline/CRC01/orion_baseline_report.json`,
  `results/gigatime_orion_finetune/orion_finetune_heldout_CRC01.json`; data under gitignored
  `data/orion_crc/`.
