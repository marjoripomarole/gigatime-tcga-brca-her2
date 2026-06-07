#!/usr/bin/env python3
"""Recompute GigaTIME per-channel virtual-vs-measured-mIF Pearson on the RELEASED test patches.

Replaces the eyeballed Fig S5 transcription (comparison-doc TODO #3) with reproducible computed
values, using GigaTIME's OWN released model (HF prov-gigatime/GigaTIME), released sample test
patches, dataset class (prov_data.HECOMETDataset_roi) and exact eval metric (8x8-box per-channel
Pearson of virtual vs measured mIF; gigatime_testing.ipynb / db_test.py).

IMPORTANT CAVEAT: these are the 50 RELEASED SAMPLE patches; the GigaTIME authors state in their
notebook that this sample is "not the actual eval set reported in the paper." So treat the output
as a faithful released-sample reproduction, not the paper's exact Fig S5 numbers.

Safety: the measured-mIF masks are gzipped pickles; their opcodes were verified (pickletools, no
execution) to reference only numpy/builtins before loading.

Run: ~/miniconda3/envs/gigatime-tcga/bin/python scripts/recompute_gigatime_mif_pearson.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

GT_SCRIPTS = Path("external/GigaTIME/scripts")
DATA_ROOT = "data/gigatime_test_patches/sample_test_data/data/"  # trailing slash => flat layout
OUT_CSV = Path("results/gigatime_mif_recompute/per_channel_pearson.csv")


# --- GigaTIME's exact box-Pearson metric (verbatim from gigatime_testing.ipynb / db_test.py) ---
def split_into_boxes(tensor, box_size):
    b, c, h, w = tensor.shape
    boxes = tensor.unfold(2, box_size, box_size).unfold(3, box_size, box_size)
    return boxes.contiguous().view(b, c, h // box_size, w // box_size, box_size, box_size)


def count_ones(boxes):
    return boxes.sum(dim=(4, 5))


def calculate_correlations(matrix1, matrix2):
    from scipy.stats import pearsonr, spearmanr
    b, c, h, w = matrix1.shape
    pear, spear = [], []
    for ch in range(c):
        pcs, scs = [], []
        for bi in range(b):
            a = matrix1[bi, ch].flatten().cpu().numpy()
            d = matrix2[bi, ch].flatten().cpu().numpy()
            valid = ~np.isnan(a) & ~np.isnan(d)
            a, d = a[valid], d[valid]
            if len(a) > 0 and len(d) > 0 and a.std() > 0 and d.std() > 0:
                pcs.append(pearsonr(a, d)[0])
                scs.append(spearmanr(a, d)[0])
            else:
                pcs.append(np.nan)
                scs.append(np.nan)
        pear.append(np.nanmean(pcs) if np.any(~np.isnan(pcs)) else np.nan)
        spear.append(np.nanmean(scs) if np.any(~np.isnan(scs)) else np.nan)
    return pear, spear


def get_box_metrics(pred, mask, box_size):
    pred_counts = count_ones(split_into_boxes(pred, box_size))
    mask_counts = count_ones(split_into_boxes(mask, box_size))
    pearson, spearman = calculate_correlations(pred_counts, mask_counts)
    return pearson, spearman


def main() -> int:
    sys.path.insert(0, str(GT_SCRIPTS))
    import warnings
    warnings.filterwarnings("ignore")
    import pandas as pd
    import torch
    import torch.nn.functional as F
    from torch.utils.data import DataLoader
    import albumentations as A
    from albumentations.augmentations import transforms as Atransforms
    from albumentations.core.composition import Compose
    from huggingface_hub import snapshot_download
    import archs
    import prov_data
    from prov_data import HECOMETDataset_roi

    channels = prov_data.common_channel_list  # 23, order == model output index
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")

    # Build metadata over the released sample patches (flat layout under DATA_ROOT).
    pairs = sorted(p[: -len("_he.png")] for p in os.listdir(DATA_ROOT) if p.endswith("_he.png"))
    meta = pd.DataFrame({"pair_name": pairs, "dir_name": DATA_ROOT})
    print(f"Patches: {len(pairs)} | device: {device}", file=sys.stderr)

    val_transform = Compose([A.Resize(512, 512), Atransforms.Normalize()])
    dataset = HECOMETDataset_roi(
        all_tile_pair=meta, tile_pair_df=meta, mask_noncell=True, transform=val_transform,
        cell_mask_label=True, dir_path=DATA_ROOT, window_size=256, split="full", standard="all",
    )
    loader = DataLoader(dataset, batch_size=4, shuffle=False, num_workers=0)

    model = archs.__dict__["gigatime"](len(channels), 3)
    local_dir = snapshot_download(repo_id="prov-gigatime/GigaTIME")
    state = torch.load(os.path.join(local_dir, "model.pth"), map_location="cpu")
    model.load_state_dict(state)
    model.eval().to(device)

    per_channel_vals: list[list[float]] = [[] for _ in channels]
    with torch.no_grad():
        for bi, (img, target_g, _info) in enumerate(loader):
            img = img.to(device)
            # Coarsen target as in GigaTIME training/eval (1/8 down then up to 512).
            ds = F.interpolate(target_g, scale_factor=1 / 8, mode="bilinear", align_corners=False)
            target = F.interpolate(ds, size=(512, 512), mode="bilinear", align_corners=False).to(device)
            b, c, h, w = img.shape
            logits = torch.zeros(b, len(channels), h, w, device=device)
            for i in range(0, h, 256):
                for j in range(0, w, 256):
                    logits[:, :, i:i + 256, j:j + 256] = model(img[:, :, i:i + 256, j:j + 256])
            pred = (torch.sigmoid(logits) > 0.5).float()
            pearson, _ = get_box_metrics(pred, target, 8)
            for ci, val in enumerate(pearson):
                if not np.isnan(val):
                    per_channel_vals[ci].append(float(val))
            print(f"  batch {bi + 1}/{len(loader)}", file=sys.stderr)

    rows = []
    for ci, ch in enumerate(channels):
        vals = per_channel_vals[ci]
        rows.append({
            "channel": ch,
            "pearson_mean": float(np.mean(vals)) if vals else float("nan"),
            "pearson_std": float(np.std(vals)) if vals else float("nan"),
            "n_patches": len(vals),
        })
    rows.sort(key=lambda r: (np.isnan(r["pearson_mean"]), -np.nan_to_num(r["pearson_mean"], nan=-9)))

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print("\nGigaTIME per-channel mIF Pearson (recomputed on the 50 released sample patches):")
    for r in rows:
        print(f"  {r['channel']:14s} {r['pearson_mean']:.3f} +/- {r['pearson_std']:.3f}  (n={r['n_patches']})")
    print(f"\nWrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
