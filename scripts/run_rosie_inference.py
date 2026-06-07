#!/usr/bin/env python3
"""ROSIE (Wu et al., Nat Commun 2025) H&E->virtual-mIF inference for the RNA-specificity audit.

Loads the public ROSIE ConvNeXt-small checkpoint (HF ericwu09/ROSIE, best_model_single.pth,
CC-BY-NC) and produces per-tile mean activations per marker on the SAME 256 px tissue-tile grid
the GigaTIME pipeline uses, mapped to the pipeline's RNA channel names, so
validate_gigatime_hest_rna can audit ROSIE with the identical stats core.

Scale handling: ROSIE's native H&E resolution is 0.3775 um/px, so its 128 px input patch covers a
48.3 um field (Wu et al. 2025). For each 256 px tissue tile we sample a grid x grid of patch
centers, extract a physical-field-matched (48.3 um) crop from the full-res H&E around each center,
downsample it to 128 px (ROSIE's native patch size) then to 224 px (its ViT/CNN input), ImageNet-
normalize, run ROSIE, and average the RAW 50-vector predictions across the grid -> per-tile channel
means. Raw outputs are used because the ROSIE README specifies raw values for quantitative analysis
(the percentile/box-blur postprocessing in evaluate.py is for human-viewable images only).

This module is imported by validate_gigatime_hest_rna.py (--model rosie); it has no CLI of its own.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# ROSIE output channel order (HF README; channel i -> name).
ROSIE_CHANNELS = [
    "DAPI", "CD45", "CD68", "CD14", "PD1", "FoxP3", "CD8", "HLA-DR", "PanCK", "CD3e",
    "CD4", "aSMA", "CD31", "Vimentin", "CD45RO", "Ki67", "CD20", "CD11c", "Podoplanin", "PDL1",
    "GranzymeB", "CD38", "CD141", "CD21", "CD163", "BCL2", "LAG3", "EpCAM", "CD44", "ICOS",
    "GATA3", "Gal3", "CD39", "CD34", "TIGIT", "ECad", "CD40", "VISTA", "HLA-A", "MPO",
    "PCNA", "ATM", "TP63", "IFNg", "Keratin8/18", "IDO1", "CD79a", "HLA-E", "CollagenIV", "CD66",
]

# Pipeline (RNA-audit) channel -> ROSIE channel name. ROSIE lacks CD16/CD138/T-bet/Tryptase.
CHANNEL_MAP = {
    "CD3": "CD3e", "CD8": "CD8", "CD4": "CD4", "CD20": "CD20", "CD68": "CD68", "CD14": "CD14",
    "CD11c": "CD11c", "PD-1": "PD1", "PD-L1": "PDL1", "CK": "PanCK", "Ki67": "Ki67", "CD34": "CD34",
}

ROSIE_NATIVE_MPP = 0.3775           # Wu et al. 2025: ROSIE native H&E resolution
ROSIE_PATCH_PX = 128                # ROSIE input patch size (native px)
ROSIE_FIELD_UM = ROSIE_PATCH_PX * ROSIE_NATIVE_MPP  # 48.32 um physical field of one patch
ROSIE_NUM_OUTPUTS = 50

IMAGENET_MEAN = np.array([0.485, 0.456, 0.406], dtype=np.float32)
IMAGENET_STD = np.array([0.229, 0.224, 0.225], dtype=np.float32)


def resolve_device(torch, device: str):
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device)


def load_rosie_model(torch, weights_path: str, device):
    import torch.nn as nn
    from torchvision import models

    model = models.convnext_small(weights=None)
    model.classifier[2] = nn.Linear(model.classifier[2].in_features, ROSIE_NUM_OUTPUTS)
    ckpt = torch.load(weights_path, map_location="cpu", weights_only=False)
    state = ckpt["model_state_dict"] if isinstance(ckpt, dict) and "model_state_dict" in ckpt else ckpt
    state = {(k[7:] if k.startswith("module.") else k): v for k, v in state.items()}
    model.load_state_dict(state)
    model.eval().to(device)
    return model


def _extract(he: np.ndarray, cx: int, cy: int, size: int, H: int, W: int) -> np.ndarray:
    """Extract a size x size RGB crop centered at (cx, cy) from the full-res H&E, zero-padded at borders."""
    half = size // 2
    y0, y1, x0, x1 = cy - half, cy - half + size, cx - half, cx - half + size
    yy0, xx0, yy1, xx1 = max(0, y0), max(0, x0), min(H, y1), min(W, x1)
    crop = he[yy0:yy1, xx0:xx1]
    if crop.shape[0] != size or crop.shape[1] != size:
        out = np.zeros((size, size, 3), dtype=np.uint8)
        out[yy0 - y0:yy0 - y0 + crop.shape[0], xx0 - x0:xx0 - x0 + crop.shape[1]] = crop
        crop = out
    return crop


def _prep_patch(crop: np.ndarray) -> np.ndarray:
    """Physical-field crop -> ROSIE input: downsample to 128 px (native) then 224 px (model), normalize, CHW."""
    import cv2

    p128 = cv2.resize(crop, (ROSIE_PATCH_PX, ROSIE_PATCH_PX), interpolation=cv2.INTER_AREA)
    p224 = cv2.resize(p128, (224, 224), interpolation=cv2.INTER_LINEAR)
    arr = p224.astype(np.float32) / 255.0
    arr = (arr - IMAGENET_MEAN) / IMAGENET_STD
    return np.transpose(arr, (2, 0, 1))


def collect_tiles_rosie(args, he_array: np.ndarray, mpp: float) -> list[dict]:
    """Tile the H&E (same tissue tiles as the GigaTIME path) and run ROSIE; return tiles with mean_<channel>."""
    import torch

    scripts_dir = str(Path(__file__).resolve().parent)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import run_gigatime_tcga_brca as gigarun
    import validate_gigatime_xenium_rna as xrna
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    device = resolve_device(torch, args.device)
    print(f"ROSIE device: {device}", file=sys.stderr)
    model = load_rosie_model(torch, str(args.rosie_weights), device)

    H, W = he_array.shape[:2]
    patch_native = max(ROSIE_PATCH_PX, int(round(ROSIE_FIELD_UM / mpp)))
    ts = args.tile_size
    g = max(1, args.rosie_grid)
    offsets = [int(round((i + 0.5) / g * ts)) for i in range(g)]
    print(f"ROSIE: mpp={mpp:.4f} -> physical patch {patch_native}px (={ROSIE_FIELD_UM:.1f}um), {g}x{g} centers/tile", file=sys.stderr)

    slide = xrna.ArraySlide(Image, he_array)
    metas = [
        {"x": int(t["x"]), "y": int(t["y"]), "tissue_fraction": float(t["tissue_fraction"])}
        for t in gigarun.iter_tissue_tiles(
            slide, 0, ts, args.tile_stride, args.tissue_threshold, args.max_tiles, args.tile_order, args.random_seed,
        )
    ]
    n = len(metas)
    if n == 0:
        return []
    mapped = [(ch, ROSIE_CHANNELS.index(rname)) for ch, rname in CHANNEL_MAP.items() if rname in ROSIE_CHANNELS]

    sums = np.zeros((n, ROSIE_NUM_OUTPUTS), dtype=np.float64)
    counts = np.zeros(n, dtype=np.float64)
    flush_at = max(8, args.rosie_batch)
    buf: list[np.ndarray] = []
    idx: list[int] = []

    def flush():
        nonlocal buf, idx
        if not buf:
            return
        x = torch.from_numpy(np.stack(buf)).to(device)
        with torch.no_grad():
            out = model(x).float().cpu().numpy()
        for j, ti in enumerate(idx):
            sums[ti] += out[j]
            counts[ti] += 1.0
        buf, idx = [], []

    for ti, m in enumerate(metas):
        for oy in offsets:
            for ox in offsets:
                crop = _extract(he_array, m["x"] + ox, m["y"] + oy, patch_native, H, W)
                buf.append(_prep_patch(crop))
                idx.append(ti)
                if len(buf) >= flush_at:
                    flush()
        if (ti + 1) % 2000 == 0:
            print(f"  ROSIE tiles {ti + 1}/{n}", file=sys.stderr)
    flush()

    tiles: list[dict] = []
    for ti, m in enumerate(metas):
        mean50 = sums[ti] / max(counts[ti], 1.0)
        d = {"x": m["x"], "y": m["y"], "tissue_fraction": m["tissue_fraction"]}
        for ch, cidx in mapped:
            d[f"mean_{ch}"] = float(mean50[cidx])
        tiles.append(d)
    return tiles
