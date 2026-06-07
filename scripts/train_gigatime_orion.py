#!/usr/bin/env python3
"""In-domain fine-tune test: does training GigaTIME on Orion CRC rescue per-channel mIF specificity?

The decisive step of the CRC-Orion controlled experiment. Trains the GigaTIME UNet++ (warm-started
from the released lung weights, with a fresh 9-channel head) on Orion CRC H&E -> cell-level CyCIF
targets, then evaluates per-channel specificity on a HELD-OUT specimen, to compare against the
released (out-of-domain) baseline. If in-domain training rescues CD68/PD-L1/CD8 specificity -> the
ceiling is domain shift (breast training would help); if not -> intrinsic.

Targets are built WITHOUT the 100+ GB raw image: per tile, each Otsu-gated marker-positive cell is
painted as a filled disk (radius from its Area) at its centroid -> per-channel binary map, then
coarsened (cv2 area-resize /8) exactly as GigaTIME coarsens its targets.

Run (gigatime-tcga env), e.g.:
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/train_gigatime_orion.py \
      --train CRC02 CRC03 CRC04 CRC05 CRC06 --heldout CRC01 \
      --max-train-tiles 800 --epochs 40
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path("external/GigaTIME/scripts").resolve()))
import validate_gigatime_xenium_rna as xrna  # noqa: E402
import validate_gigatime_orion as vorion  # noqa: E402

CHANNELS = list(vorion.ORION_TO_GIGATIME.values())  # CD3,CD8,CD4,CD20,CD68,PD-L1,PD-1,Ki67,CK
IMAGENET_MEAN = np.array([0.485, 0.456, 0.406], np.float32)
IMAGENET_STD = np.array([0.229, 0.224, 0.225], np.float32)


def _rgb(x) -> np.ndarray:
    """Coerce an openslide/PIL tile to HxWx3 uint8."""
    a = np.asarray(x)
    if a.ndim == 2:
        a = np.stack([a] * 3, axis=-1)
    return np.ascontiguousarray(a[:, :, :3]).astype(np.uint8)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--train", nargs="+", default=["CRC02", "CRC03", "CRC04", "CRC05", "CRC06"])
    p.add_argument("--heldout", default="CRC01")
    p.add_argument("--orion-dir", type=Path, default=Path("data/orion_crc"))
    p.add_argument("--tile-size", type=int, default=256)
    p.add_argument("--coarse", type=int, default=32, help="Target coarse resolution (256/8=32, GigaTIME-style).")
    p.add_argument("--max-train-tiles", type=int, default=800, help="Tissue tiles sampled per training specimen.")
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--scratch", action="store_true", help="Train from scratch instead of warm-starting released weights.")
    p.add_argument("--tissue-threshold", type=float, default=0.35)
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--block-tiles", type=int, default=4)
    p.add_argument("--min-cells-per-channel", type=int, default=50)
    p.add_argument("--eval-max-tiles", type=int, default=3000)
    p.add_argument("--out-dir", type=Path, default=Path("results/gigatime_orion_finetune"))
    return p.parse_args()


def resolve_device(torch, device):
    if device == "auto":
        return torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    return torch.device(device)


def gate_cells(cells_csv: Path):
    """Per-cell centroid, disk radius, and per-channel Otsu-positive boolean from the single-cell table."""
    import pandas as pd
    from skimage.filters import threshold_otsu

    df = pd.read_csv(cells_csv)
    cx = df["X_centroid"].to_numpy(float)
    cy = df["Y_centroid"].to_numpy(float)
    radius = np.clip(np.sqrt(df["Area"].to_numpy(float) / np.pi), 2, 12)
    pos = {}
    for orion_marker, channel in vorion.ORION_TO_GIGATIME.items():
        v = np.log1p(df[orion_marker].to_numpy(float))
        try:
            thr = float(threshold_otsu(v))
        except Exception:
            thr = float(np.quantile(v, 0.9))
        pos[channel] = v > thr
    return cx, cy, radius, pos


def build_coarse_target(cx, cy, radius, pos, x0, y0, T, coarse):
    """(C, coarse, coarse) float target: paint positive cells as disks at full res, area-resize /8."""
    import cv2

    sel = (cx >= x0) & (cx < x0 + T) & (cy >= y0) & (cy < y0 + T)
    idx = np.nonzero(sel)[0]
    out = np.zeros((len(CHANNELS), coarse, coarse), np.float32)
    if len(idx) == 0:
        return out
    lx = (cx[idx] - x0).astype(int)
    ly = (cy[idx] - y0).astype(int)
    rr = radius[idx].astype(int)
    for ci, ch in enumerate(CHANNELS):
        m = np.zeros((T, T), np.uint8)
        chan_pos = pos[ch][idx]
        for k in np.nonzero(chan_pos)[0]:
            cv2.circle(m, (lx[k], ly[k]), int(rr[k]), 1, -1)
        out[ci] = cv2.resize(m.astype(np.float32), (coarse, coarse), interpolation=cv2.INTER_AREA)
    return out


def precompute_training(args, gigarun):
    """Sample tissue tiles across training specimens; return (H&E uint8 NxTxTx3, target NxCxcoarsexcoarse)."""
    import openslide

    he_list, tgt_list = [], []
    for spec in args.train:
        sd = args.orion_dir / spec
        he_file = next(iter(sorted(sd.glob("*-registered.ome.tif"))), sd / "he.ome.tif")
        cells_file = next((c for c in sorted(sd.glob("*.csv")) if c.name != "markers.csv"), sd / "cells.csv")
        if not he_file.exists() or not cells_file.exists():
            print(f"  skip {spec}: missing files", file=sys.stderr)
            continue
        cx, cy, radius, pos = gate_cells(cells_file)
        slide = openslide.OpenSlide(str(he_file))
        n = 0
        for tile in gigarun.iter_tissue_tiles(slide, 0, args.tile_size, args.tile_size,
                                              args.tissue_threshold, args.max_train_tiles, "random", args.random_seed):
            rgb = _rgb(tile["rgb"])
            if rgb.shape[:2] != (args.tile_size, args.tile_size):
                continue
            tgt = build_coarse_target(cx, cy, radius, pos, tile["x"], tile["y"], args.tile_size, args.coarse)
            if tgt.sum() <= 0:
                continue
            he_list.append(rgb)
            tgt_list.append(tgt)
            n += 1
        slide.close()
        print(f"  {spec}: {n} training tiles", file=sys.stderr)
    return np.stack(he_list), np.stack(tgt_list)


def load_warmstart_state(torch):
    from huggingface_hub import snapshot_download

    local = snapshot_download(repo_id="prov-gigatime/GigaTIME")
    sd = torch.load(Path(local) / "model.pth", map_location="cpu")
    # Drop the 23-channel final head so the fresh 9-channel head is kept.
    return {k: v for k, v in sd.items() if not k.startswith("final.")}


def train_model(args, torch, he, tgt, device):
    import torch.nn.functional as F
    from torch.utils.data import DataLoader, TensorDataset
    import archs
    import losses

    model = archs.gigatime(num_classes=len(CHANNELS), input_channels=3)
    if not args.scratch:
        missing, unexpected = model.load_state_dict(load_warmstart_state(torch), strict=False)
        print(f"warm-start: missing {len(missing)} (incl. final head), unexpected {len(unexpected)}", file=sys.stderr)
    model.to(device).train()
    criterion = losses.BCEDiceLoss()
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs, eta_min=1e-6)

    he_t = torch.from_numpy(he)          # N,T,T,3 uint8
    tgt_t = torch.from_numpy(tgt)        # N,C,coarse,coarse float
    loader = DataLoader(TensorDataset(he_t, tgt_t), batch_size=args.batch_size, shuffle=True, drop_last=True)
    mean = torch.tensor(IMAGENET_MEAN).view(1, 3, 1, 1).to(device)
    std = torch.tensor(IMAGENET_STD).view(1, 3, 1, 1).to(device)

    for ep in range(args.epochs):
        tot = 0.0
        for hb, tb in loader:
            x = hb.to(device).permute(0, 3, 1, 2).contiguous().float() / 255.0
            x = (x - mean) / std
            target = F.interpolate(tb.to(device).float(), size=x.shape[2:], mode="bilinear", align_corners=False)
            out = model(x)  # native 256x256, matching GigaTIME inference (preprocess_tile does no resize)
            loss = criterion(out, target)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item()
        sched.step()
        if ep % 5 == 0 or ep == args.epochs - 1:
            print(f"  epoch {ep+1}/{args.epochs} loss {tot/max(1,len(loader)):.4f}", file=sys.stderr)
    model.eval()
    return model


def eval_heldout(args, torch, model, gigarun, device, spec):
    """Per-channel specificity of the trained model on a held-out specimen (reuses the audit stats)."""
    import torch.nn.functional as F
    import openslide

    sd = args.orion_dir / spec
    he_file = next(iter(sorted(sd.glob("*-registered.ome.tif"))), sd / "he.ome.tif")
    cells_file = next((c for c in sorted(sd.glob("*.csv")) if c.name != "markers.csv"), sd / "cells.csv")
    slide = openslide.OpenSlide(str(he_file))
    W, H = slide.level_dimensions[0]
    stride = args.tile_size
    n_cols, n_rows = W // stride + 2, H // stride + 2
    total_grid, gene_grids, n_cells, _gate = vorion.load_cell_grids(cells_file, stride, n_rows, n_cols, W, H)

    mean = torch.tensor(IMAGENET_MEAN).view(1, 3, 1, 1).to(device)
    std = torch.tensor(IMAGENET_STD).view(1, 3, 1, 1).to(device)
    tiles = []
    batch, meta = [], []

    def flush():
        nonlocal batch, meta
        if not batch:
            return
        x = torch.from_numpy(np.stack(batch)).to(device).permute(0, 3, 1, 2).contiguous().float() / 255.0
        x = (x - mean) / std
        with torch.no_grad():
            probs = torch.sigmoid(model(x)).mean(dim=(2, 3)).cpu().numpy()  # N x C (native 256)
        for row, p in zip(meta, probs):
            d = {"x": row["x"], "y": row["y"], "tissue_fraction": row["tissue_fraction"]}
            for ci, ch in enumerate(CHANNELS):
                d[f"mean_{ch}"] = float(p[ci])
            tiles.append(d)
        batch, meta = [], []

    for tile in gigarun.iter_tissue_tiles(slide, 0, stride, stride, args.tissue_threshold,
                                          args.eval_max_tiles, "random", args.random_seed):
        rgb = _rgb(tile["rgb"])
        if rgb.shape[:2] != (stride, stride):
            continue
        batch.append(rgb)
        meta.append({"x": int(tile["x"]), "y": int(tile["y"]), "tissue_fraction": float(tile["tissue_fraction"])})
        if len(batch) == args.batch_size:
            flush()
    flush()
    slide.close()

    tcol = np.array([t["x"] // stride for t in tiles], np.int64)
    trow = np.array([t["y"] // stride for t in tiles], np.int64)
    total_counts = total_grid[trow, tcol].astype(float)
    occ = total_counts > 0
    tiles = [t for t, k in zip(tiles, occ) if k]
    tcol, trow, total_counts = tcol[occ], trow[occ], total_counts[occ]

    import os
    import matplotlib
    os.environ.setdefault("MPLCONFIGDIR", str(args.out_dir / ".matplotlib"))
    (args.out_dir / ".matplotlib").mkdir(parents=True, exist_ok=True)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args.asset_dir = Path("docs/assets/gigatime_orion_finetune")
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    stats = xrna.require_stats()
    rng = np.random.default_rng(args.random_seed)
    vbc, dbc = {}, {}
    for ch, grid in gene_grids.items():
        if f"mean_{ch}" not in tiles[0]:
            continue
        counts = grid[trow, tcol].astype(float)
        if counts.sum() < args.min_cells_per_channel:
            continue
        vbc[ch] = np.array([float(t[f"mean_{ch}"]) for t in tiles], float)
        dbc[ch] = counts
    spec_res = xrna.compute_specificity(stats, rng, plt, args, vbc, dbc, total_counts, tcol, trow)
    return spec_res, len(tiles)


def main():
    args = parse_args()
    import torch
    import run_gigatime_tcga_brca as gigarun

    args.out_dir.mkdir(parents=True, exist_ok=True)
    device = resolve_device(torch, args.device)
    print(f"device {device}; channels {CHANNELS}", file=sys.stderr)

    print("Precomputing training tiles...", file=sys.stderr)
    he, tgt = precompute_training(args, gigarun)
    print(f"training set: {he.shape[0]} tiles", file=sys.stderr)

    tag = "scratch" if args.scratch else "finetune"
    model = train_model(args, torch, he, tgt, device)
    ckpt = args.out_dir / f"gigatime_orion_{tag}.pth"
    torch.save(model.state_dict(), ckpt)
    print(f"saved {ckpt}", file=sys.stderr)

    print(f"Evaluating on held-out {args.heldout}...", file=sys.stderr)
    spec_res, n_eval = eval_heldout(args, torch, model, gigarun, device, args.heldout)

    report = {"tag": tag, "train": args.train, "heldout": args.heldout, "channels": CHANNELS,
              "n_train_tiles": int(he.shape[0]), "n_eval_tiles": n_eval, "specificity": spec_res}
    (args.out_dir / f"orion_{tag}_heldout_{args.heldout}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\n=== {tag} model on held-out {args.heldout} (in-domain, protein-vs-protein) ===", file=sys.stderr)
    if spec_res:
        print(f"own-marker row-max {spec_res['n_own_is_row_max']}/{spec_res['n_channels']}; "
              f"partial r>0 {spec_res['n_partial_survives']}/{spec_res['n_channels']}", file=sys.stderr)
        for r in spec_res["per_channel"]:
            print(f"  {r['channel']:6s} partial={r['partial_r_control_total']:.3f} rowmax={r['own_is_row_max']}", file=sys.stderr)
    print(f"Wrote {args.out_dir / f'orion_{tag}_heldout_{args.heldout}.json'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
