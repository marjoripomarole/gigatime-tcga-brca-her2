#!/usr/bin/env python3
"""In-domain (CRC) specificity audit of a virtual-mIF model against Orion measured protein.

Tests whether GigaTIME's marker-channel specificity is better IN-DOMAIN and IN-MODALITY (protein
vs protein) than it was out-of-domain on breast vs RNA. Uses the Orion CRC dataset (Lin et al.
2023; same-section H&E + 17-plex CyCIF). Crucially, we DON'T need the 100+ GB raw CyCIF image: the
release ships a cell-segmentation single-cell table with per-cell marker mean intensities + X/Y
centroids, so we reconstruct cell-level marker maps directly. Per 256 px H&E tile we bin
marker-positive cell centroids (Otsu gate on log intensity) -> per-channel cell density, with total
cells/tile as the cellularity control -- exactly the within-slide specificity design used for the
Xenium/Visium RNA audit, reusing that stats core verbatim.

H&E is read region-wise via openslide (the full image is ~13 GB decompressed).

Run (gigatime-tcga env):
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/validate_gigatime_orion.py --specimen CRC01 --max-tiles 3000
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_gigatime_xenium_rna as xrna  # noqa: E402

# Orion marker (single-cell-table column) -> GigaTIME channel name.
ORION_TO_GIGATIME = {
    "CD3e": "CD3", "CD8a": "CD8", "CD4": "CD4", "CD20": "CD20", "CD68": "CD68",
    "PD-L1": "PD-L1", "PD-1": "PD-1", "Ki67": "Ki67", "Pan-CK": "CK",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--specimen", default="CRC01")
    p.add_argument("--orion-dir", type=Path, default=Path("data/orion_crc"))
    p.add_argument("--he-file", type=Path, default=None, help="Defaults to the *-registered.ome.tif in the specimen dir.")
    p.add_argument("--cells-file", type=Path, default=None, help="Defaults to the single-cell .csv in the specimen dir.")
    p.add_argument("--gigatime-repo", default="external/GigaTIME")
    p.add_argument("--tile-size", type=int, default=256)
    p.add_argument("--tile-stride", type=int, default=256)
    p.add_argument("--tissue-threshold", type=float, default=0.35)
    p.add_argument("--max-tiles", type=int, default=3000, help="Cap tissue tiles for a pilot (0 = all).")
    p.add_argument("--tile-order", default="random", choices=["row-major", "random"])
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--block-tiles", type=int, default=4)
    p.add_argument("--min-cells-per-channel", type=int, default=50)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--out-markdown", type=Path, default=None)
    args = p.parse_args()
    sd = args.orion_dir / args.specimen
    if args.he_file is None:
        cands = sorted(sd.glob("*-registered.ome.tif"))
        args.he_file = cands[0] if cands else sd / "he.ome.tif"
    if args.cells_file is None:
        cands = [c for c in sorted(sd.glob("*.csv")) if c.name != "markers.csv"]
        args.cells_file = cands[0] if cands else sd / "cells.csv"
    args.out_dir = args.out_dir or Path("results/gigatime_orion_baseline") / args.specimen
    args.out_markdown = args.out_markdown or Path("docs") / f"gigatime_orion_baseline_{args.specimen}.md"
    return args


def load_cell_grids(cells_csv: Path, stride: int, n_rows: int, n_cols: int, width: int, height: int):
    """Per-channel marker-positive cell-count grids + total-cell grid, from the single-cell table."""
    import pandas as pd
    from skimage.filters import threshold_otsu

    df = pd.read_csv(cells_csv)
    cx = df["X_centroid"].to_numpy(dtype=float)
    cy = df["Y_centroid"].to_numpy(dtype=float)
    total_grid = xrna.bin_grid_counts(cx, cy, stride, n_rows, n_cols, width, height)

    gene_grids: dict[str, np.ndarray] = {}
    gate_info: dict[str, dict] = {}
    for orion_marker, channel in ORION_TO_GIGATIME.items():
        if orion_marker not in df.columns:
            continue
        vals = np.log1p(df[orion_marker].to_numpy(dtype=float))
        try:
            thr = float(threshold_otsu(vals))
        except Exception:
            thr = float(np.quantile(vals, 0.9))
        pos = vals > thr
        if pos.sum() < 1:
            continue
        gene_grids[channel] = xrna.bin_grid_counts(cx[pos], cy[pos], stride, n_rows, n_cols, width, height)
        gate_info[channel] = {"orion_marker": orion_marker, "n_positive": int(pos.sum()),
                              "frac_positive": round(float(pos.mean()), 4)}
    return total_grid, gene_grids, int(len(df)), gate_info


def collect_tiles_gigatime(args, slide):
    """Tile the Orion H&E (openslide) and run GigaTIME inference -> tiles with mean_<channel>."""
    import run_gigatime_tcga_brca as gigarun

    torch, gigatime_class, snapshot_download, Image, _openslide = gigarun.import_runtime(Path(args.gigatime_repo))
    device = gigarun.resolve_device(torch, args.device)
    print(f"GigaTIME device: {device}", file=sys.stderr)
    model = gigarun.load_model(torch, gigatime_class, snapshot_download, device)
    tiles, batch, meta = [], [], []
    with torch.no_grad():
        for tile in gigarun.iter_tissue_tiles(slide, 0, args.tile_size, args.tile_stride,
                                              args.tissue_threshold, args.max_tiles, args.tile_order, args.random_seed):
            batch.append(gigarun.preprocess_tile(torch, tile["rgb"], device))
            meta.append({"x": int(tile["x"]), "y": int(tile["y"]), "tissue_fraction": float(tile["tissue_fraction"])})
            if len(batch) == args.batch_size:
                tiles.extend(gigarun.infer_batch(torch, model, batch, meta, 0.5))
                batch, meta = [], []
        if batch:
            tiles.extend(gigarun.infer_batch(torch, model, batch, meta, 0.5))
    return tiles


def main() -> int:
    args = parse_args()
    stats = xrna.require_stats()
    if not args.he_file.exists() or not args.cells_file.exists():
        raise SystemExit(f"Missing H&E ({args.he_file}) or cells ({args.cells_file}).")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)

    import openslide
    slide = openslide.OpenSlide(str(args.he_file))
    width, height = slide.level_dimensions[0]
    stride = args.tile_stride
    n_cols, n_rows = width // stride + 2, height // stride + 2
    print(f"H&E {width}x{height}; binning cells...", file=sys.stderr)
    total_grid, gene_grids, n_cells, gate = load_cell_grids(args.cells_file, stride, n_rows, n_cols, width, height)
    print(f"{n_cells:,} cells; channels: {list(gene_grids)}", file=sys.stderr)

    print("Tiling H&E + GigaTIME inference...", file=sys.stderr)
    tiles = collect_tiles_gigatime(args, slide)
    if not tiles:
        raise SystemExit("No tissue tiles found.")
    tcol = np.array([t["x"] // stride for t in tiles], dtype=np.int64)
    trow = np.array([t["y"] // stride for t in tiles], dtype=np.int64)
    total_counts = total_grid[trow, tcol].astype(float)
    # Restrict to tiles that contain >=1 segmented cell (so the protein target is defined).
    occ = total_counts > 0
    tiles = [t for t, k in zip(tiles, occ) if k]
    tcol, trow, total_counts = tcol[occ], trow[occ], total_counts[occ]
    print(f"Tiles with >=1 cell: {len(tiles)}", file=sys.stderr)

    import os
    mpl_dir = args.out_dir / ".matplotlib"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    asset_dir = Path("docs/assets") / f"gigatime_orion_baseline_{args.specimen}"
    asset_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir = asset_dir

    channel_results, virtual_by_channel, density_by_channel = [], {}, {}
    for channel, grid in gene_grids.items():
        mean_key = f"mean_{channel}"
        if mean_key not in tiles[0]:
            continue
        counts = grid[trow, tcol].astype(float)
        if counts.sum() < args.min_cells_per_channel:
            continue
        virtual = np.array([float(t[mean_key]) for t in tiles], dtype=float)
        r, pval = xrna.spearman(stats, virtual, counts)
        lo, hi = xrna.block_bootstrap_ci(stats, rng, tcol, trow, virtual, counts, args.block_tiles, args.bootstrap)
        virtual_by_channel[channel] = virtual
        density_by_channel[channel] = counts
        channel_results.append({"channel": channel, "orion_marker": gate.get(channel, {}).get("orion_marker"),
                                "n_positive_cells_on_grid": int(counts.sum()), "spearman_r": r, "spearman_p": pval,
                                "ci95_low": lo, "ci95_high": hi})
    channel_results.sort(key=lambda d: np.nan_to_num(d["spearman_r"], nan=-9), reverse=True)
    specificity = xrna.compute_specificity(stats, rng, plt, args, virtual_by_channel, density_by_channel,
                                           total_counts, tcol, trow)

    report = {
        "specimen": args.specimen, "model": "gigatime", "domain": "CRC (Orion)", "modality": "protein (CyCIF cells)",
        "he_shape": [height, width], "n_cells": n_cells, "n_tiles": len(tiles),
        "tile_size": args.tile_size, "channel_results": channel_results, "specificity": specificity, "gate": gate,
    }
    (args.out_dir / "orion_baseline_report.json").write_text(json.dumps(report, indent=2) + "\n")
    _write_md(args.out_markdown, report)
    print(f"\nReleased GigaTIME on CRC {args.specimen} (in-domain, protein-vs-protein):", file=sys.stderr)
    for r in channel_results:
        print(f"  {r['channel']:6s} raw={r['spearman_r']:.3f}", file=sys.stderr)
    if specificity:
        print(f"own-gene row-max {specificity['n_own_is_row_max']}/{specificity['n_channels']}; "
              f"partial r>0 {specificity['n_partial_survives']}/{specificity['n_channels']}", file=sys.stderr)
        for r in specificity["per_channel"]:
            print(f"  {r['channel']:6s} partial={r['partial_r_control_total']:.3f} "
                  f"CI[{r['partial_ci95_low']:.3f},{r['partial_ci95_high']:.3f}] rowmax={r['own_is_row_max']}", file=sys.stderr)
    print(f"Wrote {args.out_markdown}", file=sys.stderr)
    return 0


def _write_md(path: Path, report: dict) -> None:
    spec = report.get("specificity") or {}
    lines = [
        f"# Orion CRC in-domain specificity baseline — {report['specimen']} (released GigaTIME)",
        "",
        f"Status: in-domain (CRC), in-modality (protein-vs-protein) specificity audit of the RELEASED (lung-trained) "
        f"GigaTIME on Orion CRC. Per 256px H&E tile, virtual channel activation is correlated against the density of "
        f"marker-positive cells (Otsu-gated single-cell intensities), with total cells/tile as the cellularity control "
        f"— the same audit used for the breast RNA work, reusing its stats core.",
        "",
        f"- H&E: {report['he_shape'][0]} x {report['he_shape'][1]} px (0.325 um/px); {report['n_cells']:,} segmented cells; "
        f"{report['n_tiles']} tiles (>=1 cell) used.",
        "",
        "## Per-channel specificity (cellularity-controlled)",
        "",
        "| Channel (Orion marker) | raw Spearman | partial r \\| total cells | partial 95% CI | own-marker row-max? |",
        "|---|---:|---:|---|:--:|",
    ]
    per = {r["channel"]: r for r in spec.get("per_channel", [])}
    for r in report["channel_results"]:
        sp = per.get(r["channel"], {})
        ci = f"[{sp.get('partial_ci95_low', float('nan')):.3f}, {sp.get('partial_ci95_high', float('nan')):.3f}]"
        lines.append(f"| {r['channel']} ({r.get('orion_marker')}) | {r['spearman_r']:.3f} | "
                     f"{sp.get('partial_r_control_total', float('nan')):.3f} | {ci} | "
                     f"{'yes' if sp.get('own_is_row_max') else 'no'} |")
    if spec:
        lines += ["", f"Own-marker is the row-max for **{spec['n_own_is_row_max']}/{spec['n_channels']}** channels; "
                  f"cellularity-controlled partial r stays >0 for **{spec['n_partial_survives']}/{spec['n_channels']}**.", ""]
    lines += ["## Output Files", "", f"- `{report['specimen']}` baseline JSON in `results/gigatime_orion_baseline/`", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
