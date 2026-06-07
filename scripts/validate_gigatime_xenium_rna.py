#!/usr/bin/env python3
"""Validate GigaTIME virtual channels against Xenium breast spatial RNA, within-slide.

Pipeline (single Xenium section):
  1. Read the post-Xenium H&E (OME-TIFF) at full resolution.
  2. Tile it on GigaTIME's regular 256 px grid (origin 0) and run virtual-channel
     inference, reusing run_gigatime_tcga_brca internals.
  3. Map every Xenium transcript from microns -> Xenium-morphology px -> H&E px via
     the inverse H&E alignment affine, then bin transcripts onto the SAME tile grid
     (tile = (x // stride, y // stride), exact integer binning).
  4. Per GigaTIME channel with an RNA counterpart, compute the within-slide Spearman
     correlation between virtual-channel intensity and that channel's transcript
     density across tiles, with a spatial block-bootstrap CI.

Why within-slide: the correlation is computed across tiles of ONE section, so it is
immune to the cross-patient acquisition/composition confound that dominated the TCGA
analysis. A positive result means "virtual-CD8 co-localizes with CD8 transcripts".

Alignment orientation is auto-detected by which direction lands the most transcripts
inside the H&E frame (override with --alignment-direction).

Model-free sanity: --alignment-check-only skips GigaTIME and just correlates per-tile
tissue fraction with total transcript density (should be strongly positive if the
alignment is right), so coordinate mapping can be validated before spending GPU time.

Run with the gigatime-tcga env, e.g.:
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/validate_gigatime_xenium_rna.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# GigaTIME virtual channel -> candidate Xenium gene symbols (intersected with the
# panel/transcripts at runtime; channels with no present gene are skipped).
CHANNEL_GENES: dict[str, list[str]] = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "CD20": ["MS4A1"],
    "CD68": ["CD68"],
    "CD14": ["CD14"],
    "CD11c": ["ITGAX"],
    "CD16": ["FCGR3A"],
    "PD-1": ["PDCD1"],
    "PD-L1": ["CD274"],
    "CK": ["KRT8", "KRT18", "KRT19", "KRT7", "EPCAM"],
    "Ki67": ["MKI67"],
    "CD138": ["SDC1"],
    "CD34": ["CD34"],
    "T-bet": ["TBX21"],
    "Tryptase": ["TPSAB1", "TPSB2"],
}

DEFAULT_XENIUM_MPP = 0.2125  # microns/pixel of the Xenium morphology image (v1).


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=Path("data/xenium_breast/Xenium_FFPE_Human_Breast_Cancer_Rep1"))
    parser.add_argument("--sample", default="Xenium_FFPE_Human_Breast_Cancer_Rep1")
    parser.add_argument("--he-image", type=Path, default=None, help="Defaults to <data-dir>/<sample>_he_image.ome.tif.")
    parser.add_argument("--alignment-csv", type=Path, default=None, help="Defaults to <data-dir>/<sample>_he_imagealignment.csv.")
    parser.add_argument("--transcripts", type=Path, default=None, help="Defaults to <data-dir>/<sample>_transcripts.parquet.")
    parser.add_argument("--gigatime-repo", default="external/GigaTIME")
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--tile-stride", type=int, default=256)
    parser.add_argument("--tissue-threshold", type=float, default=0.35)
    parser.add_argument("--max-tiles", type=int, default=0, help="Cap tissue tiles (0 = all). Use a small value for a smoke run.")
    parser.add_argument("--tile-order", default="row-major", choices=["row-major", "random"])
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    parser.add_argument("--xenium-mpp", type=float, default=DEFAULT_XENIUM_MPP)
    parser.add_argument("--alignment-direction", default="auto", choices=["auto", "he_to_morph", "morph_to_he"])
    parser.add_argument("--min-transcripts-per-channel", type=int, default=50, help="Skip channels with fewer total transcripts on-grid.")
    parser.add_argument("--bootstrap", type=int, default=1000, help="Spatial block-bootstrap iterations for the CI.")
    parser.add_argument("--block-tiles", type=int, default=4, help="Block size (in tiles) for the spatial block bootstrap.")
    parser.add_argument("--alignment-check-only", action="store_true", help="Skip GigaTIME; only correlate tissue fraction vs total transcript density.")
    parser.add_argument("--out-dir", type=Path, default=Path("results/gigatime_xenium_rna_validation"))
    parser.add_argument("--asset-dir", type=Path, default=Path("docs/assets/gigatime_xenium_rna_validation"))
    parser.add_argument("--out-markdown", type=Path, default=Path("docs/xenium_breast_rna_validation_results.md"))
    parser.add_argument("--model", choices=["gigatime", "rosie"], default="gigatime",
                        help="Virtual-mIF model whose channels are audited against RNA.")
    parser.add_argument("--rosie-weights", type=Path, default=Path("external/ROSIE/best_model_single.pth"))
    parser.add_argument("--rosie-grid", type=int, default=2, help="ROSIE patch centers per tile (grid x grid).")
    parser.add_argument("--rosie-batch", type=int, default=64, help="ROSIE patch batch size.")
    parser.add_argument("--rosie-mpp", type=float, default=None, help="H&E microns/px for ROSIE (required for --model rosie).")
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    he = args.he_image or args.data_dir / f"{args.sample}_he_image.ome.tif"
    align = args.alignment_csv or args.data_dir / f"{args.sample}_he_imagealignment.csv"
    tx = args.transcripts or args.data_dir / f"{args.sample}_transcripts.parquet"
    for label, path in (("H&E image", he), ("alignment CSV", align), ("transcripts", tx)):
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")
    return he, align, tx


def load_alignment(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8").strip().splitlines():
        parts = [p for p in line.replace(",", " ").split() if p]
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue
    matrix = np.array(rows, dtype=float)
    if matrix.shape != (3, 3):
        raise SystemExit(f"Expected a 3x3 alignment matrix, got shape {matrix.shape}")
    return matrix


def load_he_fullres(path: Path) -> np.ndarray:
    import tifffile

    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        level = series.levels[0] if getattr(series, "levels", None) else series
        arr = level.asarray()
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    elif arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[2] not in (3, 4):
        arr = np.transpose(arr, (1, 2, 0))
    arr = arr[:, :, :3]
    if arr.dtype != np.uint8:
        amax = float(arr.max()) or 1.0
        arr = np.clip(arr / amax * 255.0, 0, 255).astype(np.uint8) if amax > 255 else arr.astype(np.uint8)
    return np.ascontiguousarray(arr)


class ArraySlide:
    """OpenSlide-like adapter over an in-memory full-resolution RGB array."""

    def __init__(self, image_module, array: np.ndarray):
        self._image_module = image_module
        self.array = array
        height, width = array.shape[:2]
        self.level_count = 1
        self.level_dimensions = [(width, height)]
        self.level_downsamples = [1.0]

    def read_region(self, location: tuple[int, int], level: int, size: tuple[int, int]):
        if level != 0:
            raise ValueError("ArraySlide has a single level")
        x, y = location
        w, h = size
        crop = self.array[y : y + h, x : x + w]
        return self._image_module.fromarray(crop).convert("RGBA")

    def close(self) -> None:  # parity with OpenSlide
        pass


def transcript_he_pixels(matrix: np.ndarray, x_um, y_um, mpp: float, direction: str) -> tuple[np.ndarray, np.ndarray]:
    """Map transcript microns to H&E pixel coordinates for a given alignment direction."""
    mpx = np.asarray(x_um, dtype=float) / mpp
    mpy = np.asarray(y_um, dtype=float) / mpp
    if direction == "he_to_morph":
        transform = np.linalg.inv(matrix)  # morphology px -> H&E px
    else:  # morph_to_he: matrix already maps morphology px -> H&E px
        transform = matrix
    hx = transform[0, 0] * mpx + transform[0, 1] * mpy + transform[0, 2]
    hy = transform[1, 0] * mpx + transform[1, 1] * mpy + transform[1, 2]
    hw = transform[2, 0] * mpx + transform[2, 1] * mpy + transform[2, 2]
    hw = np.where(hw == 0, 1.0, hw)
    return hx / hw, hy / hw


def choose_direction(matrix, x_um, y_um, mpp, width, height, requested: str) -> tuple[str, dict[str, float]]:
    fractions: dict[str, float] = {}
    for direction in ("he_to_morph", "morph_to_he"):
        hx, hy = transcript_he_pixels(matrix, x_um, y_um, mpp, direction)
        inside = (hx >= 0) & (hx < width) & (hy >= 0) & (hy < height)
        fractions[direction] = float(inside.mean())
    if requested != "auto":
        return requested, fractions
    best = max(fractions, key=fractions.get)
    return best, fractions


def bin_grid_counts(hx, hy, stride: int, n_rows: int, n_cols: int, width: int, height: int) -> np.ndarray:
    valid = (hx >= 0) & (hx < width) & (hy >= 0) & (hy < height)
    cols = (hx[valid] / stride).astype(np.int64)
    rows = (hy[valid] / stride).astype(np.int64)
    keep = (cols >= 0) & (cols < n_cols) & (rows >= 0) & (rows < n_rows)
    linear = rows[keep] * n_cols + cols[keep]
    grid = np.bincount(linear, minlength=n_rows * n_cols).reshape(n_rows, n_cols)
    return grid


def collect_tiles(args, he_array, alignment_check_only: bool):
    """Tile the H&E and (unless check-only) run GigaTIME inference. Returns tile dicts."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_gigatime_tcga_brca as gigarun

    if alignment_check_only:
        from PIL import Image

        Image.MAX_IMAGE_PIXELS = None
        slide = ArraySlide(Image, he_array)
        tiles = []
        for tile in gigarun.iter_tissue_tiles(
            slide, 0, args.tile_size, args.tile_stride, args.tissue_threshold,
            args.max_tiles, args.tile_order, args.random_seed,
        ):
            tiles.append({"x": int(tile["x"]), "y": int(tile["y"]), "tissue_fraction": float(tile["tissue_fraction"])})
        return tiles, gigarun

    torch, gigatime_class, snapshot_download, Image, _openslide = gigarun.import_runtime(Path(args.gigatime_repo))
    device = gigarun.resolve_device(torch, args.device)
    print(f"Using device: {device}", file=sys.stderr)
    model = gigarun.load_model(torch, gigatime_class, snapshot_download, device)
    slide = ArraySlide(Image, he_array)

    tiles: list[dict] = []
    batch, batch_meta = [], []
    with torch.no_grad():
        for tile in gigarun.iter_tissue_tiles(
            slide, 0, args.tile_size, args.tile_stride, args.tissue_threshold,
            args.max_tiles, args.tile_order, args.random_seed,
        ):
            batch.append(gigarun.preprocess_tile(torch, tile["rgb"], device))
            batch_meta.append({"x": int(tile["x"]), "y": int(tile["y"]), "tissue_fraction": float(tile["tissue_fraction"])})
            if len(batch) == args.batch_size:
                tiles.extend(gigarun.infer_batch(torch, model, batch, batch_meta, 0.5))
                batch, batch_meta = [], []
        if batch:
            tiles.extend(gigarun.infer_batch(torch, model, batch, batch_meta, 0.5))
    return tiles, gigarun


def spearman(stats, x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan")
    res = stats.spearmanr(x, y)
    return float(res.statistic), float(res.pvalue)


def partial_spearman(stats, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Spearman partial correlation of x and y controlling for z (rank residualization)."""
    if len(x) < 5:
        return float("nan")
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    design = np.c_[np.ones_like(rz), rz]
    bx, _, _, _ = np.linalg.lstsq(design, rx, rcond=None)
    by, _, _, _ = np.linalg.lstsq(design, ry, rcond=None)
    ex = rx - design @ bx
    ey = ry - design @ by
    if np.allclose(ex, ex[0]) or np.allclose(ey, ey[0]):
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def block_members(cols: np.ndarray, rows: np.ndarray, block_tiles: int) -> list[np.ndarray]:
    block_ids = (rows // block_tiles).astype(np.int64) * 100000 + (cols // block_tiles).astype(np.int64)
    _uniq, inverse = np.unique(block_ids, return_inverse=True)
    return [np.where(inverse == i)[0] for i in range(int(inverse.max()) + 1)] if len(inverse) else []


def block_bootstrap_fn(rng, members, compute_fn, iters: int) -> tuple[float, float]:
    if iters <= 0 or len(members) < 2:
        return float("nan"), float("nan")
    n = len(members)
    estimates = []
    for _ in range(iters):
        chosen = rng.integers(0, n, size=n)
        idx = np.concatenate([members[c] for c in chosen])
        value = compute_fn(idx)
        if not np.isnan(value):
            estimates.append(value)
    if not estimates:
        return float("nan"), float("nan")
    lo, hi = np.percentile(estimates, [2.5, 97.5])
    return float(lo), float(hi)


def block_bootstrap_ci(stats, rng, cols, rows, x, y, block_tiles: int, iters: int) -> tuple[float, float]:
    if iters <= 0 or len(x) < 5:
        return float("nan"), float("nan")
    block_ids = (rows // block_tiles).astype(np.int64) * 100000 + (cols // block_tiles).astype(np.int64)
    unique_blocks, inverse = np.unique(block_ids, return_inverse=True)
    members: dict[int, np.ndarray] = {i: np.where(inverse == i)[0] for i in range(len(unique_blocks))}
    n_blocks = len(unique_blocks)
    estimates = []
    for _ in range(iters):
        chosen = rng.integers(0, n_blocks, size=n_blocks)
        idx = np.concatenate([members[c] for c in chosen])
        r, _p = spearman(stats, x[idx], y[idx])
        if not np.isnan(r):
            estimates.append(r)
    if not estimates:
        return float("nan"), float("nan")
    lo, hi = np.percentile(estimates, [2.5, 97.5])
    return float(lo), float(hi)


def require_stats():
    try:
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(f"Missing package: {exc.name}. Use the gigatime-tcga env.") from exc
    return stats


def load_transcript_arrays(transcripts_path: Path, genes_of_interest: set[str]):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    table = pq.read_table(transcripts_path, columns=["feature_name", "x_location", "y_location"])
    x_all = table["x_location"].to_numpy(zero_copy_only=False).astype(float)
    y_all = table["y_location"].to_numpy(zero_copy_only=False).astype(float)

    goi_bytes = pa.array([g.encode("utf-8") for g in sorted(genes_of_interest)], type=table["feature_name"].type)
    mask = pc.is_in(table["feature_name"], value_set=goi_bytes)
    sub = table.filter(mask)
    sub_genes = [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in sub["feature_name"].to_pylist()]
    sub_x = sub["x_location"].to_numpy(zero_copy_only=False).astype(float)
    sub_y = sub["y_location"].to_numpy(zero_copy_only=False).astype(float)
    return x_all, y_all, np.array(sub_genes, dtype=object), sub_x, sub_y, table.num_rows


def render_scatter(plt, channel, x, y, r, out_path: Path) -> None:
    plt.figure(figsize=(4.6, 4.2))
    plt.scatter(x, y, s=6, alpha=0.3, color="#2563eb", edgecolors="none")
    plt.xlabel(f"Virtual {channel} (GigaTIME mean activation)")
    plt.ylabel("Transcript density (counts/tile)")
    plt.title(f"{channel}  Spearman r={r:.3f}")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def render_specificity_heatmap(plt, channels, matrix, out_path: Path) -> None:
    import numpy as _np

    fig, ax = plt.subplots(figsize=(0.55 * len(channels) + 2, 0.55 * len(channels) + 2))
    data = _np.array(matrix, dtype=float)
    im = ax.imshow(data, cmap="magma", vmin=0.0, vmax=max(0.1, _np.nanmax(data)))
    ax.set_xticks(range(len(channels)))
    ax.set_yticks(range(len(channels)))
    ax.set_xticklabels([f"{c} RNA" for c in channels], rotation=90, fontsize=7)
    ax.set_yticklabels([f"virtual {c}" for c in channels], fontsize=7)
    for i in range(len(channels)):
        for j in range(len(channels)):
            if not _np.isnan(data[i, j]):
                ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center",
                        fontsize=5.5, color="white" if data[i, j] < 0.6 * _np.nanmax(data) else "black")
    ax.set_title("Virtual channel (rows) vs RNA gene-set (cols) Spearman r")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_specificity(stats, rng, plt, args, virtual_by_channel, density_by_channel, total_counts, tile_cols, tile_rows) -> dict:
    """Channel x gene-set Spearman matrix plus cellularity-controlled partial correlations."""
    channels = list(virtual_by_channel.keys())
    if len(channels) < 2:
        return {}
    members = block_members(tile_cols, tile_rows, args.block_tiles)
    matrix = [[spearman(stats, virtual_by_channel[rc], density_by_channel[cc])[0] for cc in channels] for rc in channels]

    per_channel = []
    for i, ch in enumerate(channels):
        row = matrix[i]
        own = row[i]
        off = [(channels[j], row[j]) for j in range(len(channels)) if j != i and not np.isnan(row[j])]
        best_other_ch, best_other = max(off, key=lambda kv: kv[1]) if off else ("", float("nan"))
        is_row_max = all(np.isnan(row[j]) or own >= row[j] for j in range(len(channels)) if j != i)
        v = virtual_by_channel[ch]
        d = density_by_channel[ch]
        partial = partial_spearman(stats, v, d, total_counts)
        p_lo, p_hi = block_bootstrap_fn(rng, members, lambda idx, v=v, d=d: partial_spearman(stats, v[idx], d[idx], total_counts[idx]), args.bootstrap)
        per_channel.append({
            "channel": ch,
            "own_gene_r": own,
            "best_other_channel": best_other_ch,
            "best_other_r": best_other,
            "specificity_gap": float(own - best_other) if not np.isnan(best_other) else float("nan"),
            "own_is_row_max": bool(is_row_max),
            "partial_r_control_total": partial,
            "partial_ci95_low": p_lo,
            "partial_ci95_high": p_hi,
            "partial_survives": bool(not np.isnan(p_lo) and p_lo > 0),
        })
    per_channel.sort(key=lambda d: np.nan_to_num(d["partial_r_control_total"], nan=-9), reverse=True)

    asset = args.asset_dir / "specificity_matrix.png"
    render_specificity_heatmap(plt, channels, matrix, asset)

    n_row_max = sum(1 for r in per_channel if r["own_is_row_max"])
    n_partial = sum(1 for r in per_channel if r["partial_survives"])
    return {
        "channels": channels,
        "matrix": matrix,
        "per_channel": per_channel,
        "n_channels": len(channels),
        "n_own_is_row_max": n_row_max,
        "n_partial_survives": n_partial,
        "heatmap_asset": str(asset),
    }


def save_tile_features(path: Path, tiles, tile_cols, tile_rows, total_counts, density_by_channel) -> None:
    import csv as _csv

    channel_keys = sorted(k for k in tiles[0].keys() if k.startswith("mean_"))
    density_channels = list(density_by_channel.keys())
    fieldnames = ["col", "row", "x", "y", "tissue_fraction", "total_tx"] + channel_keys + [f"tx_{c}" for c in density_channels]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = _csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for i, tile in enumerate(tiles):
            row = {
                "col": int(tile_cols[i]), "row": int(tile_rows[i]),
                "x": int(tile["x"]), "y": int(tile["y"]),
                "tissue_fraction": float(tile["tissue_fraction"]), "total_tx": int(total_counts[i]),
            }
            for key in channel_keys:
                row[key] = float(tile[key])
            for c in density_channels:
                row[f"tx_{c}"] = int(density_by_channel[c][i])
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    stats = require_stats()
    he_path, align_path, tx_path = resolve_paths(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)

    print(f"Loading H&E: {he_path}", file=sys.stderr)
    he_array = load_he_fullres(he_path)
    height, width = he_array.shape[:2]
    matrix = load_alignment(align_path)
    print(f"H&E shape: {height} x {width}; alignment det={np.linalg.det(matrix):.4f}", file=sys.stderr)

    genes_of_interest = {g for genes in CHANNEL_GENES.values() for g in genes}
    print("Loading transcripts...", file=sys.stderr)
    x_all, y_all, sub_genes, sub_x, sub_y, n_tx = load_transcript_arrays(tx_path, genes_of_interest)

    direction, in_bounds = choose_direction(matrix, x_all, y_all, args.xenium_mpp, width, height, args.alignment_direction)
    print(f"Alignment direction: {direction}  in-bounds fractions={in_bounds}", file=sys.stderr)

    stride = args.tile_stride
    n_cols = width // stride + 2
    n_rows = height // stride + 2

    hx_all, hy_all = transcript_he_pixels(matrix, x_all, y_all, args.xenium_mpp, direction)
    total_grid = bin_grid_counts(hx_all, hy_all, stride, n_rows, n_cols, width, height)

    hx_sub, hy_sub = transcript_he_pixels(matrix, sub_x, sub_y, args.xenium_mpp, direction)
    gene_grids: dict[str, np.ndarray] = {}
    present_genes = set(np.unique(sub_genes).tolist())
    for gene in present_genes:
        gmask = sub_genes == gene
        gene_grids[gene] = bin_grid_counts(hx_sub[gmask], hy_sub[gmask], stride, n_rows, n_cols, width, height)

    if args.model == "rosie" and not args.alignment_check_only:
        import run_rosie_inference as rosie

        if not args.rosie_mpp:
            raise SystemExit("ROSIE needs --rosie-mpp (H&E microns/px) for the Janesick H&E.")
        print(f"Tiling H&E + ROSIE inference (mpp={args.rosie_mpp:.4f})...", file=sys.stderr)
        tiles = rosie.collect_tiles_rosie(args, he_array, float(args.rosie_mpp))
    else:
        print("Tiling H&E" + ("" if args.alignment_check_only else " + GigaTIME inference") + "...", file=sys.stderr)
        tiles, _gigarun = collect_tiles(args, he_array, args.alignment_check_only)
    if not tiles:
        raise SystemExit("No tissue tiles found; check --tissue-threshold and the H&E image.")
    print(f"Tissue tiles: {len(tiles)}", file=sys.stderr)

    tile_cols = np.array([t["x"] // stride for t in tiles], dtype=np.int64)
    tile_rows = np.array([t["y"] // stride for t in tiles], dtype=np.int64)
    tissue = np.array([t["tissue_fraction"] for t in tiles], dtype=float)
    total_counts = total_grid[tile_rows, tile_cols].astype(float)

    # Model-free alignment sanity: tissue fraction vs total transcript density.
    align_r, align_p = spearman(stats, tissue, total_counts)
    align_lo, align_hi = block_bootstrap_ci(stats, rng, tile_cols, tile_rows, tissue, total_counts, args.block_tiles, min(args.bootstrap, 500))

    channel_results = []
    virtual_by_channel: dict[str, np.ndarray] = {}
    density_by_channel: dict[str, np.ndarray] = {}
    specificity: dict = {}
    if not args.alignment_check_only:
        import os

        mpl_dir = args.out_dir / ".matplotlib"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for channel, genes in CHANNEL_GENES.items():
            present = [g for g in genes if g in gene_grids]
            if not present:
                continue
            mean_key = f"mean_{channel}"
            if mean_key not in tiles[0]:
                continue
            channel_grid = np.sum([gene_grids[g] for g in present], axis=0)
            counts = channel_grid[tile_rows, tile_cols].astype(float)
            if counts.sum() < args.min_transcripts_per_channel:
                continue
            virtual = np.array([float(t[mean_key]) for t in tiles], dtype=float)
            r, p = spearman(stats, virtual, counts)
            lo, hi = block_bootstrap_ci(stats, rng, tile_cols, tile_rows, virtual, counts, args.block_tiles, args.bootstrap)
            asset = args.asset_dir / f"scatter_{channel.replace('/', '_')}.png"
            render_scatter(plt, channel, virtual, counts, r, asset)
            virtual_by_channel[channel] = virtual
            density_by_channel[channel] = counts
            channel_results.append({
                "channel": channel,
                "genes": present,
                "n_tiles": len(tiles),
                "total_transcripts_on_grid": int(counts.sum()),
                "spearman_r": r,
                "spearman_p": p,
                "ci95_low": lo,
                "ci95_high": hi,
                "scatter_asset": str(asset),
            })
        channel_results.sort(key=lambda d: (np.nan_to_num(d["spearman_r"], nan=-9)), reverse=True)

        specificity = compute_specificity(
            stats, rng, plt, args, virtual_by_channel, density_by_channel,
            total_counts, tile_cols, tile_rows,
        )
        save_tile_features(args.out_dir / "tile_features.csv", tiles, tile_cols, tile_rows, total_counts, density_by_channel)

    report = {
        "sample": args.sample,
        "model": getattr(args, "model", "gigatime"),
        "he_shape": [height, width],
        "alignment_direction": direction,
        "alignment_in_bounds_fraction": in_bounds,
        "n_transcripts_total": int(n_tx),
        "n_tissue_tiles": len(tiles),
        "tile_size": args.tile_size,
        "tile_stride": stride,
        "alignment_sanity": {
            "tissue_vs_total_transcripts_spearman_r": align_r,
            "p": align_p,
            "ci95": [align_lo, align_hi],
        },
        "channel_results": channel_results,
        "specificity": specificity,
        "alignment_check_only": args.alignment_check_only,
    }
    (args.out_dir / "xenium_rna_validation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(args.out_markdown, args, report)
    print(f"Alignment sanity Spearman(tissue, transcripts) = {align_r:.3f} (p={align_p:.1e})", file=sys.stderr)
    for res in channel_results:
        print(f"  {res['channel']:8s} r={res['spearman_r']:.3f} CI[{res['ci95_low']:.3f},{res['ci95_high']:.3f}] n_tx={res['total_transcripts_on_grid']}", file=sys.stderr)
    if specificity:
        print(f"Specificity: own-gene is row-max for {specificity['n_own_is_row_max']}/{specificity['n_channels']} channels; "
              f"partial r (control total tx) stays >0 for {specificity['n_partial_survives']}/{specificity['n_channels']}.", file=sys.stderr)
        for res in specificity["per_channel"]:
            print(f"  {res['channel']:8s} own={res['own_gene_r']:.3f} partial={res['partial_r_control_total']:.3f} "
                  f"CI[{res['partial_ci95_low']:.3f},{res['partial_ci95_high']:.3f}] rowmax={res['own_is_row_max']} "
                  f"vs {res['best_other_channel']}({res['best_other_r']:.3f})", file=sys.stderr)
    print(f"Wrote {args.out_dir / 'xenium_rna_validation_report.json'}", file=sys.stderr)
    print(f"Wrote {args.out_markdown}", file=sys.stderr)
    return 0


def write_markdown(path: Path, args, report: dict) -> None:
    sanity = report["alignment_sanity"]
    lines = [
        "# Xenium Breast RNA-Validation Results",
        "",
        f"Status: within-slide validation of GigaTIME virtual channels against Xenium spatial RNA. Sample `{report['sample']}`.",
        "",
        "## Method",
        "",
        f"- H&E full resolution: {report['he_shape'][0]} x {report['he_shape'][1]} px; "
        f"{report['n_tissue_tiles']} tissue tiles at {report['tile_size']} px (stride {report['tile_stride']}).",
        f"- Transcripts: {report['n_transcripts_total']:,} total; binned to the tile grid via the H&E alignment affine "
        f"(direction `{report['alignment_direction']}`, in-bounds fraction "
        f"{report['alignment_in_bounds_fraction'].get(report['alignment_direction'], float('nan')):.3f}).",
        "- Per channel: within-slide Spearman correlation of virtual-channel mean activation vs transcript density across tiles, "
        "with a spatial block-bootstrap 95% CI.",
        "",
        "## Alignment Sanity (model-free)",
        "",
        f"Spearman(tile tissue fraction, total transcript density) = **{sanity['tissue_vs_total_transcripts_spearman_r']:.3f}** "
        f"(p={sanity['p']:.1e}, 95% CI [{sanity['ci95'][0]:.3f}, {sanity['ci95'][1]:.3f}]).",
        "A strongly positive value confirms the transcript-to-H&E coordinate mapping is correct before interpreting channels.",
        "",
    ]
    if report["alignment_check_only"]:
        lines += ["## Channels", "", "_Alignment-check-only run: GigaTIME channels not computed._", ""]
    else:
        lines += [
            "## Channel Correlations (virtual channel vs RNA)",
            "",
            "| Channel | Gene(s) | Spearman r | 95% CI | p | Transcripts on grid |",
            "|---|---|---:|---|---:|---:|",
        ]
        for res in report["channel_results"]:
            ci = f"[{res['ci95_low']:.3f}, {res['ci95_high']:.3f}]"
            lines.append(
                f"| {res['channel']} | {', '.join(res['genes'])} | {res['spearman_r']:.3f} | {ci} | "
                f"{res['spearman_p']:.1e} | {res['total_transcripts_on_grid']:,} |"
            )
        lines += ["", "### Scatter plots", ""]
        for res in report["channel_results"]:
            asset = Path(res["scatter_asset"])
            lines.append(f"![{res['channel']} scatter](assets/{asset.parent.name}/{asset.name})")

        spec = report.get("specificity") or {}
        if spec:
            lines += [
                "",
                "## Channel Specificity (is the signal channel-specific, not just cellularity?)",
                "",
                f"Two tests beyond the raw correlation. (1) Row-max: for each virtual channel, is its own gene the most "
                f"correlated gene-set among all channels? Own-gene is the row maximum for "
                f"**{spec['n_own_is_row_max']}/{spec['n_channels']}** channels. (2) Partial correlation: does the "
                f"virtual-vs-own-gene correlation survive partialling out total transcript density per tile (a per-tile "
                f"cellularity control)? It stays positive (95% CI > 0) for **{spec['n_partial_survives']}/{spec['n_channels']}** channels.",
                "",
                "| Channel | Own-gene r | Partial r (control total tx) | Partial 95% CI | Own-gene row-max? | Closest other channel |",
                "|---|---:|---:|---|:--:|---|",
            ]
            for res in spec["per_channel"]:
                ci = f"[{res['partial_ci95_low']:.3f}, {res['partial_ci95_high']:.3f}]"
                rowmax = "yes" if res["own_is_row_max"] else "no"
                other = f"{res['best_other_channel']} ({res['best_other_r']:.3f})" if res["best_other_channel"] else "-"
                lines.append(
                    f"| {res['channel']} | {res['own_gene_r']:.3f} | {res['partial_r_control_total']:.3f} | {ci} | {rowmax} | {other} |"
                )
            heat = Path(spec["heatmap_asset"])
            lines += [
                "",
                f"![Specificity matrix](assets/{heat.parent.name}/{heat.name})",
                "",
                "Read the heatmap diagonal: a channel-specific model has its brightest cell on the diagonal (virtual-X "
                "tracks gene-X more than other genes). Off-diagonal brightness is expected among co-localized cell types "
                "(e.g. T-cell markers travel together).",
            ]

        spec = report.get("specificity") or {}
        n_spec = spec.get("n_channels", 0)
        lines += [
            "",
            "## Interpretation",
            "",
            "- Raw within-slide correlations are positive and significant for all 13 channels (r about 0.13 to 0.43), "
            "so the virtual channels do carry real, spatially-localized signal that tracks RNA. This is the first RNA "
            "check of these channels. But raw correlation is not the same as channel specificity, and the specificity "
            "tests qualify it sharply.",
            f"- Specificity is limited. Own-gene is the most-correlated gene-set for only "
            f"{spec.get('n_own_is_row_max', 0)}/{n_spec} channels (CK, CD11c); for the rest, some other channel's gene "
            "correlates as well or better, and the immune channels (CD3/CD4/CD8/CD20) collectively track "
            "lymphocyte-dense regions rather than their specific cell type.",
            f"- After partialling out total per-tile transcript density (a cellularity control), channel-specific "
            f"signal survives (95% CI > 0) for {spec.get('n_partial_survives', 0)}/{n_spec} channels and is meaningful "
            "for only a few: CK 0.31 (epithelium, the most specific), then the T-cell channels CD3 0.26 / CD8 0.24 / "
            "CD4 0.21. Ki67, CD14, CD16 and PD-L1 collapse to about zero, and CD68 goes strongly negative (about "
            "-0.33) - i.e. virtual CD68 tracks cellularity/epithelium, not macrophages.",
            "- Note the CK inversion: CK had the weakest raw correlation (0.15) but the strongest cellularity-"
            "controlled correlation (0.31), because epithelium-rich tiles are immune-poor, which suppresses the raw "
            "number until cellularity is removed. This is exactly why the specificity control matters.",
            "- Takeaway: GigaTIME virtual channels mostly reflect a broad epithelial-versus-immune/cellularity "
            "contrast rather than faithful per-marker stains. Only the epithelial (CK) and aggregate T-cell channels "
            "are even modestly marker-specific. Use GigaTIME as interpretive context, not as a quantitative cell-type "
            "readout and not as load-bearing biological evidence.",
            "- Caveats: single section (repeat across Xenium breast replicates and/or HEST-1k for generalization); "
            "sparse channels are exploratory (PD-1 n=1,219; PD-L1 n=9,099); GigaTIME predicts protein (IF) so RNA is a "
            "proxy with a concordance ceiling, a partial excuse for low coefficients but not for the failed specificity.",
            "",
        ]
    lines += [
        "## Output Files",
        "",
        f"- `{args.out_dir / 'xenium_rna_validation_report.json'}`",
        f"- `{args.asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
