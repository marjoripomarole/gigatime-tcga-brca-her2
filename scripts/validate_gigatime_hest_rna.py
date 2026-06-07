#!/usr/bin/env python3
"""Validate GigaTIME virtual channels against HEST-1k breast spatial RNA, within-slide.

Generalization of scripts/validate_gigatime_xenium_rna.py to public HEST-1k breast
samples (Mahmood Lab; gated HF dataset MahmoodLab/hest). This closes the long-standing
"single Xenium section" caveat by replicating the RNA-specificity audit on independent
breast tissues / patients and on a second platform (Visium).

Design: the audited statistics core (within-slide Spearman, channel x gene-set
specificity matrix, cellularity-controlled partial correlation, spatial block-bootstrap
CIs, figures, tile-feature dump) is REUSED VERBATIM by importing the leaf helpers from
validate_gigatime_xenium_rna -- that file is not modified, so the statistical method is
byte-for-byte the same one documented for Rep1/Rep2. Only the data loader and write-up differ.

Two HEST modalities (--source):
  xenium: transcripts/<id>_transcripts.parquet stores every transcript already aligned into
          H&E pixel space in columns he_x/he_y, so binning onto GigaTIME's 256 px tile grid is
          direct -- no alignment affine, no microns->pixel, no choose_direction. is_gene drops
          negative-control / deprecated codewords.
  visium: st/<id>.h5ad (scanpy AnnData); spot centroids are full-res H&E pixels in obs columns
          pxl_col_in_fullres / pxl_row_in_fullres and X holds raw counts (whole transcriptome).
          Per-tile density = sum of spot counts binned onto the grid; because Visium spots
          (~100 um pitch) are sparser than 256 px tiles, the analysis is restricted to tiles
          that contain >=1 spot.
H&E is wsis/<id>.tif (pyramidal, level 0 = full resolution).

Run (gigatime-tcga env):
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/validate_gigatime_hest_rna.py --id TENX199
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/validate_gigatime_hest_rna.py --id TENX39 --source visium
  # cheap pre-check (no GPU): add --alignment-check-only
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Reuse the audited statistics/render/inference core unchanged.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_gigatime_xenium_rna as xrna  # noqa: E402

CHANNEL_GENES = xrna.CHANNEL_GENES


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--id", required=True, help="HEST sample id, e.g. TENX199.")
    p.add_argument("--source", choices=["xenium", "visium"], default="xenium",
                   help="HEST modality: Xenium transcripts (he_x/he_y) or Visium spots (h5ad).")
    p.add_argument("--hest-dir", type=Path, default=Path("data/hest"))
    p.add_argument("--wsi", type=Path, default=None, help="Defaults to <hest-dir>/wsis/<id>.tif.")
    p.add_argument("--transcripts", type=Path, default=None, help="Xenium; defaults to <hest-dir>/transcripts/<id>_transcripts.parquet.")
    p.add_argument("--st", type=Path, default=None, help="Visium AnnData; defaults to <hest-dir>/st/<id>.h5ad.")
    p.add_argument("--metadata", type=Path, default=None, help="Defaults to <hest-dir>/metadata/<id>.json.")
    p.add_argument("--gigatime-repo", default="external/GigaTIME")
    p.add_argument("--tile-size", type=int, default=256)
    p.add_argument("--tile-stride", type=int, default=256)
    p.add_argument("--tissue-threshold", type=float, default=0.35)
    p.add_argument("--max-tiles", type=int, default=0, help="Cap tissue tiles (0 = all); small value for a smoke run.")
    p.add_argument("--tile-order", default="row-major", choices=["row-major", "random"])
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    p.add_argument("--min-transcripts-per-channel", type=int, default=50)
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--block-tiles", type=int, default=4)
    p.add_argument("--alignment-check-only", action="store_true", help="Skip GigaTIME; only correlate tissue fraction vs total transcript density.")
    p.add_argument("--out-dir", type=Path, default=None, help="Defaults to results/gigatime_hest_rna_validation/<id>.")
    p.add_argument("--asset-dir", type=Path, default=None, help="Defaults to docs/assets/gigatime_hest_rna_validation_<id>.")
    p.add_argument("--out-markdown", type=Path, default=None, help="Defaults to docs/<model>_rna_validation_<id>.md.")
    p.add_argument("--model", choices=["gigatime", "rosie"], default="gigatime",
                   help="Virtual-mIF model whose channels are audited against RNA.")
    p.add_argument("--rosie-weights", type=Path, default=Path("external/ROSIE/best_model_single.pth"))
    p.add_argument("--rosie-grid", type=int, default=3, help="ROSIE patch centers per tile (grid x grid).")
    p.add_argument("--rosie-batch", type=int, default=64, help="ROSIE patch batch size.")
    p.add_argument("--mpp", type=float, default=None, help="H&E microns/px (ROSIE only; defaults to metadata pixel_size_um_estimated).")
    args = p.parse_args()
    args.wsi = args.wsi or args.hest_dir / "wsis" / f"{args.id}.tif"
    args.transcripts = args.transcripts or args.hest_dir / "transcripts" / f"{args.id}_transcripts.parquet"
    args.st = args.st or args.hest_dir / "st" / f"{args.id}.h5ad"
    args.metadata = args.metadata or args.hest_dir / "metadata" / f"{args.id}.json"
    tag = "gigatime" if args.model == "gigatime" else args.model
    md_stem = "hest_rna_validation" if args.model == "gigatime" else f"{tag}_rna_validation"
    args.out_dir = args.out_dir or Path(f"results/{tag}_hest_rna_validation") / args.id
    args.asset_dir = args.asset_dir or Path("docs/assets") / f"{tag}_hest_rna_validation_{args.id}"
    args.out_markdown = args.out_markdown or Path("docs") / f"{md_stem}_{args.id}.md"
    return args


def load_metadata(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def load_he_robust(path: Path) -> np.ndarray:
    """Load a full-res RGB H&E array, tolerating HEST pyramids that break tifffile series detection.

    Tries the shared series-based loader first; on failure (e.g. 'incompatible keyframe' on some
    HEST LZW pyramids) falls back to decoding page 0 (full resolution) directly, then normalizes
    to HxWx3 uint8 exactly as xrna.load_he_fullres does.
    """
    try:
        return xrna.load_he_fullres(path)
    except Exception as exc:  # noqa: BLE001 - any tifffile series/keyframe failure
        print(f"load_he_fullres failed ({type(exc).__name__}: {exc}); falling back to page 0.", file=sys.stderr)
    import tifffile

    with tifffile.TiffFile(str(path)) as tif:
        arr = np.asarray(tif.pages[0].asarray())
    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    elif arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[2] not in (3, 4):
        arr = np.transpose(arr, (1, 2, 0))
    arr = arr[:, :, :3]
    if arr.dtype != np.uint8:
        amax = float(arr.max()) or 1.0
        arr = np.clip(arr / amax * 255.0, 0, 255).astype(np.uint8) if amax > 255 else arr.astype(np.uint8)
    return np.ascontiguousarray(arr)


def bin_grid_weighted(x, y, w, stride: int, n_rows: int, n_cols: int, width: int, height: int) -> np.ndarray:
    """Sum per-point weights w into the tile grid (Visium: spot counts onto tiles)."""
    valid = (x >= 0) & (x < width) & (y >= 0) & (y < height)
    cols = (x[valid] / stride).astype(np.int64)
    rows = (y[valid] / stride).astype(np.int64)
    wv = np.asarray(w, dtype=float)[valid]
    keep = (cols >= 0) & (cols < n_cols) & (rows >= 0) & (rows < n_rows)
    linear = rows[keep] * n_cols + cols[keep]
    grid = np.bincount(linear, weights=wv[keep], minlength=n_rows * n_cols).reshape(n_rows, n_cols)
    return grid


def load_hest_xenium_grids(args, width: int, height: int):
    """Bin HEST-Xenium transcripts (already in H&E px via he_x/he_y) onto the tile grid.

    Returns (total_grid, gene_grids, n_rows, n_cols, n_gene_transcripts, n_all_rows).
    """
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    table = pq.read_table(args.transcripts, columns=["feature_name", "he_x", "he_y", "is_gene"])
    n_all_rows = table.num_rows
    table = table.filter(table["is_gene"])  # drop negative-control / deprecated codewords

    stride = args.tile_stride
    n_cols = width // stride + 2
    n_rows = height // stride + 2

    hx_all = table["he_x"].to_numpy(zero_copy_only=False).astype(float)
    hy_all = table["he_y"].to_numpy(zero_copy_only=False).astype(float)
    total_grid = xrna.bin_grid_counts(hx_all, hy_all, stride, n_rows, n_cols, width, height)

    goi = sorted({g for genes in CHANNEL_GENES.values() for g in genes})
    goi_arr = pa.array(goi, type=table["feature_name"].type)
    sub = table.filter(pc.is_in(table["feature_name"], value_set=goi_arr))
    sub_names = np.array(
        [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in sub["feature_name"].to_pylist()],
        dtype=object,
    )
    sub_x = sub["he_x"].to_numpy(zero_copy_only=False).astype(float)
    sub_y = sub["he_y"].to_numpy(zero_copy_only=False).astype(float)

    gene_grids: dict[str, np.ndarray] = {}
    for gene in np.unique(sub_names):
        m = sub_names == gene
        gene_grids[str(gene)] = xrna.bin_grid_counts(sub_x[m], sub_y[m], stride, n_rows, n_cols, width, height)
    return total_grid, gene_grids, n_rows, n_cols, int(table.num_rows), int(n_all_rows)


def load_hest_visium_grids(args, width: int, height: int):
    """Bin HEST-Visium spot counts (centroids in full-res H&E px) onto the tile grid.

    Returns (total_grid, gene_grids, n_rows, n_cols, n_total_umi, n_spots).
    """
    import anndata as ad
    import scipy.sparse as sp

    a = ad.read_h5ad(args.st)
    x = a.obs["pxl_col_in_fullres"].to_numpy().astype(float)
    y = a.obs["pxl_row_in_fullres"].to_numpy().astype(float)
    X = a.X
    total_per_spot = np.asarray(X.sum(axis=1)).ravel().astype(float)

    stride = args.tile_stride
    n_cols = width // stride + 2
    n_rows = height // stride + 2
    total_grid = bin_grid_weighted(x, y, total_per_spot, stride, n_rows, n_cols, width, height)

    var_index = {str(g): i for i, g in enumerate(a.var_names)}
    goi = sorted({g for genes in CHANNEL_GENES.values() for g in genes})
    gene_grids: dict[str, np.ndarray] = {}
    for gene in goi:
        j = var_index.get(gene)
        if j is None:
            continue
        col = X[:, j]
        col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        if col.sum() <= 0:
            continue
        gene_grids[gene] = bin_grid_weighted(x, y, col.astype(float), stride, n_rows, n_cols, width, height)
    return total_grid, gene_grids, n_rows, n_cols, int(total_per_spot.sum()), int(a.n_obs)


def run_channel_analysis(args, stats, rng, plt, tiles, gene_grids, tile_cols, tile_rows, total_counts):
    """Per-channel raw Spearman + figures, then the shared specificity/partial analysis.

    Mirrors validate_gigatime_xenium_rna.main()'s analysis block; every statistic comes
    from the imported xrna helpers so the method is identical.
    """
    channel_results = []
    virtual_by_channel: dict[str, np.ndarray] = {}
    density_by_channel: dict[str, np.ndarray] = {}
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
        r, pval = xrna.spearman(stats, virtual, counts)
        lo, hi = xrna.block_bootstrap_ci(stats, rng, tile_cols, tile_rows, virtual, counts, args.block_tiles, args.bootstrap)
        asset = args.asset_dir / f"scatter_{channel.replace('/', '_')}.png"
        xrna.render_scatter(plt, channel, virtual, counts, r, asset)
        virtual_by_channel[channel] = virtual
        density_by_channel[channel] = counts
        channel_results.append({
            "channel": channel,
            "genes": present,
            "n_tiles": len(tiles),
            "total_transcripts_on_grid": int(counts.sum()),
            "spearman_r": r,
            "spearman_p": pval,
            "ci95_low": lo,
            "ci95_high": hi,
            "scatter_asset": str(asset),
        })
    channel_results.sort(key=lambda d: np.nan_to_num(d["spearman_r"], nan=-9), reverse=True)
    specificity = xrna.compute_specificity(
        stats, rng, plt, args, virtual_by_channel, density_by_channel, total_counts, tile_cols, tile_rows,
    )
    return channel_results, virtual_by_channel, density_by_channel, specificity


def main() -> int:
    args = parse_args()
    stats = xrna.require_stats()
    data_label, data_path = (("transcripts", args.transcripts) if args.source == "xenium" else ("ST h5ad", args.st))
    for label, path in (("WSI", args.wsi), (data_label, data_path)):
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)
    meta = load_metadata(args.metadata)

    print(f"Loading H&E WSI: {args.wsi}", file=sys.stderr)
    he_array = load_he_robust(args.wsi)
    height, width = he_array.shape[:2]
    print(f"H&E shape: {height} x {width}", file=sys.stderr)

    if args.source == "xenium":
        print("Loading + binning transcripts (he_x/he_y, pre-aligned)...", file=sys.stderr)
        total_grid, gene_grids, n_rows, n_cols, n_primary, n_secondary = load_hest_xenium_grids(args, width, height)
    else:
        print("Loading + binning Visium spots (pxl_*_in_fullres)...", file=sys.stderr)
        total_grid, gene_grids, n_rows, n_cols, n_primary, n_secondary = load_hest_visium_grids(args, width, height)
    stride = args.tile_stride
    present_channels = [c for c, gs in CHANNEL_GENES.items() if any(g in gene_grids for g in gs)]
    missing_channels = [c for c in CHANNEL_GENES if c not in present_channels]
    print(f"On-grid signal units: {n_primary:,}; channels with a panel gene: {present_channels}", file=sys.stderr)

    rosie_mpp = None
    if args.model == "gigatime" or args.alignment_check_only:
        print("Tiling H&E" + ("" if args.alignment_check_only else " + GigaTIME inference") + "...", file=sys.stderr)
        tiles, _gigarun = xrna.collect_tiles(args, he_array, args.alignment_check_only)
    else:
        import run_rosie_inference as rosie

        rosie_mpp = args.mpp if args.mpp else meta.get("pixel_size_um_estimated")
        if not rosie_mpp:
            raise SystemExit("ROSIE needs an H&E mpp: pass --mpp or provide metadata pixel_size_um_estimated.")
        print(f"Tiling H&E + ROSIE inference (mpp={float(rosie_mpp):.4f})...", file=sys.stderr)
        tiles = rosie.collect_tiles_rosie(args, he_array, float(rosie_mpp))
    if not tiles:
        raise SystemExit("No tissue tiles found; check --tissue-threshold and the WSI.")

    tile_cols = np.array([t["x"] // stride for t in tiles], dtype=np.int64)
    tile_rows = np.array([t["y"] // stride for t in tiles], dtype=np.int64)
    tissue = np.array([t["tissue_fraction"] for t in tiles], dtype=float)
    total_counts = total_grid[tile_rows, tile_cols].astype(float)

    # Visium spots (~100 um pitch) are sparser than 256 px tiles -> restrict to spot-occupied tiles.
    n_occupied = None
    if args.source == "visium":
        occ = total_counts > 0
        tiles = [t for t, k in zip(tiles, occ) if k]
        tile_cols, tile_rows = tile_cols[occ], tile_rows[occ]
        tissue, total_counts = tissue[occ], total_counts[occ]
        n_occupied = int(len(tiles))
        if not tiles:
            raise SystemExit("No tissue tiles overlap Visium spots.")
        print(f"Visium: {n_occupied} of {len(occ)} tissue tiles contain >=1 spot.", file=sys.stderr)
    print(f"Tiles used: {len(tiles)}", file=sys.stderr)

    align_r, align_p = xrna.spearman(stats, tissue, total_counts)
    align_lo, align_hi = xrna.block_bootstrap_ci(stats, rng, tile_cols, tile_rows, tissue, total_counts, args.block_tiles, min(args.bootstrap, 500))

    channel_results: list = []
    specificity: dict = {}
    if not args.alignment_check_only:
        import os

        mpl_dir = args.out_dir / ".matplotlib"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        channel_results, _vbc, density_by_channel, specificity = run_channel_analysis(
            args, stats, rng, plt, tiles, gene_grids, tile_cols, tile_rows, total_counts,
        )
        xrna.save_tile_features(args.out_dir / "tile_features.csv", tiles, tile_cols, tile_rows, total_counts, density_by_channel)

    report = {
        "sample": args.id,
        "model": args.model,
        "rosie_mpp": rosie_mpp,
        "source": "HEST-1k",
        "modality": args.source,
        "st_technology": meta.get("st_technology", args.source),
        "patient": meta.get("patient"),
        "subseries": meta.get("subseries"),
        "dataset_title": meta.get("dataset_title"),
        "disease_comment": meta.get("disease_comment"),
        "oncotree_code": meta.get("oncotree_code"),
        "pixel_size_um": meta.get("pixel_size_um_estimated"),
        "he_shape": [height, width],
        "n_primary_units": n_primary,    # xenium: gene transcripts on grid; visium: total UMI
        "n_secondary_units": n_secondary,  # xenium: all parquet rows; visium: n_spots
        "n_tissue_tiles": int(len(tiles)),
        "n_occupied_tiles": n_occupied,
        "tile_size": args.tile_size,
        "tile_stride": stride,
        "channels_present": present_channels,
        "channels_missing_in_panel": missing_channels,
        "alignment_sanity": {
            "tissue_vs_total_transcripts_spearman_r": align_r,
            "p": align_p,
            "ci95": [align_lo, align_hi],
        },
        "channel_results": channel_results,
        "specificity": specificity,
        "alignment_check_only": args.alignment_check_only,
    }
    (args.out_dir / "hest_rna_validation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown_hest(args.out_markdown, args, report)

    print(f"Alignment sanity Spearman(tissue, transcripts) = {align_r:.3f} (p={align_p:.1e}) CI[{align_lo:.3f},{align_hi:.3f}]", file=sys.stderr)
    for res in channel_results:
        print(f"  {res['channel']:8s} r={res['spearman_r']:.3f} CI[{res['ci95_low']:.3f},{res['ci95_high']:.3f}] n={res['total_transcripts_on_grid']}", file=sys.stderr)
    if specificity:
        print(f"Specificity: own-gene row-max {specificity['n_own_is_row_max']}/{specificity['n_channels']}; "
              f"partial r >0 for {specificity['n_partial_survives']}/{specificity['n_channels']}.", file=sys.stderr)
        for res in specificity["per_channel"]:
            print(f"  {res['channel']:8s} own={res['own_gene_r']:.3f} partial={res['partial_r_control_total']:.3f} "
                  f"CI[{res['partial_ci95_low']:.3f},{res['partial_ci95_high']:.3f}] rowmax={res['own_is_row_max']}", file=sys.stderr)
    print(f"Wrote {args.out_dir / 'hest_rna_validation_report.json'}", file=sys.stderr)
    print(f"Wrote {args.out_markdown}", file=sys.stderr)
    return 0


def _interpretation_lines(report: dict) -> list[str]:
    """Build an honest, computed-from-results interpretation (no hard-coded numbers)."""
    spec = report.get("specificity") or {}
    per = {r["channel"]: r for r in spec.get("per_channel", [])}
    lines: list[str] = []
    if not per:
        return ["- Alignment-check-only run: GigaTIME channels not computed."]

    survivors = [r for r in spec["per_channel"] if r.get("partial_survives")]
    survivors_txt = ", ".join(f"{r['channel']} {r['partial_r_control_total']:.2f}" for r in survivors[:8]) or "none"
    negatives = [r for r in spec["per_channel"] if not np.isnan(r["partial_r_control_total"]) and r["partial_r_control_total"] < -0.05]
    neg_txt = ", ".join(f"{r['channel']} {r['partial_r_control_total']:.2f}" for r in negatives) or "none"

    lines.append(
        f"- Own-gene is the most-correlated gene-set for **{spec['n_own_is_row_max']}/{spec['n_channels']}** channels; "
        f"after partialling out total per-tile transcript density (cellularity), channel-specific signal stays positive "
        f"(95% CI > 0) for **{spec['n_partial_survives']}/{spec['n_channels']}** channels: {survivors_txt}."
    )
    if negatives:
        lines.append(f"- Channels going negative after the cellularity control (track epithelium/cellularity, not their marker): {neg_txt}.")

    checks = []
    if "CK" in per:
        ck = per["CK"]["partial_r_control_total"]
        checks.append(f"CK partial r = {ck:.2f} ({'specific/positive' if ck > 0.05 else 'not positive'})")
    tcell = [c for c in ("CD3", "CD8", "CD4") if c in per]
    if tcell:
        checks.append("T-cell " + ", ".join(f"{c} {per[c]['partial_r_control_total']:.2f}" for c in tcell))
    if "CD68" in per:
        cd68 = per["CD68"]["partial_r_control_total"]
        checks.append(f"CD68 = {cd68:.2f} ({'negative' if cd68 < 0 else 'not negative'})")
    if checks:
        lines.append("- Headline-channel check (CK epithelium; T-cell; CD68 macrophage): " + "; ".join(checks) + ".")
    return lines


def write_markdown_hest(path: Path, args, report: dict) -> None:
    sanity = report["alignment_sanity"]
    pid = report.get("patient") or "n/a"
    px = report.get("pixel_size_um")
    px_txt = f"{px:.4f}" if isinstance(px, (int, float)) else "n/a"
    modality = report.get("modality", "xenium")
    if modality == "xenium":
        data_line = (
            f"- Transcripts: {report['n_primary_units']:,} gene transcripts (of {report['n_secondary_units']:,} incl. controls), "
            f"binned onto the tile grid directly via the HEST-provided H&E pixel coordinates (`he_x`/`he_y`) — no alignment affine."
        )
    else:
        data_line = (
            f"- Visium: {report['n_secondary_units']:,} spots ({report['n_primary_units']:,} total UMI), binned onto the tile grid via "
            f"`pxl_col/row_in_fullres`. Analysis restricted to the **{report.get('n_occupied_tiles')}** tiles containing >=1 spot "
            f"(spots are ~100 um apart, sparser than 256 px tiles)."
        )
    model_name = {"gigatime": "GigaTIME", "rosie": "ROSIE"}.get(report.get("model", "gigatime"), report.get("model"))
    lines = [
        f"# HEST-1k Breast RNA-Validation Results — {report['sample']} ({model_name})",
        "",
        f"Status: within-slide validation of {model_name} virtual channels against HEST-1k spatial RNA "
        f"({report.get('st_technology')}). Same audited pipeline as the GigaTIME run, applied to a second "
        f"H&E->virtual-mIF model for a field-level specificity claim.",
        "",
        f"- Sample: `{report['sample']}` ({report.get('st_technology')}, HEST-1k); {pid}; "
        f"`{report.get('subseries') or ''}`. Dataset: {report.get('dataset_title') or 'n/a'}.",
        f"- Clinical (from HEST metadata): {report.get('oncotree_code') or 'n/a'}; {report.get('disease_comment') or 'n/a'}.",
        "",
        "## Method",
        "",
        f"- H&E full resolution: {report['he_shape'][0]} x {report['he_shape'][1]} px ({px_txt} um/px); "
        f"{report['n_tissue_tiles']} tiles used at {report['tile_size']} px (stride {report['tile_stride']}).",
        data_line,
        f"- Channels with a panel gene ({len(report['channels_present'])}/{len(report['channels_present']) + len(report['channels_missing_in_panel'])}): "
        f"{', '.join(report['channels_present'])}. Not in this panel: {', '.join(report['channels_missing_in_panel']) or 'none'}.",
        "- Statistics are computed by the same audited core as the Xenium Rep1/Rep2 run "
        "(`scripts/validate_gigatime_xenium_rna.py`, imported unchanged): within-slide Spearman, channel x gene-set "
        "specificity matrix, cellularity-controlled partial correlation, spatial block-bootstrap 95% CIs.",
        "",
        "## Alignment Sanity (model-free)",
        "",
        f"Spearman(tile tissue fraction, total transcript density) = **{sanity['tissue_vs_total_transcripts_spearman_r']:.3f}** "
        f"(p={sanity['p']:.1e}, 95% CI [{sanity['ci95'][0]:.3f}, {sanity['ci95'][1]:.3f}]). "
        "A strongly positive value confirms the transcript-to-H&E mapping before interpreting channels.",
        "",
    ]
    if report["alignment_check_only"]:
        lines += ["## Channels", "", "_Alignment-check-only run: GigaTIME channels not computed._", ""]
    else:
        lines += [
            "## Channel Correlations (virtual channel vs RNA)",
            "",
            "| Channel | Gene(s) | Spearman r | 95% CI | p | Counts on grid |",
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
                f"(1) Row-max: own-gene is the most-correlated gene-set for **{spec['n_own_is_row_max']}/{spec['n_channels']}** "
                f"channels. (2) Partial correlation controlling for total per-tile transcript density stays positive "
                f"(95% CI > 0) for **{spec['n_partial_survives']}/{spec['n_channels']}** channels.",
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
            lines += ["", f"![Specificity matrix](assets/{heat.parent.name}/{heat.name})", ""]

        lines += ["## Interpretation", ""] + _interpretation_lines(report) + [""]

    lines += [
        "## Output Files",
        "",
        f"- `{args.out_dir / 'hest_rna_validation_report.json'}`",
        f"- `{args.asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
