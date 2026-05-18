#!/usr/bin/env python3
"""Render H&E thumbnails, tile overlays, and tile montages for analyzed slides."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np


MARKER_COLORS = {
    "CD3": "#2ca02c",
    "CD8": "#1f77b4",
    "PD-L1": "#d62728",
    "CK": "#9467bd",
}


def require_runtime(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    import matplotlib.pyplot as plt
    import openslide
    import pandas as pd
    from PIL import Image, ImageDraw

    return pd, plt, openslide, Image, ImageDraw


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca_extremes/slide_scores.csv")
    parser.add_argument("--tile-scores", default="results/gigatime_tcga_brca_extremes/tile_scores.csv")
    parser.add_argument("--case-groups", default="data/tcga_brca/her2_extreme_cases.csv")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca_extremes/he_examples")
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--thumbnail-width", type=int, default=1800)
    parser.add_argument("--max-slides", type=int, default=12)
    parser.add_argument("--markers", default="CD3,CD8,PD-L1,CK")
    parser.add_argument("--montage-tiles", type=int, default=12)
    return parser.parse_args()


def scale_rect(x: float, y: float, tile_size: int, scale: float) -> tuple[int, int, int, int]:
    return (
        int(round(x * scale)),
        int(round(y * scale)),
        int(round((x + tile_size) * scale)),
        int(round((y + tile_size) * scale)),
    )


def case_label(case_id: str, case_groups) -> str:
    if case_groups is None or case_id not in case_groups.index:
        return ""
    row = case_groups.loc[case_id]
    return f"{row['her2_group']} | ERBB2 TPM {float(row['erbb2_tpm']):.1f}"


def read_thumbnail(slide, width: int):
    w, h = slide.dimensions
    height = max(1, int(round(width * h / w)))
    return slide.get_thumbnail((width, height)).convert("RGB")


def draw_tile_outline_image(thumbnail, slide, tiles, title: str, tile_size: int, out_path: Path, ImageDraw):
    canvas = thumbnail.copy()
    draw = ImageDraw.Draw(canvas, "RGBA")
    scale = canvas.width / slide.dimensions[0]
    for _, tile in tiles.iterrows():
        rect = scale_rect(tile["x"], tile["y"], tile_size, scale)
        draw.rectangle(rect, outline=(255, 215, 0, 230), width=3)
    draw.rectangle((0, 0, canvas.width, 54), fill=(255, 255, 255, 210))
    draw.text((18, 17), title, fill=(0, 0, 0, 255))
    canvas.save(out_path)


def draw_marker_overlay(thumbnail, slide, tiles, marker: str, title: str, tile_size: int, out_path: Path, plt):
    score_col = f"mean_{marker}"
    if score_col not in tiles.columns:
        return
    scale = thumbnail.width / slide.dimensions[0]
    values = tiles[score_col].astype(float)
    vmin = float(values.min())
    vmax = float(values.max())
    norm = plt.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1e-9)
    cmap = plt.get_cmap("magma")

    fig, ax = plt.subplots(figsize=(10, max(3, 10 * thumbnail.height / thumbnail.width)))
    ax.imshow(thumbnail)
    for _, tile in tiles.iterrows():
        x0, y0, x1, y1 = scale_rect(tile["x"], tile["y"], tile_size, scale)
        color = cmap(norm(float(tile[score_col])))
        rect = plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color, edgecolor="white", linewidth=0.5, alpha=0.58)
        ax.add_patch(rect)
    ax.set_title(f"{title} | virtual {marker}")
    ax.axis("off")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label(f"Mean virtual {marker}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def save_tile_montage(slide, tiles, marker: str, title: str, tile_size: int, n_tiles: int, out_path: Path, plt):
    score_col = f"mean_{marker}"
    if score_col not in tiles.columns:
        return
    selected = tiles.sort_values(score_col, ascending=False).head(n_tiles)
    if selected.empty:
        return
    ncols = 4
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 2.6, nrows * 2.9))
    axes = np.atleast_1d(axes).flatten()
    for ax, (_, tile) in zip(axes, selected.iterrows()):
        patch = slide.read_region((int(tile["x"]), int(tile["y"])), 0, (tile_size, tile_size)).convert("RGB")
        ax.imshow(patch)
        ax.set_title(f"{marker} {float(tile[score_col]):.3f}", fontsize=9)
        ax.axis("off")
    for ax in axes[len(selected) :]:
        ax.axis("off")
    fig.suptitle(f"{title} | top H&E tiles by virtual {marker}", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, openslide, _Image, ImageDraw = require_runtime(out_dir / ".matplotlib")

    slide_scores = pd.read_csv(args.slide_scores)
    tile_scores = pd.read_csv(args.tile_scores)
    case_groups = None
    if args.case_groups and Path(args.case_groups).exists():
        case_groups = pd.read_csv(args.case_groups).set_index("case_submitter_id")

    markers = [marker.strip() for marker in args.markers.split(",") if marker.strip()]
    manifest_rows = []
    for _, slide_row in slide_scores.head(args.max_slides).iterrows():
        slide_path = Path(slide_row["slide_path"])
        slide_id = slide_row["slide_id"]
        case_id = slide_row["case_submitter_id"]
        if not slide_path.exists():
            print(f"Skipping missing slide: {slide_path}")
            continue
        tiles = tile_scores[tile_scores["slide_id"] == slide_id].copy()
        if tiles.empty:
            continue
        title = f"{case_id} | {case_label(case_id, case_groups)}".strip(" |")
        slide_out = out_dir / slide_id
        slide_out.mkdir(parents=True, exist_ok=True)
        slide = openslide.OpenSlide(str(slide_path))
        thumbnail = read_thumbnail(slide, args.thumbnail_width)

        overview_path = slide_out / f"{slide_id}_he_thumbnail.png"
        thumbnail.save(overview_path)
        outlined_path = slide_out / f"{slide_id}_he_sampled_tiles.png"
        draw_tile_outline_image(thumbnail, slide, tiles, title, args.tile_size, outlined_path, ImageDraw)

        manifest_rows.append({"case_submitter_id": case_id, "slide_id": slide_id, "image_type": "he_thumbnail", "path": str(overview_path)})
        manifest_rows.append({"case_submitter_id": case_id, "slide_id": slide_id, "image_type": "sampled_tile_overlay", "path": str(outlined_path)})

        for marker in markers:
            overlay_path = slide_out / f"{slide_id}_he_virtual_{marker}_overlay.png"
            draw_marker_overlay(thumbnail, slide, tiles, marker, title, args.tile_size, overlay_path, plt)
            montage_path = slide_out / f"{slide_id}_he_top_{marker}_tiles.png"
            save_tile_montage(slide, tiles, marker, title, args.tile_size, args.montage_tiles, montage_path, plt)
            manifest_rows.append({"case_submitter_id": case_id, "slide_id": slide_id, "image_type": f"virtual_{marker}_overlay", "path": str(overlay_path)})
            manifest_rows.append({"case_submitter_id": case_id, "slide_id": slide_id, "image_type": f"top_{marker}_tile_montage", "path": str(montage_path)})
        slide.close()

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(out_dir / "he_image_manifest.csv", index=False)
    print(f"Wrote {len(manifest)} images to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
