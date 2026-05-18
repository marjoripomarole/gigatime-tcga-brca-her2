#!/usr/bin/env python3
"""Render fluorescence-style virtual mIF composites from GigaTIME map outputs."""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from run_gigatime_tcga_brca import (
    GIGATIME_CHANNELS,
    import_runtime,
    load_model,
    preprocess_tile,
    resolve_device,
)


PANELS = {
    "immune_checkpoint": [
        ("DAPI", "#2D5BFF"),
        ("CK", "#FF3B30"),
        ("CD3", "#00E676"),
        ("CD8", "#00D5FF"),
        ("PD-L1", "#FF2DCE"),
        ("PD-1", "#FFD60A"),
    ],
    "myeloid_bcell": [
        ("DAPI", "#2D5BFF"),
        ("CD68", "#FF9500"),
        ("CD11c", "#64D2FF"),
        ("CD20", "#BF5AF2"),
        ("CD14", "#FFD60A"),
        ("CD4", "#30D158"),
    ],
    "tumor_proliferation": [
        ("DAPI", "#2D5BFF"),
        ("CK", "#FF3B30"),
        ("Ki67", "#FFD60A"),
        ("Caspase3-D", "#30D158"),
        ("PHH3-B", "#FF2DCE"),
        ("CD34", "#64D2FF"),
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca_extremes/slide_scores.csv")
    parser.add_argument("--tile-scores", default="results/gigatime_tcga_brca_extremes/tile_scores.csv")
    parser.add_argument("--joined", default="results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv")
    parser.add_argument("--out-dir", default="docs/assets/virtual_mif_composites")
    parser.add_argument("--gigatime-repo", default="external/GigaTIME")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--tiles-per-panel", type=int, default=16)
    parser.add_argument("--batch-size", type=int, default=8)
    return parser.parse_args()


def require_runtime(mpl_config_dir: Path):
    cache_dir = Path(tempfile.gettempdir()) / "gigatime_tcga_mplconfig"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import pandas as pd

    return pd, plt, mpatches


def hex_to_rgb(hex_color: str) -> np.ndarray:
    clean = hex_color.lstrip("#")
    return np.array([int(clean[i : i + 2], 16) for i in (0, 2, 4)], dtype=np.float32) / 255.0


def reference_rows(joined):
    return [
        ("HER2-high", joined.sort_values("erbb2_tpm", ascending=False).iloc[0]),
        ("HER2-low", joined.sort_values("erbb2_tpm", ascending=True).iloc[0]),
    ]


def select_tiles(tile_scores, slide_id: str, panel: list[tuple[str, str]], n_tiles: int):
    slide_tiles = tile_scores[tile_scores["slide_id"] == slide_id].copy()
    marker_cols = [f"mean_{marker}" for marker, _color in panel if marker != "DAPI" and f"mean_{marker}" in slide_tiles.columns]
    if not marker_cols:
        return slide_tiles.head(n_tiles)
    slide_tiles["panel_signal"] = slide_tiles[marker_cols].astype(float).sum(axis=1)
    return slide_tiles.sort_values(["panel_signal", "tissue_fraction"], ascending=False).head(n_tiles)


def infer_tile_maps(torch, model, openslide, slide_path: Path, tile_rows, tile_size: int, batch_size: int, device):
    slide = openslide.OpenSlide(str(slide_path))
    records: list[dict[str, Any]] = []
    batch = []
    meta = []
    with torch.no_grad():
        for _idx, row in tile_rows.iterrows():
            x = int(row["x"])
            y = int(row["y"])
            region = slide.read_region((x, y), 0, (tile_size, tile_size)).convert("RGB")
            rgb = np.asarray(region)
            batch.append(preprocess_tile(torch, rgb, device))
            meta.append({"x": x, "y": y, "he_rgb": rgb})
            if len(batch) == batch_size:
                records.extend(run_batch(torch, model, batch, meta))
                batch = []
                meta = []
        if batch:
            records.extend(run_batch(torch, model, batch, meta))
    slide.close()
    return records


def run_batch(torch, model, batch, meta):
    tensor = torch.stack(batch, dim=0)
    maps = torch.sigmoid(model(tensor)).detach().cpu().numpy()
    rows = []
    for idx, item in enumerate(meta):
        rows.append({**item, "maps": maps[idx]})
    return rows


def panel_limits(records: list[dict[str, Any]], panel: list[tuple[str, str]]) -> dict[str, tuple[float, float]]:
    limits = {}
    for marker, _color in panel:
        channel_idx = GIGATIME_CHANNELS.index(marker)
        values = np.concatenate([record["maps"][channel_idx].ravel() for record in records])
        lower_quantile = 0.55 if marker == "DAPI" else 0.80
        upper_quantile = 0.995
        low = float(np.quantile(values, lower_quantile))
        high = float(np.quantile(values, upper_quantile))
        if high <= low:
            high = low + 1e-6
        limits[marker] = (low, high)
    return limits


def make_composite(channel_maps: np.ndarray, panel: list[tuple[str, str]], limits: dict[str, tuple[float, float]]) -> np.ndarray:
    height, width = channel_maps.shape[1:]
    composite = np.zeros((height, width, 3), dtype=np.float32)
    for marker, color in panel:
        channel_idx = GIGATIME_CHANNELS.index(marker)
        low, high = limits[marker]
        signal = np.clip((channel_maps[channel_idx] - low) / (high - low), 0, 1)
        signal = np.power(signal, 0.75)
        composite += signal[..., None] * hex_to_rgb(color)
    return np.clip(composite, 0, 1)


def save_montage(records, panel_name: str, panel: list[tuple[str, str]], title: str, out_path: Path, plt, mpatches) -> None:
    limits = panel_limits(records, panel)
    n_tiles = len(records)
    ncols = int(math.ceil(math.sqrt(n_tiles)))
    nrows = int(math.ceil(n_tiles / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 2.6, nrows * 2.8), facecolor="black")
    axes = np.atleast_1d(axes).flatten()
    for axis, record in zip(axes, records):
        composite = make_composite(record["maps"], panel, limits)
        axis.imshow(composite)
        axis.set_facecolor("black")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(f"x={record['x']} y={record['y']}", color="white", fontsize=8)
        for spine in axis.spines.values():
            spine.set_edgecolor("#333333")
    for axis in axes[n_tiles:]:
        axis.axis("off")
        axis.set_facecolor("black")
    handles = [mpatches.Patch(color=color, label=marker) for marker, color in panel]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), facecolor="black", labelcolor="white")
    fig.suptitle(title, color="white", fontsize=15, y=0.985)
    fig.text(
        0.5,
        0.045,
        "Fluorescence-style rendering of GigaTIME predictions from H&E tiles. This is virtual mIF, not experimental mIF.",
        color="white",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def save_he_vs_mif(records, panel: list[tuple[str, str]], title: str, out_path: Path, plt, mpatches) -> None:
    limits = panel_limits(records, panel)
    selected = records[:6]
    fig, axes = plt.subplots(len(selected), 2, figsize=(7, len(selected) * 3.0), facecolor="black")
    for row_idx, record in enumerate(selected):
        axes[row_idx, 0].imshow(record["he_rgb"])
        axes[row_idx, 0].set_title("H&E tile", color="white", fontsize=10)
        axes[row_idx, 1].imshow(make_composite(record["maps"], panel, limits))
        axes[row_idx, 1].set_title("Virtual mIF composite", color="white", fontsize=10)
        for axis in axes[row_idx]:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_facecolor("black")
    handles = [mpatches.Patch(color=color, label=marker) for marker, color in panel]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), facecolor="black", labelcolor="white")
    fig.suptitle(title, color="white", fontsize=15, y=0.985)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, mpatches = require_runtime(out_dir / ".matplotlib")
    torch, gigatime_class, snapshot_download, _image, openslide = import_runtime(Path(args.gigatime_repo))
    device = resolve_device(torch, args.device)
    model = load_model(torch, gigatime_class, snapshot_download, device)

    slide_scores = pd.read_csv(args.slide_scores)
    tile_scores = pd.read_csv(args.tile_scores)
    joined = pd.read_csv(args.joined)
    slide_lookup = slide_scores.set_index("slide_id")

    written = []
    for group_label, ref_row in reference_rows(joined):
        slide_id = ref_row["slide_id"]
        slide_path = Path(slide_lookup.loc[slide_id]["slide_path"])
        case_id = ref_row["case_submitter_id"]
        erbb2 = float(ref_row["erbb2_tpm"])
        for panel_name, panel in PANELS.items():
            selected_tiles = select_tiles(tile_scores, slide_id, panel, args.tiles_per_panel)
            records = infer_tile_maps(
                torch,
                model,
                openslide,
                slide_path,
                selected_tiles,
                args.tile_size,
                args.batch_size,
                device,
            )
            prefix = f"{group_label.lower().replace('-', '_')}_{panel_name}"
            montage_path = out_dir / f"{prefix}_virtual_mif_montage.png"
            title = f"{case_id} {group_label} | {panel_name.replace('_', ' ')} | ERBB2 TPM {erbb2:.1f}"
            save_montage(records, panel_name, panel, title, montage_path, plt, mpatches)
            written.append(montage_path)
            if panel_name == "immune_checkpoint":
                comparison_path = out_dir / f"{prefix}_he_vs_virtual_mif.png"
                save_he_vs_mif(records, panel, title, comparison_path, plt, mpatches)
                written.append(comparison_path)

    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
