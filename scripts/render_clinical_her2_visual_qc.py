#!/usr/bin/env python3
"""Render visual QC panels for the clinical HER2 GigaTIME pilot."""

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

GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
QC_MARKERS = ["CD68", "PD-L1", "CD11c"]

PANELS = {
    "immune_checkpoint": [
        ("DAPI", "#2D5BFF"),
        ("CK", "#FF3B30"),
        ("CD3", "#00E676"),
        ("CD8", "#00D5FF"),
        ("PD-L1", "#FF2DCE"),
        ("PD-1", "#FFD60A"),
    ],
    "myeloid_checkpoint": [
        ("DAPI", "#2D5BFF"),
        ("CD68", "#FF9500"),
        ("CD11c", "#64D2FF"),
        ("PD-L1", "#FF2DCE"),
        ("CD20", "#BF5AF2"),
        ("CD4", "#30D158"),
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--joined",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/joined_slide_clinical_her2_gigatime.csv",
    )
    parser.add_argument("--tile-scores", default="results/gigatime_tcga_brca_clinical_her2/tile_scores.csv")
    parser.add_argument("--out-dir", default="docs/assets/clinical_her2_visual_qc")
    parser.add_argument("--gigatime-repo", default="external/GigaTIME")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--tiles-per-case", type=int, default=6)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--thumbnail-width", type=int, default=1600)
    return parser.parse_args()


def require_runtime():
    cache_dir = Path(tempfile.gettempdir()) / "gigatime_tcga_visual_qc_mplconfig"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import pandas as pd

    return pd, plt, mpatches


def slug(value: str) -> str:
    return value.lower().replace(" ", "_").replace("-", "_")


def hex_to_rgb(hex_color: str) -> np.ndarray:
    clean = hex_color.lstrip("#")
    return np.array([int(clean[i : i + 2], 16) for i in (0, 2, 4)], dtype=np.float32) / 255.0


def add_qc_signal(frame, markers: list[str]):
    frame = frame.copy()
    score_cols = [f"mean_{marker}" for marker in markers if f"mean_{marker}" in frame.columns]
    if not score_cols:
        raise ValueError(f"None of the requested marker columns are present: {markers}")
    frame["qc_signal"] = frame[score_cols].astype(float).sum(axis=1)
    return frame


def select_cases(joined, markers: list[str]):
    joined = add_qc_signal(joined, markers)
    rows = []
    for group in GROUP_ORDER:
        group_rows = joined[joined["clinical_her2_group"] == group].copy()
        if group_rows.empty:
            continue
        rows.append(group_rows.sort_values(["qc_signal", "case_submitter_id"], ascending=[False, True]).iloc[0])
    return rows


def select_tiles(tile_scores, slide_id: str, markers: list[str], n_tiles: int):
    slide_tiles = tile_scores[tile_scores["slide_id"] == slide_id].copy()
    if slide_tiles.empty:
        raise ValueError(f"No tile scores found for slide_id={slide_id}")
    slide_tiles = add_qc_signal(slide_tiles, markers)
    return slide_tiles.sort_values(["qc_signal", "tissue_fraction"], ascending=False).head(n_tiles)


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
            meta.append(
                {
                    "x": x,
                    "y": y,
                    "he_rgb": rgb,
                    "qc_signal": float(row["qc_signal"]),
                    "mean_CD68": float(row.get("mean_CD68", float("nan"))),
                    "mean_PD-L1": float(row.get("mean_PD-L1", float("nan"))),
                    "mean_CD11c": float(row.get("mean_CD11c", float("nan"))),
                }
            )
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
    return [{**item, "maps": maps[idx]} for idx, item in enumerate(meta)]


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


def save_qc_panel(records, case_row, out_path: Path, plt, mpatches) -> None:
    panel_limits_by_name = {name: panel_limits(records, panel) for name, panel in PANELS.items()}
    nrows = len(records)
    fig, axes = plt.subplots(nrows=nrows, ncols=3, figsize=(9.6, max(8, nrows * 2.45)), facecolor="black")
    axes = np.atleast_2d(axes)
    for row_idx, record in enumerate(records):
        axes[row_idx, 0].imshow(record["he_rgb"])
        axes[row_idx, 0].set_title(
            f"H&E\nCD68 {record['mean_CD68']:.3f}  PD-L1 {record['mean_PD-L1']:.3f}  CD11c {record['mean_CD11c']:.3f}",
            color="white",
            fontsize=8,
        )
        for col_idx, (panel_name, panel) in enumerate(PANELS.items(), start=1):
            axes[row_idx, col_idx].imshow(make_composite(record["maps"], panel, panel_limits_by_name[panel_name]))
            axes[row_idx, col_idx].set_title(panel_name.replace("_", " "), color="white", fontsize=10)
        for axis in axes[row_idx]:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_facecolor("black")
            for spine in axis.spines.values():
                spine.set_edgecolor("#333333")
    handles = []
    seen = set()
    for panel in PANELS.values():
        for marker, color in panel:
            if marker in seen:
                continue
            handles.append(mpatches.Patch(color=color, label=marker))
            seen.add(marker)
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 8), facecolor="black", labelcolor="white")
    title = (
        f"{case_row['case_submitter_id']} | {case_row['clinical_her2_group']} | "
        f"top sampled tiles by CD68 + PD-L1 + CD11c"
    )
    fig.suptitle(title, color="white", fontsize=15, y=0.992)
    fig.text(
        0.5,
        0.035,
        "Virtual mIF composites are GigaTIME predictions from H&E tiles, not experimental mIF.",
        color="white",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.965))
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def scale_rect(x: float, y: float, tile_size: int, scale: float) -> tuple[int, int, int, int]:
    return (
        int(round(x * scale)),
        int(round(y * scale)),
        int(round((x + tile_size) * scale)),
        int(round((y + tile_size) * scale)),
    )


def save_sampled_tile_overlay(openslide, tile_rows, slide_path: Path, case_row, tile_size: int, thumbnail_width: int, out_path: Path, plt):
    slide = openslide.OpenSlide(str(slide_path))
    slide_width, slide_height = slide.dimensions
    thumbnail_height = max(1, int(round(thumbnail_width * slide_height / slide_width)))
    thumbnail = slide.get_thumbnail((thumbnail_width, thumbnail_height)).convert("RGB")
    slide.close()

    scale = thumbnail.width / slide_width
    values = tile_rows["qc_signal"].astype(float)
    norm = plt.Normalize(vmin=float(values.min()), vmax=float(values.max()) if values.max() > values.min() else float(values.min()) + 1e-9)
    cmap = plt.get_cmap("magma")
    fig, ax = plt.subplots(figsize=(10.5, max(3.5, 10.5 * thumbnail.height / thumbnail.width)))
    ax.imshow(thumbnail)
    for _, tile in tile_rows.iterrows():
        x0, y0, x1, y1 = scale_rect(tile["x"], tile["y"], tile_size, scale)
        color = cmap(norm(float(tile["qc_signal"])))
        rect = plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color, edgecolor="white", linewidth=0.7, alpha=0.62)
        ax.add_patch(rect)
    ax.set_title(
        f"{case_row['case_submitter_id']} | {case_row['clinical_her2_group']} | sampled tile signal overlay"
    )
    ax.axis("off")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("CD68 + PD-L1 + CD11c tile score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, mpatches = require_runtime()
    torch, gigatime_class, snapshot_download, _image, openslide = import_runtime(Path(args.gigatime_repo))
    device = resolve_device(torch, args.device)
    model = load_model(torch, gigatime_class, snapshot_download, device)

    joined = pd.read_csv(args.joined)
    tile_scores = pd.read_csv(args.tile_scores)
    selected_cases = select_cases(joined, QC_MARKERS)
    manifest_rows = []
    selected_rows = []
    for case_row in selected_cases:
        slide_id = case_row["slide_id"]
        slide_path = Path(case_row["slide_path"])
        selected_tiles = select_tiles(tile_scores, slide_id, QC_MARKERS, args.tiles_per_case)
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
        prefix = f"{slug(case_row['clinical_her2_group'])}_{case_row['case_submitter_id']}"
        qc_path = out_dir / f"{prefix}_he_vs_virtual_mif_qc.png"
        overlay_path = out_dir / f"{prefix}_sampled_tile_overlay.png"
        save_qc_panel(records, case_row, qc_path, plt, mpatches)
        save_sampled_tile_overlay(
            openslide,
            selected_tiles,
            slide_path,
            case_row,
            args.tile_size,
            args.thumbnail_width,
            overlay_path,
            plt,
        )
        manifest_rows.extend(
            [
                {
                    "case_submitter_id": case_row["case_submitter_id"],
                    "clinical_her2_group": case_row["clinical_her2_group"],
                    "image_type": "he_vs_virtual_mif_qc",
                    "path": str(qc_path),
                },
                {
                    "case_submitter_id": case_row["case_submitter_id"],
                    "clinical_her2_group": case_row["clinical_her2_group"],
                    "image_type": "sampled_tile_overlay",
                    "path": str(overlay_path),
                },
            ]
        )
        selected_rows.append(
            {
                "case_submitter_id": case_row["case_submitter_id"],
                "clinical_her2_group": case_row["clinical_her2_group"],
                "slide_id": slide_id,
                "qc_signal": float(case_row["qc_signal"]),
                "mean_CD68": float(case_row["mean_CD68"]),
                "mean_PD-L1": float(case_row["mean_PD-L1"]),
                "mean_CD11c": float(case_row["mean_CD11c"]),
                "selected_tile_count": int(len(selected_tiles)),
                "qc_panel": str(qc_path),
                "sampled_tile_overlay": str(overlay_path),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(out_dir / "clinical_her2_visual_qc_manifest.csv", index=False)
    pd.DataFrame(selected_rows).to_csv(out_dir / "clinical_her2_visual_qc_selected_cases.csv", index=False)
    print(f"Wrote visual QC assets to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
