#!/usr/bin/env python3
"""Render overview images for GigaTIME virtual mIF channel predictions."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


GIGATIME_CHANNELS = [
    "DAPI",
    "TRITC",
    "Cy5",
    "PD-1",
    "CD14",
    "CD4",
    "T-bet",
    "CD34",
    "CD68",
    "CD16",
    "CD11c",
    "CD138",
    "CD20",
    "CD3",
    "CD8",
    "PD-L1",
    "CK",
    "Ki67",
    "Tryptase",
    "Actin-D",
    "Caspase3-D",
    "PHH3-B",
    "Transgelin",
]


CHANNEL_GROUPS = {
    "DAPI": "nuclei",
    "TRITC": "fluorescence channel",
    "Cy5": "fluorescence channel",
    "PD-1": "immune checkpoint",
    "CD14": "myeloid/monocyte",
    "CD4": "T cell",
    "T-bet": "T-cell activation",
    "CD34": "vascular/stromal",
    "CD68": "macrophage",
    "CD16": "NK/myeloid",
    "CD11c": "dendritic/myeloid",
    "CD138": "plasma cell/epithelial",
    "CD20": "B cell",
    "CD3": "T cell",
    "CD8": "cytotoxic T cell",
    "PD-L1": "immune checkpoint",
    "CK": "epithelial/tumor",
    "Ki67": "proliferation",
    "Tryptase": "mast cell",
    "Actin-D": "structural",
    "Caspase3-D": "cell death",
    "PHH3-B": "mitosis",
    "Transgelin": "stromal/smooth muscle",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca_extremes/slide_scores.csv")
    parser.add_argument("--tile-scores", default="results/gigatime_tcga_brca_extremes/tile_scores.csv")
    parser.add_argument("--joined", default="results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv")
    parser.add_argument("--out-dir", default="docs/assets/virtual_mif_channels")
    parser.add_argument("--max-reference-slides", type=int, default=2)
    return parser.parse_args()


def require_runtime(mpl_config_dir: Path):
    cache_dir = Path(tempfile.gettempdir()) / "gigatime_tcga_mplconfig"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    return pd, np, plt


def available_channels(frame, prefix: str = "mean_") -> list[str]:
    return [channel for channel in GIGATIME_CHANNELS if f"{prefix}{channel}" in frame.columns]


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["figure", "path", "caption"])
        writer.writeheader()
        writer.writerows(rows)


def channel_limits(tile_scores, channels: list[str]) -> dict[str, tuple[float, float]]:
    limits = {}
    for channel in channels:
        values = tile_scores[f"mean_{channel}"].astype(float)
        low = float(values.quantile(0.02))
        high = float(values.quantile(0.98))
        if high <= low:
            high = low + 1e-9
        limits[channel] = (low, high)
    return limits


def plot_reference_slide_grid(tile_scores, slide_row, channels: list[str], limits, out_path: Path, plt) -> None:
    slide_id = slide_row["slide_id"]
    case_id = slide_row["case_submitter_id"]
    group = slide_row["her2_group"]
    erbb2 = float(slide_row["erbb2_tpm"])
    slide_tiles = tile_scores[tile_scores["slide_id"] == slide_id].copy()
    ncols = 5
    nrows = 5
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 11), constrained_layout=True)
    axes = axes.flatten()
    for axis, channel in zip(axes, channels):
        values = slide_tiles[f"mean_{channel}"].astype(float)
        low, high = limits[channel]
        axis.scatter(
            slide_tiles["x"].astype(float),
            slide_tiles["y"].astype(float),
            c=values,
            cmap="magma",
            vmin=low,
            vmax=high,
            s=32,
            edgecolors="none",
        )
        axis.invert_yaxis()
        axis.set_aspect("auto")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(f"{channel}\n{CHANNEL_GROUPS.get(channel, '')}", fontsize=10)
    for axis in axes[len(channels) :]:
        axis.axis("off")
    fig.suptitle(
        f"Virtual mIF channel maps for {case_id} ({group}, ERBB2 TPM {erbb2:.1f})\n"
        "Each dot is one sampled H&E tile; brighter colors mean higher predicted activation for that channel.",
        fontsize=15,
    )
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_group_means(slide_scores, joined, channels: list[str], out_path: Path, plt) -> None:
    merged = slide_scores.merge(
        joined[["slide_id", "her2_group", "erbb2_tpm"]],
        on="slide_id",
        how="inner",
        suffixes=("", "_joined"),
    )
    rows = []
    for channel in channels:
        score_col = f"mean_{channel}"
        for group, group_df in merged.groupby("her2_group"):
            rows.append(
                {
                    "channel": channel,
                    "group": group,
                    "mean": float(group_df[score_col].astype(float).mean()),
                }
            )
    import pandas as pd

    plot_df = pd.DataFrame(rows)
    order = (
        plot_df.pivot(index="channel", columns="group", values="mean")
        .assign(delta=lambda df: df.get("HER2-high", 0) - df.get("HER2-low", 0))
        .sort_values("delta")
        .index.tolist()
    )
    y_positions = range(len(order))
    high = plot_df[plot_df["group"] == "HER2-high"].set_index("channel").reindex(order)["mean"]
    low = plot_df[plot_df["group"] == "HER2-low"].set_index("channel").reindex(order)["mean"]

    fig, ax = plt.subplots(figsize=(9, 10))
    offset = 0.18
    ax.barh([y - offset for y in y_positions], low, height=0.34, label="HER2-low / low ERBB2", color="#4C78A8")
    ax.barh([y + offset for y in y_positions], high, height=0.34, label="HER2-high / high ERBB2", color="#F58518")
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(order)
    ax.set_xlabel("Mean virtual mIF activation across processed slides")
    ax.set_title("GigaTIME virtual mIF channel means by ERBB2-expression group")
    ax.legend(loc="lower right")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_slide_channel_matrix(slide_scores, joined, channels: list[str], out_path: Path, np, plt) -> None:
    merged = slide_scores.merge(
        joined[["slide_id", "case_submitter_id", "her2_group", "erbb2_tpm"]],
        on=["slide_id", "case_submitter_id"],
        how="inner",
    ).sort_values(["her2_group", "erbb2_tpm"], ascending=[True, False])
    matrix = merged[[f"mean_{channel}" for channel in channels]].astype(float).to_numpy()
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0)
    stds[stds == 0] = 1
    z = (matrix - means) / stds
    labels = [
        f"{row.case_submitter_id} | {row.her2_group.replace('HER2-', '')} | ERBB2 {row.erbb2_tpm:.0f}"
        for row in merged.itertuples()
    ]
    fig, ax = plt.subplots(figsize=(15, max(6, 0.45 * len(labels))))
    image = ax.imshow(z, aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
    ax.set_xticks(range(len(channels)))
    ax.set_xticklabels(channels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title("Relative virtual mIF activation by slide and channel")
    ax.set_xlabel("GigaTIME virtual mIF channel")
    ax.set_ylabel("Processed TCGA-BRCA slide")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("Within-channel z-score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, np, plt = require_runtime(out_dir / ".matplotlib")

    slide_scores = pd.read_csv(args.slide_scores)
    tile_scores = pd.read_csv(args.tile_scores)
    joined = pd.read_csv(args.joined)
    channels = available_channels(tile_scores)
    if not channels:
        raise ValueError("No GigaTIME mean_* channel columns found in tile scores.")

    manifest_rows: list[dict[str, str]] = []

    group_means_path = out_dir / "virtual_mif_all_channel_group_means.png"
    plot_group_means(slide_scores, joined, channels, group_means_path, plt)
    manifest_rows.append(
        {
            "figure": "All-channel group means",
            "path": str(group_means_path),
            "caption": "Mean GigaTIME virtual mIF activation for all 23 channels, compared between ERBB2-high and ERBB2-low processed slides.",
        }
    )

    matrix_path = out_dir / "virtual_mif_slide_channel_matrix.png"
    plot_slide_channel_matrix(slide_scores, joined, channels, matrix_path, np, plt)
    manifest_rows.append(
        {
            "figure": "Slide by channel matrix",
            "path": str(matrix_path),
            "caption": "Relative activation heatmap showing each processed slide against all GigaTIME virtual mIF channels.",
        }
    )

    limits = channel_limits(tile_scores, channels)
    references = [
        ("HER2-high reference", joined.sort_values("erbb2_tpm", ascending=False).iloc[0]),
        ("HER2-low reference", joined.sort_values("erbb2_tpm", ascending=True).iloc[0]),
    ][: args.max_reference_slides]
    for label, row in references:
        safe_label = label.lower().replace(" ", "_").replace("-", "_")
        out_path = out_dir / f"{safe_label}_all_virtual_mif_channels.png"
        plot_reference_slide_grid(tile_scores, row, channels, limits, out_path, plt)
        manifest_rows.append(
            {
                "figure": label,
                "path": str(out_path),
                "caption": "Spatial tile-level GigaTIME virtual mIF maps for all 23 channels on one representative processed slide.",
            }
        )

    write_manifest(out_dir / "virtual_mif_channel_figure_manifest.csv", manifest_rows)
    print(f"Wrote {len(manifest_rows)} virtual mIF figures to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
