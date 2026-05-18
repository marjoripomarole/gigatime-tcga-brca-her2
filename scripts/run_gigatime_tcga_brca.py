#!/usr/bin/env python3
"""Run GigaTIME inference on TCGA-BRCA diagnostic slide images."""

from __future__ import annotations

import argparse
import csv
import os
import random
import re
import sys
from pathlib import Path

import numpy as np

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

IMAGENET_MEAN = np.array([0.485, 0.456, 0.406], dtype=np.float32)
IMAGENET_STD = np.array([0.229, 0.224, 0.225], dtype=np.float32)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slides-dir", default="data/tcga_brca/slides", help="Directory containing TCGA .svs files.")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca", help="Output directory.")
    parser.add_argument("--gigatime-repo", default="external/GigaTIME", help="Path to the official GigaTIME repo clone.")
    parser.add_argument("--tile-size", type=int, default=256, help="Tile size passed to GigaTIME.")
    parser.add_argument("--level", type=int, default=0, help="OpenSlide pyramid level to tile.")
    parser.add_argument("--tile-stride", type=int, default=256, help="Stride in pixels at selected OpenSlide level.")
    parser.add_argument("--tile-limit", type=int, default=512, help="Maximum tissue tiles per slide. Use 0 for all.")
    parser.add_argument("--batch-size", type=int, default=16, help="Inference batch size.")
    parser.add_argument("--tissue-threshold", type=float, default=0.35, help="Minimum tissue fraction for a tile.")
    parser.add_argument("--activation-threshold", type=float, default=0.5, help="Threshold for per-channel activation fraction.")
    parser.add_argument("--tile-order", default="random", choices=["random", "row-major"], help="Tile traversal order before applying --tile-limit.")
    parser.add_argument("--random-seed", type=int, default=42, help="Seed used when --tile-order=random.")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"], help="Torch device.")
    parser.add_argument("--max-slides", type=int, default=0, help="Maximum slides to process. Use 0 for all.")
    parser.add_argument("--save-tile-csv", action="store_true", help="Write one row per tile.")
    parser.add_argument("--heatmap-channels", default="CD3,CD8,PD-L1,CK", help="Comma-separated channels for tile heatmaps.")
    return parser.parse_args()


def import_runtime(gigatime_repo: Path):
    scripts_dir = gigatime_repo / "scripts"
    if not scripts_dir.exists():
        raise FileNotFoundError(f"Could not find GigaTIME scripts directory: {scripts_dir}")
    sys.path.insert(0, str(scripts_dir.resolve()))
    import torch
    from archs import gigatime
    from huggingface_hub import snapshot_download
    from PIL import Image
    import openslide

    return torch, gigatime, snapshot_download, Image, openslide


def resolve_device(torch, requested: str):
    if requested == "cuda":
        return torch.device("cuda")
    if requested == "mps":
        return torch.device("mps")
    if requested == "cpu":
        return torch.device("cpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_model(torch, gigatime_class, snapshot_download, device):
    model = gigatime_class(num_classes=len(GIGATIME_CHANNELS), input_channels=3)
    local_dir = snapshot_download(repo_id="prov-gigatime/GigaTIME")
    weights_path = Path(local_dir) / "model.pth"
    state_dict = torch.load(weights_path, map_location="cpu")
    model.load_state_dict(state_dict)
    model.to(device)
    model.eval()
    return model


def find_slides(slides_dir: Path) -> list[Path]:
    suffixes = {".svs", ".tif", ".tiff"}
    return sorted(path for path in slides_dir.rglob("*") if path.suffix.lower() in suffixes)


def case_from_slide_path(path: Path) -> str:
    text = str(path)
    match = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", text)
    return match.group(1) if match else path.stem[:12]


def tissue_fraction(rgb: np.ndarray) -> float:
    arr = rgb.astype(np.float32)
    mean_intensity = arr.mean(axis=2)
    chroma = arr.max(axis=2) - arr.min(axis=2)
    tissue = (mean_intensity < 235.0) & (chroma > 8.0)
    return float(tissue.mean())


def preprocess_tile(torch, rgb: np.ndarray, device):
    arr = rgb.astype(np.float32) / 255.0
    arr = (arr - IMAGENET_MEAN) / IMAGENET_STD
    arr = np.transpose(arr, (2, 0, 1))
    return torch.from_numpy(arr).to(device)


def iter_tissue_tiles(
    slide,
    level: int,
    tile_size: int,
    stride: int,
    tissue_threshold: float,
    tile_limit: int,
    tile_order: str,
    random_seed: int,
):
    level_downsample = int(round(float(slide.level_downsamples[level])))
    width, height = slide.level_dimensions[level]
    count = 0
    xs = range(0, max(width - tile_size + 1, 1), stride)
    ys = range(0, max(height - tile_size + 1, 1), stride)
    coordinates = [(x, y) for y in ys for x in xs]
    if tile_order == "random":
        rng = random.Random(random_seed)
        rng.shuffle(coordinates)
    for x, y in coordinates:
        location = (x * level_downsample, y * level_downsample)
        region = slide.read_region(location, level, (tile_size, tile_size)).convert("RGB")
        rgb = np.asarray(region)
        frac = tissue_fraction(rgb)
        if frac < tissue_threshold:
            continue
        yield {"x": x, "y": y, "tissue_fraction": frac, "rgb": rgb}
        count += 1
        if tile_limit and count >= tile_limit:
            return


def summarize_slide(tile_rows: list[dict[str, float | int | str]], slide_path: Path) -> dict[str, float | int | str]:
    row: dict[str, float | int | str] = {
        "slide_path": str(slide_path),
        "slide_id": slide_path.stem,
        "case_submitter_id": case_from_slide_path(slide_path),
        "n_tiles": len(tile_rows),
    }
    if not tile_rows:
        return row
    row["mean_tissue_fraction"] = float(np.mean([float(tile["tissue_fraction"]) for tile in tile_rows]))
    for channel in GIGATIME_CHANNELS:
        mean_key = f"mean_{channel}"
        frac_key = f"frac_{channel}"
        row[mean_key] = float(np.mean([float(tile[mean_key]) for tile in tile_rows]))
        row[frac_key] = float(np.mean([float(tile[frac_key]) for tile in tile_rows]))
    return row


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_heatmaps(tile_rows: list[dict[str, float | int | str]], channels: list[str], out_dir: Path, slide_id: str) -> None:
    if not tile_rows:
        return
    mpl_config_dir = out_dir.parent / ".matplotlib"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    xs = [int(tile["x"]) for tile in tile_rows]
    ys = [int(tile["y"]) for tile in tile_rows]
    for channel in channels:
        key = f"mean_{channel}"
        if key not in tile_rows[0]:
            continue
        values = [float(tile[key]) for tile in tile_rows]
        plt.figure(figsize=(6, 5))
        scatter = plt.scatter(xs, ys, c=values, s=12, cmap="viridis")
        plt.gca().invert_yaxis()
        plt.axis("equal")
        plt.xlabel("Tile x")
        plt.ylabel("Tile y")
        plt.title(f"{slide_id} virtual {channel}")
        plt.colorbar(scatter, label="Mean activation")
        plt.tight_layout()
        plt.savefig(out_dir / f"{slide_id}_{channel}.png", dpi=160)
        plt.close()


def run_slide(torch, model, openslide, slide_path: Path, args: argparse.Namespace, device) -> tuple[dict[str, object], list[dict[str, object]]]:
    slide = openslide.OpenSlide(str(slide_path))
    if args.level >= slide.level_count:
        raise ValueError(f"{slide_path} has {slide.level_count} levels, cannot use level {args.level}")
    tile_records: list[dict[str, object]] = []
    batch = []
    batch_meta = []
    with torch.no_grad():
        for tile in iter_tissue_tiles(
            slide,
            args.level,
            args.tile_size,
            args.tile_stride,
            args.tissue_threshold,
            args.tile_limit,
            args.tile_order,
            args.random_seed,
        ):
            batch.append(preprocess_tile(torch, tile["rgb"], device))
            batch_meta.append({key: value for key, value in tile.items() if key != "rgb"})
            if len(batch) == args.batch_size:
                tile_records.extend(infer_batch(torch, model, batch, batch_meta, args.activation_threshold))
                batch = []
                batch_meta = []
        if batch:
            tile_records.extend(infer_batch(torch, model, batch, batch_meta, args.activation_threshold))
    slide.close()
    slide_row = summarize_slide(tile_records, slide_path)
    return slide_row, tile_records


def infer_batch(torch, model, batch, batch_meta, activation_threshold: float) -> list[dict[str, object]]:
    tensor = torch.stack(batch, dim=0)
    output = torch.sigmoid(model(tensor)).detach().cpu().numpy()
    rows: list[dict[str, object]] = []
    for idx, meta in enumerate(batch_meta):
        row: dict[str, object] = dict(meta)
        for channel_idx, channel in enumerate(GIGATIME_CHANNELS):
            channel_map = output[idx, channel_idx]
            row[f"mean_{channel}"] = float(channel_map.mean())
            row[f"frac_{channel}"] = float((channel_map > activation_threshold).mean())
        rows.append(row)
    return rows


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    torch, gigatime_class, snapshot_download, _image, openslide = import_runtime(Path(args.gigatime_repo))
    device = resolve_device(torch, args.device)
    if device.type != "cpu" and not os.environ.get("HF_TOKEN"):
        print("HF_TOKEN is not set; model download may fail if terms are not accepted.", file=sys.stderr)
    print(f"Using device: {device}", file=sys.stderr)
    model = load_model(torch, gigatime_class, snapshot_download, device)

    slides = find_slides(Path(args.slides_dir))
    if args.max_slides:
        slides = slides[: args.max_slides]
    if not slides:
        raise FileNotFoundError(f"No .svs/.tif/.tiff slides found under {args.slides_dir}")

    heatmap_channels = [channel.strip() for channel in args.heatmap_channels.split(",") if channel.strip()]
    slide_rows: list[dict[str, object]] = []
    all_tile_rows: list[dict[str, object]] = []
    for index, slide_path in enumerate(slides, start=1):
        print(f"[{index}/{len(slides)}] Processing {slide_path}", file=sys.stderr)
        slide_row, tile_rows = run_slide(torch, model, openslide, slide_path, args, device)
        slide_rows.append(slide_row)
        for tile_row in tile_rows:
            tile_row["slide_id"] = slide_row["slide_id"]
            tile_row["case_submitter_id"] = slide_row["case_submitter_id"]
        all_tile_rows.extend(tile_rows)
        save_heatmaps(tile_rows, heatmap_channels, out_dir / "heatmaps", str(slide_row["slide_id"]))
        write_csv(out_dir / "slide_scores.csv", slide_rows)
        if args.save_tile_csv:
            write_csv(out_dir / "tile_scores.csv", all_tile_rows)

    print(f"Done. Wrote {out_dir / 'slide_scores.csv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
