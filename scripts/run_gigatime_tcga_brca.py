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
from typing import Any

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
SLIDE_SUFFIXES = {".svs", ".tif", ".tiff", ".ndpi", ".mrxs", ".jpg", ".jpeg", ".png"}
PIL_SUFFIXES = {".jpg", ".jpeg", ".png"}
DEFAULT_METADATA_COLUMNS = (
    "patient_id,clinical_her2_group,her2_ihc,her2_status,grade,ER,PR,"
    "ki67,molecular_subtype,aln_status,tumor_size,age"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slides-dir", default="data/tcga_brca/slides", help="Directory containing TCGA .svs files.")
    parser.add_argument("--slide-table", default=None, help="Optional CSV listing slides to process.")
    parser.add_argument("--slide-path-column", default="slide_local_path", help="Column in --slide-table containing slide paths.")
    parser.add_argument(
        "--metadata-columns",
        default=DEFAULT_METADATA_COLUMNS,
        help="Comma-separated slide-table columns to copy into slide_scores.csv when --slide-table is used.",
    )
    parser.add_argument(
        "--missing-slide-policy",
        default="error",
        choices=["error", "skip"],
        help="What to do when --slide-table includes paths not present locally.",
    )
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
    parser.add_argument(
        "--slide-backend",
        default="auto",
        choices=["auto", "openslide", "pil"],
        help="Slide reader. Use pil for BCNB full-slide .jpg files; auto chooses by suffix and falls back when possible.",
    )
    parser.add_argument("--max-slides", type=int, default=0, help="Maximum slides to process. Use 0 for all.")
    parser.add_argument("--save-tile-csv", action="store_true", help="Write one row per tile.")
    parser.add_argument("--heatmap-channels", default="CD3,CD8,PD-L1,CK", help="Comma-separated channels for tile heatmaps.")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="If output CSVs already exist, keep existing rows and skip slides already present in slide_scores.csv.",
    )
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

    Image.MAX_IMAGE_PIXELS = None
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
    return sorted(path for path in slides_dir.rglob("*") if path.suffix.lower() in SLIDE_SUFFIXES)


def parse_metadata_columns(raw: str) -> list[str]:
    return [column.strip() for column in raw.split(",") if column.strip()]


def find_slide_records_from_table(
    path: Path,
    slide_path_column: str,
    metadata_columns: list[str],
    missing_policy: str,
) -> list[tuple[Path, dict[str, str]]]:
    slide_records: list[tuple[Path, dict[str, str]]] = []
    seen: set[Path] = set()
    missing: list[Path] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No columns found in slide table: {path}")
        if slide_path_column not in reader.fieldnames:
            raise ValueError(
                f"Slide table {path} does not contain column {slide_path_column!r}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        for row in reader:
            raw_value = (row.get(slide_path_column) or "").strip()
            if not raw_value:
                continue
            slide_path = Path(raw_value)
            if slide_path.exists():
                if slide_path in seen:
                    continue
                seen.add(slide_path)
                metadata = {
                    column: row.get(column, "")
                    for column in metadata_columns
                    if column in reader.fieldnames and column != slide_path_column
                }
                slide_records.append((slide_path, metadata))
            else:
                missing.append(slide_path)
    if missing:
        message = f"{len(missing)} slide paths from {path} are missing locally."
        if missing_policy == "error":
            example = "\n".join(f"- {slide}" for slide in missing[:10])
            raise FileNotFoundError(f"{message}\n{example}")
        print(f"{message} Skipping missing slides.", file=sys.stderr)
    return slide_records


def case_from_slide_path(path: Path) -> str:
    text = str(path)
    match = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", text)
    if match:
        return match.group(1)
    stem = path.stem
    bcnb_match = re.match(r"(?:patient|case|sample|slide|bcnb|p)?[_-]?0*(\d{1,4})(?:[_\-.]|$)", stem, flags=re.IGNORECASE)
    if bcnb_match:
        return str(int(bcnb_match.group(1)))
    return stem[:64]


class RasterImageSlide:
    """Tiny OpenSlide-like adapter for large flat RGB images such as BCNB JPG WSIs."""

    def __init__(self, image_module: Any, path: Path):
        self.path = path
        self.image = image_module.open(path)
        self.level_count = 1
        self.level_dimensions = [self.image.size]
        self.level_downsamples = [1.0]

    def read_region(self, location: tuple[int, int], level: int, size: tuple[int, int]):
        if level != 0:
            raise ValueError(f"{self.path} has one raster level, cannot use level {level}")
        x, y = location
        width, height = size
        return self.image.crop((x, y, x + width, y + height)).convert("RGBA")

    def close(self) -> None:
        self.image.close()


def open_slide_any(openslide, Image, slide_path: Path, backend: str):
    suffix = slide_path.suffix.lower()
    if backend == "pil" or (backend == "auto" and suffix in PIL_SUFFIXES):
        return RasterImageSlide(Image, slide_path), "pil"
    if backend in {"auto", "openslide"}:
        try:
            return openslide.OpenSlide(str(slide_path)), "openslide"
        except Exception:
            if backend == "openslide":
                raise
            return RasterImageSlide(Image, slide_path), "pil"
    raise ValueError(f"Unsupported slide backend: {backend}")


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


def read_existing_csv(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


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


def run_slide(torch, model, openslide, Image, slide_path: Path, args: argparse.Namespace, device) -> tuple[dict[str, object], list[dict[str, object]]]:
    slide, backend = open_slide_any(openslide, Image, slide_path, args.slide_backend)
    if args.level >= slide.level_count:
        raise ValueError(f"{slide_path} has {slide.level_count} levels, cannot use level {args.level}")
    slide_width = int(slide.level_dimensions[args.level][0])
    slide_height = int(slide.level_dimensions[args.level][1])
    tile_records: list[dict[str, object]] = []
    batch = []
    batch_meta = []
    try:
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
    finally:
        slide.close()
    slide_row = summarize_slide(tile_records, slide_path)
    slide_row["slide_backend"] = backend
    slide_row["slide_width"] = slide_width
    slide_row["slide_height"] = slide_height
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
    torch, gigatime_class, snapshot_download, Image, openslide = import_runtime(Path(args.gigatime_repo))
    device = resolve_device(torch, args.device)
    if device.type != "cpu" and not os.environ.get("HF_TOKEN"):
        print("HF_TOKEN is not set; model download may fail if terms are not accepted.", file=sys.stderr)
    print(f"Using device: {device}", file=sys.stderr)
    model = load_model(torch, gigatime_class, snapshot_download, device)

    if args.slide_table:
        slide_records = find_slide_records_from_table(
            Path(args.slide_table),
            args.slide_path_column,
            parse_metadata_columns(args.metadata_columns),
            args.missing_slide_policy,
        )
    else:
        slide_records = [(slide_path, {}) for slide_path in find_slides(Path(args.slides_dir))]
    if args.max_slides:
        slide_records = slide_records[: args.max_slides]
    if not slide_records:
        source = args.slide_table or args.slides_dir
        raise FileNotFoundError(f"No supported slide files found from {source}")

    heatmap_channels = [channel.strip() for channel in args.heatmap_channels.split(",") if channel.strip()]
    slide_scores_path = out_dir / "slide_scores.csv"
    tile_scores_path = out_dir / "tile_scores.csv"
    slide_rows: list[dict[str, object]] = read_existing_csv(slide_scores_path) if args.resume else []
    all_tile_rows: list[dict[str, object]] = read_existing_csv(tile_scores_path) if args.resume and args.save_tile_csv else []
    processed_slide_ids = {str(row.get("slide_id", "")) for row in slide_rows if row.get("slide_id")}
    for index, (slide_path, metadata) in enumerate(slide_records, start=1):
        if args.resume and slide_path.stem in processed_slide_ids:
            print(f"[{index}/{len(slide_records)}] Skipping existing {slide_path}", file=sys.stderr)
            continue
        print(f"[{index}/{len(slide_records)}] Processing {slide_path}", file=sys.stderr)
        slide_row, tile_rows = run_slide(torch, model, openslide, Image, slide_path, args, device)
        for key, value in metadata.items():
            if key not in slide_row:
                slide_row[key] = value
        slide_rows.append(slide_row)
        for tile_row in tile_rows:
            tile_row["slide_id"] = slide_row["slide_id"]
            tile_row["case_submitter_id"] = slide_row["case_submitter_id"]
        all_tile_rows.extend(tile_rows)
        save_heatmaps(tile_rows, heatmap_channels, out_dir / "heatmaps", str(slide_row["slide_id"]))
        write_csv(slide_scores_path, slide_rows)
        if args.save_tile_csv:
            write_csv(tile_scores_path, all_tile_rows)

    print(f"Done. Wrote {slide_scores_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
