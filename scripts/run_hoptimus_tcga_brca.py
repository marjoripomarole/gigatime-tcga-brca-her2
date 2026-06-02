#!/usr/bin/env python3
"""Extract H-Optimus/H0-mini tile embeddings from TCGA-BRCA H&E slides."""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import os
import random
import re
import sys
from pathlib import Path

import numpy as np


HOPTIMUS0_MEAN = np.array([0.707223, 0.578729, 0.703617], dtype=np.float32)
HOPTIMUS0_STD = np.array([0.211883, 0.230117, 0.177517], dtype=np.float32)

MODEL_PRESETS = {
    "hoptimus0": {
        "model_id": "hf-hub:bioptimus/H-optimus-0",
        "embedding_mode": "direct",
        "model_kwargs": {"init_values": 1e-5, "dynamic_img_size": False},
        "manual_normalize": True,
    },
    "h0-mini": {
        "model_id": "hf-hub:bioptimus/H0-mini",
        "embedding_mode": "cls",
        "model_kwargs": "h0-mini",
        "manual_normalize": False,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slides-dir", default="data/tcga_brca/slides")
    parser.add_argument(
        "--slide-table",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
        help="CSV listing slides to process. Defaults to the strict high-trust 171-slide list.",
    )
    parser.add_argument("--slide-path-column", default="slide_local_path")
    parser.add_argument("--missing-slide-policy", choices=["error", "skip"], default="error")
    parser.add_argument("--out-dir", default="results/hoptimus_tcga_brca_high_trust_tile224")
    parser.add_argument("--model-preset", choices=sorted(MODEL_PRESETS), default="h0-mini")
    parser.add_argument("--model-id", default=None, help="Override the preset Hugging Face timm model id.")
    parser.add_argument(
        "--embedding-mode",
        choices=["auto", "direct", "cls", "concat_cls_patch_mean", "patch_mean"],
        default="auto",
        help="How to convert model output to one vector per tile.",
    )
    parser.add_argument("--input-size", type=int, default=224, help="Model input size in pixels.")
    parser.add_argument("--target-mpp", type=float, default=0.5, help="Target microns per pixel for extracted tiles.")
    parser.add_argument(
        "--fallback-mpp",
        type=float,
        default=0.5,
        help="MPP to assume when a slide lacks OpenSlide mpp-x/mpp-y metadata.",
    )
    parser.add_argument("--tile-stride", type=int, default=224, help="Stride in target-MPP pixels.")
    parser.add_argument("--tile-limit", type=int, default=64, help="Maximum tissue tiles per slide. Use 0 for all.")
    parser.add_argument("--tissue-threshold", type=float, default=0.35)
    parser.add_argument("--tile-order", choices=["random", "row-major"], default="random")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    parser.add_argument(
        "--precision",
        choices=["auto", "float32", "float16"],
        default="auto",
        help="Mixed precision mode. Auto uses float16 only on CUDA.",
    )
    parser.add_argument("--max-slides", type=int, default=0, help="Maximum slides to process. Use 0 for all.")
    parser.add_argument("--save-tile-csv", action="store_true", help="Write per-tile embeddings. Can be large.")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check slide discovery and runtime imports without loading model weights.",
    )
    parser.add_argument(
        "--inspect-slide",
        action="store_true",
        help="With --dry-run, open one slide and report the H-Optimus extraction geometry without model inference.",
    )
    return parser.parse_args()


def import_runtime():
    try:
        import torch
        import timm
        from PIL import Image
        import openslide
        from timm.data import resolve_data_config
        from timm.data.transforms_factory import create_transform
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Install the project environment with "
            "`conda env create -f envs/gigatime-tcga.yml`, or add the missing package to "
            "the active `gigatime-tcga` environment."
        ) from exc
    return torch, timm, Image, openslide, resolve_data_config, create_transform


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


def model_kwargs_for_preset(torch, timm, preset: str) -> dict[str, object]:
    raw = MODEL_PRESETS[preset]["model_kwargs"]
    if raw == "h0-mini":
        return {"mlp_layer": timm.layers.SwiGLUPacked, "act_layer": torch.nn.SiLU}
    return dict(raw)


def load_model(torch, timm, preset: str, model_id: str | None, device):
    selected_id = model_id or str(MODEL_PRESETS[preset]["model_id"])
    kwargs = model_kwargs_for_preset(torch, timm, preset)
    try:
        model = timm.create_model(selected_id, pretrained=True, **kwargs)
    except Exception as exc:
        raise SystemExit(
            f"Could not load {selected_id}. H-Optimus models are gated on Hugging Face; "
            "accept the model terms and set HF_TOKEN/HUGGING_FACE_HUB_TOKEN before running."
        ) from exc
    model.to(device)
    model.eval()
    return model, selected_id


def make_transform(preset: str, model, input_size: int, Image, resolve_data_config, create_transform):
    if bool(MODEL_PRESETS[preset]["manual_normalize"]):

        def transform(rgb: np.ndarray):
            image = Image.fromarray(rgb)
            if image.size != (input_size, input_size):
                image = image.resize((input_size, input_size), Image.BILINEAR)
            arr = np.asarray(image).astype(np.float32) / 255.0
            arr = (arr - HOPTIMUS0_MEAN) / HOPTIMUS0_STD
            return np.transpose(arr, (2, 0, 1))

        return transform

    timm_transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))

    def transform(rgb: np.ndarray):
        image = Image.fromarray(rgb)
        if image.size != (input_size, input_size):
            image = image.resize((input_size, input_size), Image.BILINEAR)
        return timm_transform(image).numpy()

    return transform


def find_slides(slides_dir: Path) -> list[Path]:
    suffixes = {".svs", ".tif", ".tiff"}
    return sorted(path for path in slides_dir.rglob("*") if path.suffix.lower() in suffixes)


def find_slides_from_table(path: Path, slide_path_column: str, missing_policy: str) -> list[Path]:
    slides: list[Path] = []
    missing: list[Path] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No columns found in slide table: {path}")
        if slide_path_column not in reader.fieldnames:
            raise ValueError(
                f"Slide table {path} does not contain {slide_path_column!r}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        for row in reader:
            raw_value = (row.get(slide_path_column) or "").strip()
            if not raw_value:
                continue
            slide_path = Path(raw_value)
            if slide_path.exists():
                slides.append(slide_path)
            else:
                missing.append(slide_path)
    if missing:
        message = f"{len(missing)} slide paths from {path} are missing locally."
        if missing_policy == "error":
            example = "\n".join(f"- {slide}" for slide in missing[:10])
            raise FileNotFoundError(f"{message}\n{example}")
        print(f"{message} Skipping missing slides.", file=sys.stderr)
    return list(dict.fromkeys(slides))


def case_from_slide_path(path: Path) -> str:
    match = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", str(path))
    return match.group(1) if match else path.stem[:12]


def tissue_fraction(rgb: np.ndarray) -> float:
    arr = rgb.astype(np.float32)
    mean_intensity = arr.mean(axis=2)
    chroma = arr.max(axis=2) - arr.min(axis=2)
    tissue = (mean_intensity < 235.0) & (chroma > 8.0)
    return float(tissue.mean())


def parse_slide_mpp(slide, fallback_mpp: float) -> tuple[float, float, str]:
    props = slide.properties

    def parse_value(key: str) -> float | None:
        raw = props.get(key)
        if raw is None:
            return None
        try:
            return float(raw)
        except ValueError:
            return None

    mpp_x = parse_value("openslide.mpp-x")
    mpp_y = parse_value("openslide.mpp-y")
    if mpp_x and mpp_y:
        return mpp_x, mpp_y, "openslide"
    return fallback_mpp, fallback_mpp, "fallback"


def iter_tissue_tiles(
    slide,
    input_size: int,
    target_mpp: float,
    fallback_mpp: float,
    tile_stride: int,
    tissue_threshold: float,
    tile_limit: int,
    tile_order: str,
    random_seed: int,
):
    mpp_x, mpp_y, mpp_source = parse_slide_mpp(slide, fallback_mpp)
    read_width = max(1, int(round(input_size * target_mpp / mpp_x)))
    read_height = max(1, int(round(input_size * target_mpp / mpp_y)))
    stride_x = max(1, int(round(tile_stride * target_mpp / mpp_x)))
    stride_y = max(1, int(round(tile_stride * target_mpp / mpp_y)))
    width, height = slide.dimensions
    xs = range(0, max(width - read_width + 1, 1), stride_x)
    ys = range(0, max(height - read_height + 1, 1), stride_y)
    coordinates = [(x, y) for y in ys for x in xs]
    if tile_order == "random":
        rng = random.Random(random_seed)
        rng.shuffle(coordinates)
    count = 0
    for x, y in coordinates:
        region = slide.read_region((x, y), 0, (read_width, read_height)).convert("RGB")
        rgb = np.asarray(region.resize((input_size, input_size)))
        frac = tissue_fraction(rgb)
        if frac < tissue_threshold:
            continue
        yield {
            "x": x,
            "y": y,
            "read_width": read_width,
            "read_height": read_height,
            "mpp_x": mpp_x,
            "mpp_y": mpp_y,
            "mpp_source": mpp_source,
            "tissue_fraction": frac,
            "rgb": rgb,
        }
        count += 1
        if tile_limit and count >= tile_limit:
            return


def slide_extraction_geometry(slide, input_size: int, target_mpp: float, fallback_mpp: float, tile_stride: int):
    mpp_x, mpp_y, mpp_source = parse_slide_mpp(slide, fallback_mpp)
    read_width = max(1, int(round(input_size * target_mpp / mpp_x)))
    read_height = max(1, int(round(input_size * target_mpp / mpp_y)))
    stride_x = max(1, int(round(tile_stride * target_mpp / mpp_x)))
    stride_y = max(1, int(round(tile_stride * target_mpp / mpp_y)))
    width, height = slide.dimensions
    grid_x = len(range(0, max(width - read_width + 1, 1), stride_x))
    grid_y = len(range(0, max(height - read_height + 1, 1), stride_y))
    return {
        "slide_width": width,
        "slide_height": height,
        "mpp_x": mpp_x,
        "mpp_y": mpp_y,
        "mpp_source": mpp_source,
        "read_width_level0": read_width,
        "read_height_level0": read_height,
        "stride_x_level0": stride_x,
        "stride_y_level0": stride_y,
        "grid_x": grid_x,
        "grid_y": grid_y,
        "candidate_tile_positions": grid_x * grid_y,
    }


def selected_embedding_mode(preset: str, requested: str) -> str:
    if requested != "auto":
        return requested
    return str(MODEL_PRESETS[preset]["embedding_mode"])


def features_from_output(torch, output, mode: str):
    if isinstance(output, (tuple, list)):
        output = output[0]
    if output.ndim == 2:
        if mode not in {"auto", "direct"}:
            raise ValueError(f"Embedding mode {mode!r} requires token output, got shape {tuple(output.shape)}")
        return output
    if output.ndim != 3:
        raise ValueError(f"Unexpected model output shape: {tuple(output.shape)}")
    cls = output[:, 0]
    patch_tokens = output[:, 1:]
    if mode == "cls":
        return cls
    if mode == "patch_mean":
        return patch_tokens.mean(dim=1)
    if mode == "concat_cls_patch_mean":
        return torch.cat([cls, patch_tokens.mean(dim=1)], dim=-1)
    if mode == "direct":
        raise ValueError(f"Embedding mode 'direct' cannot handle token output shape {tuple(output.shape)}")
    return cls


def autocast_context(torch, device, precision: str):
    if precision == "float32":
        return contextlib.nullcontext()
    if precision == "float16" or (precision == "auto" and device.type == "cuda"):
        return torch.autocast(device_type=device.type, dtype=torch.float16)
    return contextlib.nullcontext()


def infer_batch(torch, model, batch, batch_meta, embedding_mode: str, device, precision: str) -> list[dict[str, object]]:
    tensor = torch.stack(batch, dim=0).to(device)
    with autocast_context(torch, device, precision):
        output = model(tensor)
        features = features_from_output(torch, output, embedding_mode)
    features_np = features.detach().cpu().numpy()
    rows: list[dict[str, object]] = []
    for idx, meta in enumerate(batch_meta):
        row: dict[str, object] = dict(meta)
        for feature_idx, value in enumerate(features_np[idx]):
            row[f"embedding_{feature_idx:04d}"] = float(value)
        rows.append(row)
    return rows


def run_slide(torch, model, openslide, transform, slide_path: Path, args: argparse.Namespace, device) -> tuple[dict[str, object], list[dict[str, object]]]:
    slide = openslide.OpenSlide(str(slide_path))
    tile_rows: list[dict[str, object]] = []
    batch = []
    batch_meta = []
    embedding_mode = selected_embedding_mode(args.model_preset, args.embedding_mode)
    with torch.inference_mode():
        for tile in iter_tissue_tiles(
            slide,
            args.input_size,
            args.target_mpp,
            args.fallback_mpp,
            args.tile_stride,
            args.tissue_threshold,
            args.tile_limit,
            args.tile_order,
            args.random_seed,
        ):
            batch.append(torch.from_numpy(transform(tile["rgb"])))
            batch_meta.append({key: value for key, value in tile.items() if key != "rgb"})
            if len(batch) == args.batch_size:
                tile_rows.extend(infer_batch(torch, model, batch, batch_meta, embedding_mode, device, args.precision))
                batch = []
                batch_meta = []
        if batch:
            tile_rows.extend(infer_batch(torch, model, batch, batch_meta, embedding_mode, device, args.precision))
    slide.close()
    slide_row = summarize_slide(tile_rows, slide_path)
    return slide_row, tile_rows


def summarize_slide(tile_rows: list[dict[str, object]], slide_path: Path) -> dict[str, object]:
    row: dict[str, object] = {
        "slide_path": str(slide_path),
        "slide_id": slide_path.stem,
        "case_submitter_id": case_from_slide_path(slide_path),
        "n_tiles": len(tile_rows),
    }
    if not tile_rows:
        return row
    row["mean_tissue_fraction"] = float(np.mean([float(tile["tissue_fraction"]) for tile in tile_rows]))
    row["mpp_x"] = tile_rows[0].get("mpp_x", "")
    row["mpp_y"] = tile_rows[0].get("mpp_y", "")
    row["mpp_source"] = tile_rows[0].get("mpp_source", "")
    embedding_keys = [key for key in tile_rows[0] if key.startswith("embedding_")]
    for key in embedding_keys:
        row[key] = float(np.mean([float(tile[key]) for tile in tile_rows]))
    return row


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
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
        return [dict(row) for row in csv.DictReader(handle)]


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    torch, timm, Image, openslide, resolve_data_config, create_transform = import_runtime()
    device = resolve_device(torch, args.device)

    if args.slide_table:
        slides = find_slides_from_table(Path(args.slide_table), args.slide_path_column, args.missing_slide_policy)
    else:
        slides = find_slides(Path(args.slides_dir))
    if args.max_slides:
        slides = slides[: args.max_slides]
    if not slides:
        raise FileNotFoundError("No local slides found to process.")

    if args.dry_run:
        print(f"Found {len(slides)} slides. Runtime imports succeeded. Dry run did not load model weights.")
        if args.inspect_slide:
            slide_path = slides[0]
            slide = openslide.OpenSlide(str(slide_path))
            geometry = slide_extraction_geometry(
                slide,
                args.input_size,
                args.target_mpp,
                args.fallback_mpp,
                args.tile_stride,
            )
            slide.close()
            print(f"Inspected slide: {slide_path}")
            for key, value in geometry.items():
                print(f"{key}: {value}")
        return 0

    if not os.environ.get("HF_TOKEN") and not os.environ.get("HUGGING_FACE_HUB_TOKEN"):
        print("HF_TOKEN/HUGGING_FACE_HUB_TOKEN is not set; gated H-Optimus download may fail.", file=sys.stderr)
    print(f"Using device: {device}", file=sys.stderr)
    model, selected_id = load_model(torch, timm, args.model_preset, args.model_id, device)
    transform = make_transform(args.model_preset, model, args.input_size, Image, resolve_data_config, create_transform)
    embedding_mode = selected_embedding_mode(args.model_preset, args.embedding_mode)

    slide_embeddings_path = out_dir / "slide_embeddings.csv"
    tile_embeddings_path = out_dir / "tile_embeddings.csv"
    slide_rows = read_existing_csv(slide_embeddings_path) if args.resume else []
    all_tile_rows = read_existing_csv(tile_embeddings_path) if args.resume and args.save_tile_csv else []
    processed_slide_ids = {str(row.get("slide_id", "")) for row in slide_rows if row.get("slide_id")}

    for index, slide_path in enumerate(slides, start=1):
        if args.resume and slide_path.stem in processed_slide_ids:
            print(f"[{index}/{len(slides)}] Skipping existing {slide_path}", file=sys.stderr)
            continue
        print(f"[{index}/{len(slides)}] Processing {slide_path}", file=sys.stderr)
        slide_row, tile_rows = run_slide(torch, model, openslide, transform, slide_path, args, device)
        slide_rows.append(slide_row)
        for tile_row in tile_rows:
            tile_row["slide_id"] = slide_row["slide_id"]
            tile_row["case_submitter_id"] = slide_row["case_submitter_id"]
        all_tile_rows.extend(tile_rows)
        write_csv(slide_embeddings_path, slide_rows)
        if args.save_tile_csv:
            write_csv(tile_embeddings_path, all_tile_rows)

    summary = {
        "model_id": selected_id,
        "model_preset": args.model_preset,
        "embedding_mode": embedding_mode,
        "input_size": args.input_size,
        "target_mpp": args.target_mpp,
        "tile_limit": args.tile_limit,
        "n_slides": len(slide_rows),
        "embedding_dimensions": len([key for key in slide_rows[-1] if key.startswith("embedding_")]) if slide_rows else 0,
        "slide_embeddings": str(slide_embeddings_path),
        "tile_embeddings": str(tile_embeddings_path) if args.save_tile_csv else "",
    }
    (out_dir / "hoptimus_embedding_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Done. Wrote {slide_embeddings_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
