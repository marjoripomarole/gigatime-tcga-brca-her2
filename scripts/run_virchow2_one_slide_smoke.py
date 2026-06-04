#!/usr/bin/env python3
"""Run a bounded Virchow2 one-slide embedding smoke test."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from run_hoptimus_tcga_brca import (
    autocast_context,
    case_from_slide_path,
    find_slides_from_table,
    import_runtime,
    iter_tissue_tiles,
    resolve_device,
    slide_extraction_geometry,
    write_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-table",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument("--slide-path-column", default="slide_local_path")
    parser.add_argument("--missing-slide-policy", choices=["error", "skip"], default="error")
    parser.add_argument("--out-dir", default="results/virchow2_one_slide_smoke")
    parser.add_argument("--model-id", default="hf-hub:paige-ai/Virchow2")
    parser.add_argument("--max-slides", type=int, default=1)
    parser.add_argument("--tile-limit", type=int, default=1)
    parser.add_argument("--input-size", type=int, default=224)
    parser.add_argument("--target-mpp", type=float, default=0.5)
    parser.add_argument("--fallback-mpp", type=float, default=0.5)
    parser.add_argument("--tile-stride", type=int, default=224)
    parser.add_argument("--tissue-threshold", type=float, default=0.35)
    parser.add_argument("--tile-order", choices=["random", "row-major"], default="random")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    parser.add_argument("--precision", choices=["auto", "float32", "float16"], default="auto")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def load_virchow2(torch, timm, model_id: str, device):
    model = timm.create_model(
        model_id,
        pretrained=True,
        mlp_layer=timm.layers.SwiGLUPacked,
        act_layer=torch.nn.SiLU,
    )
    model.to(device)
    model.eval()
    return model


def make_transform(model, input_size: int, Image, resolve_data_config, create_transform):
    timm_transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))

    def transform(rgb: np.ndarray):
        image = Image.fromarray(rgb)
        if image.size != (input_size, input_size):
            image = image.resize((input_size, input_size), Image.BILINEAR)
        return timm_transform(image)

    return transform


def virchow2_embedding(torch, output):
    """Official Virchow2 tile embedding: class token + mean patch tokens.

    Virchow2 returns 261 tokens: class token, 4 register tokens, and 256 patch
    tokens. The model card recommends skipping register tokens and concatenating
    the class token with the mean pooled patch tokens, giving 2560 dimensions.
    """
    if isinstance(output, (tuple, list)):
        output = output[0]
    if output.ndim == 2:
        return output
    if output.ndim != 3:
        raise ValueError(f"Unexpected Virchow2 output shape: {tuple(output.shape)}")
    if output.shape[1] < 6:
        raise ValueError(f"Virchow2 token output is too short: {tuple(output.shape)}")
    class_token = output[:, 0]
    patch_tokens = output[:, 5:]
    return torch.cat([class_token, patch_tokens.mean(dim=1)], dim=-1)


def write_one_vector_csv(path: Path, metadata: dict[str, object], vector: np.ndarray) -> None:
    row = dict(metadata)
    for idx, value in enumerate(vector):
        row[f"embedding_{idx:04d}"] = float(value)
    write_csv(path, [row])


def main() -> int:
    args = parse_args()
    torch, timm, Image, openslide, resolve_data_config, create_transform = import_runtime()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    slides = find_slides_from_table(Path(args.slide_table), args.slide_path_column, args.missing_slide_policy)
    if args.max_slides:
        slides = slides[: args.max_slides]
    if not slides:
        raise SystemExit("No slides found for Virchow2 smoke test.")
    slide_path = slides[0]

    slide = openslide.OpenSlide(str(slide_path))
    geometry = slide_extraction_geometry(slide, args.input_size, args.target_mpp, args.fallback_mpp, args.tile_stride)
    first_tile = None
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
        first_tile = tile
        break
    slide.close()
    if first_tile is None:
        raise SystemExit(f"No tissue tile found in {slide_path}.")

    tile_metadata = {key: value for key, value in first_tile.items() if key != "rgb"}
    tile_metadata.update({"slide_id": slide_path.stem, "case_submitter_id": case_from_slide_path(slide_path)})
    Image.fromarray(first_tile["rgb"]).save(out_dir / "tile.png")

    if args.dry_run:
        summary = {
            "model": "Virchow2 one-slide smoke dry run",
            "model_id": args.model_id,
            "slide_path": str(slide_path),
            "geometry": geometry,
            "tile": tile_metadata,
            "outputs": {"tile_png": str(out_dir / "tile.png")},
        }
        (out_dir / "virchow2_one_slide_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 0

    device = resolve_device(torch, args.device)
    model = load_virchow2(torch, timm, args.model_id, device)
    transform = make_transform(model, args.input_size, Image, resolve_data_config, create_transform)
    tensor = transform(first_tile["rgb"]).unsqueeze(0).to(device)
    with torch.inference_mode(), autocast_context(torch, device, args.precision):
        output = model(tensor)
        embedding = virchow2_embedding(torch, output)
    output_shape = list(output.shape) if hasattr(output, "shape") else []
    vector = embedding.detach().cpu().float().numpy()[0]

    write_one_vector_csv(out_dir / "tile_embedding.csv", tile_metadata, vector)
    write_one_vector_csv(out_dir / "slide_embedding.csv", {**tile_metadata, "n_tiles": 1}, vector)

    summary = {
        "model": "Virchow2 one-slide smoke",
        "model_id": args.model_id,
        "slide_path": str(slide_path),
        "device": str(device),
        "precision": args.precision,
        "geometry": geometry,
        "tile": tile_metadata,
        "raw_output_shape": output_shape,
        "embedding_recipe": "concat(class_token, mean(patch_tokens_after_4_register_tokens))",
        "embedding_dimensions": int(vector.size),
        "embedding_stats": {
            "mean": float(np.mean(vector)),
            "std": float(np.std(vector)),
            "min": float(np.min(vector)),
            "max": float(np.max(vector)),
        },
        "embedding_first_12": [float(value) for value in vector[:12]],
        "outputs": {
            "tile_png": str(out_dir / "tile.png"),
            "tile_embedding_csv": str(out_dir / "tile_embedding.csv"),
            "slide_embedding_csv": str(out_dir / "slide_embedding.csv"),
            "summary_json": str(out_dir / "virchow2_one_slide_summary.json"),
        },
    }
    (out_dir / "virchow2_one_slide_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
