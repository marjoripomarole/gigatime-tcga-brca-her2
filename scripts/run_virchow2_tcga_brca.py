#!/usr/bin/env python3
"""Extract Virchow2 mean-pooled slide embeddings from TCGA-BRCA H&E slides.

This is the cohort-scale companion to ``run_virchow2_one_slide_smoke.py``. It
reuses the H-Optimus runner's slide tiling, tissue filtering, and slide-mean
pooling, but loads Virchow2 and uses the official Virchow2 tile embedding
(class token concatenated with the mean of the patch tokens after the four
register tokens, giving 2560 dimensions per tile). Output matches the H-Optimus
runner format (`slide_embeddings.csv` with one mean-pooled row per slide), so
the same `analyze_hoptimus_embedding_control.py` control consumes it directly.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from run_hoptimus_tcga_brca import (
    autocast_context,
    find_slides,
    find_slides_from_table,
    import_runtime,
    iter_tissue_tiles,
    read_existing_csv,
    resolve_device,
    summarize_slide,
    write_csv,
)
from run_virchow2_one_slide_smoke import load_virchow2, make_transform, virchow2_embedding


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
    parser.add_argument("--out-dir", default="results/virchow2_tcga_brca_high_trust_tile128")
    parser.add_argument("--model-id", default="hf-hub:paige-ai/Virchow2")
    parser.add_argument("--input-size", type=int, default=224)
    parser.add_argument("--target-mpp", type=float, default=0.5)
    parser.add_argument("--fallback-mpp", type=float, default=0.5)
    parser.add_argument("--tile-stride", type=int, default=224)
    parser.add_argument("--tile-limit", type=int, default=128, help="Maximum tissue tiles per slide. Use 0 for all.")
    parser.add_argument("--tissue-threshold", type=float, default=0.35)
    parser.add_argument("--tile-order", choices=["random", "row-major"], default="random")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    parser.add_argument("--precision", choices=["auto", "float32", "float16"], default="auto")
    parser.add_argument("--max-slides", type=int, default=0, help="Maximum slides to process. Use 0 for all.")
    parser.add_argument("--save-tile-csv", action="store_true", help="Write per-tile embeddings. Can be large.")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def infer_batch(torch, model, transform, batch_rgb, batch_meta, device, precision):
    tensors = [transform(rgb) for rgb in batch_rgb]
    tensor = torch.stack(tensors, dim=0).to(device)
    with autocast_context(torch, device, precision):
        output = model(tensor)
        features = virchow2_embedding(torch, output)
    features_np = features.detach().cpu().float().numpy()
    rows = []
    for idx, meta in enumerate(batch_meta):
        row = dict(meta)
        for feature_idx, value in enumerate(features_np[idx]):
            row[f"embedding_{feature_idx:04d}"] = float(value)
        rows.append(row)
    return rows


def run_slide(torch, model, openslide, transform, slide_path: Path, args, device):
    slide = openslide.OpenSlide(str(slide_path))
    tile_rows = []
    batch_rgb = []
    batch_meta = []
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
            batch_rgb.append(tile["rgb"])
            batch_meta.append({key: value for key, value in tile.items() if key != "rgb"})
            if len(batch_rgb) == args.batch_size:
                tile_rows.extend(infer_batch(torch, model, transform, batch_rgb, batch_meta, device, args.precision))
                batch_rgb = []
                batch_meta = []
        if batch_rgb:
            tile_rows.extend(infer_batch(torch, model, transform, batch_rgb, batch_meta, device, args.precision))
    slide.close()
    return summarize_slide(tile_rows, slide_path), tile_rows


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

    print(f"Using device: {device}", file=sys.stderr)
    model = load_virchow2(torch, timm, args.model_id, device)
    transform = make_transform(model, args.input_size, Image, resolve_data_config, create_transform)

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
        "model_id": args.model_id,
        "embedding_recipe": "concat(class_token, mean(patch_tokens_after_4_register_tokens))",
        "input_size": args.input_size,
        "target_mpp": args.target_mpp,
        "tile_limit": args.tile_limit,
        "n_slides": len(slide_rows),
        "embedding_dimensions": len([key for key in slide_rows[-1] if key.startswith("embedding_")]) if slide_rows else 0,
        "slide_embeddings": str(slide_embeddings_path),
        "tile_embeddings": str(tile_embeddings_path) if args.save_tile_csv else "",
    }
    (out_dir / "virchow2_embedding_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Done. Wrote {slide_embeddings_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
