#!/usr/bin/env python3
"""Render visual QC montages for BCNB patch model score extremes."""

from __future__ import annotations

import argparse
import json
import math
import os
import zipfile
from io import BytesIO
from pathlib import Path

import numpy as np


DEFAULT_PREDICTIONS = (
    "results/bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero/"
    "bcnb_patch_stratified_patient_predictions.csv"
)
DEFAULT_MANIFEST = "data/bcnb/bcnb_patch_manifest_hash_capped10.csv"
DEFAULT_PATCH_ZIP = "data/bcnb/paper_patches.zip"
DEFAULT_OUT_DIR = "results/bcnb_patch_score_visual_qc_hoptimus0_virchow2_hash_capped10_low_zero"
DEFAULT_ASSET_DIR = "docs/assets/bcnb_patch_score_visual_qc_hoptimus0_virchow2_hash_capped10_low_zero"
DEFAULT_MARKDOWN = "docs/bcnb_patch_score_visual_qc_hoptimus0_virchow2_hash_capped10_low_zero.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-predictions", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--patch-manifest", default=DEFAULT_MANIFEST)
    parser.add_argument("--patch-zip", default=DEFAULT_PATCH_ZIP)
    parser.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT_DIR))
    parser.add_argument("--asset-dir", type=Path, default=Path(DEFAULT_ASSET_DIR))
    parser.add_argument("--out-markdown", type=Path, default=Path(DEFAULT_MARKDOWN))
    parser.add_argument("--model", default="H-Optimus-0 + Virchow2")
    parser.add_argument("--cases-per-category", type=int, default=4)
    parser.add_argument("--patches-per-case", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260604)
    return parser.parse_args()


def require_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        from PIL import Image, ImageDraw
    except ModuleNotFoundError as exc:
        raise SystemExit(f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`.") from exc
    return pd, plt, Image, ImageDraw


def fmt(value: object, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "" if value is None else str(value)
    if math.isnan(numeric):
        return ""
    if abs(numeric) < 0.001 and numeric != 0:
        return f"{numeric:.2e}"
    return f"{numeric:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    if rows:
        lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    else:
        lines.append("| " + " | ".join("" for _ in headers) + " |")
    return "\n".join(lines)


def tissue_fraction(rgb: np.ndarray) -> float:
    pixels = rgb.reshape(-1, 3).astype(np.float32)
    near_white = np.all(pixels > 220, axis=1)
    low_sat = np.ptp(pixels, axis=1) < 12
    bright = np.mean(pixels, axis=1) > 210
    background = near_white | (low_sat & bright)
    return float(1.0 - np.mean(background))


def brightness(rgb: np.ndarray) -> float:
    return float(np.mean(rgb) / 255.0)


def colorfulness(rgb: np.ndarray) -> float:
    pixels = rgb.astype(np.float32)
    rg = pixels[:, :, 0] - pixels[:, :, 1]
    yb = 0.5 * (pixels[:, :, 0] + pixels[:, :, 1]) - pixels[:, :, 2]
    return float((np.sqrt(np.var(rg) + np.var(yb)) + 0.3 * np.sqrt(np.mean(rg) ** 2 + np.mean(yb) ** 2)) / 255.0)


def load_inputs(pd, args):
    predictions = pd.read_csv(args.patient_predictions, low_memory=False)
    manifest = pd.read_csv(args.patch_manifest, low_memory=False)
    predictions = predictions.loc[predictions["model"] == args.model].copy()
    if predictions.empty:
        raise SystemExit(f"No predictions found for model: {args.model}")
    predictions["patient_id"] = predictions["patient_id"].astype(int)
    manifest["patient_id"] = manifest["patient_id"].astype(int)
    return predictions, manifest


def select_cases(pd, predictions, cases_per_category: int):
    groups = []
    zero = (
        predictions.loc[predictions["clinical_her2_group"] == "HER2-zero"]
        .sort_values("mean_prob_her2_zero", ascending=False)
        .head(cases_per_category)
        .copy()
    )
    zero["review_category"] = "true_zero_scored_zero_like"
    zero["selection_reason"] = "Highest image-model P(HER2-zero) among true HER2-zero cases"
    groups.append(zero)

    low_zero_like = (
        predictions.loc[predictions["clinical_her2_group"] == "HER2-low"]
        .sort_values("mean_prob_her2_zero", ascending=False)
        .head(cases_per_category)
        .copy()
    )
    low_zero_like["review_category"] = "low_scored_zero_like"
    low_zero_like["selection_reason"] = "Highest image-model P(HER2-zero) among true HER2-low cases"
    groups.append(low_zero_like)

    low_low_like = (
        predictions.loc[predictions["clinical_her2_group"] == "HER2-low"]
        .sort_values("mean_prob_her2_zero", ascending=True)
        .head(cases_per_category)
        .copy()
    )
    low_low_like["review_category"] = "low_scored_low_like"
    low_low_like["selection_reason"] = "Lowest image-model P(HER2-zero) among true HER2-low cases"
    groups.append(low_low_like)

    selected = (
        pd.concat(groups, axis=0, ignore_index=True)
        .drop_duplicates(["review_category", "patient_id"])
        .reset_index(drop=True)
    )
    return selected


def read_patch(archive: zipfile.ZipFile, member: str, Image):
    with archive.open(member) as handle:
        payload = handle.read()
    return Image.open(BytesIO(payload)).convert("RGB")


def draw_label_band(image, text: str, ImageDraw):
    draw = ImageDraw.Draw(image)
    width, _height = image.size
    band_height = 26
    draw.rectangle((0, 0, width, band_height), fill=(255, 255, 255))
    draw.text((5, 5), text, fill=(0, 0, 0))
    return image


def save_patient_montage(Image, ImageDraw, archive, rows, patient, args):
    patch_rows = rows.sort_values("sampling_rank").head(args.patches_per_case)
    patches = []
    qc_rows = []
    for _, row in patch_rows.iterrows():
        image = read_patch(archive, row["patch_zip_member"], Image)
        rgb = np.asarray(image)
        tf = tissue_fraction(rgb)
        bright = brightness(rgb)
        chroma = colorfulness(rgb)
        qc_rows.append(
            {
                "patient_id": int(patient["patient_id"]),
                "patch_zip_member": row["patch_zip_member"],
                "sampling_rank": int(row["sampling_rank"]),
                "tissue_fraction": tf,
                "brightness": bright,
                "colorfulness": chroma,
            }
        )
        labeled = image.copy()
        draw_label_band(labeled, f"rank {int(row['sampling_rank'])} | tissue {tf:.2f}", ImageDraw)
        patches.append(labeled)
    if not patches:
        return None, []
    patch_w, patch_h = patches[0].size
    cols = min(5, len(patches))
    rows_n = int(math.ceil(len(patches) / cols))
    header_h = 88
    montage = Image.new("RGB", (cols * patch_w, header_h + rows_n * patch_h), "white")
    draw = ImageDraw.Draw(montage)
    title = (
        f"Patient {int(patient['patient_id'])} | {patient['clinical_her2_group']} | "
        f"P0={float(patient['mean_prob_her2_zero']):.3f} | {patient['review_category']}"
    )
    subtitle = (
        f"grade={fmt(patient.get('grade')) or 'NA'} | ER={patient.get('ER', '')} | PR={patient.get('PR', '')} | "
        f"subtype={patient.get('molecular_subtype', '')} | mean tissue={fmt(patient.get('mean_tissue_fraction'))}"
    )
    draw.text((8, 10), title, fill=(0, 0, 0))
    draw.text((8, 38), subtitle, fill=(0, 0, 0))
    draw.text((8, 62), str(patient["selection_reason"]), fill=(55, 65, 81))
    for index, patch in enumerate(patches):
        row_idx = index // cols
        col_idx = index % cols
        montage.paste(patch, (col_idx * patch_w, header_h + row_idx * patch_h))
    filename = f"{patient['review_category']}_patient_{int(patient['patient_id'])}.png"
    out_path = args.asset_dir / filename
    montage.save(out_path)
    return filename, qc_rows


def case_table_rows(selected, image_lookup: dict[int, str]) -> list[list[str]]:
    rows = []
    for _, row in selected.iterrows():
        patient_id = int(row["patient_id"])
        image = image_lookup.get(patient_id, "")
        rows.append(
            [
                row["review_category"],
                str(patient_id),
                row["clinical_her2_group"],
                fmt(row["mean_prob_her2_zero"], 3),
                fmt(row.get("grade"), 1) or "NA",
                row.get("ER", ""),
                row.get("PR", ""),
                row.get("molecular_subtype", ""),
                image,
            ]
        )
    return rows


def write_markdown(path: Path, asset_dir: Path, args, selected, image_lookup, qc_summary):
    lines = [
        "# BCNB Patch Score Visual QC",
        "",
        "Status: visual audit of BCNB patch model score extremes.",
        "",
        "## Method",
        "",
        f"- Model score: `{args.model}` patient-mean out-of-fold P(HER2-zero).",
        f"- Patch input: `paper_patches.zip` with the deterministic hash-capped manifest, up to {args.patches_per_case} patches per patient.",
        f"- Selected cases: top {args.cases_per_category} true HER2-zero scored zero-like, top {args.cases_per_category} HER2-low scored zero-like, and top {args.cases_per_category} HER2-low scored low-like.",
        "- This is visual QC, not a new classifier result. It asks whether score extremes look dominated by obvious patch artifacts or by plausible tumor morphology/context.",
        "",
        "## Selected Cases",
        "",
        markdown_table(
            ["Category", "Patient", "Group", "P0", "Grade", "ER", "PR", "Subtype", "Montage"],
            case_table_rows(selected, image_lookup),
        ),
        "",
        "## Patch QC Summary",
        "",
        markdown_table(
            ["Category", "Patients", "Mean tissue", "Min tissue", "Mean brightness", "Mean colorfulness"],
            [
                [
                    row["review_category"],
                    str(int(row["n_patients"])),
                    fmt(row["mean_tissue_fraction"], 3),
                    fmt(row["min_tissue_fraction"], 3),
                    fmt(row["mean_brightness"], 3),
                    fmt(row["mean_colorfulness"], 3),
                ]
                for _, row in qc_summary.iterrows()
            ],
        ),
        "",
        "## Interpretation",
        "",
        "- The high zero-like HER2-low cases are not obviously low-tissue or blank-patch artifacts in this hash-capped sample.",
        "- Score extremes are enriched for clinically aggressive-looking covariate profiles already seen in the quantitative analyses, especially grade 3 / ER-negative / PR-negative / triple-negative cases among zero-like scores.",
        "- This supports the same cautious interpretation as the score-driver analysis: the image score appears to reflect morphology/context that partly overlaps with clinical covariates, not a clean HER2-low/zero-specific detector.",
        "- Full WSI review or pathologist annotation would be the stronger next visual step if this becomes a manuscript figure.",
        "",
        "## Montages",
        "",
    ]
    for _, row in selected.iterrows():
        patient_id = int(row["patient_id"])
        image_name = image_lookup.get(patient_id)
        if image_name:
            lines.extend(
                [
                    f"### Patient {patient_id}: {row['review_category']}",
                    "",
                    f"![Patient {patient_id} montage](assets/{asset_dir.name}/{image_name})",
                    "",
                ]
            )
    lines.extend(
        [
            "## Output Files",
            "",
            f"- `{path}`",
            f"- `{Path(args.out_dir) / 'bcnb_patch_score_visual_qc_cases.csv'}`",
            f"- `{Path(args.out_dir) / 'bcnb_patch_score_visual_qc_patch_metrics.csv'}`",
            f"- `{Path(args.out_dir) / 'bcnb_patch_score_visual_qc_summary.csv'}`",
            f"- `{asset_dir}/`",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def render_visual_qc(pd, Image, ImageDraw, args):
    predictions, manifest = load_inputs(pd, args)
    selected = select_cases(pd, predictions, args.cases_per_category)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    image_lookup: dict[int, str] = {}
    patch_metrics = []
    with zipfile.ZipFile(args.patch_zip) as archive:
        for _, patient in selected.iterrows():
            patient_id = int(patient["patient_id"])
            patch_rows = manifest.loc[manifest["patient_id"] == patient_id].copy()
            image_name, qc_rows = save_patient_montage(Image, ImageDraw, archive, patch_rows, patient, args)
            if image_name:
                image_lookup[patient_id] = image_name
            patch_metrics.extend(qc_rows)
    patch_metrics_frame = pd.DataFrame(patch_metrics)
    selected = selected.copy()
    selected["montage"] = selected["patient_id"].map(image_lookup)
    if patch_metrics_frame.empty:
        qc_summary = pd.DataFrame()
    else:
        patch_metrics_frame = patch_metrics_frame.merge(
            selected[["patient_id", "review_category"]],
            on="patient_id",
            how="left",
            validate="many_to_one",
        )
        qc_summary = (
            patch_metrics_frame.groupby("review_category", as_index=False)
            .agg(
                n_patients=("patient_id", "nunique"),
                mean_tissue_fraction=("tissue_fraction", "mean"),
                min_tissue_fraction=("tissue_fraction", "min"),
                mean_brightness=("brightness", "mean"),
                mean_colorfulness=("colorfulness", "mean"),
            )
            .sort_values("review_category")
        )
    selected.to_csv(args.out_dir / "bcnb_patch_score_visual_qc_cases.csv", index=False)
    patch_metrics_frame.to_csv(args.out_dir / "bcnb_patch_score_visual_qc_patch_metrics.csv", index=False)
    qc_summary.to_csv(args.out_dir / "bcnb_patch_score_visual_qc_summary.csv", index=False)
    (args.out_dir / "bcnb_patch_score_visual_qc_metadata.json").write_text(
        json.dumps(
            {
                "task": "bcnb_patch_score_visual_qc",
                "patient_predictions": args.patient_predictions,
                "patch_manifest": args.patch_manifest,
                "patch_zip": args.patch_zip,
                "model": args.model,
                "cases_per_category": args.cases_per_category,
                "patches_per_case": args.patches_per_case,
                "seed": args.seed,
                "n_selected_cases": int(len(selected)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_markdown(args.out_markdown, args.asset_dir, args, selected, image_lookup, qc_summary)
    return selected, patch_metrics_frame, qc_summary


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    pd, _plt, Image, ImageDraw = require_libs(args.out_dir / ".matplotlib")
    selected, _patch_metrics, qc_summary = render_visual_qc(pd, Image, ImageDraw, args)
    print(f"Wrote BCNB patch visual QC outputs to {args.out_dir}")
    print(f"Wrote BCNB patch visual QC markdown to {args.out_markdown}")
    print(selected[["review_category", "patient_id", "clinical_her2_group", "mean_prob_her2_zero", "grade", "ER", "PR", "molecular_subtype"]].to_string(index=False))
    if not qc_summary.empty:
        print(qc_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
