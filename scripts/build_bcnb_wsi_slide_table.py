#!/usr/bin/env python3
"""Build a patient-linked BCNB full-WSI slide table for GigaTIME runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from audit_bcnb_image_inputs import infer_patient_id, load_labels


DEFAULT_LABELS = Path("data/bcnb/bcnb_her2_labels.csv")
DEFAULT_WSI_DIR = Path("data/bcnb/WSIs")
DEFAULT_OUTPUT = Path("data/bcnb/bcnb_wsi_slide_table.csv")
WSI_SUFFIXES = {".svs", ".tif", ".tiff", ".ndpi", ".mrxs", ".jpg", ".jpeg", ".png"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--wsi-dir", type=Path, default=DEFAULT_WSI_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--groups",
        default="HER2-zero,HER2-low,HER2-positive",
        help="Comma-separated clinical HER2 groups to include.",
    )
    parser.add_argument("--require-all-patients", action="store_true", help="Fail unless every selected patient has a mapped WSI.")
    return parser.parse_args()


def parse_groups(raw: str) -> set[str]:
    return {item.strip() for item in raw.split(",") if item.strip()}


def read_label_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No columns found in {path}")
        if "patient_id" not in reader.fieldnames or "clinical_her2_group" not in reader.fieldnames:
            raise ValueError(f"{path} must contain patient_id and clinical_her2_group")
        return {row["patient_id"]: dict(row) for row in reader if row.get("patient_id")}


def find_wsi_files(path: Path) -> list[Path]:
    if not path.exists():
        raise FileNotFoundError(f"BCNB WSI directory is not present: {path}")
    return sorted(file for file in path.rglob("*") if file.is_file() and file.suffix.lower() in WSI_SUFFIXES)


def main() -> int:
    args = parse_args()
    label_groups = load_labels(args.labels)
    label_rows = read_label_rows(args.labels)
    groups = parse_groups(args.groups)
    selected_ids = {pid for pid, group in label_groups.items() if group in groups}
    known_ids = set(label_groups)

    rows: list[dict[str, str]] = []
    unmatched: list[str] = []
    ambiguous: list[str] = []
    seen_patient_ids: set[str] = set()
    for slide_path in find_wsi_files(args.wsi_dir):
        relative = str(slide_path.relative_to(args.wsi_dir))
        match = infer_patient_id(relative, known_ids)
        if match is None:
            unmatched.append(relative)
            continue
        if match.confidence == "ambiguous" or not match.patient_id:
            ambiguous.append(relative)
            continue
        patient_id = match.patient_id
        if patient_id not in selected_ids:
            continue
        label_row = label_rows[patient_id]
        row = dict(label_row)
        row.update(
            {
                "slide_id": slide_path.stem,
                "slide_local_path": str(slide_path),
                "slide_filename": slide_path.name,
                "slide_suffix": slide_path.suffix.lower(),
                "slide_file_size": str(slide_path.stat().st_size),
                "patient_match_confidence": match.confidence,
                "patient_match_reason": match.reason,
            }
        )
        rows.append(row)
        seen_patient_ids.add(patient_id)

    missing_patient_ids = sorted(selected_ids - seen_patient_ids, key=lambda value: int(value))
    if args.require_all_patients and missing_patient_ids:
        examples = ", ".join(missing_patient_ids[:20])
        raise FileNotFoundError(f"Missing WSIs for {len(missing_patient_ids)} selected patients. Examples: {examples}")

    rows.sort(key=lambda row: (row["clinical_her2_group"], int(row["patient_id"]), row["slide_local_path"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    group_counts: dict[str, int] = {}
    for row in rows:
        group_counts[row["clinical_her2_group"]] = group_counts.get(row["clinical_her2_group"], 0) + 1
    print(f"Wrote {len(rows)} BCNB WSI rows to {args.output}")
    print("Mapped groups: " + ", ".join(f"{group}={group_counts[group]}" for group in sorted(group_counts)))
    print(f"Missing selected patients: {len(missing_patient_ids)}")
    print(f"Unmatched image files: {len(unmatched)}")
    print(f"Ambiguous image files: {len(ambiguous)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
