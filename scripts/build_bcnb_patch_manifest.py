#!/usr/bin/env python3
"""Build a patient-linked manifest for BCNB precomputed paper patches."""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
from collections import Counter
from pathlib import Path


DEFAULT_LABELS = Path("data/bcnb/bcnb_her2_labels.csv")
DEFAULT_PATCH_ZIP = Path("data/bcnb/paper_patches.zip")
DEFAULT_OUTPUT = Path("data/bcnb/bcnb_patch_manifest.csv")

OUTPUT_COLUMNS = [
    "patient_id",
    "clinical_her2_group",
    "her2_status",
    "her2_ihc",
    "grade",
    "ER",
    "PR",
    "ki67",
    "molecular_subtype",
    "aln_status",
    "patch_zip_member",
    "patch_filename",
    "patch_file_size",
    "filename_num_1",
    "filename_num_2",
    "filename_num_3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS, help="BCNB label CSV.")
    parser.add_argument("--patch-zip", type=Path, default=DEFAULT_PATCH_ZIP, help="BCNB paper_patches.zip path.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output patch manifest CSV.")
    parser.add_argument(
        "--max-patches-per-patient",
        type=int,
        default=0,
        help="Optional deterministic cap per patient. Use 0 to include every patch.",
    )
    return parser.parse_args()


def load_labels(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"BCNB label CSV not found: {path}")
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "patient_id" not in reader.fieldnames:
            raise ValueError(f"{path} must contain a patient_id column")
        return {row["patient_id"]: row for row in reader if row.get("patient_id")}


def numeric_filename_tokens(filename: str) -> tuple[str, str, str]:
    tokens = re.findall(r"\d+", Path(filename).stem)
    padded = (tokens + ["", "", ""])[:4]
    # Token 0 is the repeated patient id in BCNB patch filenames.
    return padded[1], padded[2], padded[3]


def iter_manifest_rows(labels: dict[str, dict[str, str]], patch_zip: Path, max_patches_per_patient: int):
    if not patch_zip.exists():
        raise FileNotFoundError(f"BCNB patch zip not found: {patch_zip}")
    per_patient_counts: Counter[str] = Counter()
    with zipfile.ZipFile(patch_zip) as archive:
        for info in sorted(archive.infolist(), key=lambda item: item.filename):
            if info.is_dir() or not info.filename.lower().endswith(".jpg"):
                continue
            parts = info.filename.split("/")
            if len(parts) < 3 or parts[0] != "patches" or not parts[1].isdigit():
                continue
            patient_id = str(int(parts[1]))
            if patient_id not in labels:
                raise ValueError(f"Patch member {info.filename} references unknown patient_id={patient_id}")
            if max_patches_per_patient and per_patient_counts[patient_id] >= max_patches_per_patient:
                continue
            per_patient_counts[patient_id] += 1

            label = labels[patient_id]
            filename = Path(info.filename).name
            filename_num_1, filename_num_2, filename_num_3 = numeric_filename_tokens(filename)
            yield {
                "patient_id": patient_id,
                "clinical_her2_group": label.get("clinical_her2_group", ""),
                "her2_status": label.get("her2_status", ""),
                "her2_ihc": label.get("her2_ihc", ""),
                "grade": label.get("grade", ""),
                "ER": label.get("ER", ""),
                "PR": label.get("PR", ""),
                "ki67": label.get("ki67", ""),
                "molecular_subtype": label.get("molecular_subtype", ""),
                "aln_status": label.get("aln_status", ""),
                "patch_zip_member": info.filename,
                "patch_filename": filename,
                "patch_file_size": str(info.file_size),
                "filename_num_1": filename_num_1,
                "filename_num_2": filename_num_2,
                "filename_num_3": filename_num_3,
            }


def write_manifest(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]], labels: dict[str, dict[str, str]], max_patches_per_patient: int) -> str:
    patient_counts = Counter(row["patient_id"] for row in rows)
    group_patch_counts = Counter(row["clinical_her2_group"] for row in rows)
    group_patient_counts = Counter(labels[patient_id]["clinical_her2_group"] for patient_id in patient_counts)
    lines = [
        f"Wrote {len(rows)} patch rows for {len(patient_counts)} patients.",
        f"Patch cap per patient: {max_patches_per_patient or 'none'}",
        "Patients by group: "
        + ", ".join(
            f"{group}={group_patient_counts[group]}" for group in ["HER2-zero", "HER2-low", "HER2-positive"]
        ),
        "Patches by group: "
        + ", ".join(f"{group}={group_patch_counts[group]}" for group in ["HER2-zero", "HER2-low", "HER2-positive"]),
    ]
    missing_patients = sorted(set(labels) - set(patient_counts), key=int)
    if missing_patients:
        lines.append(f"Patients without patches in manifest: {len(missing_patients)}")
        lines.append("First missing patient IDs: " + ", ".join(missing_patients[:10]))
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    labels = load_labels(args.labels)
    rows = list(iter_manifest_rows(labels, args.patch_zip, args.max_patches_per_patient))
    write_manifest(rows, args.output)
    print(summarize(rows, labels, args.max_patches_per_patient))


if __name__ == "__main__":
    main()
