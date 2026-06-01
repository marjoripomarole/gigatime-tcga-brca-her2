#!/usr/bin/env python3
"""Build clinical HER2 labels for TCGA-BRCA from GDC clinical supplement data."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import requests

GDC_API = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-BRCA"
PATIENT_BIOTAB_NAME = "nationwidechildrens.org_clinical_patient_brca.txt"

MISSING_VALUES = {
    "",
    "[Not Applicable]",
    "[Not Available]",
    "[Not Evaluated]",
    "[Unknown]",
    "Not Reported",
    "Unknown",
}

OUTPUT_COLUMNS = [
    "case_submitter_id",
    "bcr_patient_uuid",
    "clinical_her2_group",
    "clinical_her2_group_rule",
    "clinical_her2_group_confidence",
    "her2_ihc_receptor_status",
    "her2_ihc_percent_category",
    "her2_ihc_score",
    "her2_ihc_other_scale",
    "her2_ihc_calculation_method",
    "her2_ish_status",
    "her2_copy_number",
    "her2_centromere17_copy_number",
    "her2_cep17_ratio",
    "her2_ish_other_scale",
    "her2_ish_calculation_method",
    "er_status",
    "pr_status",
    "source_file_id",
    "source_file_name",
]


def eq_filter(field: str, value: str | list[str]) -> dict[str, Any]:
    return {"op": "=", "content": {"field": field, "value": value}}


def and_filter(*parts: dict[str, Any]) -> dict[str, Any]:
    return {"op": "and", "content": list(parts)}


def clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def is_missing(value: str) -> bool:
    return clean(value) in MISSING_VALUES


def normalize_status(value: str) -> str:
    value = clean(value)
    if is_missing(value):
        return value
    normalized = value.lower()
    if normalized in {"positive", "negative", "equivocal", "indeterminate"}:
        return normalized.title()
    return value


def normalize_ihc_score(value: str) -> str:
    value = clean(value)
    if is_missing(value):
        return value
    value = value.replace(" ", "")
    if value in {"0", "1+", "2+", "3+"}:
        return value
    return clean(value)


def classify_her2(row: dict[str, str]) -> tuple[str, str, str]:
    ihc_status = normalize_status(row.get("lab_proc_her2_neu_immunohistochemistry_receptor_status", ""))
    ihc_score = normalize_ihc_score(row.get("her2_immunohistochemistry_level_result", ""))
    ish_status = normalize_status(row.get("lab_procedure_her2_neu_in_situ_hybrid_outcome_type", ""))

    if ihc_score == "3+":
        return "HER2-positive", "IHC score 3+", "direct"
    if ish_status == "Positive":
        return "HER2-positive", "ISH positive", "direct"

    if ihc_score == "1+":
        return "HER2-low", "IHC score 1+ with no positive ISH", "direct"
    if ihc_score == "2+" and ish_status == "Negative":
        return "HER2-low", "IHC score 2+ and ISH negative", "direct"
    if ihc_score == "0":
        return "HER2-zero", "IHC score 0 with no positive ISH", "direct"

    if ihc_status == "Positive" and is_missing(ihc_score) and is_missing(ish_status):
        return "HER2-positive", "IHC receptor status positive; detailed score/ISH missing", "inferred"

    if ihc_score == "2+":
        return "HER2-unknown", "IHC score 2+ without negative or positive ISH", "unknown"
    if ihc_status == "Negative" and is_missing(ihc_score):
        return "HER2-unknown", "IHC receptor status negative but score missing; cannot split zero versus low", "unknown"
    if ihc_status == "Equivocal":
        return "HER2-unknown", "IHC receptor status equivocal without definitive ISH", "unknown"
    if ihc_status == "Indeterminate" or ish_status == "Indeterminate":
        return "HER2-unknown", "HER2 status indeterminate", "unknown"
    if is_missing(ihc_status) and is_missing(ihc_score) and is_missing(ish_status):
        return "HER2-unknown", "HER2 IHC/ISH fields missing or not evaluated", "unknown"
    return "HER2-unknown", "HER2 fields incomplete or not classifiable by rule", "unknown"


def query_patient_biotab(project_id: str) -> dict[str, Any]:
    filters = and_filter(
        eq_filter("cases.project.project_id", project_id),
        eq_filter("data_type", "Clinical Supplement"),
        eq_filter("data_format", "BCR Biotab"),
    )
    payload = {
        "filters": filters,
        "fields": ",".join(["file_id", "file_name", "data_type", "data_format", "data_category"]),
        "format": "JSON",
        "size": 2000,
        "sort": "file_name:asc",
    }
    response = requests.post(f"{GDC_API}/files", json=payload, timeout=120)
    response.raise_for_status()
    hits = response.json()["data"]["hits"]
    for hit in hits:
        if hit.get("file_name") == PATIENT_BIOTAB_NAME:
            return hit
    names = ", ".join(hit.get("file_name", "") for hit in hits)
    raise ValueError(f"Could not find {PATIENT_BIOTAB_NAME}. Available BCR Biotabs: {names}")


def download_text_file(file_id: str, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(f"{GDC_API}/data/{file_id}", timeout=120)
    response.raise_for_status()
    text = response.text
    destination.write_text(text, encoding="utf-8")
    return text


def parse_patient_biotab(text: str) -> list[dict[str, str]]:
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 4:
        raise ValueError("Clinical patient Biotab is unexpectedly short.")
    header = lines[1].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[3:]:
        values = line.split("\t")
        row = dict(zip(header, values))
        if row.get("bcr_patient_barcode"):
            rows.append(row)
    return rows


def build_label_row(row: dict[str, str], source: dict[str, Any]) -> dict[str, str]:
    group, rule, confidence = classify_her2(row)
    return {
        "case_submitter_id": clean(row.get("bcr_patient_barcode", "")),
        "bcr_patient_uuid": clean(row.get("bcr_patient_uuid", "")),
        "clinical_her2_group": group,
        "clinical_her2_group_rule": rule,
        "clinical_her2_group_confidence": confidence,
        "her2_ihc_receptor_status": clean(row.get("lab_proc_her2_neu_immunohistochemistry_receptor_status", "")),
        "her2_ihc_percent_category": clean(row.get("her2_erbb_pos_finding_cell_percent_category", "")),
        "her2_ihc_score": clean(row.get("her2_immunohistochemistry_level_result", "")),
        "her2_ihc_other_scale": clean(row.get("pos_finding_her2_erbb2_other_measurement_scale_text", "")),
        "her2_ihc_calculation_method": clean(row.get("her2_erbb_method_calculation_method_text", "")),
        "her2_ish_status": clean(row.get("lab_procedure_her2_neu_in_situ_hybrid_outcome_type", "")),
        "her2_copy_number": clean(row.get("her2_neu_breast_carcinoma_copy_analysis_input_total_number", "")),
        "her2_centromere17_copy_number": clean(
            row.get("her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count", "")
        ),
        "her2_cep17_ratio": clean(row.get("her2_neu_chromosone_17_signal_ratio_value", "")),
        "her2_ish_other_scale": clean(row.get("her2_and_centromere_17_positive_finding_other_measurement_scale_text", "")),
        "her2_ish_calculation_method": clean(
            row.get("her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text", "")
        ),
        "er_status": clean(row.get("breast_carcinoma_estrogen_receptor_status", "")),
        "pr_status": clean(row.get("breast_carcinoma_progesterone_receptor_status", "")),
        "source_file_id": clean(source.get("file_id", "")),
        "source_file_name": clean(source.get("file_name", "")),
    }


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def read_case_set(path: Path, column: str = "case_submitter_id") -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or column not in reader.fieldnames:
            return set()
        return {row[column] for row in reader if row.get(column)}


def overlap_counts(label_rows: list[dict[str, str]], case_ids: set[str]) -> Counter[str]:
    labels_by_case = {row["case_submitter_id"]: row["clinical_her2_group"] for row in label_rows}
    return Counter(labels_by_case[case_id] for case_id in case_ids if case_id in labels_by_case)


def print_counts(title: str, counts: Counter[str]) -> None:
    ordered_groups = ["HER2-positive", "HER2-low", "HER2-zero", "HER2-unknown"]
    print(title)
    for group in ordered_groups:
        print(f"  {group}: {counts.get(group, 0)}")


def write_metadata(path: Path, source: dict[str, Any], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "project_id": PROJECT_ID,
        "gdc_api": GDC_API,
        "source_file": source,
        "n_rows": len(rows),
        "clinical_her2_group_counts": dict(Counter(row["clinical_her2_group"] for row in rows)),
        "label_rules": {
            "HER2-positive": "IHC 3+, ISH positive, or positive IHC receptor status when detailed score/ISH are missing.",
            "HER2-low": "IHC 1+ with no positive ISH, or IHC 2+ with ISH negative.",
            "HER2-zero": "IHC 0 with no positive ISH.",
            "HER2-unknown": "Missing, not evaluated, equivocal without definitive ISH, or otherwise incomplete fields.",
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", default=PROJECT_ID)
    parser.add_argument("--out", default="data/tcga_brca/clinical_her2_labels.csv")
    parser.add_argument("--raw-out", default=f"data/tcga_brca/clinical/{PATIENT_BIOTAB_NAME}")
    parser.add_argument("--metadata-out", default="data/tcga_brca/clinical_her2_labels_metadata.json")
    parser.add_argument("--slides", default="data/tcga_brca/tcga_brca_diagnostic_slides_files.csv")
    parser.add_argument("--expression", default="data/tcga_brca/erbb2_expression.csv")
    parser.add_argument("--selected-cases", default="data/tcga_brca/her2_extreme_cases.csv")
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca_extremes/slide_scores.csv")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = query_patient_biotab(args.project_id)
    text = download_text_file(source["file_id"], Path(args.raw_out))
    patient_rows = parse_patient_biotab(text)
    label_rows = [build_label_row(row, source) for row in patient_rows]
    label_rows = sorted(label_rows, key=lambda row: row["case_submitter_id"])

    write_csv(Path(args.out), label_rows, OUTPUT_COLUMNS)
    write_metadata(Path(args.metadata_out), source, label_rows)

    print(f"Wrote {len(label_rows)} clinical HER2 labels to {args.out}")
    print(f"Wrote raw clinical Biotab to {args.raw_out}")
    print(f"Wrote metadata to {args.metadata_out}")
    print_counts("Clinical HER2 group counts:", Counter(row["clinical_her2_group"] for row in label_rows))

    overlap_inputs = [
        ("Overlap with available slide metadata:", Path(args.slides)),
        ("Overlap with ERBB2 expression table:", Path(args.expression)),
        ("Overlap with selected ERBB2-extreme cases:", Path(args.selected_cases)),
        ("Overlap with processed GigaTIME slides:", Path(args.slide_scores)),
    ]
    for title, path in overlap_inputs:
        case_ids = read_case_set(path)
        if case_ids:
            print_counts(title, overlap_counts(label_rows, case_ids))
        else:
            print(f"{title}\n  skipped: no readable case_submitter_id values at {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
