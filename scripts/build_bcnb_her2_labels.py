#!/usr/bin/env python3
"""Build derived BCNB clinical HER2 labels from the gated clinical workbook."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_INPUT = Path("data/bcnb/patient-clinical-data.xlsx")
DEFAULT_OUTPUT = Path("data/bcnb/bcnb_her2_labels.csv")

REQUIRED_COLUMNS = [
    "Patient ID",
    "Age(years)",
    "Tumour Size(cm)",
    "Tumour Type",
    "ER",
    "PR",
    "HER2",
    "HER2 Expression",
    "Histological grading",
    "Ki67",
    "Molecular subtype",
    "ALN status",
]

OUTPUT_COLUMNS = [
    "patient_id",
    "age",
    "tumor_size",
    "tumor_type",
    "ER",
    "PR",
    "her2_status",
    "her2_ihc",
    "clinical_her2_group",
    "grade",
    "ki67",
    "molecular_subtype",
    "aln_status",
]


def clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_ihc(value: Any) -> str:
    value_str = clean(value).replace(" ", "")
    if value_str in {"0.0", "0"}:
        return "0"
    if value_str in {"1+", "2+", "3+"}:
        return value_str
    return value_str


def classify_her2(ihc_value: Any, her2_status_value: Any) -> str:
    ihc = normalize_ihc(ihc_value)
    her2_status = clean(her2_status_value).title()

    if ihc == "0":
        return "HER2-zero"
    if ihc == "1+":
        return "HER2-low"
    if ihc == "2+" and her2_status == "Negative":
        return "HER2-low"
    if ihc == "2+" and her2_status == "Positive":
        return "HER2-positive"
    if ihc == "3+":
        return "HER2-positive"
    return "HER2-unknown"


def require_columns(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"BCNB clinical workbook is missing required columns: {joined}")


def build_labels(clinical: pd.DataFrame) -> pd.DataFrame:
    require_columns(clinical)
    output = pd.DataFrame(
        {
            "patient_id": clinical["Patient ID"],
            "age": clinical["Age(years)"],
            "tumor_size": clinical["Tumour Size(cm)"],
            "tumor_type": clinical["Tumour Type"],
            "ER": clinical["ER"],
            "PR": clinical["PR"],
            "her2_status": clinical["HER2"],
            "her2_ihc": clinical["HER2 Expression"].map(normalize_ihc),
            "clinical_her2_group": [
                classify_her2(ihc, status)
                for ihc, status in zip(clinical["HER2 Expression"], clinical["HER2"], strict=True)
            ],
            "grade": clinical["Histological grading"],
            "ki67": clinical["Ki67"],
            "molecular_subtype": clinical["Molecular subtype"],
            "aln_status": clinical["ALN status"],
        },
        columns=OUTPUT_COLUMNS,
    )
    unknown = output[output["clinical_her2_group"] == "HER2-unknown"]
    if not unknown.empty:
        examples = unknown[["patient_id", "her2_status", "her2_ihc"]].head(10).to_dict(orient="records")
        raise ValueError(f"Could not classify {len(unknown)} BCNB rows by HER2 rule. Examples: {examples}")
    return output


def summarize(labels: pd.DataFrame) -> str:
    ihc_counts = Counter(labels["her2_ihc"])
    group_counts = Counter(labels["clinical_her2_group"])

    two_plus = labels[labels["her2_ihc"] == "2+"]
    two_plus_status_counts = Counter(two_plus["her2_status"])

    low_zero = labels[labels["clinical_her2_group"].isin(["HER2-zero", "HER2-low"])].copy()
    low_zero_graded = low_zero[low_zero["grade"].notna()]
    grade3 = low_zero_graded[low_zero_graded["grade"] == 3.0]
    grade3_counts = Counter(grade3["clinical_her2_group"])
    graded_counts = Counter(low_zero_graded["clinical_her2_group"])

    lines = [
        f"Wrote {len(labels)} BCNB label rows.",
        "HER2 IHC counts: "
        + ", ".join(f"{key}={ihc_counts[key]}" for key in ["0", "1+", "2+", "3+"] if key in ihc_counts),
        "Derived group counts: "
        + ", ".join(
            f"{key}={group_counts[key]}" for key in ["HER2-zero", "HER2-low", "HER2-positive"] if key in group_counts
        ),
        "2+ split by binary HER2 status: "
        + ", ".join(f"{key}={two_plus_status_counts[key]}" for key in sorted(two_plus_status_counts)),
    ]
    for group in ["HER2-zero", "HER2-low"]:
        if graded_counts[group]:
            pct = 100 * grade3_counts[group] / graded_counts[group]
            lines.append(f"{group} grade 3 among graded: {grade3_counts[group]}/{graded_counts[group]} ({pct:.1f}%)")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Path to BCNB patient-clinical-data.xlsx.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output CSV path for derived labels.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    clinical = pd.read_excel(args.input)
    labels = build_labels(clinical)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    labels.to_csv(args.output, index=False)
    print(summarize(labels))


if __name__ == "__main__":
    main()
