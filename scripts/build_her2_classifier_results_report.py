#!/usr/bin/env python3
"""Build a classifier-only notebook and HTML report for clinical HER2 results."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path


TASK_ORDER = [
    "her2_low_vs_zero",
    "her2_positive_vs_negative",
    "her2_three_class",
]

FEATURE_LABELS = {
    "gigatime_mean_channels": "GigaTIME mean channels",
    "gigatime_mean_and_fraction_channels": "GigaTIME mean + fraction channels",
    "interpretable_marker_means": "Interpretable marker means",
    "virtual_programs": "Virtual programs",
    "erbb2_rna_reference_not_h_e": "ERBB2 RNA reference",
}

CLASSIFIER_ASSETS = "../docs/assets/clinical_her2_classifier_baseline"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metrics",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/"
            "classifier_metrics.csv"
        ),
    )
    parser.add_argument(
        "--predictions",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/"
            "classifier_crossval_predictions.csv"
        ),
    )
    parser.add_argument(
        "--feature-sets",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/"
            "classifier_feature_sets.json"
        ),
    )
    parser.add_argument("--out-notebook", default="notebooks/clinical_her2_classifier_results.ipynb")
    parser.add_argument("--out-html", default="notebooks/clinical_her2_classifier_results.html")
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    return f"{value:.{digits}f}"


def feature_label(name: str) -> str:
    return FEATURE_LABELS.get(name, name)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def html_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def sort_key(row: dict[str, str]) -> tuple[int, float]:
    task = row.get("task", "")
    task_idx = TASK_ORDER.index(task) if task in TASK_ORDER else len(TASK_ORDER)
    return task_idx, -as_float(row, "balanced_accuracy")


def best_rows(metrics: list[dict[str, str]], h_e_only: bool) -> list[dict[str, str]]:
    candidates = [row for row in metrics if row.get("model") == "regularized_logistic"]
    if h_e_only:
        candidates = [row for row in candidates if row.get("feature_set") != "erbb2_rna_reference_not_h_e"]
    else:
        candidates = [row for row in candidates if row.get("feature_set") == "erbb2_rna_reference_not_h_e"]
    rows = []
    for task in TASK_ORDER:
        task_rows = [row for row in candidates if row.get("task") == task]
        if task_rows:
            rows.append(max(task_rows, key=lambda row: as_float(row, "balanced_accuracy")))
    return rows


def best_summary_rows(rows: list[dict[str, str]], include_feature: bool = True) -> list[list[str]]:
    table_rows = []
    for row in rows:
        item = [row["task_label"]]
        if include_feature:
            item.append(feature_label(row["feature_set"]))
        item.extend(
            [
                row["n_cases"],
                fmt(as_float(row, "accuracy")),
                fmt(as_float(row, "balanced_accuracy")),
                fmt(as_float(row, "macro_auc_ovr")),
                fmt(as_float(row, "sensitivity")),
                fmt(as_float(row, "specificity")),
            ]
        )
        table_rows.append(item)
    return table_rows


def all_logistic_rows(metrics: list[dict[str, str]]) -> list[list[str]]:
    rows = [row for row in metrics if row.get("model") == "regularized_logistic"]
    rows = sorted(rows, key=sort_key)
    return [
        [
            row["task_label"],
            feature_label(row["feature_set"]),
            row["n_cases"],
            fmt(as_float(row, "accuracy")),
            fmt(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr")),
        ]
        for row in rows
    ]


def mistake_rows(predictions: list[dict[str, str]], best_h_e: list[dict[str, str]], limit: int = 18) -> list[list[str]]:
    rows = []
    for best in best_h_e:
        subset = [
            row
            for row in predictions
            if row.get("task") == best["task"]
            and row.get("feature_set") == best["feature_set"]
            and row.get("model") == best["model"]
            and row.get("correct") == "False"
        ]
        for row in subset:
            probabilities = [
                as_float(row, key)
                for key in row
                if key.startswith("prob_") and as_float(row, key) == as_float(row, key)
            ]
            confidence = max(probabilities) if probabilities else float("nan")
            rows.append(
                [
                    row["task_label"],
                    row["case_submitter_id"],
                    row["true_label"],
                    row["predicted_label"],
                    fmt(confidence),
                ]
            )
    return rows[:limit]


def feature_set_rows(feature_sets: dict[str, list[str]]) -> list[list[str]]:
    return [[feature_label(name), str(len(cols)), ", ".join(cols[:8]) + ("..." if len(cols) > 8 else "")] for name, cols in feature_sets.items()]


def notebook_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source.strip().splitlines()],
    }


def build_notebook(
    path: Path,
    best_h_e_table: list[list[str]],
    reference_table: list[list[str]],
    all_table: list[list[str]],
    mistake_table: list[list[str]],
    features_table: list[list[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        notebook_cell(
            """
# Clinical HER2 Classifier Results

Classifier-only report for the TCGA-BRCA GigaTIME HER2 pilot.

**Plain-language bottom line:** the first real classifier is working, but it is not clinically reliable yet. The strongest image-derived result is HER2-low versus HER2-zero. HER2-positive classification from GigaTIME/H&E features is weak, and full three-class prediction is still at chance.
            """
        ),
        notebook_cell(
            """
## What Was Trained

- Input: slide-level GigaTIME virtual mIF features from the 256-tile run.
- Labels: clinical HER2 groups from TCGA IHC/ISH fields.
- Models: regularized logistic classifier and nearest-centroid baseline.
- Evaluation: leave-one-out cross-validation.

Every reported prediction is held out: each slide is predicted after that slide is removed from the training set.
            """
        ),
        notebook_cell(
            "## Feature Sets\n\n"
            + markdown_table(["Feature set", "Number of features", "Example features"], features_table)
        ),
        notebook_cell(
            "## Best GigaTIME/H&E Classifier Results\n\n"
            + markdown_table(
                [
                    "Task",
                    "Best feature set",
                    "N",
                    "Accuracy",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                best_h_e_table,
            )
            + "\n\n**Interpretation:** the HER2-low versus HER2-zero classifier is the only promising image-derived result so far."
        ),
        notebook_cell(
            "## Performance Plot\n\n"
            f"![Classifier balanced accuracy]({CLASSIFIER_ASSETS}/classifier_balanced_accuracy.png)"
        ),
        notebook_cell(
            "## Confusion Matrices\n\n"
            "These matrices show the best GigaTIME/H&E regularized logistic model for each task. They do not use ERBB2 RNA.\n\n"
            "### HER2-low versus HER2-zero\n\n"
            f"![HER2-low versus HER2-zero confusion matrix]({CLASSIFIER_ASSETS}/confusion_her2_low_vs_zero.png)\n\n"
            "### HER2-positive versus HER2-negative\n\n"
            f"![HER2-positive versus HER2-negative confusion matrix]({CLASSIFIER_ASSETS}/confusion_her2_positive_vs_negative.png)\n\n"
            "### Three-class HER2 group\n\n"
            f"![Three-class HER2 confusion matrix]({CLASSIFIER_ASSETS}/confusion_her2_three_class.png)"
        ),
        notebook_cell(
            "## ERBB2 RNA Reference\n\n"
            "ERBB2 RNA is included only as a non-H&E reference. It is not an image classifier.\n\n"
            + markdown_table(
                ["Task", "N", "Accuracy", "Balanced accuracy", "Macro AUC", "Sensitivity", "Specificity"],
                reference_table,
            )
            + "\n\nThis matters because ERBB2 RNA performs much better for HER2-positive versus HER2-negative than the current GigaTIME/H&E features. So the labels contain HER2 signal, but the current image-derived model is not capturing that HER2-positive signal reliably yet."
        ),
        notebook_cell(
            "## All Regularized Logistic Results\n\n"
            + markdown_table(["Task", "Feature set", "N", "Accuracy", "Balanced accuracy", "Macro AUC"], all_table)
        ),
        notebook_cell(
            "## Example Held-Out Mistakes\n\n"
            + markdown_table(["Task", "Case", "True label", "Predicted label", "Model confidence"], mistake_table)
            + "\n\nMistakes are useful here. They tell us what to inspect next: the H&E tiles and virtual channels for cases the classifier confuses."
        ),
        notebook_cell(
            """
## Scientific Interpretation

What we gained:

- We now have a real classifier baseline, not only group-average comparisons.
- The HER2-low versus HER2-zero result matches the earlier group-average pattern.
- The current slide-level GigaTIME features are not enough for reliable HER2-positive diagnosis.
- The full three-class problem is not solved.

What this means clinically:

- This does not support clinical HER2 diagnosis yet.
- It does support a sharper research hypothesis: GigaTIME may capture a microenvironment signal that separates HER2-low from HER2-zero better than it captures HER2-positive status.
            """
        ),
        notebook_cell(
            """
## Next Classifier Step

The next model should improve the input, not just swap the classifier:

1. Restrict features to tumor-rich tiles instead of all tissue tiles.
2. Add tile distribution features: percentiles, maximum signal, and heterogeneity.
3. Add H&E tile embeddings if available.
4. Move toward multiple-instance learning after tile-level features are organized.
5. Evaluate with nested cross-validation or a separate held-out test set before making stronger claims.
            """
        ),
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path.write_text(json.dumps(notebook, indent=2) + "\n", encoding="utf-8")


def section(title: str, body: str) -> str:
    return f"<section><h2>{html.escape(title)}</h2>{body}</section>"


def build_html(
    path: Path,
    best_h_e_table: list[list[str]],
    reference_table: list[list[str]],
    all_table: list[list[str]],
    mistake_table: list[list[str]],
    features_table: list[list[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    css = """
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #1f2933; background: #f6f7f9; }
main { max-width: 1120px; margin: 0 auto; padding: 42px 28px 70px; }
.hero { background: #111827; color: white; padding: 34px; border-radius: 8px; }
.hero h1 { margin: 0 0 10px; font-size: 34px; letter-spacing: 0; }
.hero p { font-size: 18px; line-height: 1.5; max-width: 900px; }
section { background: white; margin-top: 18px; padding: 26px; border-radius: 8px; border: 1px solid #e5e7eb; }
h2 { margin: 0 0 14px; font-size: 24px; }
h3 { margin: 20px 0 10px; }
p, li { font-size: 16px; line-height: 1.55; }
.callout { background: #eef6ff; border-left: 5px solid #2563eb; padding: 14px 16px; margin: 16px 0; }
.warning { background: #fff7ed; border-left: 5px solid #f97316; padding: 14px 16px; margin: 16px 0; }
table { width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 14px; }
th, td { text-align: left; border-bottom: 1px solid #e5e7eb; padding: 8px 9px; vertical-align: top; }
th { background: #f3f4f6; font-weight: 650; }
img { display: block; width: 100%; max-width: 980px; height: auto; margin: 14px auto; border: 1px solid #e5e7eb; border-radius: 6px; background: white; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }
.small { color: #52606d; font-size: 14px; }
"""
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Clinical HER2 Classifier Results</title>",
        f"<style>{css}</style></head><body><main>",
        """
<div class="hero">
  <h1>Clinical HER2 Classifier Results</h1>
  <p>Classifier-only summary for the TCGA-BRCA GigaTIME HER2 pilot: model inputs, cross-validated metrics, confusion matrices, ERBB2 RNA reference, and next modeling steps.</p>
</div>
""",
        section(
            "Bottom Line",
            """
<div class="callout">
  <p><strong>Best image-derived signal:</strong> HER2-low versus HER2-zero reached balanced accuracy 0.800 and macro AUC 0.870 using GigaTIME mean + fraction features.</p>
</div>
<div class="warning">
  <p><strong>Important caution:</strong> HER2-positive classification from GigaTIME/H&amp;E features was weak, and full three-class HER2 prediction was at chance. This is not a clinical diagnostic model.</p>
</div>
""",
        ),
        section(
            "What Was Trained",
            """
<ul>
  <li>Input: slide-level GigaTIME virtual mIF features from the 256-tile run.</li>
  <li>Labels: clinical HER2 groups from TCGA IHC/ISH fields.</li>
  <li>Models: regularized logistic classifier and nearest-centroid baseline.</li>
  <li>Evaluation: leave-one-out cross-validation.</li>
</ul>
<p>Every reported prediction is held out: each slide is predicted after that slide is removed from the training set.</p>
""",
        ),
        section("Feature Sets", html_table(["Feature set", "Number of features", "Example features"], features_table)),
        section(
            "Best GigaTIME/H&E Results",
            html_table(
                [
                    "Task",
                    "Best feature set",
                    "N",
                    "Accuracy",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                best_h_e_table,
            )
            + "<p class='small'>The HER2-low versus HER2-zero classifier is the only promising image-derived result so far.</p>",
        ),
        section(
            "Performance Plot",
            f"<img src='{CLASSIFIER_ASSETS}/classifier_balanced_accuracy.png' alt='Classifier balanced accuracy'>",
        ),
        section(
            "Confusion Matrices",
            f"""
<p>These matrices show the best GigaTIME/H&amp;E regularized logistic model for each task. They do not use ERBB2 RNA.</p>
<div class="grid">
  <div><h3>HER2-low vs HER2-zero</h3><img src="{CLASSIFIER_ASSETS}/confusion_her2_low_vs_zero.png" alt="HER2-low versus HER2-zero confusion matrix"></div>
  <div><h3>HER2-positive vs HER2-negative</h3><img src="{CLASSIFIER_ASSETS}/confusion_her2_positive_vs_negative.png" alt="HER2-positive versus HER2-negative confusion matrix"></div>
  <div><h3>Three-class HER2</h3><img src="{CLASSIFIER_ASSETS}/confusion_her2_three_class.png" alt="Three-class HER2 confusion matrix"></div>
</div>
""",
        ),
        section(
            "ERBB2 RNA Reference",
            "<p>ERBB2 RNA is included only as a non-H&amp;E reference. It is not an image classifier.</p>"
            + html_table(
                ["Task", "N", "Accuracy", "Balanced accuracy", "Macro AUC", "Sensitivity", "Specificity"],
                reference_table,
            )
            + "<p class='small'>ERBB2 RNA performs much better for HER2-positive versus HER2-negative than the current GigaTIME/H&amp;E features, showing that the labels contain HER2 signal that the current image-derived model is not capturing reliably yet.</p>",
        ),
        section(
            "All Regularized Logistic Results",
            html_table(["Task", "Feature set", "N", "Accuracy", "Balanced accuracy", "Macro AUC"], all_table),
        ),
        section(
            "Example Held-Out Mistakes",
            html_table(["Task", "Case", "True label", "Predicted label", "Model confidence"], mistake_table)
            + "<p class='small'>Mistakes are useful here. They tell us what to inspect next: the H&amp;E tiles and virtual channels for cases the classifier confuses.</p>",
        ),
        section(
            "Scientific Interpretation",
            """
<ul>
  <li>We now have a real classifier baseline, not only group-average comparisons.</li>
  <li>The HER2-low versus HER2-zero result matches the earlier group-average pattern.</li>
  <li>The current slide-level GigaTIME features are not enough for reliable HER2-positive diagnosis.</li>
  <li>The full three-class problem is not solved.</li>
  <li>The result supports a sharper research hypothesis: GigaTIME may capture a microenvironment signal that separates HER2-low from HER2-zero better than it captures HER2-positive status.</li>
</ul>
""",
        ),
        section(
            "Next Classifier Step",
            """
<ol>
  <li>Restrict features to tumor-rich tiles instead of all tissue tiles.</li>
  <li>Add tile distribution features: percentiles, maximum signal, and heterogeneity.</li>
  <li>Add H&amp;E tile embeddings if available.</li>
  <li>Move toward multiple-instance learning after tile-level features are organized.</li>
  <li>Evaluate with nested cross-validation or a separate held-out test set before making stronger claims.</li>
</ol>
""",
        ),
        "</main></body></html>",
    ]
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    args = parse_args()
    metrics = read_rows(Path(args.metrics))
    predictions = read_rows(Path(args.predictions))
    feature_sets = json.loads(Path(args.feature_sets).read_text(encoding="utf-8"))

    best_h_e = best_rows(metrics, h_e_only=True)
    reference = best_rows(metrics, h_e_only=False)
    best_h_e_table = best_summary_rows(best_h_e, include_feature=True)
    reference_table = best_summary_rows(reference, include_feature=False)
    all_table = all_logistic_rows(metrics)
    mistake_table = mistake_rows(predictions, best_h_e)
    features_table = feature_set_rows(feature_sets)

    build_notebook(Path(args.out_notebook), best_h_e_table, reference_table, all_table, mistake_table, features_table)
    build_html(Path(args.out_html), best_h_e_table, reference_table, all_table, mistake_table, features_table)
    print(f"Wrote {args.out_notebook}")
    print(f"Wrote {args.out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
