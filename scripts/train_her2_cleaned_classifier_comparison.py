#!/usr/bin/env python3
"""Compare HER2 classifiers across cleaned GigaTIME feature views."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from train_her2_classifier_baseline import (
    TASKS,
    leave_one_out_predictions,
    metrics_for_predictions,
    require_libs,
)


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
INTERPRETABLE_CHANNELS = ["CD3", "CD8", "CD4", "CD20", "CD68", "CD11c", "PD-1", "PD-L1", "CK", "Ki67"]
FILTER_ORDER = ["all_sampled_tissue", "qc_cellular_tissue", "ck_enriched_top50", "ck_enriched_top25"]
FILTER_LABELS = {
    "all_sampled_tissue": "All sampled tissue",
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_enriched_top50": "CK-enriched top 50%",
    "ck_enriched_top25": "CK-enriched top 25%",
}
FEATURE_LABELS = {
    "gigatime_mean_channels": "Mean channels",
    "gigatime_mean_and_fraction_channels": "Mean + fraction channels",
    "interpretable_marker_means": "Interpretable means",
    "interpretable_distribution_features": "Interpretable distribution features",
    "virtual_programs": "Virtual programs",
    "erbb2_rna_reference_not_h_e": "ERBB2 RNA reference",
}
TASK_ORDER = ["her2_low_vs_zero", "her2_positive_vs_negative", "her2_three_class"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cleaned-features",
        default="results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleaned_slide_features.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison",
    )
    parser.add_argument("--asset-dir", default="docs/assets/clinical_her2_cleaned_classifier")
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def existing_columns(rows, columns: list[str]) -> list[str]:
    return [col for col in columns if col in rows.columns and rows[col].dtype.kind in "biufc"]


def feature_sets_for(rows) -> dict[str, list[str]]:
    mean_cols = existing_columns(rows, [f"mean_{channel}" for channel in GIGATIME_CHANNELS])
    frac_cols = existing_columns(rows, [f"frac_{channel}" for channel in GIGATIME_CHANNELS])
    interpretable_means = existing_columns(rows, [f"mean_{channel}" for channel in INTERPRETABLE_CHANNELS])
    distribution_cols = []
    for channel in INTERPRETABLE_CHANNELS:
        distribution_cols.extend(
            [
                f"mean_{channel}",
                f"median_{channel}",
                f"p90_{channel}",
                f"top10_mean_{channel}",
                f"std_{channel}",
            ]
        )
    distribution_cols = existing_columns(rows, distribution_cols)
    program_cols = existing_columns(
        rows,
        [
            "virtual_myeloid_checkpoint_score",
            "virtual_t_cell_score",
            "mean_CK",
            "mean_Ki67",
            "p90_CK",
            "p90_Ki67",
            "top10_mean_CK",
            "top10_mean_Ki67",
        ],
    )
    feature_sets = {
        "gigatime_mean_channels": mean_cols,
        "gigatime_mean_and_fraction_channels": mean_cols + frac_cols,
        "interpretable_marker_means": interpretable_means,
        "interpretable_distribution_features": distribution_cols,
        "virtual_programs": program_cols,
    }
    if "erbb2_tpm" in rows.columns:
        feature_sets["erbb2_rna_reference_not_h_e"] = ["erbb2_tpm"]
    return {name: cols for name, cols in feature_sets.items() if cols}


def run_view(pd, optimize, stats, rows, feature_view: str, l2_penalty: float):
    view_rows = rows.loc[rows["feature_view"] == feature_view].copy()
    feature_sets = feature_sets_for(view_rows)
    prediction_tables = []
    for task in TASKS:
        for feature_name, feature_cols in feature_sets.items():
            predictions = leave_one_out_predictions(pd, optimize, view_rows, task, feature_name, feature_cols, l2_penalty)
            predictions["feature_view"] = feature_view
            predictions["feature_view_label"] = FILTER_LABELS.get(feature_view, feature_view)
            prediction_tables.append(predictions)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    metrics, confusion = metrics_for_predictions(pd, stats, predictions, {task.name: task for task in TASKS})
    metrics["feature_view"] = feature_view
    metrics["feature_view_label"] = FILTER_LABELS.get(feature_view, feature_view)
    confusion["feature_view"] = feature_view
    confusion["feature_view_label"] = FILTER_LABELS.get(feature_view, feature_view)
    return predictions, metrics, confusion, feature_sets


def best_h_e_rows(metrics):
    h_e = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["feature_set"] != "erbb2_rna_reference_not_h_e")
    ].copy()
    h_e["task_order"] = h_e["task"].map({task: idx for idx, task in enumerate(TASK_ORDER)})
    h_e["view_order"] = h_e["feature_view"].map({view: idx for idx, view in enumerate(FILTER_ORDER)})
    return (
        h_e.sort_values(["view_order", "task_order", "balanced_accuracy"], ascending=[True, True, False])
        .groupby(["feature_view", "task"], as_index=False)
        .head(1)
        .sort_values(["view_order", "task_order"])
    )


def best_reference_rows(metrics):
    ref = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["feature_set"] == "erbb2_rna_reference_not_h_e")
    ].copy()
    ref["task_order"] = ref["task"].map({task: idx for idx, task in enumerate(TASK_ORDER)})
    ref["view_order"] = ref["feature_view"].map({view: idx for idx, view in enumerate(FILTER_ORDER)})
    return (
        ref.sort_values(["view_order", "task_order", "balanced_accuracy"], ascending=[True, True, False])
        .groupby(["feature_view", "task"], as_index=False)
        .head(1)
        .sort_values(["view_order", "task_order"])
    )


def plot_best_by_view(plt, sns, best_rows, asset_dir: Path) -> None:
    plot_df = best_rows.copy()
    plot_df["feature_view_label"] = plot_df["feature_view"].map(FILTER_LABELS)
    plt.figure(figsize=(10.5, 5.6))
    sns.barplot(
        data=plot_df,
        x="feature_view_label",
        y="balanced_accuracy",
        hue="task_label",
        order=[FILTER_LABELS[view] for view in FILTER_ORDER],
    )
    plt.axhline(0.5, color="#7a7a7a", linestyle="--", linewidth=1, label="Binary chance")
    plt.axhline(1 / 3, color="#b0b0b0", linestyle=":", linewidth=1, label="Three-class chance")
    plt.ylim(0, 1)
    plt.xlabel("GigaTIME cleanup view")
    plt.ylabel("Best leave-one-out balanced accuracy")
    plt.title("Best GigaTIME/H&E Classifier by Cleanup View")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "cleaned_classifier_best_by_view.png", dpi=180)
    plt.close()


def plot_low_zero_feature_sets(plt, sns, metrics, asset_dir: Path) -> None:
    plot_df = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["task"] == "her2_low_vs_zero")
        & (metrics["feature_set"] != "erbb2_rna_reference_not_h_e")
    ].copy()
    plot_df["feature_view_label"] = plot_df["feature_view"].map(FILTER_LABELS)
    plot_df["feature_set_label"] = plot_df["feature_set"].map(FEATURE_LABELS).fillna(plot_df["feature_set"])
    plt.figure(figsize=(11.5, 5.8))
    sns.barplot(
        data=plot_df,
        x="feature_view_label",
        y="balanced_accuracy",
        hue="feature_set_label",
        order=[FILTER_LABELS[view] for view in FILTER_ORDER],
    )
    plt.axhline(0.5, color="#7a7a7a", linestyle="--", linewidth=1, label="Binary chance")
    plt.ylim(0, 1)
    plt.xlabel("GigaTIME cleanup view")
    plt.ylabel("HER2-low vs HER2-zero balanced accuracy")
    plt.title("HER2-Low vs HER2-Zero Classifier Across Cleaned Feature Sets")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "cleaned_classifier_low_zero_feature_sets.png", dpi=180)
    plt.close()


def plot_low_zero_confusions(plt, sns, pd, predictions, best_rows, asset_dir: Path) -> None:
    low_zero = best_rows.loc[best_rows["task"] == "her2_low_vs_zero"].copy()
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 7.6))
    axes = axes.ravel()
    task = next(task for task in TASKS if task.name == "her2_low_vs_zero")
    for axis, view in zip(axes, FILTER_ORDER):
        best = low_zero.loc[low_zero["feature_view"] == view]
        if best.empty:
            axis.axis("off")
            continue
        best_row = best.iloc[0]
        subset = predictions.loc[
            (predictions["feature_view"] == view)
            & (predictions["task"] == "her2_low_vs_zero")
            & (predictions["feature_set"] == best_row["feature_set"])
            & (predictions["model"] == "regularized_logistic")
        ]
        matrix = pd.DataFrame(0, index=task.class_order, columns=task.class_order)
        for _, row in subset.iterrows():
            matrix.loc[row["true_label"], row["predicted_label"]] += 1
        sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=0.4, ax=axis)
        axis.set_title(f"{FILTER_LABELS[view]}\n{FEATURE_LABELS.get(best_row['feature_set'], best_row['feature_set'])}")
        axis.set_xlabel("Predicted")
        axis.set_ylabel("True")
    fig.suptitle("HER2-Low vs HER2-Zero Confusion Matrices by Cleanup View")
    fig.tight_layout()
    fig.savefig(asset_dir / "cleaned_classifier_low_zero_confusions.png", dpi=180)
    plt.close(fig)


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def write_summary(path: Path, metrics, best_h_e, best_ref, feature_sets_by_view: dict[str, dict[str, list[str]]]) -> None:
    best_rows = []
    for _, row in best_h_e.iterrows():
        best_rows.append(
            [
                row["feature_view_label"],
                row["task_label"],
                FEATURE_LABELS.get(row["feature_set"], row["feature_set"]),
                row["n_cases"],
                fmt(float(row["accuracy"])),
                fmt(float(row["balanced_accuracy"])),
                fmt(float(row["macro_auc_ovr"])),
                fmt(float(row.get("sensitivity", np.nan))),
                fmt(float(row.get("specificity", np.nan))),
            ]
        )

    low_zero_rows = []
    for _, row in best_h_e.loc[best_h_e["task"] == "her2_low_vs_zero"].iterrows():
        low_zero_rows.append(
            [
                row["feature_view_label"],
                FEATURE_LABELS.get(row["feature_set"], row["feature_set"]),
                fmt(float(row["accuracy"])),
                fmt(float(row["balanced_accuracy"])),
                fmt(float(row["macro_auc_ovr"])),
                fmt(float(row.get("sensitivity", np.nan))),
                fmt(float(row.get("specificity", np.nan))),
            ]
        )

    ref_rows = []
    for _, row in best_ref.iterrows():
        ref_rows.append(
            [
                row["feature_view_label"],
                row["task_label"],
                fmt(float(row["accuracy"])),
                fmt(float(row["balanced_accuracy"])),
                fmt(float(row["macro_auc_ovr"])),
            ]
        )

    def metric_lookup(view: str, task: str):
        matches = best_h_e.loc[(best_h_e["feature_view"] == view) & (best_h_e["task"] == task)]
        return matches.iloc[0] if not matches.empty else None

    all_low_zero = metric_lookup("all_sampled_tissue", "her2_low_vs_zero")
    qc_low_zero = metric_lookup("qc_cellular_tissue", "her2_low_vs_zero")
    ck50_low_zero = metric_lookup("ck_enriched_top50", "her2_low_vs_zero")
    ck25_low_zero = metric_lookup("ck_enriched_top25", "her2_low_vs_zero")
    ck25_pos = metric_lookup("ck_enriched_top25", "her2_positive_vs_negative")

    feature_rows = []
    first_view = next(iter(feature_sets_by_view))
    for name, cols in feature_sets_by_view[first_view].items():
        feature_rows.append([FEATURE_LABELS.get(name, name), str(len(cols))])

    lines = [
        "# Cleaned GigaTIME HER2 Classifier Comparison",
        "",
        "This analysis reruns the slide-level HER2 classifier after cleaning the GigaTIME tile inputs. It compares all sampled tissue against cellular and virtual CK-enriched feature views.",
        "",
        "Every prediction is leave-one-out cross-validated. This remains a 30-case pilot, not a clinical model.",
        "",
        "## Feature Views",
        "",
        "- All sampled tissue: the original 256-tile slide aggregation.",
        "- QC cellular tissue: tissue fraction >= 0.70 and virtual DAPI mean >= 0.05.",
        "- CK-enriched top 50%: top half of virtual CK tiles within each slide after QC.",
        "- CK-enriched top 25%: top quarter of virtual CK tiles within each slide after QC.",
        "",
        "Virtual CK and DAPI are GigaTIME predictions, not real stains or pathologist tumor annotations.",
        "",
        "## Feature Sets",
        "",
        markdown_table(["Feature set", "Number of features"], feature_rows),
        "",
        "## Best GigaTIME/H&E Result Per View and Task",
        "",
        markdown_table(
            [
                "Cleanup view",
                "Task",
                "Best feature set",
                "N",
                "Accuracy",
                "Balanced accuracy",
                "Macro AUC",
                "Sensitivity",
                "Specificity",
            ],
            best_rows,
        ),
        "",
        "![Best classifier by cleanup view](assets/clinical_her2_cleaned_classifier/cleaned_classifier_best_by_view.png)",
        "",
        "## HER2-Low Versus HER2-Zero Focus",
        "",
        markdown_table(
            ["Cleanup view", "Best feature set", "Accuracy", "Balanced accuracy", "Macro AUC", "Sensitivity", "Specificity"],
            low_zero_rows,
        ),
        "",
        "![HER2-low versus HER2-zero feature-set comparison](assets/clinical_her2_cleaned_classifier/cleaned_classifier_low_zero_feature_sets.png)",
        "",
        "![HER2-low versus HER2-zero confusion matrices](assets/clinical_her2_cleaned_classifier/cleaned_classifier_low_zero_confusions.png)",
        "",
        "## Main Result",
        "",
        (
            f"- All sampled tissue HER2-low versus HER2-zero balanced accuracy: {fmt(float(all_low_zero['balanced_accuracy']))}, "
            f"macro AUC: {fmt(float(all_low_zero['macro_auc_ovr']))}."
            if all_low_zero is not None
            else "- All sampled tissue result was not available."
        ),
        (
            f"- QC cellular tissue preserved the HER2-low versus HER2-zero result: balanced accuracy {fmt(float(qc_low_zero['balanced_accuracy']))}, "
            f"macro AUC {fmt(float(qc_low_zero['macro_auc_ovr']))}."
            if qc_low_zero is not None
            else "- QC cellular tissue result was not available."
        ),
        (
            f"- CK-enriched top 50% reduced HER2-low versus HER2-zero performance to balanced accuracy {fmt(float(ck50_low_zero['balanced_accuracy']))}."
            if ck50_low_zero is not None
            else "- CK-enriched top 50% result was not available."
        ),
        (
            f"- CK-enriched top 25% also reduced HER2-low versus HER2-zero performance to balanced accuracy {fmt(float(ck25_low_zero['balanced_accuracy']))}."
            if ck25_low_zero is not None
            else "- CK-enriched top 25% result was not available."
        ),
        (
            f"- CK-enriched top 25% modestly improved HER2-positive versus HER2-negative balanced accuracy to {fmt(float(ck25_pos['balanced_accuracy']))}, "
            f"but sensitivity remained low at {fmt(float(ck25_pos.get('sensitivity', np.nan)))}."
            if ck25_pos is not None
            else "- CK-enriched HER2-positive result was not available."
        ),
        "",
        "## ERBB2 RNA Reference",
        "",
        "ERBB2 RNA is repeated as a non-H&E reference. It does not depend on the cleanup view and should not be interpreted as image-derived performance.",
        "",
        markdown_table(["Cleanup view", "Task", "Accuracy", "Balanced accuracy", "Macro AUC"], ref_rows),
        "",
        "## Interpretation",
        "",
        "The cleaned-view comparison asks whether the classifier signal is stronger in tumor-enriched tile views or in broader tissue context.",
        "",
        "In this run, cellular-tissue cleanup preserved the HER2-low versus HER2-zero classifier signal, which argues against the result being only blank/background artifact. However, the signal weakened when restricted to the most CK-enriched tiles. That suggests the current GigaTIME signal may be capturing broader tissue or microenvironment context more than a purely epithelial tumor-cell HER2 phenotype.",
        "",
        "The practical next step is not to claim diagnosis. It is to inspect the cases that change classification between all-tissue/QC-cellular and CK-enriched views, because those flips can reveal whether GigaTIME is using tumor regions, stromal context, immune infiltrates, or tile-selection artifacts.",
        "",
        "## Outputs",
        "",
        "- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_predictions.csv`",
        "- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_metrics.csv`",
        "- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_confusion_matrices.csv`",
        "- `docs/assets/clinical_her2_cleaned_classifier/`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_libs(out_dir / ".matplotlib")
    rows = pd.read_csv(args.cleaned_features)

    prediction_tables = []
    metric_tables = []
    confusion_tables = []
    feature_sets_by_view = {}
    for view in FILTER_ORDER:
        predictions, metrics, confusion, feature_sets = run_view(pd, optimize, stats, rows, view, args.l2_penalty)
        prediction_tables.append(predictions)
        metric_tables.append(metrics)
        confusion_tables.append(confusion)
        feature_sets_by_view[view] = feature_sets

    predictions = pd.concat(prediction_tables, ignore_index=True)
    metrics = pd.concat(metric_tables, ignore_index=True)
    confusion = pd.concat(confusion_tables, ignore_index=True)
    best_h_e = best_h_e_rows(metrics)
    best_ref = best_reference_rows(metrics)

    predictions.to_csv(out_dir / "cleaned_classifier_predictions.csv", index=False)
    metrics.to_csv(out_dir / "cleaned_classifier_metrics.csv", index=False)
    confusion.to_csv(out_dir / "cleaned_classifier_confusion_matrices.csv", index=False)
    best_h_e.to_csv(out_dir / "cleaned_classifier_best_h_e_metrics.csv", index=False)
    best_ref.to_csv(out_dir / "cleaned_classifier_erbb2_reference_metrics.csv", index=False)
    (out_dir / "cleaned_classifier_feature_sets.json").write_text(
        json.dumps(feature_sets_by_view, indent=2) + "\n",
        encoding="utf-8",
    )

    plot_best_by_view(plt, sns, best_h_e, asset_dir)
    plot_low_zero_feature_sets(plt, sns, metrics, asset_dir)
    plot_low_zero_confusions(plt, sns, pd, predictions, best_h_e, asset_dir)
    write_summary(out_dir / "cleaned_classifier_summary.md", metrics, best_h_e, best_ref, feature_sets_by_view)
    write_summary(
        Path("docs/clinical_her2_cleaned_classifier_comparison.md"),
        metrics,
        best_h_e,
        best_ref,
        feature_sets_by_view,
    )
    print(f"Wrote cleaned classifier outputs to {out_dir}")
    print(f"Wrote cleaned classifier figures to {asset_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
