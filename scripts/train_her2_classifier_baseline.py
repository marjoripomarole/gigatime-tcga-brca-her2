#!/usr/bin/env python3
"""Train baseline slide-level HER2 classifiers from GigaTIME features."""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np


GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
VIRTUAL_PROGRAMS = {
    "virtual_myeloid_checkpoint": ["mean_CD68", "mean_CD11c", "mean_PD-L1"],
    "virtual_t_cell_checkpoint": ["mean_CD3", "mean_CD4", "mean_CD8", "mean_PD-1"],
    "virtual_all_immune_checkpoint": [
        "mean_CD3",
        "mean_CD4",
        "mean_CD8",
        "mean_CD20",
        "mean_CD68",
        "mean_CD11c",
        "mean_PD-1",
        "mean_PD-L1",
    ],
    "virtual_proliferation": ["mean_Ki67"],
    "virtual_epithelial": ["mean_CK"],
}


@dataclass(frozen=True)
class TaskSpec:
    name: str
    label: str
    class_order: list[str]
    positive_class: str | None = None


TASKS = [
    TaskSpec(
        name="her2_positive_vs_negative",
        label="HER2-positive vs HER2-negative",
        class_order=["HER2-negative", "HER2-positive"],
        positive_class="HER2-positive",
    ),
    TaskSpec(
        name="her2_low_vs_zero",
        label="HER2-low vs HER2-zero",
        class_order=["HER2-low", "HER2-zero"],
        positive_class="HER2-zero",
    ),
    TaskSpec(
        name="her2_three_class",
        label="HER2-positive vs HER2-low vs HER2-zero",
        class_order=GROUP_ORDER,
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--joined",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/"
            "joined_slide_clinical_her2_gigatime.csv"
        ),
        help="Joined clinical HER2 + GigaTIME slide score table.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline",
    )
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def require_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import optimize, stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd, plt, sns, optimize, stats


def feature_sets_for(joined) -> dict[str, list[str]]:
    mean_cols = [
        col
        for col in joined.columns
        if col.startswith("mean_") and col != "mean_tissue_fraction" and joined[col].dtype.kind in "biufc"
    ]
    frac_cols = [col for col in joined.columns if col.startswith("frac_") and joined[col].dtype.kind in "biufc"]
    interpretable = [
        col
        for col in [
            "mean_CD3",
            "mean_CD8",
            "mean_CD4",
            "mean_CD20",
            "mean_CD68",
            "mean_CD11c",
            "mean_PD-1",
            "mean_PD-L1",
            "mean_CK",
            "mean_Ki67",
        ]
        if col in joined.columns
    ]
    feature_sets = {
        "gigatime_mean_channels": mean_cols,
        "gigatime_mean_and_fraction_channels": mean_cols + frac_cols,
        "interpretable_marker_means": interpretable,
        "virtual_programs": list(VIRTUAL_PROGRAMS),
    }
    if "erbb2_tpm" in joined.columns:
        feature_sets["erbb2_rna_reference_not_h_e"] = ["erbb2_tpm"]
    return feature_sets


def add_virtual_programs(joined):
    joined = joined.copy()
    for program, cols in VIRTUAL_PROGRAMS.items():
        available = [col for col in cols if col in joined.columns]
        joined[program] = joined[available].mean(axis=1) if available else np.nan
    return joined


def labels_for_task(rows, task: TaskSpec):
    if task.name == "her2_positive_vs_negative":
        return rows["clinical_her2_group"].map(
            {
                "HER2-positive": "HER2-positive",
                "HER2-low": "HER2-negative",
                "HER2-zero": "HER2-negative",
            }
        )
    if task.name == "her2_low_vs_zero":
        return rows["clinical_her2_group"].where(rows["clinical_her2_group"].isin(["HER2-low", "HER2-zero"]))
    if task.name == "her2_three_class":
        return rows["clinical_her2_group"].where(rows["clinical_her2_group"].isin(GROUP_ORDER))
    raise ValueError(task.name)


def standardize_train_test(x_train: np.ndarray, x_test: np.ndarray):
    center = np.nanmean(x_train, axis=0)
    scale = np.nanstd(x_train, axis=0, ddof=0)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    x_train = np.where(np.isfinite(x_train), x_train, center)
    x_test = np.where(np.isfinite(x_test), x_test, center)
    return (x_train - center) / scale, (x_test - center) / scale


def softmax(scores: np.ndarray) -> np.ndarray:
    scores = scores - np.max(scores, axis=1, keepdims=True)
    exp_scores = np.exp(scores)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def fit_predict_nearest_centroid(x_train, y_train, x_test, classes):
    centroids = []
    for class_idx in range(len(classes)):
        centroids.append(x_train[y_train == class_idx].mean(axis=0))
    centroids = np.vstack(centroids)
    distances = ((x_test[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    scores = -distances
    return softmax(scores)


def fit_predict_logistic(optimize, x_train, y_train, x_test, classes, l2_penalty: float):
    n_features = x_train.shape[1]
    n_classes = len(classes)
    x_aug = np.c_[np.ones(x_train.shape[0]), x_train]
    x_test_aug = np.c_[np.ones(x_test.shape[0]), x_test]

    def unpack(params):
        return params.reshape(n_features + 1, n_classes)

    def objective(params):
        weights = unpack(params)
        logits = x_aug @ weights
        probs = softmax(logits)
        nll = -np.log(probs[np.arange(len(y_train)), y_train] + 1e-12).mean()
        penalty = 0.5 * l2_penalty * np.sum(weights[1:] ** 2) / len(y_train)
        grad_logits = probs
        grad_logits[np.arange(len(y_train)), y_train] -= 1.0
        grad = x_aug.T @ grad_logits / len(y_train)
        grad[1:] += l2_penalty * weights[1:] / len(y_train)
        return nll + penalty, grad.ravel()

    initial = np.zeros((n_features + 1, n_classes), dtype=float).ravel()
    result = optimize.minimize(objective, initial, jac=True, method="L-BFGS-B", options={"maxiter": 500})
    weights = unpack(result.x)
    return softmax(x_test_aug @ weights)


def leave_one_out_predictions(pd, optimize, rows, task: TaskSpec, feature_name: str, feature_cols: list[str], l2_penalty: float):
    labels = labels_for_task(rows, task)
    task_df = rows.loc[labels.notna()].copy()
    task_df["task_label"] = labels.loc[labels.notna()].values
    task_df = task_df.loc[task_df["task_label"].isin(task.class_order)].copy()
    class_to_idx = {label: idx for idx, label in enumerate(task.class_order)}
    y = task_df["task_label"].map(class_to_idx).to_numpy(dtype=int)
    x = task_df[feature_cols].to_numpy(dtype=float)

    prediction_rows = []
    for test_idx in range(len(task_df)):
        train_idx = np.array([idx for idx in range(len(task_df)) if idx != test_idx])
        x_train, x_test = standardize_train_test(x[train_idx], x[[test_idx]])
        y_train = y[train_idx]
        if len(set(y_train.tolist())) != len(task.class_order):
            continue
        model_probs = {
            "nearest_centroid": fit_predict_nearest_centroid(x_train, y_train, x_test, task.class_order)[0],
            "regularized_logistic": fit_predict_logistic(
                optimize, x_train, y_train, x_test, task.class_order, l2_penalty
            )[0],
        }
        row = task_df.iloc[test_idx]
        for model_name, probs in model_probs.items():
            pred_idx = int(np.argmax(probs))
            pred_label = task.class_order[pred_idx]
            pred_row = {
                "task": task.name,
                "task_label": task.label,
                "feature_set": feature_name,
                "model": model_name,
                "case_submitter_id": row["case_submitter_id"],
                "clinical_her2_group": row["clinical_her2_group"],
                "true_label": row["task_label"],
                "predicted_label": pred_label,
                "correct": pred_label == row["task_label"],
            }
            for class_idx, class_label in enumerate(task.class_order):
                pred_row[f"prob_{class_label}"] = float(probs[class_idx])
            prediction_rows.append(pred_row)
    return pd.DataFrame(prediction_rows)


def binary_auc(stats, labels: np.ndarray, scores: np.ndarray) -> float:
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = stats.rankdata(np.concatenate([pos, neg]))
    pos_ranks = ranks[: len(pos)]
    return float((pos_ranks.sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def metrics_for_predictions(pd, stats, predictions, task_lookup: dict[str, TaskSpec]):
    metric_rows = []
    confusion_rows = []
    grouped = predictions.groupby(["task", "feature_set", "model"], sort=False)
    for (task_name, feature_set, model), group in grouped:
        task = task_lookup[task_name]
        class_order = task.class_order
        true_labels = group["true_label"].tolist()
        pred_labels = group["predicted_label"].tolist()
        n = len(group)
        accuracy = float(np.mean([true == pred for true, pred in zip(true_labels, pred_labels)]))
        recalls = []
        specificities = []
        aucs = []
        for class_label in class_order:
            tp = sum(t == class_label and p == class_label for t, p in zip(true_labels, pred_labels))
            fn = sum(t == class_label and p != class_label for t, p in zip(true_labels, pred_labels))
            tn = sum(t != class_label and p != class_label for t, p in zip(true_labels, pred_labels))
            fp = sum(t != class_label and p == class_label for t, p in zip(true_labels, pred_labels))
            recall = tp / (tp + fn) if (tp + fn) else float("nan")
            specificity = tn / (tn + fp) if (tn + fp) else float("nan")
            recalls.append(recall)
            specificities.append(specificity)
            binary_labels = np.array([1 if value == class_label else 0 for value in true_labels], dtype=int)
            auc = binary_auc(stats, binary_labels, group[f"prob_{class_label}"].to_numpy(dtype=float))
            aucs.append(auc)
            for predicted_class in class_order:
                confusion_rows.append(
                    {
                        "task": task_name,
                        "feature_set": feature_set,
                        "model": model,
                        "true_label": class_label,
                        "predicted_label": predicted_class,
                        "count": sum(
                            t == class_label and p == predicted_class for t, p in zip(true_labels, pred_labels)
                        ),
                    }
                )
        metric_row = {
            "task": task_name,
            "task_label": task.label,
            "feature_set": feature_set,
            "model": model,
            "n_cases": n,
            "n_classes": len(class_order),
            "accuracy": accuracy,
            "balanced_accuracy": float(np.nanmean(recalls)),
            "macro_specificity": float(np.nanmean(specificities)),
            "macro_auc_ovr": float(np.nanmean(aucs)),
        }
        if task.positive_class:
            positive = task.positive_class
            metric_row["positive_class"] = positive
            metric_row["sensitivity"] = recalls[class_order.index(positive)]
            metric_row["specificity"] = specificities[class_order.index(positive)]
            binary_labels = np.array([1 if value == positive else 0 for value in true_labels], dtype=int)
            metric_row["positive_auc"] = binary_auc(
                stats, binary_labels, group[f"prob_{positive}"].to_numpy(dtype=float)
            )
        metric_rows.append(metric_row)
    return pd.DataFrame(metric_rows), pd.DataFrame(confusion_rows)


def plot_performance(plt, sns, metrics, out_dir: Path) -> None:
    plot_df = metrics.loc[metrics["model"] == "regularized_logistic"].copy()
    if plot_df.empty:
        return
    label_map = {
        "gigatime_mean_channels": "GigaTIME\nmean channels",
        "gigatime_mean_and_fraction_channels": "GigaTIME\nmean + fraction",
        "interpretable_marker_means": "Interpretable\nmarker means",
        "virtual_programs": "Virtual\nprograms",
        "erbb2_rna_reference_not_h_e": "ERBB2 RNA\nreference",
    }
    plot_df["feature_set_display"] = plot_df["feature_set"].map(label_map).fillna(plot_df["feature_set"])
    display_order = [label_map.get(name, name) for name in plot_df["feature_set"].drop_duplicates()]
    plt.figure(figsize=(13.5, 6.0))
    sns.barplot(data=plot_df, x="feature_set_display", y="balanced_accuracy", hue="task_label", order=display_order)
    plt.axhline(0.5, color="#7a7a7a", linestyle="--", linewidth=1, label="Binary chance")
    plt.axhline(1 / 3, color="#b0b0b0", linestyle=":", linewidth=1, label="Three-class chance")
    plt.ylim(0, 1)
    plt.xticks(rotation=0, ha="center")
    plt.ylabel("Leave-one-out balanced accuracy")
    plt.xlabel("Feature set")
    plt.title("Cross-validated HER2 classifier performance")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(out_dir / "classifier_balanced_accuracy.png", dpi=180)
    plt.close()


def plot_confusion_matrices(plt, sns, pd, predictions, metrics, task_lookup: dict[str, TaskSpec], out_dir: Path) -> None:
    best_rows = (
        metrics.loc[metrics["model"] == "regularized_logistic"]
        .loc[metrics["feature_set"] != "erbb2_rna_reference_not_h_e"]
        .sort_values(["task", "balanced_accuracy"], ascending=[True, False])
        .groupby("task", as_index=False)
        .head(1)
    )
    for _, best in best_rows.iterrows():
        task = task_lookup[best["task"]]
        subset = predictions.loc[
            (predictions["task"] == best["task"])
            & (predictions["feature_set"] == best["feature_set"])
            & (predictions["model"] == best["model"])
        ]
        matrix = pd.DataFrame(0, index=task.class_order, columns=task.class_order)
        for _, row in subset.iterrows():
            matrix.loc[row["true_label"], row["predicted_label"]] += 1
        plt.figure(figsize=(4.8, 4.0))
        sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=0.4)
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.title(f"{task.label}\n{best['feature_set']}")
        plt.tight_layout()
        plt.savefig(out_dir / f"confusion_{best['task']}.png", dpi=180)
        plt.close()


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_markdown(path: Path, metrics, feature_metadata: dict[str, list[str]]) -> None:
    logistic = metrics.loc[metrics["model"] == "regularized_logistic"].copy()
    best = (
        logistic.sort_values(["task", "balanced_accuracy"], ascending=[True, False])
        .groupby("task", as_index=False)
        .head(1)
    )
    all_rows = logistic.sort_values(["task", "balanced_accuracy"], ascending=[True, False])
    lines = [
        "# Clinical HER2 Classifier Baseline",
        "",
        "This is the first diagnostic-model style analysis in the project. It uses slide-level GigaTIME features to predict clinical HER2 labels with leave-one-out cross-validation.",
        "",
        "Important: this is a tiny 30-case pilot. These results are useful for feasibility and failure-mode analysis, not for clinical diagnosis.",
        "",
        "## Tasks",
        "",
        "- HER2-positive vs HER2-negative, where HER2-negative combines HER2-low and HER2-zero.",
        "- HER2-low vs HER2-zero.",
        "- Three-class HER2-positive vs HER2-low vs HER2-zero.",
        "",
        "## Feature Sets",
        "",
    ]
    for name, cols in feature_metadata.items():
        note = "not H&E; RNA reference only" if name == "erbb2_rna_reference_not_h_e" else "GigaTIME/H&E-derived"
        lines.append(f"- `{name}`: {len(cols)} features ({note}).")
    lines.extend(
        [
            "",
            "## Best Cross-Validated Result Per Task",
            "",
            markdown_table(
                ["Task", "Best feature set", "Accuracy", "Balanced accuracy", "Macro AUC", "Sensitivity", "Specificity"],
                [
                    [
                        row["task_label"],
                        row["feature_set"],
                        fmt(row["accuracy"]),
                        fmt(row["balanced_accuracy"]),
                        fmt(row["macro_auc_ovr"]),
                        fmt(row.get("sensitivity", float("nan"))),
                        fmt(row.get("specificity", float("nan"))),
                    ]
                    for _, row in best.iterrows()
                ],
            ),
            "",
            "## All Regularized Logistic Results",
            "",
            markdown_table(
                ["Task", "Feature set", "Accuracy", "Balanced accuracy", "Macro AUC"],
                [
                    [
                        row["task_label"],
                        row["feature_set"],
                        fmt(row["accuracy"]),
                        fmt(row["balanced_accuracy"]),
                        fmt(row["macro_auc_ovr"]),
                    ]
                    for _, row in all_rows.iterrows()
                ],
            ),
            "",
            "![Classifier balanced accuracy](assets/clinical_her2_classifier_baseline/classifier_balanced_accuracy.png)",
            "",
            "## Interpretation Guardrails",
            "",
            "- Every prediction is cross-validated: the tested slide was left out of training.",
            "- The cohort is only 30 cases, so estimates can move a lot when cases are added.",
            "- HER2-low versus HER2-zero is especially difficult because the clinical boundary between IHC 0 and 1+ is subtle and noisy.",
            "- ERBB2 RNA is included only as a non-H&E reference. It should not be confused with an image-based diagnostic model.",
            "- If the GigaTIME-only feature sets do not perform well, the next step is not to claim failure; it is to improve tumor-region selection, add tile-level models, and validate on more cases.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_libs(out_dir / ".matplotlib")
    joined = pd.read_csv(args.joined)
    joined = add_virtual_programs(joined)
    feature_sets = feature_sets_for(joined)
    task_lookup = {task.name: task for task in TASKS}

    prediction_tables = []
    for task in TASKS:
        for feature_name, feature_cols in feature_sets.items():
            predictions = leave_one_out_predictions(
                pd, optimize, joined, task, feature_name, feature_cols, args.l2_penalty
            )
            prediction_tables.append(predictions)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    metrics, confusion = metrics_for_predictions(pd, stats, predictions, task_lookup)

    predictions.to_csv(out_dir / "classifier_crossval_predictions.csv", index=False)
    metrics.to_csv(out_dir / "classifier_metrics.csv", index=False)
    confusion.to_csv(out_dir / "classifier_confusion_matrices.csv", index=False)
    (out_dir / "classifier_feature_sets.json").write_text(json.dumps(feature_sets, indent=2) + "\n", encoding="utf-8")

    plot_performance(plt, sns, metrics, out_dir)
    plot_confusion_matrices(plt, sns, pd, predictions, metrics, task_lookup, out_dir)
    write_markdown(out_dir / "classifier_baseline_summary.md", metrics, feature_sets)
    print(f"Wrote classifier baseline outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
