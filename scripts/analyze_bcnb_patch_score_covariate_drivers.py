#!/usr/bin/env python3
"""Analyze clinical and patch-QC drivers of BCNB patch model scores."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from analyze_bcnb_patch_embedding_control import NEGATIVE_CLASS, POSITIVE_CLASS
from analyze_classifier_permutation_sanity import fmt, make_repeated_stratified_folds, markdown_table
from train_her2_classifier_baseline import standardize_train_test


IMAGE_MODELS = [
    "H-Optimus-0 embedding",
    "Virchow2 embedding",
    "H-Optimus-0 + Virchow2",
    "Average probability ensemble",
]
PRIMARY_IMAGE_MODELS = [
    "H-Optimus-0 embedding",
    "Virchow2 embedding",
    "H-Optimus-0 + Virchow2",
]
PREDICTOR_LABELS = {
    "label_only": "HER2 label only",
    "grade_only": "Grade only",
    "er_pr_only": "ER/PR only",
    "subtype_only": "Subtype only",
    "ki67_only": "Ki67 only",
    "patch_qc": "Patch QC only",
    "clinical_covariates": "Clinical covariates",
    "clinical_plus_patch_qc": "Clinical + patch QC",
    "clinical_plus_patch_qc_plus_label": "Clinical + patch QC + HER2 label",
}
PREDICTOR_ORDER = list(PREDICTOR_LABELS)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--patient-predictions",
        default=(
            "results/bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero/"
            "bcnb_patch_stratified_patient_predictions.csv"
        ),
        help="Patient-level mean out-of-fold prediction table from the stratified BCNB analysis.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/bcnb_patch_score_covariate_drivers_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=Path("docs/assets/bcnb_patch_score_covariate_drivers_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--out-markdown",
        type=Path,
        default=Path("docs/bcnb_patch_score_covariate_drivers_hoptimus0_virchow2_hash_capped10_low_zero.md"),
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260604)
    parser.add_argument("--ridge-alpha", type=float, default=1.0)
    return parser.parse_args()


def require_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`.") from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return pd, plt, sns, stats


def load_predictions(pd, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Patient prediction table not found: {path}")
    rows = pd.read_csv(path, low_memory=False)
    rows["grade"] = pd.to_numeric(rows["grade"], errors="coerce")
    rows["ki67"] = pd.to_numeric(rows["ki67"], errors="coerce")
    for col in ["ER", "PR", "molecular_subtype", "aln_status"]:
        rows[col] = rows[col].fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown"})
    for col in ["n_manifest_patches", "n_used_patches", "mean_tissue_fraction", "min_tissue_fraction"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows["label_zero"] = (rows["clinical_her2_group"] == POSITIVE_CLASS).astype(float)
    return rows


def numeric_piece(pd, rows, cols: list[str], add_missing: bool = True):
    pieces = []
    for col in cols:
        if col not in rows.columns:
            continue
        values = pd.to_numeric(rows[col], errors="coerce")
        if add_missing:
            pieces.append(pd.DataFrame({f"{col}_missing": values.isna().astype(float)}))
        fill = float(values.median()) if values.notna().any() else 0.0
        pieces.append(pd.DataFrame({col: values.fillna(fill).astype(float)}))
    if not pieces:
        return pd.DataFrame(index=rows.index)
    return pd.concat(pieces, axis=1)


def categorical_piece(pd, rows, cols: list[str]):
    available = [col for col in cols if col in rows.columns]
    if not available:
        return pd.DataFrame(index=rows.index)
    return pd.get_dummies(rows[available].fillna("Unknown").astype(str), prefix=available, dummy_na=False)


def design_matrix(pd, rows, predictor_set: str):
    pieces = []
    if predictor_set in {"label_only", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(rows[["label_zero"]].astype(float).reset_index(drop=True))
    if predictor_set in {"grade_only", "clinical_covariates", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(numeric_piece(pd, rows, ["grade"]).reset_index(drop=True))
    if predictor_set in {"er_pr_only", "clinical_covariates", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(categorical_piece(pd, rows, ["ER", "PR"]).reset_index(drop=True))
    if predictor_set in {"subtype_only", "clinical_covariates", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(categorical_piece(pd, rows, ["molecular_subtype"]).reset_index(drop=True))
    if predictor_set in {"clinical_covariates", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(categorical_piece(pd, rows, ["aln_status"]).reset_index(drop=True))
    if predictor_set in {"ki67_only", "clinical_covariates", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(numeric_piece(pd, rows, ["ki67"]).reset_index(drop=True))
    if predictor_set in {"patch_qc", "clinical_plus_patch_qc", "clinical_plus_patch_qc_plus_label"}:
        pieces.append(
            numeric_piece(
                pd,
                rows,
                ["n_manifest_patches", "n_used_patches", "mean_tissue_fraction", "min_tissue_fraction"],
                add_missing=False,
            ).reset_index(drop=True)
        )
    pieces = [piece for piece in pieces if not piece.empty]
    if not pieces:
        raise ValueError(f"No design columns for predictor set: {predictor_set}")
    design = pd.concat(pieces, axis=1)
    return design.loc[:, ~design.columns.duplicated()].astype(float)


def fit_predict_ridge(x_train, y_train, x_test, alpha: float):
    x_aug = np.c_[np.ones(x_train.shape[0]), x_train]
    x_test_aug = np.c_[np.ones(x_test.shape[0]), x_test]
    penalty = np.eye(x_aug.shape[1], dtype=float) * alpha
    penalty[0, 0] = 0.0
    weights = np.linalg.solve(x_aug.T @ x_aug + penalty, x_aug.T @ y_train)
    return x_test_aug @ weights


def r2_score(y_true, y_pred):
    denom = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if denom <= 0:
        return float("nan")
    return 1.0 - float(np.sum((y_true - y_pred) ** 2)) / denom


def cv_predict_score(pd, rows, predictor_set: str, folds, alpha: float):
    design = design_matrix(pd, rows, predictor_set)
    x = design.to_numpy(dtype=float)
    y = rows["mean_prob_her2_zero"].to_numpy(dtype=float)
    prediction_rows = []
    for fold_id, (train_idx, test_idx) in enumerate(folds):
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        pred = fit_predict_ridge(x_train, y[train_idx], x_test, alpha)
        for local_idx, row_idx in enumerate(test_idx):
            prediction_rows.append(
                {
                    "fold_id": int(fold_id),
                    "patient_id": rows.iloc[row_idx]["patient_id"],
                    "observed_score": float(y[row_idx]),
                    "predicted_score": float(pred[local_idx]),
                }
            )
    pred_frame = pd.DataFrame(prediction_rows)
    patient_pred = (
        pred_frame.groupby("patient_id", as_index=False)
        .agg(observed_score=("observed_score", "first"), predicted_score=("predicted_score", "mean"))
        .merge(rows[["patient_id", "clinical_her2_group", "label_zero"]], on="patient_id", how="left", validate="one_to_one")
    )
    y_true = patient_pred["observed_score"].to_numpy(dtype=float)
    y_pred = patient_pred["predicted_score"].to_numpy(dtype=float)
    return {
        "n_predictors": int(design.shape[1]),
        "r2": r2_score(y_true, y_pred),
        "mae": float(np.mean(np.abs(y_true - y_pred))),
        "pearson_r": float(np.corrcoef(y_true, y_pred)[0, 1]) if len(y_true) > 1 else float("nan"),
        "patient_predictions": patient_pred,
    }


def auc_from_scores(stats, y_true, scores):
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    positive = scores[y_true == 1]
    negative = scores[y_true == 0]
    if len(positive) == 0 or len(negative) == 0:
        return float("nan"), float("nan")
    test = stats.mannwhitneyu(positive, negative, alternative="two-sided")
    auc = float(test.statistic) / float(len(positive) * len(negative))
    return auc, float(test.pvalue)


def residual_label_test(pd, stats, rows, residual_predictions):
    merged = rows.merge(
        residual_predictions[["patient_id", "predicted_score"]],
        on="patient_id",
        how="left",
        validate="one_to_one",
    )
    merged["residual_score"] = merged["mean_prob_her2_zero"] - merged["predicted_score"]
    y = merged["label_zero"].to_numpy(dtype=int)
    raw_auc, raw_p = auc_from_scores(stats, y, merged["mean_prob_her2_zero"].to_numpy(dtype=float))
    residual_auc, residual_p = auc_from_scores(stats, y, merged["residual_score"].to_numpy(dtype=float))
    low = merged.loc[merged["clinical_her2_group"] == NEGATIVE_CLASS]
    zero = merged.loc[merged["clinical_her2_group"] == POSITIVE_CLASS]
    return {
        "n_low": int(len(low)),
        "n_zero": int(len(zero)),
        "raw_auc": raw_auc,
        "raw_mannwhitney_p_value": raw_p,
        "raw_delta_zero_minus_low": float(zero["mean_prob_her2_zero"].mean() - low["mean_prob_her2_zero"].mean()),
        "residual_auc_after_clinical_patch_qc": residual_auc,
        "residual_mannwhitney_p_value": residual_p,
        "residual_delta_zero_minus_low": float(zero["residual_score"].mean() - low["residual_score"].mean()),
        "mean_absolute_residual": float(np.mean(np.abs(merged["residual_score"]))),
    }


def score_covariate_associations(pd, stats, rows):
    records = []
    numeric_cols = ["grade", "ki67", "n_manifest_patches", "mean_tissue_fraction", "min_tissue_fraction", "clinical_covariate_score"]
    for col in numeric_cols:
        if col not in rows.columns:
            continue
        pair = rows[[col, "mean_prob_her2_zero"]].dropna()
        if len(pair) < 4:
            continue
        if pair[col].nunique(dropna=True) < 2:
            continue
        corr = stats.spearmanr(pair[col], pair["mean_prob_her2_zero"])
        records.append(
            {
                "covariate": col,
                "association_type": "spearman",
                "n": int(len(pair)),
                "statistic": float(corr.statistic),
                "p_value": float(corr.pvalue),
            }
        )
    for col in ["ER", "PR", "molecular_subtype", "aln_status"]:
        groups = [
            group["mean_prob_her2_zero"].dropna().to_numpy(dtype=float)
            for _level, group in rows.groupby(col, dropna=False)
            if len(group["mean_prob_her2_zero"].dropna()) >= 3
        ]
        if len(groups) < 2:
            continue
        test = stats.kruskal(*groups)
        total = rows["mean_prob_her2_zero"].dropna().to_numpy(dtype=float)
        total_mean = np.mean(total)
        ss_total = float(np.sum((total - total_mean) ** 2))
        ss_between = 0.0
        for group_values in groups:
            ss_between += float(len(group_values) * (np.mean(group_values) - total_mean) ** 2)
        records.append(
            {
                "covariate": col,
                "association_type": "kruskal_eta2",
                "n": int(len(total)),
                "statistic": float(ss_between / ss_total) if ss_total > 0 else float("nan"),
                "p_value": float(test.pvalue),
            }
        )
    return pd.DataFrame(records)


def add_clinical_score(pd, model_rows, all_predictions):
    clinical = all_predictions.loc[all_predictions["model"] == "Clinical covariates", ["patient_id", "mean_prob_her2_zero"]].rename(
        columns={"mean_prob_her2_zero": "clinical_covariate_score"}
    )
    return model_rows.merge(clinical, on="patient_id", how="left", validate="one_to_one")


def analyze_models(pd, stats, predictions, args):
    rng = np.random.default_rng(args.seed)
    metric_rows = []
    residual_rows = []
    association_frames = []
    residual_prediction_frames = []
    for model in IMAGE_MODELS:
        rows = predictions.loc[predictions["model"] == model].copy().reset_index(drop=True)
        rows = add_clinical_score(pd, rows, predictions)
        y_label = rows["label_zero"].to_numpy(dtype=int)
        folds = make_repeated_stratified_folds(y_label, args.folds, args.repeats, rng)
        clinical_patch_predictions = None
        for predictor_set in PREDICTOR_ORDER:
            result = cv_predict_score(pd, rows, predictor_set, folds, args.ridge_alpha)
            metric_rows.append(
                {
                    "model": model,
                    "predictor_set": predictor_set,
                    "predictor_label": PREDICTOR_LABELS[predictor_set],
                    "n_predictors": result["n_predictors"],
                    "r2": result["r2"],
                    "mae": result["mae"],
                    "pearson_r": result["pearson_r"],
                }
            )
            if predictor_set == "clinical_plus_patch_qc":
                clinical_patch_predictions = result["patient_predictions"]
        if clinical_patch_predictions is None:
            raise RuntimeError("clinical_plus_patch_qc predictions were not produced")
        residual = residual_label_test(pd, stats, rows, clinical_patch_predictions)
        residual_rows.append({"model": model, **residual})
        residual_prediction = clinical_patch_predictions.copy()
        residual_prediction["model"] = model
        residual_prediction_frames.append(residual_prediction)
        associations = score_covariate_associations(pd, stats, rows)
        associations.insert(0, "model", model)
        association_frames.append(associations)
    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(residual_rows),
        pd.concat(association_frames, axis=0, ignore_index=True),
        pd.concat(residual_prediction_frames, axis=0, ignore_index=True),
    )


def plot_r2(plt, sns, metrics, asset_dir: Path):
    selected = metrics.loc[metrics["model"].isin(PRIMARY_IMAGE_MODELS)].copy()
    selected["predictor_label"] = selected["predictor_set"].map(PREDICTOR_LABELS)
    selected["predictor_label"] = selected["predictor_label"].astype(str)
    fig, axis = plt.subplots(figsize=(11.2, 5.8))
    sns.barplot(
        data=selected,
        x="predictor_label",
        y="r2",
        hue="model",
        order=[PREDICTOR_LABELS[key] for key in PREDICTOR_ORDER],
        ax=axis,
    )
    axis.axhline(0, color="#374151", linewidth=1)
    axis.set_xlabel("Predictors used to explain image-model score")
    axis.set_ylabel("Cross-validated R2 for image-model P(HER2-zero)")
    axis.set_title("Clinical and patch-QC drivers of BCNB image-model scores")
    axis.tick_params(axis="x", rotation=30)
    for label in axis.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_score_driver_cv_r2.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_residual_auc(plt, sns, residuals, asset_dir: Path):
    records = []
    for _, row in residuals.iterrows():
        if row["model"] not in PRIMARY_IMAGE_MODELS:
            continue
        records.append({"model": row["model"], "score_type": "Raw image score", "auc": row["raw_auc"]})
        records.append(
            {
                "model": row["model"],
                "score_type": "Residual after clinical + patch QC",
                "auc": row["residual_auc_after_clinical_patch_qc"],
            }
        )
    if not records:
        return
    import pandas as pd

    plot_df = pd.DataFrame(records)
    fig, axis = plt.subplots(figsize=(9.6, 5.0))
    sns.barplot(data=plot_df, x="model", y="auc", hue="score_type", ax=axis)
    axis.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    axis.set_ylim(0.35, 0.75)
    axis.set_xlabel("Image model")
    axis.set_ylabel("AUC for HER2-zero vs HER2-low")
    axis.set_title("Low/zero signal before and after clinical + patch-QC residualization")
    axis.tick_params(axis="x", rotation=18)
    for label in axis.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_score_residual_auc.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_clinical_score_scatter(plt, sns, predictions, asset_dir: Path):
    clinical = predictions.loc[predictions["model"] == "Clinical covariates", ["patient_id", "mean_prob_her2_zero"]].rename(
        columns={"mean_prob_her2_zero": "clinical_covariate_score"}
    )
    dual = predictions.loc[predictions["model"] == "H-Optimus-0 + Virchow2"].merge(
        clinical, on="patient_id", how="left", validate="one_to_one"
    )
    plot_df = dual.rename(columns={"clinical_her2_group": "HER2 group"})
    fig, axis = plt.subplots(figsize=(6.2, 5.6))
    sns.scatterplot(
        data=plot_df,
        x="clinical_covariate_score",
        y="mean_prob_her2_zero",
        hue="HER2 group",
        alpha=0.75,
        linewidth=0,
        ax=axis,
    )
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Clinical covariate model P(HER2-zero)")
    axis.set_ylabel("Dual image model P(HER2-zero)")
    axis.set_title("BCNB clinical-score versus dual-image-score agreement")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_clinical_vs_dual_image_score.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def table_rows(metrics, predictor_sets: list[str]) -> list[list[str]]:
    selected = metrics.loc[
        metrics["model"].isin(PRIMARY_IMAGE_MODELS) & metrics["predictor_set"].isin(predictor_sets)
    ].copy()
    model_order = {model: idx for idx, model in enumerate(PRIMARY_IMAGE_MODELS)}
    pred_order = {name: idx for idx, name in enumerate(predictor_sets)}
    selected["_model_order"] = selected["model"].map(model_order)
    selected["_pred_order"] = selected["predictor_set"].map(pred_order)
    selected = selected.sort_values(["_model_order", "_pred_order"])
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["model"],
                row["predictor_label"],
                str(int(row["n_predictors"])),
                fmt(row["r2"], 3),
                fmt(row["pearson_r"], 3),
                fmt(row["mae"], 3),
            ]
        )
    return rows


def residual_table_rows(residuals) -> list[list[str]]:
    selected = residuals.loc[residuals["model"].isin(PRIMARY_IMAGE_MODELS)].copy()
    model_order = {model: idx for idx, model in enumerate(PRIMARY_IMAGE_MODELS)}
    selected["_model_order"] = selected["model"].map(model_order)
    selected = selected.sort_values("_model_order")
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["model"],
                fmt(row["raw_auc"], 3),
                fmt(row["raw_delta_zero_minus_low"], 3),
                fmt(row["residual_auc_after_clinical_patch_qc"], 3),
                fmt(row["residual_delta_zero_minus_low"], 3),
                fmt(row["residual_mannwhitney_p_value"], 3),
            ]
        )
    return rows


def association_table_rows(associations) -> list[list[str]]:
    selected = associations.loc[
        associations["model"].isin(PRIMARY_IMAGE_MODELS)
        & associations["covariate"].isin(["clinical_covariate_score", "grade", "ki67", "ER", "PR", "molecular_subtype"])
    ].copy()
    model_order = {model: idx for idx, model in enumerate(PRIMARY_IMAGE_MODELS)}
    cov_order = {"clinical_covariate_score": 0, "grade": 1, "ki67": 2, "ER": 3, "PR": 4, "molecular_subtype": 5}
    selected["_model_order"] = selected["model"].map(model_order)
    selected["_cov_order"] = selected["covariate"].map(cov_order)
    selected = selected.sort_values(["_model_order", "_cov_order"])
    rows = []
    for _, row in selected.iterrows():
        stat_label = "rho" if row["association_type"] == "spearman" else "eta2"
        rows.append(
            [
                row["model"],
                row["covariate"],
                stat_label,
                fmt(row["statistic"], 3),
                fmt(row["p_value"], 3),
            ]
        )
    return rows


def write_markdown(path: Path, asset_dir: Path, args, metrics, residuals, associations):
    lines = [
        "# BCNB Patch Model Score Covariate Drivers",
        "",
        "Status: clinical and patch-QC driver analysis for BCNB image-model patient scores.",
        "",
        "## Method",
        "",
        "- Input scores: patient-mean out-of-fold P(HER2-zero) from the BCNB stratified analysis.",
        "- Target: explain each image model score using clinical covariates, patch-QC variables, and the true HER2-low/zero label.",
        f"- Score-explainer model: ridge linear regression with repeated stratified {args.folds}-fold CV ({args.repeats} repeats), evaluated by patient-mean out-of-fold R2.",
        "- Residualization: predict the image score from clinical covariates + patch QC, subtract that prediction, then test whether residual score still separates HER2-zero from HER2-low.",
        "",
        "## Cross-Validated Score Explainability",
        "",
        markdown_table(
            ["Image model", "Predictor set", "Features", "CV R2", "Pred/obs r", "MAE"],
            table_rows(
                metrics,
                [
                    "label_only",
                    "grade_only",
                    "er_pr_only",
                    "clinical_covariates",
                    "patch_qc",
                    "clinical_plus_patch_qc",
                    "clinical_plus_patch_qc_plus_label",
                ],
            ),
        ),
        "",
        f"![BCNB score driver R2](assets/{asset_dir.name}/bcnb_score_driver_cv_r2.png)",
        "",
        "## Residual Low/Zero Signal After Clinical + Patch QC",
        "",
        markdown_table(
            ["Image model", "Raw AUC", "Raw delta", "Residual AUC", "Residual delta", "Residual p"],
            residual_table_rows(residuals),
        ),
        "",
        f"![BCNB score residual AUC](assets/{asset_dir.name}/bcnb_score_residual_auc.png)",
        "",
        "## Direct Score-Covariate Associations",
        "",
        markdown_table(
            ["Image model", "Covariate", "Statistic", "Value", "p"],
            association_table_rows(associations),
        ),
        "",
        f"![BCNB clinical versus dual image score](assets/{asset_dir.name}/bcnb_clinical_vs_dual_image_score.png)",
        "",
        "## Interpretation",
        "",
        "- The true HER2-low/zero label alone explains very little of the image-model score, which is expected for a weak classifier.",
        "- Clinical covariates and patch QC explain a modest but non-trivial fraction of the image scores, especially for the dual embedding.",
        "- After clinical + patch-QC residualization, residual AUC remains above 0.5 but is still modest, so the image signal is not fully reducible to the measured covariates and is not strong enough for classifier-grade claims.",
        "- This supports the current paper framing: BCNB contains weak image-readable morphology/covariate signal around the low/zero boundary, but the available patch-level evidence does not justify a standalone HER2-low/zero detector.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_score_driver_r2.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_score_residual_label_tests.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_score_covariate_associations.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_score_residual_predictions.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def run_analysis(pd, plt, sns, stats, args):
    predictions = load_predictions(pd, Path(args.patient_predictions))
    metrics, residuals, associations, residual_predictions = analyze_models(pd, stats, predictions, args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.out_dir / "bcnb_patch_score_driver_r2.csv", index=False)
    residuals.to_csv(args.out_dir / "bcnb_patch_score_residual_label_tests.csv", index=False)
    associations.to_csv(args.out_dir / "bcnb_patch_score_covariate_associations.csv", index=False)
    residual_predictions.to_csv(args.out_dir / "bcnb_patch_score_residual_predictions.csv", index=False)
    (args.out_dir / "bcnb_patch_score_covariate_driver_metadata.json").write_text(
        json.dumps(
            {
                "task": "bcnb_patch_score_covariate_drivers",
                "patient_predictions": args.patient_predictions,
                "folds": args.folds,
                "repeats": args.repeats,
                "seed": args.seed,
                "ridge_alpha": args.ridge_alpha,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    plot_r2(plt, sns, metrics, args.asset_dir)
    plot_residual_auc(plt, sns, residuals, args.asset_dir)
    plot_clinical_score_scatter(plt, sns, predictions, args.asset_dir)
    write_markdown(args.out_markdown, args.asset_dir, args, metrics, residuals, associations)
    return metrics, residuals, associations


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_libs(args.out_dir / ".matplotlib")
    metrics, residuals, _associations = run_analysis(pd, plt, sns, stats, args)
    print(f"Wrote BCNB score-driver outputs to {args.out_dir}")
    print(f"Wrote BCNB score-driver markdown to {args.out_markdown}")
    print(
        metrics.loc[
            metrics["model"].isin(PRIMARY_IMAGE_MODELS)
            & metrics["predictor_set"].isin(["label_only", "clinical_covariates", "clinical_plus_patch_qc"])
        ][["model", "predictor_label", "r2", "pearson_r", "mae"]].to_string(index=False)
    )
    print(residuals.loc[residuals["model"].isin(PRIMARY_IMAGE_MODELS)].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
