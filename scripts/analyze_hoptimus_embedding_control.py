#!/usr/bin/env python3
"""Generic H&E embedding control for HER2-low versus HER2-zero.

This asks a single disciplined question: does a generic pathology foundation-model
embedding (H-Optimus-0), which has no virtual-immune interpretation, separate
HER2-low from HER2-zero as well as the GigaTIME virtual-mIF channels do, and does
it collapse under leave-source-site-out validation the same way?

If a generic embedding also separates the groups and also collapses under
source-site holdout while slide-size covariates stay strong, then the GigaTIME
"virtual immune/checkpoint" framing is not required to explain the signal: the
low-versus-zero axis is better described as generic morphology/tissue-composition
that tracks TCGA acquisition structure.

The embedding is one mean-pooled vector per slide. It is compared on identical
cross-validation folds against slide-size, source-site, and GigaTIME baselines,
with PCA fit inside each training fold to avoid leakage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_classifier_permutation_sanity import (
    LOW_ZERO_GROUPS,
    NEGATIVE_CLASS,
    POSITIVE_CLASS,
    benjamini_hochberg,
    fmt,
    make_repeated_stratified_folds,
    markdown_table,
    metric_dict,
)
from analyze_clinical_covariate_sensitivity import (
    add_covariates,
    design_matrix_for_set,
    require_analysis_libs,
)
from analyze_source_site_generalization import leave_source_site_out_folds
from train_her2_classifier_baseline import fit_predict_logistic, standardize_train_test


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
HOPTIMUS_DIR = Path("results/hoptimus_tcga_brca_high_trust_tile128")

# Feature sets compared on identical folds. Order controls table/plot order.
FEATURE_SETS = [
    ("slide_size_only", "Slide-size covariates"),
    ("source_site_only", "Source-site covariates"),
    ("gigatime_mean_channels", "GigaTIME mean channels"),
    ("hoptimus_embedding", "H-Optimus-0 embedding (PCA)"),
    ("hoptimus_plus_slide_size", "H-Optimus-0 + slide-size"),
]
EMBEDDING_SETS = {"hoptimus_embedding", "hoptimus_plus_slide_size"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-features",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_slide_features.csv"),
        help="GigaTIME tumor-proxy slide feature table (carries GigaTIME channels and covariate joins).",
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--embeddings",
        default=str(HOPTIMUS_DIR / "slide_embeddings.csv"),
        help="H-Optimus-0 mean-pooled slide embeddings.",
    )
    parser.add_argument(
        "--reference-view",
        default="qc_cellular_tissue",
        help="GigaTIME feature view used as the whole-slide comparator for the embedding.",
    )
    parser.add_argument("--out-dir", default=str(HOPTIMUS_DIR / "embedding_low_zero_control"))
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_hoptimus_embedding_control",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_hoptimus_embedding_control.md",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--pca-components", type=int, default=20)
    parser.add_argument(
        "--pca-grid",
        default="10,20,30",
        help="Comma-separated PCA component counts for the embedding robustness mini-table.",
    )
    parser.add_argument("--seed", type=int, default=20260604)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    parser.add_argument("--min-site-count", type=int, default=5)
    return parser.parse_args()


def embedding_columns(frame) -> list[str]:
    return sorted(col for col in frame.columns if col.startswith("embedding_"))


def pca_fit_transform(x_train: np.ndarray, x_test: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Fit PCA on the (already standardized) training fold only, project both."""
    k = int(min(k, x_train.shape[0], x_train.shape[1]))
    mean = x_train.mean(axis=0)
    x_train_c = x_train - mean
    x_test_c = x_test - mean
    _u, _s, vt = np.linalg.svd(x_train_c, full_matrices=False)
    components = vt[:k]
    return x_train_c @ components.T, x_test_c @ components.T


def design_for_feature_set(pd, rows, feature_set: str, embed_cols: list[str]):
    if feature_set == "hoptimus_embedding":
        return rows[embed_cols].astype(float)
    if feature_set == "hoptimus_plus_slide_size":
        size = design_matrix_for_set(pd, rows, "slide_size_only")
        return pd.concat([rows[embed_cols].astype(float).reset_index(drop=True), size.reset_index(drop=True)], axis=1)
    return design_matrix_for_set(pd, rows, feature_set)


def evaluate(optimize, stats, x: np.ndarray, y: np.ndarray, folds, feature_set: str, l2_penalty: float, pca_k: int) -> dict[str, float]:
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    true_labels: list[int] = []
    pred_labels: list[int] = []
    positive_probs: list[float] = []
    for train_idx, test_idx in folds:
        if len(np.unique(y[train_idx])) < 2:
            continue
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        if feature_set in EMBEDDING_SETS and pca_k and x_train.shape[1] > pca_k:
            x_train, x_test = pca_fit_transform(x_train, x_test, pca_k)
        probs = fit_predict_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(np.argmax(probs, axis=1).astype(int).tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
    metrics = metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )
    metrics["n_cv_predictions"] = len(true_labels)
    return metrics


def run_analysis(pd, optimize, stats, args):
    rng = np.random.default_rng(args.seed)
    features = pd.read_csv(args.slide_features)
    metadata = pd.read_csv(args.high_trust_slides)
    features = add_covariates(pd, features, metadata, args.min_site_count)

    rows = features.loc[
        (features["feature_view"] == args.reference_view)
        & (features["clinical_her2_group"].isin(LOW_ZERO_GROUPS))
    ].copy()

    embeddings = pd.read_csv(args.embeddings)
    embed_cols = embedding_columns(embeddings)
    if not embed_cols:
        raise SystemExit(f"No embedding_* columns found in {args.embeddings}.")
    rows = rows.merge(embeddings[["slide_id", *embed_cols]], on="slide_id", how="inner").reset_index(drop=True)
    rows["tss_code"] = rows["case_submitter_id"].astype(str).str.split("-").str[1]

    n_low = int((rows["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    n_zero = int((rows["clinical_her2_group"] == POSITIVE_CLASS).sum())
    if n_low < 2 or n_zero < 2:
        raise SystemExit(
            f"Too few matched low/zero slides with embeddings (low={n_low}, zero={n_zero}). "
            "Is the H-Optimus run finished?"
        )

    y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
    random_folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
    site_folds = leave_source_site_out_folds(rows, y)

    metrics_rows = []
    for feature_set, feature_set_label in FEATURE_SETS:
        design = design_for_feature_set(pd, rows, feature_set, embed_cols)
        if design.empty:
            continue
        x = design.to_numpy(dtype=float)
        for scheme, scheme_label, folds in [
            ("repeated_stratified_cv", "Repeated stratified CV", random_folds),
            ("leave_source_site_out", "Leave source site out", site_folds),
        ]:
            metrics = evaluate(optimize, stats, x, y, folds, feature_set, args.l2_penalty, args.pca_components)
            metrics_rows.append(
                {
                    "feature_set": feature_set,
                    "feature_set_label": feature_set_label,
                    "validation_scheme": scheme,
                    "validation_scheme_label": scheme_label,
                    "n_cases": int(len(rows)),
                    "n_low": n_low,
                    "n_zero": n_zero,
                    "n_features": int(design.shape[1]),
                    "pca_components": args.pca_components if feature_set in EMBEDDING_SETS else "",
                    "n_folds": len(folds),
                    **metrics,
                }
            )
    metrics = pd.DataFrame(metrics_rows)

    # PCA-k robustness for the embedding under repeated CV.
    pca_grid = [int(k) for k in str(args.pca_grid).split(",") if k.strip()]
    embed_x = rows[embed_cols].astype(float).to_numpy(dtype=float)
    pca_rows = []
    for k in pca_grid:
        m = evaluate(optimize, stats, embed_x, y, random_folds, "hoptimus_embedding", args.l2_penalty, k)
        pca_rows.append({"pca_components": k, **m})
    pca_robustness = pd.DataFrame(pca_rows)

    # Shuffled-label permutation null for the embedding under repeated CV.
    observed = float(
        metrics.loc[
            (metrics["feature_set"] == "hoptimus_embedding")
            & (metrics["validation_scheme"] == "repeated_stratified_cv"),
            "balanced_accuracy",
        ].iloc[0]
    )
    null_ba = []
    for _ in range(args.permutations):
        y_perm = y.copy()
        rng.shuffle(y_perm)
        m = evaluate(optimize, stats, embed_x, y_perm, random_folds, "hoptimus_embedding", args.l2_penalty, args.pca_components)
        null_ba.append(m["balanced_accuracy"])
    null_arr = np.array(null_ba, dtype=float)
    permutation = {
        "feature_set": "hoptimus_embedding",
        "pca_components": args.pca_components,
        "observed_repeated_cv_balanced_accuracy": observed,
        "n_permutations": args.permutations,
        "null_balanced_accuracy_mean": float(np.nanmean(null_arr)),
        "null_balanced_accuracy_sd": float(np.nanstd(null_arr, ddof=1)),
        "null_balanced_accuracy_p95": float(np.nanquantile(null_arr, 0.95)),
        "empirical_p_balanced_accuracy": float((1 + np.sum(null_arr >= observed)) / (1 + args.permutations)),
    }
    return rows, metrics, pca_robustness, permutation, pd.DataFrame({"null_balanced_accuracy": null_arr})


def metric_table_rows(metrics) -> list[list[str]]:
    order = {name: idx for idx, (name, _) in enumerate(FEATURE_SETS)}
    scheme_order = {"repeated_stratified_cv": 0, "leave_source_site_out": 1}
    selected = metrics.copy()
    selected["_f"] = selected["feature_set"].map(order)
    selected["_s"] = selected["validation_scheme"].map(scheme_order)
    selected = selected.sort_values(["_f", "_s"])
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["feature_set_label"],
                row["validation_scheme_label"],
                str(int(row["n_features"])),
                fmt(row["balanced_accuracy"], 3),
                fmt(row["macro_auc_ovr"], 3),
                fmt(row["sensitivity"], 3),
                fmt(row["specificity"], 3),
            ]
        )
    return rows


def plot_metric_comparison(plt, sns, metrics, asset_dir: Path) -> None:
    order = [label for _, label in FEATURE_SETS]
    plt.figure(figsize=(11.0, 5.6))
    sns.barplot(
        data=metrics,
        x="feature_set_label",
        y="balanced_accuracy",
        hue="validation_scheme_label",
        order=order,
    )
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Feature set")
    plt.ylabel("Balanced accuracy (HER2-low vs HER2-zero)")
    plt.title("Generic H-Optimus-0 Embedding vs GigaTIME and Confound Baselines")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "embedding_control_balanced_accuracy.png", dpi=180)
    plt.close()


def pick(metrics, feature_set: str, scheme: str, column: str) -> float:
    subset = metrics.loc[
        (metrics["feature_set"] == feature_set) & (metrics["validation_scheme"] == scheme), column
    ]
    return float(subset.iloc[0]) if len(subset) else float("nan")


def build_verdict(metrics, permutation) -> list[str]:
    emb_cv = pick(metrics, "hoptimus_embedding", "repeated_stratified_cv", "balanced_accuracy")
    emb_site = pick(metrics, "hoptimus_embedding", "leave_source_site_out", "balanced_accuracy")
    giga_cv = pick(metrics, "gigatime_mean_channels", "repeated_stratified_cv", "balanced_accuracy")
    giga_site = pick(metrics, "gigatime_mean_channels", "leave_source_site_out", "balanced_accuracy")
    size_cv = pick(metrics, "slide_size_only", "repeated_stratified_cv", "balanced_accuracy")
    size_site = pick(metrics, "slide_size_only", "leave_source_site_out", "balanced_accuracy")
    emp_p = permutation["empirical_p_balanced_accuracy"]

    separates = emb_cv >= 0.60 and emp_p < 0.05
    collapses = (emb_cv - emb_site) >= 0.05
    size_competitive = size_cv >= emb_cv - 0.02
    size_robust_under_holdout = size_site >= emb_site

    lines = [
        f"- The generic H-Optimus-0 embedding reaches balanced accuracy {fmt(emb_cv, 3)} under repeated stratified CV "
        f"(shuffled-label null mean {fmt(permutation['null_balanced_accuracy_mean'], 3)}, empirical p {fmt(emp_p, 4)}), "
        f"versus {fmt(giga_cv, 3)} for GigaTIME mean channels and {fmt(size_cv, 3)} for slide-size covariates.",
        f"- Under leave-source-site-out validation the embedding moves to {fmt(emb_site, 3)} "
        f"(GigaTIME {fmt(giga_site, 3)}, slide-size {fmt(size_site, 3)}).",
    ]
    if separates and collapses and (size_competitive or size_robust_under_holdout):
        lines.append(
            "- Read together, a generic morphology embedding with no immune interpretation separates HER2-low from "
            "HER2-zero about as well as GigaTIME and also loses ground under source-site holdout while slide-size "
            "covariates stay strong. This supports the confound reading: the low-versus-zero axis is better explained "
            "as generic tissue/morphology that tracks TCGA acquisition structure than as GigaTIME-specific virtual "
            "immune biology."
        )
    elif separates and not collapses:
        lines.append(
            "- The generic embedding separates the groups and, unlike GigaTIME, holds up better under source-site "
            "holdout. This is a more interesting result than expected and warrants a closer look before concluding "
            "the signal is purely a confound."
        )
    elif not separates:
        lines.append(
            "- The generic embedding does not clearly separate HER2-low from HER2-zero above its shuffled-label null. "
            "That is itself informative: it weakens the idea that simple generic morphology drives the GigaTIME "
            "signal, and keeps the GigaTIME-specific result worth a second look."
        )
    else:
        lines.append(
            "- The pattern is mixed; see the table above and treat this as a control to discuss, not a settled verdict."
        )
    return lines


def asset_link(asset_dir: Path, filename: str) -> str:
    return str(asset_dir / filename).replace("docs/", "")


def write_markdown(path: Path, asset_dir: Path, args, rows, metrics, pca_robustness, permutation) -> None:
    n_low = int((rows["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    n_zero = int((rows["clinical_her2_group"] == POSITIVE_CLASS).sum())
    pca_rows = [
        [str(int(r["pca_components"])), fmt(r["balanced_accuracy"], 3), fmt(r["macro_auc_ovr"], 3)]
        for _, r in pca_robustness.iterrows()
    ]
    lines = [
        "# Generic H&E Embedding Control (H-Optimus-0)",
        "",
        "Status: control experiment for the HER2-low versus HER2-zero result.",
        "",
        "## Why This Control",
        "",
        "The primary result uses GigaTIME virtual mIF channels and is presented as immune/myeloid/checkpoint biology. "
        "This control asks whether that interpretation is necessary. H-Optimus-0 is a generic pathology "
        "foundation-model embedding with no immune-channel meaning. If a generic embedding separates HER2-low from "
        "HER2-zero about as well as GigaTIME, and if it also collapses under source-site holdout while slide-size "
        "covariates stay strong, then the low-versus-zero axis is better described as generic morphology/tissue "
        "composition tracking TCGA acquisition structure than as GigaTIME-specific virtual immune biology.",
        "",
        "## Method",
        "",
        f"- Cohort: {len(rows)} strict high-trust slides with H-Optimus-0 embeddings ({n_low} HER2-low, {n_zero} HER2-zero).",
        "- Embedding: `bioptimus/H-optimus-0`, 1536-d, mean-pooled over 128 random tissue tiles per slide (same slide list and 128-tile sampling as the GigaTIME primary run).",
        f"- Classifier: regularized logistic regression on identical folds as the GigaTIME analyses; PCA ({args.pca_components} components) fit inside each training fold only.",
        f"- Comparators: slide-size covariates, source-site covariates, and GigaTIME mean channels (`{args.reference_view}` whole-slide view).",
        f"- Validation: repeated stratified {args.folds}-fold CV ({args.repeats} repeats) and leave-one-source-site-out.",
        f"- Sanity: {args.permutations} shuffled-label permutations for the embedding under repeated CV.",
        "",
        "## Head-to-Head Results",
        "",
        markdown_table(
            ["Feature set", "Validation", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            metric_table_rows(metrics),
        ),
        "",
        f"![Embedding vs GigaTIME and confound baselines]({asset_link(asset_dir, 'embedding_control_balanced_accuracy.png')})",
        "",
        "## Embedding PCA-Component Robustness (Repeated CV)",
        "",
        markdown_table(["PCA components", "Balanced accuracy", "AUC"], pca_rows),
        "",
        "## Shuffled-Label Sanity (Embedding, Repeated CV)",
        "",
        markdown_table(
            ["Observed bal acc", "Null mean", "Null 95%", "Empirical p"],
            [[
                fmt(permutation["observed_repeated_cv_balanced_accuracy"], 3),
                fmt(permutation["null_balanced_accuracy_mean"], 3),
                fmt(permutation["null_balanced_accuracy_p95"], 3),
                fmt(permutation["empirical_p_balanced_accuracy"], 4),
            ]],
        ),
        "",
        "## Interpretation",
        "",
        *build_verdict(metrics, permutation),
        "",
        "- This is an internal control on TCGA, not external validation. It cannot prove the GigaTIME signal is "
        "biologically meaningless; it tests whether the GigaTIME-specific virtual-immune framing is required to "
        "reproduce the low-versus-zero separation.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'embedding_control_metrics.csv'}`",
        f"- `{Path(args.out_dir) / 'embedding_control_pca_robustness.csv'}`",
        f"- `{Path(args.out_dir) / 'embedding_control_permutation.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_analysis_libs(out_dir / ".matplotlib")

    rows, metrics, pca_robustness, permutation, null_df = run_analysis(pd, optimize, stats, args)

    metrics.to_csv(out_dir / "embedding_control_metrics.csv", index=False)
    pca_robustness.to_csv(out_dir / "embedding_control_pca_robustness.csv", index=False)
    pd.DataFrame([permutation]).to_csv(out_dir / "embedding_control_permutation.csv", index=False)
    null_df.to_csv(out_dir / "embedding_control_permutation_null.csv", index=False)
    (out_dir / "embedding_control_metadata.json").write_text(
        json.dumps(
            {
                "task": "her2_low_vs_zero",
                "embedding_model": "bioptimus/H-optimus-0",
                "reference_view": args.reference_view,
                "pca_components": args.pca_components,
                "folds": args.folds,
                "repeats": args.repeats,
                "permutations": args.permutations,
                "seed": args.seed,
                "l2_penalty": args.l2_penalty,
                "n_cases": int(len(rows)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    plot_metric_comparison(plt, sns, metrics, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, args, rows, metrics, pca_robustness, permutation)
    print(f"Wrote embedding control outputs to {out_dir}")
    print(f"Wrote embedding control figure to {asset_dir}")
    print(f"Wrote embedding control markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
