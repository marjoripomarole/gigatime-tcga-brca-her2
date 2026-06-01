#!/usr/bin/env python3
"""Build cleaned GigaTIME tile feature views before classifier training."""

from __future__ import annotations

import argparse
import math
import os
from itertools import combinations
from pathlib import Path


GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
KEY_CHANNELS = ["CD68", "PD-L1", "CD11c", "CD3", "CD8", "CD4", "CD20", "CK", "Ki67"]
PLOT_CHANNELS = ["CD68", "PD-L1", "CD11c", "CK", "Ki67"]
FILTER_LABELS = {
    "all_sampled_tissue": "All sampled tissue",
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_enriched_top50": "CK-enriched top 50%",
    "ck_enriched_top25": "CK-enriched top 25%",
}
FILTER_ORDER = list(FILTER_LABELS)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tile-scores",
        default="results/gigatime_tcga_brca_clinical_her2_tile256/tile_scores.csv",
        help="Tile-level GigaTIME scores from the 256-tile clinical HER2 run.",
    )
    parser.add_argument(
        "--cohort",
        default="data/tcga_brca/clinical_her2_cohort_cases.csv",
        help="Clinical HER2 cohort table.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup",
        help="Directory for cleaned feature tables.",
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_gigatime_cleanup",
        help="Directory for tracked cleanup figures.",
    )
    parser.add_argument("--min-tissue-fraction", type=float, default=0.70)
    parser.add_argument("--min-dapi", type=float, default=0.05)
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd, plt, sns, stats


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    if not indexed:
        return [float("nan")] * len(p_values)
    ranked = sorted(indexed, key=lambda item: item[1])
    m = len(ranked)
    adjusted = [float("nan")] * len(p_values)
    prev = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(ranked, start=1))):
        q_value = min(prev, p_value * m / rank)
        adjusted[idx] = q_value
        prev = q_value
    return adjusted


def safe_mean(values) -> float:
    values = list(values)
    return float(sum(values) / len(values)) if values else float("nan")


def cliffs_delta(values_a, values_b) -> float:
    a = list(values_a)
    b = list(values_b)
    if not a or not b:
        return float("nan")
    greater = 0
    less = 0
    for value_a in a:
        for value_b in b:
            if value_a > value_b:
                greater += 1
            elif value_a < value_b:
                less += 1
    return (greater - less) / (len(a) * len(b))


def load_inputs(pd, tile_scores: Path, cohort_path: Path):
    tiles = pd.read_csv(tile_scores)
    cohort = pd.read_csv(cohort_path)
    keep_cols = [
        "case_submitter_id",
        "clinical_her2_group",
        "clinical_her2_group_rule",
        "clinical_her2_group_confidence",
        "her2_ihc_score",
        "her2_ish_status",
        "erbb2_tpm",
        "er_status",
        "pr_status",
    ]
    cohort = cohort[[col for col in keep_cols if col in cohort.columns]].drop_duplicates("case_submitter_id")
    tiles = tiles.merge(cohort, on="case_submitter_id", how="left", validate="many_to_one")
    return tiles, cohort


def add_cleanup_flags(tiles, min_tissue_fraction: float, min_dapi: float):
    tiles = tiles.copy()
    tiles["qc_cellular_tissue"] = (tiles["tissue_fraction"] >= min_tissue_fraction) & (tiles["mean_DAPI"] >= min_dapi)
    tiles["ck_percentile_within_qc_slide"] = float("nan")
    qc = tiles.loc[tiles["qc_cellular_tissue"]].copy()
    if not qc.empty:
        percentiles = qc.groupby("slide_id")["mean_CK"].rank(method="average", pct=True)
        tiles.loc[qc.index, "ck_percentile_within_qc_slide"] = percentiles
    tiles["ck_enriched_top50"] = tiles["qc_cellular_tissue"] & (tiles["ck_percentile_within_qc_slide"] >= 0.50)
    tiles["ck_enriched_top25"] = tiles["qc_cellular_tissue"] & (tiles["ck_percentile_within_qc_slide"] >= 0.75)
    tiles["all_sampled_tissue"] = True
    tiles["virtual_myeloid_checkpoint_score"] = tiles[["mean_CD68", "mean_CD11c", "mean_PD-L1"]].mean(axis=1)
    tiles["virtual_t_cell_score"] = tiles[["mean_CD3", "mean_CD4", "mean_CD8"]].mean(axis=1)
    return tiles


def top_fraction_mean(values, fraction: float = 0.10) -> float:
    if len(values) == 0:
        return float("nan")
    ordered = sorted([float(value) for value in values if not math.isnan(float(value))], reverse=True)
    if not ordered:
        return float("nan")
    n_top = max(1, math.ceil(len(ordered) * fraction))
    return safe_mean(ordered[:n_top])


def aggregate_view(pd, tiles, view: str):
    rows = []
    mean_cols = [col for col in tiles.columns if col.startswith("mean_")]
    frac_cols = [col for col in tiles.columns if col.startswith("frac_")]
    total_counts = tiles.groupby("slide_id").size().to_dict()
    selected = tiles.loc[tiles[view]].copy()
    for slide_id, slide_tiles in selected.groupby("slide_id", sort=False):
        first = slide_tiles.iloc[0]
        row = {
            "feature_view": view,
            "feature_view_label": FILTER_LABELS[view],
            "slide_id": slide_id,
            "case_submitter_id": first["case_submitter_id"],
            "clinical_her2_group": first["clinical_her2_group"],
            "n_tiles_total": int(total_counts.get(slide_id, 0)),
            "n_tiles_retained": int(len(slide_tiles)),
            "retained_fraction": float(len(slide_tiles) / total_counts.get(slide_id, len(slide_tiles))),
            "mean_tissue_fraction": float(slide_tiles["tissue_fraction"].mean()),
            "median_tissue_fraction": float(slide_tiles["tissue_fraction"].median()),
            "mean_ck_percentile_within_qc_slide": float(slide_tiles["ck_percentile_within_qc_slide"].mean()),
        }
        for optional_col in [
            "clinical_her2_group_rule",
            "clinical_her2_group_confidence",
            "her2_ihc_score",
            "her2_ish_status",
            "erbb2_tpm",
            "er_status",
            "pr_status",
        ]:
            if optional_col in slide_tiles.columns:
                row[optional_col] = first[optional_col]
        for col in mean_cols:
            channel = col.replace("mean_", "", 1)
            values = slide_tiles[col].dropna()
            row[col] = float(values.mean())
            row[f"median_{channel}"] = float(values.median())
            row[f"p90_{channel}"] = float(values.quantile(0.90))
            row[f"std_{channel}"] = float(values.std(ddof=1)) if len(values) > 1 else float("nan")
            row[f"top10_mean_{channel}"] = top_fraction_mean(values)
        for col in frac_cols:
            row[col] = float(slide_tiles[col].mean())
        row["virtual_myeloid_checkpoint_score"] = float(slide_tiles["virtual_myeloid_checkpoint_score"].mean())
        row["virtual_t_cell_score"] = float(slide_tiles["virtual_t_cell_score"].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_all_views(pd, tiles):
    return pd.concat([aggregate_view(pd, tiles, view) for view in FILTER_ORDER], ignore_index=True)


def build_retention_summary(slide_features):
    return slide_features[
        [
            "feature_view",
            "feature_view_label",
            "slide_id",
            "case_submitter_id",
            "clinical_her2_group",
            "n_tiles_total",
            "n_tiles_retained",
            "retained_fraction",
            "mean_tissue_fraction",
            "mean_DAPI",
            "mean_CK",
        ]
    ].copy()


def build_channel_summary(stats, slide_features, channels: list[str]):
    rows = []
    for view in FILTER_ORDER:
        view_df = slide_features.loc[slide_features["feature_view"] == view]
        for channel in channels:
            col = f"mean_{channel}"
            if col not in view_df.columns:
                continue
            by_group = {
                group: view_df.loc[view_df["clinical_her2_group"] == group, col].dropna()
                for group in GROUP_ORDER
            }
            nonempty = [values for values in by_group.values() if len(values)]
            if len(nonempty) >= 2:
                test = stats.kruskal(*nonempty)
                p_value = float(test.pvalue)
                statistic = float(test.statistic)
            else:
                p_value = float("nan")
                statistic = float("nan")
            means = {group: float(values.mean()) if len(values) else float("nan") for group, values in by_group.items()}
            available = {group: value for group, value in means.items() if not math.isnan(value)}
            highest = max(available, key=available.get) if available else ""
            lowest = min(available, key=available.get) if available else ""
            row = {
                "feature_view": view,
                "feature_view_label": FILTER_LABELS[view],
                "channel": channel,
                "kruskal_h": statistic,
                "kruskal_p_value": p_value,
                "highest_mean_group": highest,
                "lowest_mean_group": lowest,
                "max_minus_min_mean": max(available.values()) - min(available.values()) if available else float("nan"),
            }
            for group in GROUP_ORDER:
                prefix = group.lower().replace("-", "_")
                values = by_group[group]
                row[f"{prefix}_n"] = int(len(values))
                row[f"{prefix}_mean"] = float(values.mean()) if len(values) else float("nan")
                row[f"{prefix}_median"] = float(values.median()) if len(values) else float("nan")
            rows.append(row)
    for view in FILTER_ORDER:
        mask = [row["feature_view"] == view for row in rows]
        p_values = [row["kruskal_p_value"] for row, include in zip(rows, mask) if include]
        q_values = benjamini_hochberg(p_values)
        q_iter = iter(q_values)
        for row, include in zip(rows, mask):
            if include:
                row["kruskal_q_value_bh_within_view"] = next(q_iter)
    return rows


def build_pairwise_tests(stats, slide_features, channels: list[str]):
    rows = []
    for view in FILTER_ORDER:
        view_df = slide_features.loc[slide_features["feature_view"] == view]
        for channel in channels:
            col = f"mean_{channel}"
            if col not in view_df.columns:
                continue
            for group_a, group_b in combinations(GROUP_ORDER, 2):
                values_a = view_df.loc[view_df["clinical_her2_group"] == group_a, col].dropna()
                values_b = view_df.loc[view_df["clinical_her2_group"] == group_b, col].dropna()
                if len(values_a) and len(values_b):
                    test = stats.mannwhitneyu(values_a, values_b, alternative="two-sided")
                    p_value = float(test.pvalue)
                    statistic = float(test.statistic)
                    delta = float(values_a.mean() - values_b.mean())
                    cliff = cliffs_delta(values_a, values_b)
                else:
                    p_value = float("nan")
                    statistic = float("nan")
                    delta = float("nan")
                    cliff = float("nan")
                rows.append(
                    {
                        "feature_view": view,
                        "feature_view_label": FILTER_LABELS[view],
                        "channel": channel,
                        "group_a": group_a,
                        "group_b": group_b,
                        "n_a": int(len(values_a)),
                        "n_b": int(len(values_b)),
                        "mean_a": float(values_a.mean()) if len(values_a) else float("nan"),
                        "mean_b": float(values_b.mean()) if len(values_b) else float("nan"),
                        "delta_mean_a_minus_b": delta,
                        "mannwhitney_u": statistic,
                        "mannwhitney_p_value": p_value,
                        "cliffs_delta": cliff,
                    }
                )
    for view in FILTER_ORDER:
        mask = [row["feature_view"] == view for row in rows]
        p_values = [row["mannwhitney_p_value"] for row, include in zip(rows, mask) if include]
        q_values = benjamini_hochberg(p_values)
        q_iter = iter(q_values)
        for row, include in zip(rows, mask):
            if include:
                row["mannwhitney_q_value_bh_within_view"] = next(q_iter)
    return rows


def plot_retention(plt, sns, retention, asset_dir: Path) -> None:
    plt.figure(figsize=(9.5, 5.0))
    sns.boxplot(
        data=retention,
        x="feature_view_label",
        y="n_tiles_retained",
        hue="clinical_her2_group",
        hue_order=GROUP_ORDER,
        order=[FILTER_LABELS[view] for view in FILTER_ORDER],
    )
    sns.stripplot(
        data=retention,
        x="feature_view_label",
        y="n_tiles_retained",
        hue="clinical_her2_group",
        hue_order=GROUP_ORDER,
        order=[FILTER_LABELS[view] for view in FILTER_ORDER],
        dodge=True,
        palette={group: "black" for group in GROUP_ORDER},
        size=2.8,
        alpha=0.45,
        legend=False,
    )
    plt.xlabel("Cleanup view")
    plt.ylabel("Retained tiles per slide")
    plt.title("Tile Retention After GigaTIME Cleanup Views")
    plt.xticks(rotation=20, ha="right")
    legend = plt.gca().get_legend()
    if legend:
        legend.set_title("Clinical HER2 group")
    plt.tight_layout()
    plt.savefig(asset_dir / "cleanup_retained_tiles_by_filter.png", dpi=180)
    plt.close()


def plot_ck_dapi_distribution(plt, sns, tiles, asset_dir: Path) -> None:
    plot_df = tiles.copy()
    plot_df["cleanup_category"] = "Excluded"
    plot_df.loc[plot_df["qc_cellular_tissue"], "cleanup_category"] = "QC cellular tissue"
    plot_df.loc[plot_df["ck_enriched_top50"], "cleanup_category"] = "CK top 50%"
    plot_df.loc[plot_df["ck_enriched_top25"], "cleanup_category"] = "CK top 25%"
    plt.figure(figsize=(7.2, 5.5))
    sns.scatterplot(
        data=plot_df,
        x="mean_DAPI",
        y="mean_CK",
        hue="cleanup_category",
        hue_order=["Excluded", "QC cellular tissue", "CK top 50%", "CK top 25%"],
        s=12,
        alpha=0.38,
        linewidth=0,
    )
    plt.xlabel("Virtual DAPI mean")
    plt.ylabel("Virtual CK mean")
    plt.title("Tile-Level Cellularity and CK Enrichment")
    legend = plt.gca().get_legend()
    if legend:
        legend.set_title("Cleanup category")
    plt.tight_layout()
    plt.savefig(asset_dir / "cleanup_ck_dapi_distribution.png", dpi=180)
    plt.close()


def plot_key_channel_heatmap(pd, plt, sns, channel_summary, asset_dir: Path) -> None:
    rows = []
    row_order = []
    for view in FILTER_ORDER:
        for channel in ["CD68", "PD-L1", "CD11c", "CK", "Ki67"]:
            matches = [row for row in channel_summary if row["feature_view"] == view and row["channel"] == channel]
            if not matches:
                continue
            row = matches[0]
            row_label = f"{row['feature_view_label']} | {channel}"
            row_order.append(row_label)
            for group in GROUP_ORDER:
                prefix = group.lower().replace("-", "_")
                rows.append(
                    {
                        "cleanup_channel": row_label,
                        "clinical_her2_group": group,
                        "mean_activation": row[f"{prefix}_mean"],
                    }
                )
    plot_df = pd.DataFrame(rows)
    heatmap = plot_df.pivot(index="cleanup_channel", columns="clinical_her2_group", values="mean_activation")
    heatmap = heatmap.reindex(row_order)
    heatmap = heatmap[GROUP_ORDER]
    plt.figure(figsize=(7.2, max(6.0, 0.36 * len(heatmap))))
    sns.heatmap(heatmap, cmap="viridis", annot=True, fmt=".3f", linewidths=0.2)
    plt.xlabel("Clinical HER2 group")
    plt.ylabel("Cleanup view and channel")
    plt.title("Cleaned GigaTIME Group Means")
    plt.tight_layout()
    plt.savefig(asset_dir / "cleanup_key_channel_heatmap.png", dpi=180)
    plt.close()


def plot_key_channel_boxplots(pd, plt, sns, slide_features, asset_dir: Path) -> None:
    records = []
    for _, row in slide_features.iterrows():
        for channel in PLOT_CHANNELS:
            records.append(
                {
                    "feature_view_label": row["feature_view_label"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "channel": channel,
                    "mean_activation": row[f"mean_{channel}"],
                }
            )
    plot_df = pd.DataFrame(records)
    plot_df["feature_view_label"] = pd.Categorical(
        plot_df["feature_view_label"],
        categories=[FILTER_LABELS[view] for view in FILTER_ORDER],
        ordered=True,
    )
    grid = sns.catplot(
        data=plot_df,
        x="channel",
        y="mean_activation",
        hue="clinical_her2_group",
        col="feature_view_label",
        col_wrap=2,
        order=PLOT_CHANNELS,
        hue_order=GROUP_ORDER,
        kind="box",
        sharey=False,
        height=3.8,
        aspect=1.25,
        legend=True,
    )
    grid.set_axis_labels("GigaTIME channel", "Mean activation")
    grid.set_titles("{col_name}")
    for axis in grid.axes.flat:
        axis.tick_params(axis="x", rotation=30)
    if grid.legend:
        grid.legend.set_title("Clinical HER2 group")
    grid.fig.subplots_adjust(top=0.90)
    grid.fig.suptitle("Key GigaTIME Channels After Cleanup")
    grid.savefig(asset_dir / "cleanup_key_channel_boxplots.png", dpi=180)
    plt.close(grid.fig)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    return f"{value:.{digits}f}"


def write_summary(path: Path, retention, channel_summary, pairwise_tests, args: argparse.Namespace) -> None:
    retention_rows = []
    for view in FILTER_ORDER:
        subset = retention.loc[retention["feature_view"] == view]
        retention_rows.append(
            [
                FILTER_LABELS[view],
                fmt(float(subset["n_tiles_retained"].median()), 1),
                fmt(float(subset["retained_fraction"].median()), 3),
                fmt(float(subset["mean_DAPI"].median()), 3),
                fmt(float(subset["mean_CK"].median()), 3),
            ]
        )

    top_channel_rows = []
    for view in FILTER_ORDER:
        view_rows = [row for row in channel_summary if row["feature_view"] == view]
        for row in sorted(view_rows, key=lambda item: item["kruskal_p_value"])[:5]:
            top_channel_rows.append(
                [
                    row["feature_view_label"],
                    row["channel"],
                    fmt(row["kruskal_p_value"], 4),
                    fmt(row["kruskal_q_value_bh_within_view"], 4),
                    row["highest_mean_group"],
                    row["lowest_mean_group"],
                    fmt(row["max_minus_min_mean"], 4),
                ]
            )

    focus_rows = []
    for view in FILTER_ORDER:
        for channel in ["CD68", "PD-L1", "CD11c"]:
            matches = [
                row
                for row in pairwise_tests
                if row["feature_view"] == view
                and row["channel"] == channel
                and row["group_a"] == "HER2-low"
                and row["group_b"] == "HER2-zero"
            ]
            if not matches:
                continue
            row = matches[0]
            focus_rows.append(
                [
                    FILTER_LABELS[view],
                    channel,
                    fmt(row["delta_mean_a_minus_b"], 4),
                    fmt(row["mannwhitney_p_value"], 4),
                    fmt(row["mannwhitney_q_value_bh_within_view"], 4),
                ]
            )

    lines = [
        "# Clinical HER2 GigaTIME Data Cleanup",
        "",
        "This cleanup step goes back before classifier training. It asks whether the GigaTIME input features should be aggregated from all sampled tissue tiles, or from cleaner tile subsets that are more cellular and more epithelial/tumor enriched.",
        "",
        "## Cleanup Views",
        "",
        f"- `all_sampled_tissue`: all tissue tiles selected by the original GigaTIME run.",
        f"- `qc_cellular_tissue`: tiles with H&E tissue fraction >= {args.min_tissue_fraction:.2f} and virtual DAPI mean >= {args.min_dapi:.2f}.",
        "- `ck_enriched_top50`: the top 50% virtual CK tiles within each slide after cellular-tissue QC.",
        "- `ck_enriched_top25`: the top 25% virtual CK tiles within each slide after cellular-tissue QC.",
        "",
        "Important: CK and DAPI are still GigaTIME virtual predictions from H&E, not laboratory stains. These filters create tumor-enriched research feature views, not confirmed tumor masks.",
        "",
        "## Tile Retention",
        "",
        markdown_table(
            ["Cleanup view", "Median retained tiles", "Median retained fraction", "Median DAPI", "Median CK"],
            retention_rows,
        ),
        "",
        "![Tile retention by cleanup view](assets/clinical_her2_gigatime_cleanup/cleanup_retained_tiles_by_filter.png)",
        "",
        "## Tile-Level Cleanup Map",
        "",
        "This plot shows how the cleanup rules move from broad tissue tiles toward cellular, CK-enriched tiles.",
        "",
        "![CK and DAPI tile distribution](assets/clinical_her2_gigatime_cleanup/cleanup_ck_dapi_distribution.png)",
        "",
        "## Group Means After Cleanup",
        "",
        "![Cleaned GigaTIME group mean heatmap](assets/clinical_her2_gigatime_cleanup/cleanup_key_channel_heatmap.png)",
        "",
        "## Top Three-Group Signals Across Cleanup Views",
        "",
        markdown_table(
            ["Cleanup view", "Channel", "Kruskal p", "BH q within view", "Highest group", "Lowest group", "Max-min mean"],
            top_channel_rows,
        ),
        "",
        "## HER2-Low Versus HER2-Zero Focus",
        "",
        "Negative delta means HER2-low is lower than HER2-zero.",
        "",
        markdown_table(
            ["Cleanup view", "Channel", "HER2-low minus HER2-zero", "Mann-Whitney p", "BH q within view"],
            focus_rows,
        ),
        "",
        "![Key channel boxplots after cleanup](assets/clinical_her2_gigatime_cleanup/cleanup_key_channel_boxplots.png)",
        "",
        "## Interpretation",
        "",
        "This cleanup does not validate the virtual markers, but it makes the next classifier input more biologically defensible. The original baseline averaged all sampled tissue tiles. The cleaned tables let us rerun summaries or classifiers using cellular tissue tiles and CK-enriched tumor-context tiles.",
        "",
        "If HER2-low versus HER2-zero signal remains after CK enrichment, it is less likely to be explained only by blank tissue or broad non-cellular sampling. If the signal disappears, then the original classifier may have been leaning on non-tumor tissue composition rather than tumor-region biology.",
        "",
        "## Outputs",
        "",
        f"- `{args.out_dir}/tile_qc_scores.csv`",
        f"- `{args.out_dir}/cleaned_slide_features.csv`",
        f"- `{args.out_dir}/filter_retention_summary.csv`",
        f"- `{args.out_dir}/cleanup_channel_summary.csv`",
        f"- `{args.out_dir}/cleanup_pairwise_tests.csv`",
        "",
        "## Next Step",
        "",
        "Use `cleaned_slide_features.csv` to rerun the classifier separately for each cleanup view, especially `qc_cellular_tissue`, `ck_enriched_top50`, and `ck_enriched_top25`. The comparison should show whether tumor-enriched GigaTIME features improve HER2 prediction or expose the current model as tissue-composition driven.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")
    tiles, _cohort = load_inputs(pd, Path(args.tile_scores), Path(args.cohort))
    tiles = add_cleanup_flags(tiles, args.min_tissue_fraction, args.min_dapi)
    slide_features = aggregate_all_views(pd, tiles)
    retention = build_retention_summary(slide_features)
    channel_summary = build_channel_summary(stats, slide_features, KEY_CHANNELS)
    pairwise_tests = build_pairwise_tests(stats, slide_features, KEY_CHANNELS)

    tiles.to_csv(out_dir / "tile_qc_scores.csv", index=False)
    slide_features.to_csv(out_dir / "cleaned_slide_features.csv", index=False)
    retention.to_csv(out_dir / "filter_retention_summary.csv", index=False)
    pd.DataFrame(channel_summary).to_csv(out_dir / "cleanup_channel_summary.csv", index=False)
    pd.DataFrame(pairwise_tests).to_csv(out_dir / "cleanup_pairwise_tests.csv", index=False)

    plot_retention(plt, sns, retention, asset_dir)
    plot_ck_dapi_distribution(plt, sns, tiles, asset_dir)
    plot_key_channel_heatmap(pd, plt, sns, channel_summary, asset_dir)
    plot_key_channel_boxplots(pd, plt, sns, slide_features, asset_dir)
    write_summary(out_dir / "cleanup_summary.md", retention, channel_summary, pairwise_tests, args)
    write_summary(Path("docs/clinical_her2_gigatime_data_cleanup.md"), retention, channel_summary, pairwise_tests, args)
    print(f"Wrote cleanup outputs to {out_dir}")
    print(f"Wrote cleanup figures to {asset_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
