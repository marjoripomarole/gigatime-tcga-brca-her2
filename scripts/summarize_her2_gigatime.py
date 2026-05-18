#!/usr/bin/env python3
"""Join GigaTIME slide scores to ERBB2 expression and create advisor summaries."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

KEY_CHANNELS = ["CD3", "CD8", "CD4", "CD20", "CD68", "CD11c", "PD-1", "PD-L1", "CK", "Ki67"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca/slide_scores.csv")
    parser.add_argument("--expression", default="data/tcga_brca/erbb2_expression.csv")
    parser.add_argument("--case-groups", default=None, help="Optional CSV with case_submitter_id and her2_group columns.")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca/advisor_summary")
    parser.add_argument("--channels", default=",".join(KEY_CHANNELS), help="Comma-separated GigaTIME channels to summarize.")
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Create/activate the environment with "
            "`conda env create -f envs/gigatime-tcga.yml && conda activate gigatime-tcga`."
        ) from exc
    return pd, plt, sns


def write_markdown(path, joined, group_summary, channels: list[str], grouping_label: str) -> None:
    n_cases = joined["case_submitter_id"].nunique()
    n_slides = len(joined)
    top = group_summary.sort_values("absolute_delta_mean_activation", ascending=False).head(5)
    if top.empty:
        top_table = "No overlapping requested channels were found in the slide score table."
    else:
        columns = [
            "channel",
            "delta_high_minus_low",
            "cohens_d",
            "welch_p_value",
            "spearman_erbb2_activation",
            "n_high_slides",
            "n_low_slides",
        ]
        table_lines = [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join(["---"] * len(columns)) + " |",
        ]
        for _, row in top[columns].iterrows():
            table_lines.append(
                "| "
                + " | ".join(
                    [
                        str(row["channel"]),
                        f"{row['delta_high_minus_low']:.4g}",
                        f"{row['cohens_d']:.4g}",
                        f"{row['welch_p_value']:.4g}",
                        f"{row['spearman_erbb2_activation']:.4g}",
                        str(int(row["n_high_slides"])),
                        str(int(row["n_low_slides"])),
                    ]
                )
                + " |"
            )
        top_table = "\n".join(table_lines)
    lines = [
        "# TCGA-BRCA GigaTIME HER2 Summary",
        "",
        f"- Joined slides: {n_slides}",
        f"- Unique TCGA cases: {n_cases}",
        f"- HER2 grouping: {grouping_label}",
        f"- Channels summarized: {', '.join(channels)}",
        "",
        "## Largest HER2-High Versus HER2-Low Virtual mIF Differences",
        "",
        top_table,
        "",
        "## Interpretation Guardrails",
        "",
        "- ERBB2 RNA expression is a proxy for HER2 biology, not clinical IHC/FISH status.",
        "- GigaTIME outputs are virtual mIF research features and require validation before biological claims.",
        "- Review tissue QC and slide sampling before scaling beyond the pilot subset.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_erbb2_distribution(joined, out_dir: Path, plt, sns) -> None:
    plt.figure(figsize=(6, 4))
    sns.histplot(joined.drop_duplicates("case_submitter_id")["erbb2_tpm"], bins=24)
    plt.xlabel("ERBB2 TPM")
    plt.ylabel("Cases")
    plt.title("TCGA-BRCA ERBB2 Expression")
    plt.tight_layout()
    plt.savefig(out_dir / "erbb2_tpm_distribution.png", dpi=180)
    plt.close()


def plot_group_deltas(group_summary, out_dir: Path, plt, sns) -> None:
    plot_df = group_summary.sort_values("delta_high_minus_low")
    plt.figure(figsize=(7, max(4, 0.35 * len(plot_df))))
    sns.barplot(data=plot_df, y="channel", x="delta_high_minus_low", color="#4C78A8")
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Mean activation difference: HER2-high minus HER2-low")
    plt.ylabel("GigaTIME channel")
    plt.title("Virtual mIF Differences by ERBB2 Group")
    plt.tight_layout()
    plt.savefig(out_dir / "her2_group_channel_deltas.png", dpi=180)
    plt.close()


def plot_group_boxplots(joined, channels: list[str], out_dir: Path, plt, sns) -> None:
    available = [channel for channel in channels if f"mean_{channel}" in joined.columns]
    if not available:
        return
    records = []
    for channel in available:
        score_col = f"mean_{channel}"
        for _, row in joined[["case_submitter_id", "her2_group", score_col]].iterrows():
            records.append(
                {
                    "case_submitter_id": row["case_submitter_id"],
                    "her2_group": row["her2_group"],
                    "channel": channel,
                    "mean_activation": row[score_col],
                }
            )
    import pandas as pd

    plot_df = pd.DataFrame(records)
    plt.figure(figsize=(max(8, 0.8 * len(available)), 5))
    sns.boxplot(data=plot_df, x="channel", y="mean_activation", hue="her2_group")
    groups = sorted(plot_df["her2_group"].dropna().unique())
    sns.stripplot(
        data=plot_df,
        x="channel",
        y="mean_activation",
        hue="her2_group",
        dodge=True,
        palette={group: "black" for group in groups},
        size=3,
        alpha=0.5,
        legend=False,
    )
    plt.xlabel("GigaTIME channel")
    plt.ylabel("Mean virtual activation")
    plt.title("Virtual mIF Distributions by HER2 Group")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / "her2_group_channel_boxplots.png", dpi=180)
    plt.close()


def plot_scatter_panels(joined, channels: list[str], out_dir: Path, plt, sns) -> None:
    available = [channel for channel in channels if f"mean_{channel}" in joined.columns]
    if not available:
        return
    ncols = 2
    nrows = (len(available) + 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, max(4, 3.3 * nrows)))
    axes_list = list(axes.flatten()) if hasattr(axes, "flatten") else [axes]
    for axis, channel in zip(axes_list, available):
        sns.scatterplot(data=joined, x="erbb2_tpm", y=f"mean_{channel}", hue="her2_group", ax=axis, s=28)
        axis.set_title(channel)
        axis.set_xlabel("ERBB2 TPM")
        axis.set_ylabel("Mean virtual activation")
        axis.legend_.remove() if axis.legend_ else None
    for axis in axes_list[len(available) :]:
        axis.axis("off")
    handles, labels = axes_list[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
    fig.tight_layout(rect=(0, 0, 0.95, 1))
    fig.savefig(out_dir / "erbb2_vs_virtual_mif_scatter.png", dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns = require_analysis_libs(out_dir / ".matplotlib")

    slide_scores = pd.read_csv(args.slide_scores)
    expression = pd.read_csv(args.expression)
    expression = expression.rename(columns={"tpm_unstranded": "erbb2_tpm"})
    expression["erbb2_tpm"] = pd.to_numeric(expression["erbb2_tpm"], errors="coerce")
    expression = expression.dropna(subset=["case_submitter_id", "erbb2_tpm"])
    expression = expression.sort_values("sample_type").drop_duplicates("case_submitter_id")

    joined = slide_scores.merge(expression[["case_submitter_id", "sample_submitter_id", "sample_type", "erbb2_tpm"]], on="case_submitter_id", how="inner")
    if joined.empty:
        raise ValueError("No overlapping case_submitter_id values between slide scores and ERBB2 expression.")

    if args.case_groups:
        case_groups = pd.read_csv(args.case_groups)
        required = {"case_submitter_id", "her2_group"}
        missing = required - set(case_groups.columns)
        if missing:
            raise ValueError(f"--case-groups is missing required columns: {', '.join(sorted(missing))}")
        joined = joined.merge(
            case_groups[["case_submitter_id", "her2_group"]].drop_duplicates("case_submitter_id"),
            on="case_submitter_id",
            how="inner",
        )
        grouping_label = f"explicit labels from {args.case_groups}"
    else:
        median_tpm = joined.drop_duplicates("case_submitter_id")["erbb2_tpm"].median()
        joined["her2_group"] = joined["erbb2_tpm"].map(lambda value: "HER2-high" if value >= median_tpm else "HER2-low")
        grouping_label = f"ERBB2 TPM median split at {median_tpm:.4g}"
    if joined.empty:
        raise ValueError("No slide scores remained after applying HER2 case groups.")
    joined.to_csv(out_dir / "joined_slide_her2_gigatime.csv", index=False)

    channels = [channel.strip() for channel in args.channels.split(",") if channel.strip()]
    rows = []
    for channel in channels:
        score_col = f"mean_{channel}"
        if score_col not in joined.columns:
            continue
        grouped = joined.groupby("her2_group")[score_col].agg(["mean", "median", "std", "count"])
        high = grouped.loc["HER2-high"] if "HER2-high" in grouped.index else None
        low = grouped.loc["HER2-low"] if "HER2-low" in grouped.index else None
        if high is None or low is None:
            continue
        high_values = joined.loc[joined["her2_group"] == "HER2-high", score_col].dropna()
        low_values = joined.loc[joined["her2_group"] == "HER2-low", score_col].dropna()
        corr = joined[["erbb2_tpm", score_col]].corr(method="spearman").iloc[0, 1]
        delta = float(high["mean"] - low["mean"])
        if len(high_values) > 1 and len(low_values) > 1:
            from scipy import stats

            welch_p = float(stats.ttest_ind(high_values, low_values, equal_var=False, nan_policy="omit").pvalue)
            mannwhitney_p = float(stats.mannwhitneyu(high_values, low_values, alternative="two-sided").pvalue)
        else:
            welch_p = float("nan")
            mannwhitney_p = float("nan")
        pooled_sd = ((high_values.var(ddof=1) + low_values.var(ddof=1)) / 2) ** 0.5 if len(high_values) > 1 and len(low_values) > 1 else float("nan")
        cohens_d = float(delta / pooled_sd) if pooled_sd and pooled_sd == pooled_sd else float("nan")
        rows.append(
            {
                "channel": channel,
                "her2_high_mean_activation": float(high["mean"]),
                "her2_high_median_activation": float(high["median"]),
                "her2_high_sd_activation": float(high["std"]) if high["std"] == high["std"] else float("nan"),
                "her2_low_mean_activation": float(low["mean"]),
                "her2_low_median_activation": float(low["median"]),
                "her2_low_sd_activation": float(low["std"]) if low["std"] == low["std"] else float("nan"),
                "delta_high_minus_low": delta,
                "absolute_delta_mean_activation": abs(delta),
                "cohens_d": cohens_d,
                "welch_p_value": welch_p,
                "mannwhitney_p_value": mannwhitney_p,
                "spearman_erbb2_activation": float(corr),
                "n_high_slides": int(high["count"]),
                "n_low_slides": int(low["count"]),
            }
        )
    group_summary = pd.DataFrame(rows)
    group_summary.to_csv(out_dir / "her2_group_channel_summary.csv", index=False)

    plot_erbb2_distribution(joined, out_dir, plt, sns)
    if not group_summary.empty:
        plot_group_deltas(group_summary, out_dir, plt, sns)
        plot_group_boxplots(joined, channels, out_dir, plt, sns)
    plot_scatter_panels(joined, channels, out_dir, plt, sns)
    write_markdown(out_dir / "advisor_summary.md", joined, group_summary, channels, grouping_label)
    print(f"Wrote advisor summary to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
