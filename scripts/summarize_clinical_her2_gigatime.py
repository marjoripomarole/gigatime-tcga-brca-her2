#!/usr/bin/env python3
"""Summarize GigaTIME outputs across clinical HER2-positive/low/zero groups."""

from __future__ import annotations

import argparse
import math
import os
from itertools import combinations
from pathlib import Path

GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
KEY_CHANNELS = ["CD3", "CD8", "CD4", "CD20", "CD68", "CD11c", "PD-1", "PD-L1", "CK", "Ki67"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-scores", default="results/gigatime_tcga_brca_clinical_her2/slide_scores.csv")
    parser.add_argument("--cohort", default="data/tcga_brca/clinical_her2_cohort_cases.csv")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca_clinical_her2/clinical_summary")
    parser.add_argument("--channels", default=",".join(KEY_CHANNELS))
    parser.add_argument("--group-column", default="cohort_group")
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


def safe_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def build_channel_summary(joined, channels: list[str], stats):
    rows = []
    for channel in channels:
        score_col = f"mean_{channel}"
        if score_col not in joined.columns:
            continue
        values_by_group = {
            group: joined.loc[joined["clinical_her2_group"] == group, score_col].dropna()
            for group in GROUP_ORDER
        }
        nonempty = [values for values in values_by_group.values() if len(values) > 0]
        if len(nonempty) >= 2:
            kruskal = stats.kruskal(*nonempty)
            kruskal_h = float(kruskal.statistic)
            kruskal_p = float(kruskal.pvalue)
            n_total = sum(len(values) for values in nonempty)
            k_groups = len(nonempty)
            epsilon_squared = (
                max(0.0, (kruskal_h - k_groups + 1) / (n_total - k_groups)) if n_total > k_groups else float("nan")
            )
        else:
            kruskal_h = float("nan")
            kruskal_p = float("nan")
            epsilon_squared = float("nan")

        means = {group: safe_float(values.mean()) if len(values) else float("nan") for group, values in values_by_group.items()}
        available_means = {group: value for group, value in means.items() if not math.isnan(value)}
        highest_group = max(available_means, key=available_means.get) if available_means else ""
        lowest_group = min(available_means, key=available_means.get) if available_means else ""
        row = {
            "channel": channel,
            "kruskal_h": kruskal_h,
            "kruskal_p_value": kruskal_p,
            "epsilon_squared": epsilon_squared,
            "highest_mean_group": highest_group,
            "lowest_mean_group": lowest_group,
            "max_minus_min_mean": max(available_means.values()) - min(available_means.values()) if available_means else float("nan"),
        }
        for group in GROUP_ORDER:
            values = values_by_group[group]
            prefix = group.lower().replace("-", "_")
            row[f"{prefix}_n"] = int(len(values))
            row[f"{prefix}_mean"] = safe_float(values.mean()) if len(values) else float("nan")
            row[f"{prefix}_median"] = safe_float(values.median()) if len(values) else float("nan")
            row[f"{prefix}_sd"] = safe_float(values.std(ddof=1)) if len(values) > 1 else float("nan")
        rows.append(row)
    return rows


def build_pairwise_tests(joined, channels: list[str], stats):
    rows = []
    for channel in channels:
        score_col = f"mean_{channel}"
        if score_col not in joined.columns:
            continue
        for group_a, group_b in combinations(GROUP_ORDER, 2):
            values_a = joined.loc[joined["clinical_her2_group"] == group_a, score_col].dropna()
            values_b = joined.loc[joined["clinical_her2_group"] == group_b, score_col].dropna()
            if len(values_a) and len(values_b):
                test = stats.mannwhitneyu(values_a, values_b, alternative="two-sided")
                p_value = float(test.pvalue)
                statistic = float(test.statistic)
                delta_mean = safe_float(values_a.mean() - values_b.mean())
                cliff = cliffs_delta(values_a, values_b)
            else:
                p_value = float("nan")
                statistic = float("nan")
                delta_mean = float("nan")
                cliff = float("nan")
            rows.append(
                {
                    "channel": channel,
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_a": int(len(values_a)),
                    "n_b": int(len(values_b)),
                    "mean_a": safe_float(values_a.mean()) if len(values_a) else float("nan"),
                    "mean_b": safe_float(values_b.mean()) if len(values_b) else float("nan"),
                    "delta_mean_a_minus_b": delta_mean,
                    "mannwhitney_u": statistic,
                    "mannwhitney_p_value": p_value,
                    "cliffs_delta": cliff,
                }
            )
    q_values = benjamini_hochberg([row["mannwhitney_p_value"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["mannwhitney_q_value_bh"] = q_value
    return rows


def plot_group_boxplots(pd, plt, sns, joined, channels: list[str], out_dir: Path) -> None:
    records = []
    for channel in channels:
        score_col = f"mean_{channel}"
        if score_col not in joined.columns:
            continue
        for _, row in joined[["case_submitter_id", "clinical_her2_group", score_col]].iterrows():
            records.append(
                {
                    "case_submitter_id": row["case_submitter_id"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "channel": channel,
                    "mean_activation": row[score_col],
                }
            )
    if not records:
        return
    plot_df = pd.DataFrame(records)
    plt.figure(figsize=(max(9, 0.85 * len(channels)), 5.5))
    sns.boxplot(data=plot_df, x="channel", y="mean_activation", hue="clinical_her2_group", hue_order=GROUP_ORDER)
    sns.stripplot(
        data=plot_df,
        x="channel",
        y="mean_activation",
        hue="clinical_her2_group",
        hue_order=GROUP_ORDER,
        dodge=True,
        palette={group: "black" for group in GROUP_ORDER},
        size=3,
        alpha=0.45,
        legend=False,
    )
    plt.xlabel("GigaTIME channel")
    plt.ylabel("Mean virtual activation")
    plt.title("Virtual mIF Features by Clinical HER2 Group")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / "clinical_her2_channel_boxplots.png", dpi=180)
    plt.close()


def plot_group_mean_heatmap(pd, plt, sns, channel_summary, out_dir: Path) -> None:
    if channel_summary.empty:
        return
    records = []
    for _, row in channel_summary.iterrows():
        for group in GROUP_ORDER:
            prefix = group.lower().replace("-", "_")
            records.append({"channel": row["channel"], "clinical_her2_group": group, "mean_activation": row[f"{prefix}_mean"]})
    heatmap_df = pd.DataFrame(records).pivot(index="channel", columns="clinical_her2_group", values="mean_activation")
    heatmap_df = heatmap_df[GROUP_ORDER]
    plt.figure(figsize=(6.5, max(4, 0.45 * len(heatmap_df))))
    sns.heatmap(heatmap_df, cmap="viridis", annot=True, fmt=".3f", linewidths=0.3)
    plt.xlabel("Clinical HER2 group")
    plt.ylabel("GigaTIME channel")
    plt.title("Mean Virtual mIF Activation")
    plt.tight_layout()
    plt.savefig(out_dir / "clinical_her2_group_mean_heatmap.png", dpi=180)
    plt.close()


def plot_erbb2_by_group(plt, sns, joined, out_dir: Path) -> None:
    plt.figure(figsize=(6, 4.5))
    sns.boxplot(data=joined.drop_duplicates("case_submitter_id"), x="clinical_her2_group", y="erbb2_tpm", order=GROUP_ORDER)
    sns.stripplot(
        data=joined.drop_duplicates("case_submitter_id"),
        x="clinical_her2_group",
        y="erbb2_tpm",
        order=GROUP_ORDER,
        color="black",
        size=4,
        alpha=0.55,
    )
    plt.xlabel("Clinical HER2 group")
    plt.ylabel("ERBB2 TPM")
    plt.title("ERBB2 RNA Expression by Clinical HER2 Group")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / "erbb2_tpm_by_clinical_her2_group.png", dpi=180)
    plt.close()


def write_markdown(path: Path, joined, channel_summary, pairwise_tests, channels: list[str]) -> None:
    group_counts = (
        joined.groupby("clinical_her2_group", observed=False)["case_submitter_id"]
        .nunique()
        .reindex(GROUP_ORDER, fill_value=0)
    )
    if channel_summary.empty:
        top_table = "No requested GigaTIME channels were available in the slide score table."
    else:
        top = channel_summary.sort_values("kruskal_p_value", na_position="last").head(8)
        lines = [
            "| Channel | Kruskal p | Highest mean group | Lowest mean group | Max-min mean |",
            "|---|---:|---|---|---:|",
        ]
        for _, row in top.iterrows():
            lines.append(
                f"| {row['channel']} | {row['kruskal_p_value']:.4g} | {row['highest_mean_group']} | "
                f"{row['lowest_mean_group']} | {row['max_minus_min_mean']:.4g} |"
            )
        top_table = "\n".join(lines)

    pairwise_note = "No pairwise tests were available."
    if not pairwise_tests.empty:
        top_pairwise = pairwise_tests.sort_values("mannwhitney_q_value_bh", na_position="last").head(8)
        rows = [
            "| Channel | Comparison | Delta mean | Mann-Whitney p | BH q |",
            "|---|---|---:|---:|---:|",
        ]
        for _, row in top_pairwise.iterrows():
            rows.append(
                f"| {row['channel']} | {row['group_a']} vs {row['group_b']} | "
                f"{row['delta_mean_a_minus_b']:.4g} | {row['mannwhitney_p_value']:.4g} | "
                f"{row['mannwhitney_q_value_bh']:.4g} |"
            )
        pairwise_note = "\n".join(rows)

    lines = [
        "# Clinical HER2 GigaTIME Summary",
        "",
        f"- Joined slides: {len(joined)}",
        f"- Unique TCGA cases: {joined['case_submitter_id'].nunique()}",
        f"- Channels summarized: {', '.join(channels)}",
        "",
        "## Clinical HER2 Group Counts",
        "",
        "| Group | Cases/slides |",
        "|---|---:|",
    ]
    for group, count in group_counts.items():
        lines.append(f"| {group} | {int(count)} |")
    lines.extend(
        [
            "",
            "## Top Three-Group Channel Differences",
            "",
            top_table,
            "",
            "## Top Pairwise Channel Tests",
            "",
            pairwise_note,
            "",
            "## Interpretation Guardrails",
            "",
            "- Clinical HER2 labels are based on TCGA IHC/ISH clinical supplement fields.",
            "- GigaTIME outputs are virtual mIF research features, not experimental mIF measurements.",
            "- Small per-group counts should be interpreted as pilot evidence only.",
            "- The group counts above reflect the current joined slide set; missing slide downloads can limit a planned cohort.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")

    slide_scores = pd.read_csv(args.slide_scores)
    cohort = pd.read_csv(args.cohort)
    if args.group_column not in cohort.columns:
        raise ValueError(f"{args.cohort} does not contain group column {args.group_column!r}.")
    cohort["analysis_her2_group"] = cohort[args.group_column]
    cohort["erbb2_tpm"] = pd.to_numeric(cohort["erbb2_tpm"], errors="coerce")
    cohort = cohort.drop_duplicates("case_submitter_id")

    joined = slide_scores.merge(
        cohort[
            [
                "case_submitter_id",
                "analysis_her2_group",
                "clinical_her2_group_rule",
                "clinical_her2_group_confidence",
                "her2_ihc_receptor_status",
                "her2_ihc_score",
                "her2_ish_status",
                "erbb2_tpm",
                "er_status",
                "pr_status",
            ]
        ],
        on="case_submitter_id",
        how="inner",
    )
    if joined.empty:
        raise ValueError("No overlapping case_submitter_id values between slide scores and clinical HER2 cohort.")
    joined = joined.rename(columns={"analysis_her2_group": "clinical_her2_group"})
    joined["clinical_her2_group"] = pd.Categorical(joined["clinical_her2_group"], categories=GROUP_ORDER, ordered=True)
    joined = joined.sort_values(["clinical_her2_group", "case_submitter_id"])
    joined.to_csv(out_dir / "joined_slide_clinical_her2_gigatime.csv", index=False)

    channels = [channel.strip() for channel in args.channels.split(",") if channel.strip()]
    channel_summary = pd.DataFrame(build_channel_summary(joined, channels, stats))
    if not channel_summary.empty:
        channel_summary["kruskal_q_value_bh"] = benjamini_hochberg(channel_summary["kruskal_p_value"].tolist())
    channel_summary.to_csv(out_dir / "clinical_her2_channel_summary.csv", index=False)

    pairwise_tests = pd.DataFrame(build_pairwise_tests(joined, channels, stats))
    pairwise_tests.to_csv(out_dir / "clinical_her2_pairwise_tests.csv", index=False)

    plot_group_boxplots(pd, plt, sns, joined, channels, out_dir)
    plot_group_mean_heatmap(pd, plt, sns, channel_summary, out_dir)
    plot_erbb2_by_group(plt, sns, joined, out_dir)
    write_markdown(out_dir / "clinical_her2_summary.md", joined, channel_summary, pairwise_tests, channels)
    print(f"Wrote clinical HER2 summary to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
