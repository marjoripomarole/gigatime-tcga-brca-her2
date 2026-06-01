#!/usr/bin/env python3
"""Build a simple notebook and HTML report for the clinical HER2 findings."""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
from pathlib import Path


ASSET_COPIES = [
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_boxplots.png",
        "docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_group_mean_heatmap.png",
        "docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/erbb2_tpm_by_clinical_her2_group.png",
        "docs/assets/clinical_her2_findings/erbb2_tpm_by_clinical_her2_group.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_channel_boxplots.png",
        "docs/assets/clinical_her2_tile256/clinical_her2_channel_boxplots.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_group_mean_heatmap.png",
        "docs/assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/gigatime_rna_correlation_heatmap.png",
        "docs/assets/clinical_her2_tile256/gigatime_rna_correlation_heatmap.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/top_gigatime_rna_signature_scatter.png",
        "docs/assets/clinical_her2_tile256/top_gigatime_rna_signature_scatter.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_rna_program_correlation_heatmap.png",
        "docs/assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/rna_programs_by_her2_group.png",
        "docs/assets/clinical_her2_rna_program_validation/rna_programs_by_her2_group.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_programs_by_her2_group.png",
        "docs/assets/clinical_her2_rna_program_validation/virtual_programs_by_her2_group.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_balanced_accuracy.png",
        "docs/assets/clinical_her2_classifier_baseline/classifier_balanced_accuracy.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/confusion_her2_low_vs_zero.png",
        "docs/assets/clinical_her2_classifier_baseline/confusion_her2_low_vs_zero.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/confusion_her2_positive_vs_negative.png",
        "docs/assets/clinical_her2_classifier_baseline/confusion_her2_positive_vs_negative.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/confusion_her2_three_class.png",
        "docs/assets/clinical_her2_classifier_baseline/confusion_her2_three_class.png",
    ),
]

ROBUSTNESS_CHANNELS = ["CD68", "PD-L1", "CD11c", "CD3", "CD4", "Ki67"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-notebook", default="notebooks/clinical_her2_findings_simple.ipynb")
    parser.add_argument("--out-html", default="notebooks/clinical_her2_findings_simple.html")
    parser.add_argument(
        "--channel-summary",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv",
    )
    parser.add_argument(
        "--pairwise-tests",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv",
    )
    parser.add_argument(
        "--rna-correlations",
        default="results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv",
    )
    parser.add_argument(
        "--visual-qc-cases",
        default="docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv",
    )
    parser.add_argument(
        "--tile256-channel-summary",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/"
            "clinical_her2_channel_summary.csv"
        ),
    )
    parser.add_argument(
        "--tile256-pairwise-tests",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/"
            "clinical_her2_pairwise_tests.csv"
        ),
    )
    parser.add_argument(
        "--tile256-rna-correlations",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/"
            "gigatime_rna_signature_correlations.csv"
        ),
    )
    parser.add_argument(
        "--tile256-visual-qc-cases",
        default="docs/assets/clinical_her2_visual_qc_tile256/clinical_her2_visual_qc_selected_cases.csv",
    )
    parser.add_argument(
        "--rna-program-correlations",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/"
            "virtual_rna_program_correlations.csv"
        ),
    )
    parser.add_argument(
        "--rna-program-group-summary",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/"
            "rna_program_group_summary.csv"
        ),
    )
    parser.add_argument(
        "--virtual-program-group-summary",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/"
            "virtual_program_group_summary.csv"
        ),
    )
    parser.add_argument(
        "--classifier-metrics",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/"
            "classifier_metrics.csv"
        ),
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_optional_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    return read_rows(path)


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    if abs(value) < 0.001 and value != 0:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def html_table(headers: list[str], rows: list[list[str]]) -> str:
    header_html = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def copy_assets() -> None:
    for source, destination in ASSET_COPIES:
        src = Path(source)
        dst = Path(destination)
        if not src.exists():
            raise FileNotFoundError(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def channel_summary_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    top_rows = sorted(rows, key=lambda row: as_float(row, "kruskal_p_value"))[:limit]
    return [
        [
            row["channel"],
            fmt(as_float(row, "kruskal_p_value"), 4),
            fmt(as_float(row, "kruskal_q_value_bh"), 4),
            row["highest_mean_group"],
            row["lowest_mean_group"],
            fmt(as_float(row, "max_minus_min_mean"), 4),
        ]
        for row in top_rows
    ]


def pairwise_summary_rows(rows: list[dict[str, str]], limit: int = 6) -> list[list[str]]:
    top_rows = sorted(rows, key=lambda row: as_float(row, "mannwhitney_p_value"))[:limit]
    return [
        [
            row["channel"],
            f"{row['group_a']} vs {row['group_b']}",
            fmt(as_float(row, "delta_mean_a_minus_b"), 4),
            fmt(as_float(row, "mannwhitney_p_value"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh"), 4),
        ]
        for row in top_rows
    ]


def rna_summary_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    sorted_rows = sorted(rows, key=lambda row: as_float(row, "spearman_rho"), reverse=True)
    return [
        [
            row["channel"],
            fmt(as_float(row, "spearman_rho"), 3),
            fmt(as_float(row, "spearman_p_value"), 4),
            fmt(as_float(row, "spearman_q_value_bh"), 4),
        ]
        for row in sorted_rows
    ]


def qc_summary_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    return [
        [
            row["clinical_her2_group"],
            row["case_submitter_id"],
            fmt(as_float(row, "qc_signal"), 3),
            fmt(as_float(row, "mean_CD68"), 3),
            fmt(as_float(row, "mean_PD-L1"), 3),
            fmt(as_float(row, "mean_CD11c"), 3),
        ]
        for row in rows
    ]


def robustness_comparison_rows(
    baseline_rows: list[dict[str, str]],
    tile256_rows: list[dict[str, str]],
) -> list[list[str]]:
    baseline_by_channel = {row["channel"]: row for row in baseline_rows}
    tile256_by_channel = {row["channel"]: row for row in tile256_rows}
    rows = []
    for channel in ROBUSTNESS_CHANNELS:
        baseline = baseline_by_channel.get(channel)
        tile256 = tile256_by_channel.get(channel)
        if not baseline or not tile256:
            continue
        rows.append(
            [
                channel,
                fmt(as_float(baseline, "kruskal_p_value"), 4),
                fmt(as_float(tile256, "kruskal_p_value"), 4),
                fmt(as_float(baseline, "max_minus_min_mean"), 4),
                fmt(as_float(tile256, "max_minus_min_mean"), 4),
                tile256["highest_mean_group"],
                tile256["lowest_mean_group"],
            ]
        )
    return rows


def program_correlation_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            as_float(row, "spearman_q_value_bh"),
            as_float(row, "spearman_p_value"),
            -abs(as_float(row, "spearman_rho")),
        ),
    )[:limit]
    return [
        [
            row["virtual_program_label"],
            row["rna_program_label"],
            fmt(as_float(row, "spearman_rho"), 3),
            fmt(as_float(row, "spearman_p_value"), 4),
            fmt(as_float(row, "spearman_q_value_bh"), 4),
        ]
        for row in sorted_rows
    ]


def program_group_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    sorted_rows = sorted(rows, key=lambda row: as_float(row, "kruskal_p_value"))[:limit]
    return [
        [
            row["program_label"],
            fmt(as_float(row, "kruskal_p_value"), 4),
            fmt(as_float(row, "kruskal_q_value_bh"), 4),
            row["highest_mean_group"],
            row["lowest_mean_group"],
            fmt(as_float(row, "max_minus_min_mean"), 3),
        ]
        for row in sorted_rows
    ]


def classifier_rows(rows: list[dict[str, str]]) -> tuple[list[list[str]], list[list[str]]]:
    logistic = [row for row in rows if row.get("model") == "regularized_logistic"]
    h_e_rows = [row for row in logistic if row.get("feature_set") != "erbb2_rna_reference_not_h_e"]
    task_order = ["her2_positive_vs_negative", "her2_low_vs_zero", "her2_three_class"]

    def best_for(task: str, candidates: list[dict[str, str]]) -> dict[str, str] | None:
        task_rows = [row for row in candidates if row.get("task") == task]
        if not task_rows:
            return None
        return max(task_rows, key=lambda row: as_float(row, "balanced_accuracy"))

    best_h_e = []
    for task in task_order:
        row = best_for(task, h_e_rows)
        if not row:
            continue
        best_h_e.append(
            [
                row["task_label"],
                row["feature_set"],
                fmt(as_float(row, "accuracy"), 3),
                fmt(as_float(row, "balanced_accuracy"), 3),
                fmt(as_float(row, "macro_auc_ovr"), 3),
                fmt(as_float(row, "sensitivity"), 3),
                fmt(as_float(row, "specificity"), 3),
            ]
        )

    reference_rows = []
    for task in task_order:
        row = best_for(task, [item for item in logistic if item.get("feature_set") == "erbb2_rna_reference_not_h_e"])
        if not row:
            continue
        reference_rows.append(
            [
                row["task_label"],
                fmt(as_float(row, "accuracy"), 3),
                fmt(as_float(row, "balanced_accuracy"), 3),
                fmt(as_float(row, "macro_auc_ovr"), 3),
            ]
        )
    return best_h_e, reference_rows


def build_content(args: argparse.Namespace):
    copy_assets()
    channels = read_rows(Path(args.channel_summary))
    pairwise = read_rows(Path(args.pairwise_tests))
    rna = read_rows(Path(args.rna_correlations))
    qc_cases = read_rows(Path(args.visual_qc_cases))
    tile256_channels = read_optional_rows(Path(args.tile256_channel_summary))
    tile256_pairwise = read_optional_rows(Path(args.tile256_pairwise_tests))
    tile256_rna = read_optional_rows(Path(args.tile256_rna_correlations))
    tile256_qc_cases = read_optional_rows(Path(args.tile256_visual_qc_cases))
    rna_program_correlations = read_optional_rows(Path(args.rna_program_correlations))
    rna_program_groups = read_optional_rows(Path(args.rna_program_group_summary))
    virtual_program_groups = read_optional_rows(Path(args.virtual_program_group_summary))
    classifier_metrics = read_optional_rows(Path(args.classifier_metrics))
    classifier_h_e_rows, classifier_reference_rows = classifier_rows(classifier_metrics)

    return {
        "channel_rows": channel_summary_rows(channels),
        "pair_rows": pairwise_summary_rows(pairwise),
        "rna_rows": rna_summary_rows(rna),
        "qc_rows": qc_summary_rows(qc_cases),
        "tile256_channel_rows": channel_summary_rows(tile256_channels),
        "tile256_pair_rows": pairwise_summary_rows(tile256_pairwise),
        "tile256_rna_rows": rna_summary_rows(tile256_rna),
        "tile256_qc_rows": qc_summary_rows(tile256_qc_cases),
        "tile256_compare_rows": robustness_comparison_rows(channels, tile256_channels),
        "program_correlation_rows": program_correlation_rows(rna_program_correlations),
        "rna_program_group_rows": program_group_rows(rna_program_groups),
        "virtual_program_group_rows": program_group_rows(virtual_program_groups),
        "classifier_h_e_rows": classifier_h_e_rows,
        "classifier_reference_rows": classifier_reference_rows,
    }


def notebook_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source.strip().splitlines()],
    }


def build_notebook(args: argparse.Namespace, content: dict[str, list[list[str]]]) -> None:
    nb_path = Path(args.out_notebook)
    nb_path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        notebook_cell(
            """
# Clinical HER2 GigaTIME Findings

**Simple display notebook**  
Updated clinical HER2 pilot summary for TCGA-BRCA.

**Core message:** we completed a balanced 30-slide pilot across HER2-positive, HER2-low, and HER2-zero cases. The first 64-tile run suggested higher immune/checkpoint-like signal in HER2-zero than HER2-low, and the denser 256-tile robustness run reproduced and strengthened the same CD68, PD-L1, and CD11c pattern. Marker-level and broader RNA-program validation remain weak. A first cross-validated classifier suggests possible HER2-low versus HER2-zero signal, but not reliable clinical HER2 diagnosis.
            """
        ),
        notebook_cell(
            """
## Study Design

- Dataset: TCGA-BRCA breast cancer.
- Images: diagnostic H&E whole-slide images.
- Model: released GigaTIME model.
- Groups: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero cases.
- Sampling: first 64 random tissue tiles per slide, then a 256-tile robustness rerun on the same 30 slides.
- Main outputs: virtual mIF channel scores, RNA validation, and visual QC panels.

**Plain-language translation:** we asked whether an AI model can see immune-like tissue patterns on ordinary H&E slides that differ across clinical HER2 groups.
            """
        ),
        notebook_cell(
            """
## Main Finding

The strongest pilot pattern was:

> HER2-zero had higher GigaTIME-predicted immune/checkpoint signals than HER2-low, especially CD68, PD-L1, and CD11c.

HER2-positive was usually between HER2-zero and HER2-low for these channels. The 256-tile rerun made the same top signal more stable: CD68, PD-L1, and CD11c again ranked as the strongest three-group differences, with HER2-zero highest and HER2-low lowest.

This should be described as a **hypothesis-generating pilot signal**, not a final biological claim.
            """
        ),
        notebook_cell(
            "## Top Three-Group Differences\n\n"
            + markdown_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["channel_rows"],
            )
        ),
        notebook_cell(
            "## Group-Level Virtual mIF Summary\n\n"
            "This heatmap shows mean GigaTIME virtual-channel activation by clinical HER2 group.\n\n"
            "![Clinical HER2 group mean heatmap](../docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png)"
        ),
        notebook_cell(
            "## Channel Distributions\n\n"
            "Each dot is one TCGA case. This is useful for seeing how noisy the pilot still is.\n\n"
            "![Clinical HER2 channel boxplots](../docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png)"
        ),
        notebook_cell(
            "## 256-Tile Robustness Check\n\n"
            "We reran the same 30 selected slides with up to 256 random tissue tiles per slide. This tests whether the original 64-tile result was just a sampling accident.\n\n"
            + markdown_table(
                [
                    "Channel",
                    "64-tile p",
                    "256-tile p",
                    "64 max-min",
                    "256 max-min",
                    "256 highest",
                    "256 lowest",
                ],
                content["tile256_compare_rows"],
            )
            + "\n\n**Interpretation:** the same HER2-zero > HER2-low immune/checkpoint pattern persisted. For CD68, PD-L1, and CD11c, the three-group p values became smaller and the group mean gaps became larger after denser sampling."
        ),
        notebook_cell(
            "## 256-Tile Group Summary\n\n"
            + markdown_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["tile256_channel_rows"],
            )
            + "\n\n![256-tile clinical HER2 heatmap](../docs/assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png)\n\n"
            "![256-tile clinical HER2 boxplots](../docs/assets/clinical_her2_tile256/clinical_her2_channel_boxplots.png)"
        ),
        notebook_cell(
            "## Top Pairwise Tests\n\n"
            + markdown_table(
                ["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"],
                content["pair_rows"],
            )
            + "\n\nThe strongest pairwise tests were mostly HER2-low versus HER2-zero, but they did not survive multiple-testing correction."
        ),
        notebook_cell(
            "## 256-Tile Top Pairwise Tests\n\n"
            + markdown_table(
                ["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"],
                content["tile256_pair_rows"],
            )
            + "\n\nThe top 256-tile pairwise comparisons again focused on HER2-low versus HER2-zero. The best BH q values improved to about 0.113 for CD68, PD-L1, and CD11c, but they still did not cross the usual 0.05 FDR threshold."
        ),
        notebook_cell(
            "## RNA Validation Check\n\n"
            "We compared GigaTIME virtual channels with matched RNA-seq marker signatures from the same 30 cases.\n\n"
            + markdown_table(["Channel", "Spearman rho", "p", "BH q"], content["rna_rows"])
            + "\n\n**Interpretation:** RNA validation did not strongly confirm the virtual immune-channel signal. Ki67 had the strongest positive trend, but no channel was FDR-significant."
        ),
        notebook_cell(
            "## 256-Tile RNA Validation Check\n\n"
            "The 256-tile rerun did not solve the RNA-discordance problem. No virtual channel had an FDR-significant correlation with its matched RNA marker signature.\n\n"
            + markdown_table(["Channel", "Spearman rho", "p", "BH q"], content["tile256_rna_rows"])
            + "\n\n![256-tile GigaTIME RNA correlation heatmap](../docs/assets/clinical_her2_tile256/gigatime_rna_correlation_heatmap.png)"
        ),
        notebook_cell(
            "## Broader RNA Program Validation\n\n"
            "We also tested broader RNA programs, such as T-cell/cytotoxic, checkpoint/IFNG, myeloid/macrophage, B-cell, proliferation, epithelial, stromal, and endothelial signatures. This is a stronger trustworthiness check than single-marker genes.\n\n"
            + markdown_table(
                ["Virtual program", "RNA program", "Spearman rho", "p", "BH q"],
                content["program_correlation_rows"],
            )
            + "\n\n**Interpretation:** this still did not positively validate the virtual immune/checkpoint signal. The strongest FDR-significant associations were negative correlations between virtual immune/checkpoint programs and endothelial RNA signal, which should be treated as a warning sign and a reason for pathologist review.\n\n"
            "![Virtual vs RNA program heatmap](../docs/assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png)"
        ),
        notebook_cell(
            "## RNA Programs Across HER2 Groups\n\n"
            + markdown_table(
                ["RNA program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["rna_program_group_rows"],
            )
            + "\n\nNo broad RNA immune program showed an FDR-significant HER2-group difference in this 30-case pilot.\n\n"
            "![RNA programs by HER2 group](../docs/assets/clinical_her2_rna_program_validation/rna_programs_by_her2_group.png)"
        ),
        notebook_cell(
            "## Virtual Programs Across HER2 Groups\n\n"
            + markdown_table(
                ["Virtual program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["virtual_program_group_rows"],
            )
            + "\n\nThe virtual myeloid/checkpoint composite kept the HER2-zero > HER2-low direction, but remained short of FDR significance.\n\n"
            "![Virtual programs by HER2 group](../docs/assets/clinical_her2_rna_program_validation/virtual_programs_by_her2_group.png)"
        ),
        notebook_cell(
            "## First HER2 Classifier Baseline\n\n"
            "We trained a first slide-level classifier using GigaTIME features and leave-one-out cross-validation. This is the first diagnostic-model style check, but it is still only a 30-case pilot.\n\n"
            + markdown_table(
                ["Task", "Best GigaTIME/H&E feature set", "Accuracy", "Balanced accuracy", "Macro AUC", "Sensitivity", "Specificity"],
                content["classifier_h_e_rows"],
            )
            + "\n\n**Interpretation:** GigaTIME features were promising for HER2-low versus HER2-zero, but did not reliably classify HER2-positive versus HER2-negative or the full three-class HER2 grouping.\n\n"
            "![Classifier balanced accuracy](../docs/assets/clinical_her2_classifier_baseline/classifier_balanced_accuracy.png)"
        ),
        notebook_cell(
            "## GigaTIME/H&E Classifier Confusion Matrices\n\n"
            "These confusion matrices show the best GigaTIME/H&E regularized logistic model for each task. They do not use the ERBB2 RNA reference feature.\n\n"
            "### HER2-low versus HER2-zero\n\n"
            "![HER2-low versus HER2-zero confusion matrix](../docs/assets/clinical_her2_classifier_baseline/confusion_her2_low_vs_zero.png)\n\n"
            "### HER2-positive versus HER2-negative\n\n"
            "![HER2-positive versus HER2-negative confusion matrix](../docs/assets/clinical_her2_classifier_baseline/confusion_her2_positive_vs_negative.png)\n\n"
            "### Three-class HER2 group\n\n"
            "![Three-class HER2 confusion matrix](../docs/assets/clinical_her2_classifier_baseline/confusion_her2_three_class.png)"
        ),
        notebook_cell(
            "## ERBB2 RNA Reference Classifier\n\n"
            "ERBB2 RNA was included only as a non-H&E reference. It is not an image-based model.\n\n"
            + markdown_table(
                ["Task", "Accuracy", "Balanced accuracy", "Macro AUC"],
                content["classifier_reference_rows"],
            )
            + "\n\nThis comparison matters because ERBB2 RNA predicted HER2-positive versus HER2-negative better than the current GigaTIME/H&E features. That tells us the clinical labels have molecular signal, but the current image-derived classifier is not capturing it reliably yet."
        ),
        notebook_cell(
            "## RNA Correlation Heatmap\n\n"
            "![GigaTIME RNA correlation heatmap](../docs/assets/clinical_her2_rna_validation/gigatime_rna_correlation_heatmap.png)"
        ),
        notebook_cell(
            "## Visual QC Check\n\n"
            "We selected the top case from each HER2 group by combined CD68 + PD-L1 + CD11c virtual signal.\n\n"
            + markdown_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["qc_rows"])
            + "\n\nVisual QC showed that high-scoring tiles were tissue-containing and cellular, not obvious blank background. This supports continuing the analysis, but it does not validate the virtual markers."
        ),
        notebook_cell(
            "## 256-Tile Visual QC Check\n\n"
            "The 256-tile visual QC selected the same representative cases, but the HER2-zero combined signal increased clearly.\n\n"
            + markdown_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["tile256_qc_rows"])
            + "\n\n![256-tile HER2-zero visual QC](../docs/assets/clinical_her2_visual_qc_tile256/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png)"
        ),
        notebook_cell(
            "## Example H&E vs Virtual mIF Panels\n\n"
            "### HER2-zero top case\n\n"
            "![HER2-zero visual QC](../docs/assets/clinical_her2_visual_qc/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png)\n\n"
            "### HER2-low top case\n\n"
            "![HER2-low visual QC](../docs/assets/clinical_her2_visual_qc/her2_low_TCGA-A2-A04Q_he_vs_virtual_mif_qc.png)\n\n"
            "### HER2-positive top case\n\n"
            "![HER2-positive visual QC](../docs/assets/clinical_her2_visual_qc/her2_positive_TCGA-A2-A0EQ_he_vs_virtual_mif_qc.png)"
        ),
        notebook_cell(
            """
## What We Can Say

- The full 30-slide clinical HER2 pilot is complete.
- The most interesting signal is HER2-zero > HER2-low for virtual CD68, PD-L1, and CD11c.
- Visual QC suggests the signal is not simply blank background.
- The 256-tile robustness run reproduced and strengthened the same CD68, PD-L1, and CD11c direction.
- Single-marker and broader RNA-program validation still did not strongly confirm the virtual immune-channel signal.
- A first classifier showed possible HER2-low versus HER2-zero signal, but not reliable HER2-positive or three-class diagnosis.

## What We Should Not Say Yet

- Do not say GigaTIME validated real mIF in TCGA.
- Do not say HER2-zero definitively has more immune infiltration.
- Do not say the model can classify HER2 status.
- Do not overinterpret p values from 10 cases per group.
            """
        ),
        notebook_cell(
            """
## Next Step

The 256-tile robustness check is now complete. The clean next step is validation:

1. Ask an advisor/pathologist to review the high-signal H&E tiles and virtual mIF panels.
2. Improve classifier inputs by restricting to tumor-rich tiles and adding tile distribution features.
3. Add tumor purity or immune deconvolution covariates if available.
4. Consider a 512-tile or more exhaustive run if compute time allows.
5. Find an external dataset with paired H&E and real mIF, if possible.
6. Expand beyond 10 cases per group only after the QC and validation logic is clear.

**One-sentence proposal framing:** GigaTIME produced a reproducible but still unvalidated virtual immune/checkpoint signal separating HER2-zero from HER2-low in a balanced TCGA-BRCA pilot, motivating pathologist review and orthogonal validation.
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
    nb_path.write_text(json.dumps(notebook, indent=2) + "\n", encoding="utf-8")


def section(title: str, body: str) -> str:
    return f"<section><h2>{html.escape(title)}</h2>{body}</section>"


def build_html(args: argparse.Namespace, content: dict[str, list[list[str]]]) -> None:
    html_path = Path(args.out_html)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    css = """
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #1f2933; background: #f6f7f9; }
main { max-width: 1080px; margin: 0 auto; padding: 42px 28px 70px; }
.hero { background: #111827; color: white; padding: 34px; border-radius: 8px; }
.hero h1 { margin: 0 0 10px; font-size: 34px; letter-spacing: 0; }
.hero p { font-size: 18px; line-height: 1.5; max-width: 880px; }
section { background: white; margin-top: 18px; padding: 26px; border-radius: 8px; border: 1px solid #e5e7eb; }
h2 { margin: 0 0 14px; font-size: 24px; }
h3 { margin-top: 24px; }
p, li { font-size: 16px; line-height: 1.55; }
.callout { background: #eef6ff; border-left: 5px solid #2563eb; padding: 14px 16px; margin: 16px 0; }
.warning { background: #fff7ed; border-left: 5px solid #f97316; padding: 14px 16px; margin: 16px 0; }
table { width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 14px; }
th, td { text-align: left; border-bottom: 1px solid #e5e7eb; padding: 8px 9px; vertical-align: top; }
th { background: #f3f4f6; font-weight: 650; }
img { display: block; width: 100%; max-width: 980px; height: auto; margin: 14px auto; border: 1px solid #e5e7eb; border-radius: 6px; background: white; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }
.small { color: #52606d; font-size: 14px; }
"""
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Clinical HER2 GigaTIME Findings</title>",
        f"<style>{css}</style></head><body><main>",
        """
<div class="hero">
  <h1>Clinical HER2 GigaTIME Findings</h1>
  <p>Simple summary of the TCGA-BRCA pilot so far: 30 slides, three clinical HER2 groups, GigaTIME virtual mIF outputs, RNA validation, RNA program validation, visual QC, a 256-tile robustness check, and a first classifier baseline.</p>
</div>
""",
        section(
            "Bottom Line",
            """
<div class="callout">
  <p><strong>Main signal:</strong> HER2-zero had higher GigaTIME-predicted immune/checkpoint-like signal than HER2-low, especially CD68, PD-L1, and CD11c. The same pattern persisted when the same 30 slides were rerun with up to 256 tissue tiles per slide.</p>
</div>
<div class="warning">
  <p><strong>Important caution:</strong> Marker-level and broader RNA-program validation were still weak, and the pairwise tests still did not survive multiple-testing correction. This is a proposal-ready hypothesis, not a validated biological claim.</p>
</div>
""",
        ),
        section(
            "Study Design",
            """
<ul>
  <li>Dataset: TCGA-BRCA breast cancer.</li>
  <li>Groups: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero cases.</li>
  <li>Input: diagnostic H&E whole-slide images.</li>
  <li>Model: released GigaTIME virtual mIF model.</li>
  <li>Sampling: first 64 random tissue tiles per slide, then a 256-tile robustness rerun on the same 30 slides.</li>
</ul>
""",
        ),
        section(
            "Top Three-Group Differences",
            html_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["channel_rows"],
            ),
        ),
        section(
            "Group-Level Virtual mIF Summary",
            "<p>This heatmap summarizes mean virtual-channel activation by clinical HER2 group.</p>"
            "<img src='../docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png' alt='Clinical HER2 group mean heatmap'>",
        ),
        section(
            "Channel Distributions",
            "<p>Each dot is one TCGA case. The spread shows why this remains a pilot.</p>"
            "<img src='../docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png' alt='Clinical HER2 channel boxplots'>",
        ),
        section(
            "256-Tile Robustness Check",
            "<p>We reran the same 30 selected slides with up to 256 random tissue tiles per slide. This checks whether the 64-tile result was sensitive to sparse sampling.</p>"
            + html_table(
                [
                    "Channel",
                    "64-tile p",
                    "256-tile p",
                    "64 max-min",
                    "256 max-min",
                    "256 highest",
                    "256 lowest",
                ],
                content["tile256_compare_rows"],
            )
            + "<p class='small'>CD68, PD-L1, and CD11c kept the same direction: HER2-zero highest and HER2-low lowest. Their p values became smaller and their group mean gaps became larger.</p>",
        ),
        section(
            "256-Tile Group-Level Summary",
            html_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["tile256_channel_rows"],
            )
            + "<img src='../docs/assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png' alt='256-tile clinical HER2 heatmap'>"
            + "<img src='../docs/assets/clinical_her2_tile256/clinical_her2_channel_boxplots.png' alt='256-tile clinical HER2 channel boxplots'>",
        ),
        section(
            "Top Pairwise Tests",
            html_table(["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"], content["pair_rows"])
            + "<p class='small'>The strongest pairwise tests were mostly HER2-low versus HER2-zero, but none were FDR-significant.</p>",
        ),
        section(
            "256-Tile Top Pairwise Tests",
            html_table(
                ["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"],
                content["tile256_pair_rows"],
            )
            + "<p class='small'>The top 256-tile pairwise tests again focused on HER2-low versus HER2-zero. The leading BH q values improved to about 0.113 for CD68, PD-L1, and CD11c, but still did not cross 0.05.</p>",
        ),
        section(
            "RNA Validation",
            "<p>We compared virtual channels with matched RNA-seq marker signatures. No channel was FDR-significant.</p>"
            + html_table(["Channel", "Spearman rho", "p", "BH q"], content["rna_rows"])
            + "<img src='../docs/assets/clinical_her2_rna_validation/gigatime_rna_correlation_heatmap.png' alt='GigaTIME RNA correlation heatmap'>",
        ),
        section(
            "256-Tile RNA Validation",
            "<p>The 256-tile rerun did not solve the RNA-discordance problem. No channel was FDR-significant against its matched RNA marker signature.</p>"
            + html_table(["Channel", "Spearman rho", "p", "BH q"], content["tile256_rna_rows"])
            + "<img src='../docs/assets/clinical_her2_tile256/gigatime_rna_correlation_heatmap.png' alt='256-tile GigaTIME RNA correlation heatmap'>",
        ),
        section(
            "Broader RNA Program Validation",
            "<p>We also tested broader RNA programs, including T-cell/cytotoxic, checkpoint/IFNG, myeloid/macrophage, B-cell, proliferation, epithelial, stromal, and endothelial signatures.</p>"
            + html_table(
                ["Virtual program", "RNA program", "Spearman rho", "p", "BH q"],
                content["program_correlation_rows"],
            )
            + "<p class='small'>This still did not positively validate the virtual immune/checkpoint signal. The strongest FDR-significant associations were negative correlations between virtual immune/checkpoint programs and endothelial RNA signal.</p>"
            + "<img src='../docs/assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png' alt='Virtual vs RNA program correlation heatmap'>",
        ),
        section(
            "RNA Programs Across HER2 Groups",
            html_table(
                ["RNA program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["rna_program_group_rows"],
            )
            + "<p class='small'>No broad RNA immune program showed an FDR-significant HER2-group difference in this 30-case pilot.</p>"
            + "<img src='../docs/assets/clinical_her2_rna_program_validation/rna_programs_by_her2_group.png' alt='RNA programs by HER2 group'>",
        ),
        section(
            "Virtual Programs Across HER2 Groups",
            html_table(
                ["Virtual program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["virtual_program_group_rows"],
            )
            + "<p class='small'>The virtual myeloid/checkpoint composite kept the HER2-zero greater than HER2-low direction, but remained short of FDR significance.</p>"
            + "<img src='../docs/assets/clinical_her2_rna_program_validation/virtual_programs_by_her2_group.png' alt='Virtual programs by HER2 group'>",
        ),
        section(
            "First HER2 Classifier Baseline",
            "<p>We trained a first slide-level classifier using GigaTIME features and leave-one-out cross-validation. This is a diagnostic-model style check, but still only a 30-case pilot.</p>"
            + html_table(
                [
                    "Task",
                    "Best GigaTIME/H&E feature set",
                    "Accuracy",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                content["classifier_h_e_rows"],
            )
            + "<p class='small'>GigaTIME features were promising for HER2-low versus HER2-zero, but did not reliably classify HER2-positive versus HER2-negative or the full three-class HER2 grouping.</p>"
            + "<img src='../docs/assets/clinical_her2_classifier_baseline/classifier_balanced_accuracy.png' alt='HER2 classifier balanced accuracy'>",
        ),
        section(
            "GigaTIME/H&E Classifier Confusion Matrices",
            """
<p>These confusion matrices show the best GigaTIME/H&E regularized logistic model for each task. They do not use the ERBB2 RNA reference feature.</p>
<div class="grid">
  <div><h3>HER2-low vs HER2-zero</h3><img src="../docs/assets/clinical_her2_classifier_baseline/confusion_her2_low_vs_zero.png" alt="HER2-low versus HER2-zero confusion matrix"></div>
  <div><h3>HER2-positive vs HER2-negative</h3><img src="../docs/assets/clinical_her2_classifier_baseline/confusion_her2_positive_vs_negative.png" alt="HER2-positive versus HER2-negative confusion matrix"></div>
  <div><h3>Three-class HER2</h3><img src="../docs/assets/clinical_her2_classifier_baseline/confusion_her2_three_class.png" alt="Three-class HER2 confusion matrix"></div>
</div>
""",
        ),
        section(
            "ERBB2 RNA Reference Classifier",
            "<p>ERBB2 RNA was included only as a non-H&E reference. It is not an image-based model.</p>"
            + html_table(
                ["Task", "Accuracy", "Balanced accuracy", "Macro AUC"],
                content["classifier_reference_rows"],
            )
            + "<p class='small'>ERBB2 RNA predicted HER2-positive versus HER2-negative better than the current GigaTIME/H&E features, which means the image-derived classifier is not capturing that diagnostic signal reliably yet.</p>",
        ),
        section(
            "Visual QC",
            "<p>Top cases were selected by combined CD68 + PD-L1 + CD11c virtual signal.</p>"
            + html_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["qc_rows"])
            + "<p>The high-scoring tiles contain tissue and cells rather than obvious blank background. This supports follow-up, but it is not biological validation.</p>",
        ),
        section(
            "256-Tile Visual QC",
            "<p>The 256-tile visual QC selected the same representative cases. The HER2-zero representative showed a clearer combined CD68 + PD-L1 + CD11c signal.</p>"
            + html_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["tile256_qc_rows"])
            + "<img src='../docs/assets/clinical_her2_visual_qc_tile256/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png' alt='256-tile HER2-zero visual QC'>",
        ),
        section(
            "Example QC Panels",
            """
<div class="grid">
  <div><h3>HER2-zero</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png" alt="HER2-zero visual QC"></div>
  <div><h3>HER2-low</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_low_TCGA-A2-A04Q_he_vs_virtual_mif_qc.png" alt="HER2-low visual QC"></div>
  <div><h3>HER2-positive</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_positive_TCGA-A2-A0EQ_he_vs_virtual_mif_qc.png" alt="HER2-positive visual QC"></div>
</div>
""",
        ),
        section(
            "What To Say",
            """
<ul>
  <li>We completed the balanced 30-case clinical HER2 pilot.</li>
  <li>The leading signal is HER2-zero greater than HER2-low for virtual CD68, PD-L1, and CD11c.</li>
  <li>The 256-tile robustness run reproduced and strengthened that same direction.</li>
  <li>Visual QC makes the signal look plausible enough to follow.</li>
  <li>Marker-level and RNA-program validation did not confirm the signal strongly, so the claim must stay cautious.</li>
  <li>A first classifier showed possible HER2-low versus HER2-zero signal, but not reliable HER2-positive or three-class diagnosis.</li>
</ul>
""",
        ),
        section(
            "Next Step",
            """
<p>The 256-tile robustness check, richer RNA-program validation, and first classifier baseline are complete. The next step is trustworthiness review and better classifier inputs: pathologist review of high-signal tiles, tumor-rich tile selection, tumor-purity or immune-deconvolution adjustment, and ideally an external dataset with paired H&E and real mIF.</p>
<p><strong>Proposal framing:</strong> GigaTIME produced a reproducible but unvalidated virtual immune/checkpoint signal separating HER2-zero from HER2-low in a balanced TCGA-BRCA pilot, motivating pathologist review and orthogonal validation.</p>
""",
        ),
        "</main></body></html>",
    ]
    html_path.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    args = parse_args()
    content = build_content(args)
    build_notebook(args, content)
    build_html(args, content)
    print(f"Wrote {args.out_notebook}")
    print(f"Wrote {args.out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
