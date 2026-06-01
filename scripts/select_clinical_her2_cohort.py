#!/usr/bin/env python3
"""Select a balanced TCGA-BRCA cohort by clinical HER2 group."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

HER2_GROUPS = ["HER2-positive", "HER2-low", "HER2-zero"]


def require_pandas():
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing Python package: pandas. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd


def local_slide_path(slides_dir: Path, case_id: str, file_name: str) -> Path:
    return slides_dir / case_id / file_name


def read_manifest(pd, path: Path):
    if not path.exists():
        return None
    manifest = pd.read_csv(path, sep="\t")
    return manifest.rename(columns={"id": "file_id", "filename": "file_name"})


def choose_one_slide_per_case(slides, slides_dir: Path):
    slides = slides.copy()
    slides["file_size"] = slides["file_size"].apply(pd_to_numeric)
    slides["slide_local_path"] = slides.apply(
        lambda row: str(local_slide_path(slides_dir, row["case_submitter_id"], row["file_name"])),
        axis=1,
    )
    slides["slide_local_exists"] = slides["slide_local_path"].map(lambda value: Path(value).exists())
    slides["primary_tumor_rank"] = slides["sample_type"].map(lambda value: 0 if value == "Primary Tumor" else 1)
    slides["local_exists_rank"] = slides["slide_local_exists"].map(lambda value: 0 if value else 1)
    slides["file_size_rank"] = slides["file_size"].fillna(float("inf"))
    slides = slides.sort_values(
        [
            "case_submitter_id",
            "primary_tumor_rank",
            "local_exists_rank",
            "file_size_rank",
            "file_name",
        ]
    )
    return slides.drop_duplicates("case_submitter_id", keep="first")


def pd_to_numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def select_cohort(pd, candidates, groups: list[str], per_group: int):
    confidence_rank = {"direct": 0, "inferred": 1, "unknown": 2}
    candidates = candidates.copy()
    candidates["confidence_rank"] = candidates["clinical_her2_group_confidence"].map(confidence_rank).fillna(9)
    candidates["local_exists_rank"] = candidates["slide_local_exists"].map(lambda value: 0 if value else 1)
    candidates["slide_file_size_rank"] = candidates["slide_file_size"].fillna(float("inf"))
    candidates = candidates.sort_values(
        [
            "clinical_her2_group",
            "confidence_rank",
            "local_exists_rank",
            "slide_file_size_rank",
            "case_submitter_id",
        ]
    )
    selected = []
    shortages: dict[str, int] = {}
    for group in groups:
        group_rows = candidates[candidates["clinical_her2_group"] == group].head(per_group).copy()
        group_rows["cohort_group"] = group
        group_rows["selection_rank"] = range(1, len(group_rows) + 1)
        selected.append(group_rows)
        if len(group_rows) < per_group:
            shortages[group] = per_group - len(group_rows)
    if not selected:
        return candidates.iloc[0:0].copy(), shortages
    return pd.concat(selected, ignore_index=True), shortages


def write_manifest(path: Path, selected) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest = selected[
        [
            "slide_file_id",
            "slide_file_name",
            "slide_md5",
            "slide_file_size",
            "slide_state",
        ]
    ].rename(
        columns={
            "slide_file_id": "id",
            "slide_file_name": "filename",
            "slide_md5": "md5",
            "slide_file_size": "size",
            "slide_state": "state",
        }
    )
    manifest.to_csv(path, sep="\t", index=False)


def write_markdown_summary(path: Path, selected, summary: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Clinical HER2 Cohort Selection",
        "",
        "This file summarizes the balanced clinical HER2-positive / HER2-low / HER2-zero pilot cohort selected for the next GigaTIME run.",
        "",
        "## Counts",
        "",
        "| Cohort group | Candidate cases | Selected cases | Already-downloaded slides |",
        "|---|---:|---:|---:|",
    ]
    downloaded_counts = selected.groupby("cohort_group")["slide_local_exists"].sum().to_dict()
    for group in summary["groups"]:
        lines.append(
            f"| {group} | {summary['candidate_counts'].get(group, 0)} | "
            f"{summary['selected_counts'].get(group, 0)} | {int(downloaded_counts.get(group, 0))} |"
        )
    lines.extend(
        [
            "",
            "## Selection Priority",
            "",
            *[f"- {item}." for item in summary["selection_priority"]],
            "",
            "## Selected Cases",
            "",
            "| Group | Rank | Case | HER2 rule | IHC score | ISH status | ER | PR | ERBB2 TPM | Slide downloaded |",
            "|---|---:|---|---|---|---|---|---|---:|---|",
        ]
    )
    for _, row in selected.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["cohort_group"]),
                    str(int(row["selection_rank"])),
                    str(row["case_submitter_id"]),
                    str(row["clinical_her2_group_rule"]),
                    str(row["her2_ihc_score"]),
                    str(row["her2_ish_status"]),
                    str(row["er_status"]),
                    str(row["pr_status"]),
                    f"{float(row['erbb2_tpm']):.4g}",
                    "yes" if bool(row["slide_local_exists"]) else "no",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Local Outputs",
            "",
            f"- Cases CSV: `{summary['outputs']['cases']}`",
            f"- Slide table: `{summary['outputs']['slides']}`",
            f"- Slide manifest: `{summary['outputs']['manifest']}`",
            "",
            "These CSV/TSV files are under `data/`, so they are local reproducible outputs rather than tracked Git files.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", default="data/tcga_brca/clinical_her2_labels.csv")
    parser.add_argument("--expression", default="data/tcga_brca/erbb2_expression.csv")
    parser.add_argument("--slides", default="data/tcga_brca/tcga_brca_diagnostic_slides_files.csv")
    parser.add_argument("--source-manifest", default="data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv")
    parser.add_argument("--slides-dir", default="data/tcga_brca/slides")
    parser.add_argument("--out-cases", default="data/tcga_brca/clinical_her2_cohort_cases.csv")
    parser.add_argument("--out-slides", default="data/tcga_brca/clinical_her2_cohort_slides_files.csv")
    parser.add_argument("--out-manifest", default="data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv")
    parser.add_argument("--out-summary", default="data/tcga_brca/clinical_her2_cohort_summary.json")
    parser.add_argument("--out-markdown", default="docs/clinical_her2_cohort_selection.md")
    parser.add_argument("--per-group", type=int, default=10)
    parser.add_argument("--groups", default=",".join(HER2_GROUPS))
    parser.add_argument(
        "--allow-missing-expression",
        action="store_true",
        help="Select cases with clinical HER2 labels and slides even when ERBB2 RNA has not been downloaded yet.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pd = require_pandas()
    groups = [group.strip() for group in args.groups.split(",") if group.strip()]
    labels = pd.read_csv(args.labels)
    expression = pd.read_csv(args.expression).rename(
        columns={
            "sample_submitter_id": "expression_sample_submitter_id",
            "sample_type": "expression_sample_type",
            "file_id": "expression_file_id",
            "file_name": "expression_file_name",
            "tpm_unstranded": "erbb2_tpm",
        }
    )
    expression["erbb2_tpm"] = pd.to_numeric(expression["erbb2_tpm"], errors="coerce")
    expression = expression.sort_values("expression_sample_type").drop_duplicates("case_submitter_id")

    slides = pd.read_csv(args.slides)
    manifest = read_manifest(pd, Path(args.source_manifest))
    if manifest is not None:
        slides = slides.merge(
            manifest[["file_id", "md5", "state"]].rename(columns={"md5": "slide_md5", "state": "slide_state"}),
            on="file_id",
            how="left",
        )
    else:
        slides["slide_md5"] = ""
        slides["slide_state"] = "released"
    slides = choose_one_slide_per_case(slides, Path(args.slides_dir)).rename(
        columns={
            "file_id": "slide_file_id",
            "file_name": "slide_file_name",
            "file_size": "slide_file_size",
            "data_type": "slide_data_type",
            "data_format": "slide_data_format",
            "experimental_strategy": "slide_experimental_strategy",
            "sample_submitter_id": "slide_sample_submitter_id",
            "sample_type": "slide_sample_type",
        }
    )

    expression_join = "left" if args.allow_missing_expression else "inner"
    candidates = labels.merge(expression, on="case_submitter_id", how=expression_join).merge(
        slides,
        on="case_submitter_id",
        how="inner",
    )
    candidates = candidates[candidates["clinical_her2_group"].isin(groups)].copy()
    selected, shortages = select_cohort(pd, candidates, groups, args.per_group)
    selected = selected.sort_values(["cohort_group", "selection_rank"])
    selected["slide_file_size"] = selected["slide_file_size"].round().astype("Int64")

    case_columns = [
        "cohort_group",
        "selection_rank",
        "case_submitter_id",
        "clinical_her2_group",
        "clinical_her2_group_rule",
        "clinical_her2_group_confidence",
        "her2_ihc_receptor_status",
        "her2_ihc_score",
        "her2_ish_status",
        "her2_cep17_ratio",
        "er_status",
        "pr_status",
        "erbb2_tpm",
        "expression_sample_submitter_id",
        "expression_sample_type",
        "expression_file_id",
        "expression_file_name",
        "slide_file_id",
        "slide_file_name",
        "slide_file_size",
        "slide_experimental_strategy",
        "slide_sample_submitter_id",
        "slide_sample_type",
        "slide_local_path",
        "slide_local_exists",
    ]
    slide_columns = [
        "cohort_group",
        "selection_rank",
        "case_submitter_id",
        "slide_file_id",
        "slide_file_name",
        "slide_file_size",
        "slide_md5",
        "slide_state",
        "slide_data_type",
        "slide_data_format",
        "slide_experimental_strategy",
        "slide_sample_submitter_id",
        "slide_sample_type",
        "slide_local_path",
        "slide_local_exists",
    ]

    out_cases = Path(args.out_cases)
    out_cases.parent.mkdir(parents=True, exist_ok=True)
    selected[case_columns].to_csv(out_cases, index=False)
    selected[slide_columns].to_csv(args.out_slides, index=False)
    write_manifest(Path(args.out_manifest), selected)

    summary = {
        "per_group_requested": args.per_group,
        "groups": groups,
        "candidate_counts": candidates.groupby("clinical_her2_group")["case_submitter_id"].nunique().to_dict(),
        "selected_counts": selected.groupby("cohort_group")["case_submitter_id"].nunique().to_dict(),
        "shortages": shortages,
        "total_selected_cases": int(selected["case_submitter_id"].nunique()),
        "selection_priority": [
            "Clinical HER2 group in requested groups",
            "Direct clinical label before inferred label",
            "Already-downloaded slide before not-yet-downloaded slide",
            "Smaller slide file before larger slide file",
            "Case submitter ID for deterministic tie-breaking",
        ],
        "outputs": {
            "cases": args.out_cases,
            "slides": args.out_slides,
            "manifest": args.out_manifest,
            "markdown": args.out_markdown,
        },
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown_summary(Path(args.out_markdown), selected, summary)

    print(f"Wrote selected cases to {args.out_cases}")
    print(f"Wrote selected slide table to {args.out_slides}")
    print(f"Wrote selected slide manifest to {args.out_manifest}")
    print(f"Wrote summary to {args.out_summary}")
    print(f"Wrote markdown summary to {args.out_markdown}")
    print("\nCandidate cases by clinical HER2 group:")
    print(candidates.groupby("clinical_her2_group")["case_submitter_id"].nunique().to_string())
    print("\nSelected cases by cohort group:")
    print(selected.groupby("cohort_group")["case_submitter_id"].nunique().to_string())
    if shortages:
        print("\nShortages:")
        for group, shortage in shortages.items():
            print(f"  {group}: {shortage}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
