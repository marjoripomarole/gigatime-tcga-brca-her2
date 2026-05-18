#!/usr/bin/env python3
"""Select balanced ERBB2-high and ERBB2-low TCGA-BRCA cases."""

from __future__ import annotations

import argparse
from pathlib import Path


def require_pandas():
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing Python package: pandas. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", default="data/tcga_brca/erbb2_expression.csv")
    parser.add_argument("--slides", default="data/tcga_brca/tcga_brca_diagnostic_slides_files.csv")
    parser.add_argument("--out", default="data/tcga_brca/her2_extreme_cases.csv")
    parser.add_argument("--per-group", type=int, default=20, help="Number of HER2-high and HER2-low cases to select.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pd = require_pandas()
    expression = pd.read_csv(args.expression)
    slides = pd.read_csv(args.slides)
    expression = expression.rename(columns={"tpm_unstranded": "erbb2_tpm"})
    expression["erbb2_tpm"] = pd.to_numeric(expression["erbb2_tpm"], errors="coerce")
    expression = expression.dropna(subset=["case_submitter_id", "erbb2_tpm"])
    expression = expression.sort_values("sample_type").drop_duplicates("case_submitter_id")

    slide_cases = set(slides["case_submitter_id"].dropna())
    expression = expression[expression["case_submitter_id"].isin(slide_cases)].copy()
    if len(expression) < args.per_group * 2:
        raise ValueError(
            f"Need at least {args.per_group * 2} cases with both ERBB2 and slide metadata; found {len(expression)}."
        )

    low = expression.nsmallest(args.per_group, "erbb2_tpm").copy()
    high = expression.nlargest(args.per_group, "erbb2_tpm").copy()
    low["her2_group"] = "HER2-low"
    high["her2_group"] = "HER2-high"
    selected = pd.concat([high, low], ignore_index=True)
    selected["selection_rank"] = selected.groupby("her2_group")["erbb2_tpm"].rank(
        method="first", ascending=False
    )
    selected = selected.sort_values(["her2_group", "selection_rank"])
    columns = [
        "case_submitter_id",
        "sample_submitter_id",
        "sample_type",
        "erbb2_tpm",
        "her2_group",
        "selection_rank",
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    selected[columns].to_csv(out, index=False)
    print(f"Wrote {len(selected)} cases to {out}")
    print(selected.groupby("her2_group")["erbb2_tpm"].agg(["count", "min", "median", "max"]).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
