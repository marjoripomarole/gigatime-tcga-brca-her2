#!/usr/bin/env python3
"""Draft comparison figure: GigaTIME held-out mIF agreement vs our RNA specificity.

Joins two per-channel quantities:
  - GigaTIME's PUBLISHED per-channel agreement with held-out measured mIF
    (Pearson of activation density, virtual vs measured), transcribed from
    Valanarasu et al., Cell 2026, Figure S5 (and Figure 2C). These are read off
    the published scatter-panel annotations, so treat as +/- 0.01 and verify
    against the source figure before publication.
  - OUR within-slide Xenium RNA validation (raw Spearman and the cellularity-
    controlled partial correlation), read live from the validation report JSON
    produced by scripts/validate_gigatime_xenium_rna.py.

The contribution in one figure: high agreement with the training modality (mIF)
does NOT imply marker specificity against an independent modality (RNA) once
per-tile cellularity is controlled.

Run with the gigatime-tcga env python, e.g.:
  ~/miniconda3/envs/gigatime-tcga/bin/python scripts/make_gigatime_vs_rna_specificity_figure.py
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

# GigaTIME per-channel Pearson (virtual vs measured mIF; GigaTIME's own 8x8-box activation-count
# metric) is RECOMPUTED on the released 50-patch sample test set with the released model
# (scripts/recompute_gigatime_mif_pearson.py), replacing the earlier eyeballed Fig S5 transcription.
# NOTE: the released sample is NOT the paper's full test set; the paper reports higher full-set values
# in Fig S5 (e.g. CK ~0.96, CD3 ~0.89). These recomputed values are a reproducible released-sample
# measurement used as the figure's x-axis.
GIGATIME_SOURCE = ("recomputed on 50 released test patches (prov-gigatime/GigaTIME, GigaTIME 8x8-box "
                   "Pearson); released sample, not paper full-set Fig S5")


# Recomputed values, committed for repo-only reproducibility (results/ and the released test data are
# gitignored). Regenerate via scripts/recompute_gigatime_mif_pearson.py; the CSV overrides this fallback.
RECOMPUTED_MIF_PEARSON: dict[str, float] = {
    "DAPI": 0.720, "CD11c": 0.558, "PD-L1": 0.550, "CK": 0.530, "CD16": 0.517,
    "CD68": 0.516, "CD138": 0.490, "CD4": 0.473, "CD3": 0.453, "CD34": 0.365,
    "Ki67": 0.357, "T-bet": 0.334, "CD14": 0.302, "CD8": 0.231, "PD-1": 0.172,
    "Tryptase": 0.135,
}


def load_gigatime_mif(csv_path: Path) -> dict[str, float]:
    """Per-channel GigaTIME-vs-measured-mIF Pearson: prefer the recompute CSV, else committed values."""
    import csv as _csv

    out: dict[str, float] = {}
    try:
        with csv_path.open(encoding="utf-8") as fh:
            for r in _csv.DictReader(fh):
                try:
                    v = float(r["pearson_mean"])
                except (TypeError, ValueError, KeyError):
                    continue
                if v == v:  # skip NaN
                    out[r["channel"]] = v
    except FileNotFoundError:
        pass
    return out or dict(RECOMPUTED_MIF_PEARSON)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rna-report", type=Path,
                   default=Path("results/gigatime_xenium_rna_validation/xenium_rna_validation_report.json"))
    p.add_argument("--asset-dir", type=Path, default=Path("docs/assets/gigatime_vs_rna_specificity"))
    p.add_argument("--out-csv", type=Path,
                   default=Path("results/gigatime_xenium_rna_validation/gigatime_vs_rna_specificity.csv"))
    p.add_argument("--gigatime-mif-csv", type=Path,
                   default=Path("results/gigatime_mif_recompute/per_channel_pearson.csv"))
    return p.parse_args()


def load_rna(report_path: Path) -> dict[str, dict]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    raw = {r["channel"]: r for r in report.get("channel_results", [])}
    spec = {r["channel"]: r for r in (report.get("specificity") or {}).get("per_channel", [])}
    rows: dict[str, dict] = {}
    for channel in set(raw) | set(spec):
        rr = raw.get(channel, {})
        sp = spec.get(channel, {})
        rows[channel] = {
            "rna_raw_spearman": rr.get("spearman_r"),
            "rna_partial_r": sp.get("partial_r_control_total"),
            "rna_partial_ci_low": sp.get("partial_ci95_low"),
            "rna_partial_ci_high": sp.get("partial_ci95_high"),
            "own_is_row_max": sp.get("own_is_row_max"),
            "n_transcripts": rr.get("total_transcripts_on_grid"),
        }
    return rows


def build_table(rna: dict[str, dict], mif: dict[str, float]) -> list[dict]:
    rows = []
    for channel, rec in rna.items():
        rows.append({
            "channel": channel,
            "gigatime_mif_pearson": mif.get(channel),
            **rec,
        })
    # Sort by GigaTIME mIF agreement descending, missing last.
    rows.sort(key=lambda d: (d["gigatime_mif_pearson"] is None, -(d["gigatime_mif_pearson"] or 0)))
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["channel", "gigatime_mif_pearson", "rna_raw_spearman", "rna_partial_r",
              "rna_partial_ci_low", "rna_partial_ci_high", "own_is_row_max", "n_transcripts"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def point_color(row: dict) -> str:
    pr = row.get("rna_partial_r")
    if pr is None:
        return "#9ca3af"
    if row.get("own_is_row_max") and pr > 0:
        return "#15803d"  # specific (row-max + positive partial)
    if pr <= 0:
        return "#b91c1c"  # fails specificity
    return "#2563eb"      # partial-positive but not row-max


def make_scatter(plt, rows: list[dict], out_path: Path) -> None:
    pts = [r for r in rows if r["gigatime_mif_pearson"] is not None and r["rna_partial_r"] is not None]
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.axhline(0.0, color="#374151", lw=1, ls="--")
    ax.axvspan(0.5, 1.0, ymin=0, ymax=0.5, color="#fee2e2", alpha=0.5, zorder=0)
    ax.text(0.98, ax.get_ylim()[0] if False else -0.30, "high mIF agreement,\nfails RNA specificity",
            ha="right", va="bottom", fontsize=8, color="#b91c1c")
    for r in pts:
        x, y = r["gigatime_mif_pearson"], r["rna_partial_r"]
        ax.scatter(x, y, s=70, color=point_color(r), edgecolors="white", linewidths=0.8, zorder=3)
        ax.annotate(r["channel"], (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)
    ax.set_xlabel("GigaTIME agreement with held-out mIF (Pearson, Fig S5)")
    ax.set_ylabel("Our RNA specificity: partial r | per-tile cellularity (Xenium)")
    ax.set_xlim(0.0, 1.02)
    ax.set_title("Held-out mIF agreement does not imply RNA channel specificity")
    # Legend
    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#15803d", markersize=9, label="own-gene row-max & partial>0 (specific)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2563eb", markersize=9, label="partial>0 but not row-max"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#b91c1c", markersize=9, label="partial<=0 (fails specificity)"),
    ]
    ax.legend(handles=legend, fontsize=7.5, loc="upper left", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def make_bars(plt, rows: list[dict], out_path: Path) -> None:
    import numpy as np

    pts = [r for r in rows if r["gigatime_mif_pearson"] is not None and r["rna_partial_r"] is not None]
    channels = [r["channel"] for r in pts]
    x = np.arange(len(channels))
    w = 0.27
    mif = [r["gigatime_mif_pearson"] for r in pts]
    raw = [r["rna_raw_spearman"] for r in pts]
    par = [r["rna_partial_r"] for r in pts]
    fig, ax = plt.subplots(figsize=(11.5, 5.0))
    ax.axhline(0.0, color="#374151", lw=0.8)
    ax.bar(x - w, mif, w, label="GigaTIME vs held-out mIF (Pearson)", color="#6b7280")
    ax.bar(x, raw, w, label="Our RNA raw Spearman", color="#93c5fd")
    ax.bar(x + w, par, w, label="Our RNA partial r | cellularity", color="#1d4ed8")
    ax.set_xticks(x)
    ax.set_xticklabels(channels, rotation=40, ha="right")
    ax.set_ylabel("Correlation")
    ax.set_title("Per-channel: training-modality agreement vs independent RNA specificity")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.asset_dir / ".matplotlib"))
    (args.asset_dir / ".matplotlib").mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rna = load_rna(args.rna_report)
    mif = load_gigatime_mif(args.gigatime_mif_csv)
    rows = build_table(rna, mif)
    write_csv(args.out_csv, rows)
    make_scatter(plt, rows, args.asset_dir / "mif_agreement_vs_rna_specificity_scatter.png")
    make_bars(plt, rows, args.asset_dir / "per_channel_bars.png")

    print(f"GigaTIME mIF source: {GIGATIME_SOURCE}")
    print(f"Channels joined: {sum(1 for r in rows if r['gigatime_mif_pearson'] is not None)}")
    for r in rows:
        if r["gigatime_mif_pearson"] is None:
            continue
        print(f"  {r['channel']:9s} mIF={r['gigatime_mif_pearson']:.2f}  "
              f"RNA_raw={r['rna_raw_spearman']:.3f}  RNA_partial={r['rna_partial_r']:.3f}  "
              f"rowmax={r['own_is_row_max']}")
    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.asset_dir}/mif_agreement_vs_rna_specificity_scatter.png and per_channel_bars.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
