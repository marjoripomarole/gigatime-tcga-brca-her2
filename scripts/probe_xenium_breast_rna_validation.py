#!/usr/bin/env python3
"""Feasibility probe for validating GigaTIME virtual channels against Xenium breast RNA.

Downloads only the minimal artifacts needed to confirm that the 10x Xenium Human
Breast dataset can serve as a within-slide RNA ground truth for GigaTIME virtual
immune/tumor channels, then reports:

  (a) channel-gene coverage   - are CD3/CD8/PD-L1/CK/Ki67/myeloid marker genes in
      the panel and present in the transcript table?
  (b) alignment transform     - does the H&E<->Xenium affine load, invert, and map
      the H&E frame onto the transcript coordinate extent?
  (c) transcript geometry      - per-channel-gene transcript counts and the x/y
      micron extent the virtual-channel grid would be correlated against.

The full ~9.86 GB outs bundle is intentionally NOT downloaded; this probe pulls
gene_panel.json (~150 KB), he_imagealignment.csv (126 B), and transcripts.parquet
(~670 MB). The 1.43 GB H&E OME-TIFF is metadata-only and optional (--include-he).

Restricted/large data stays under the gitignored data/ tree.

Default sample: Xenium_FFPE_Human_Breast_Cancer_Rep1 (Janesick et al., Nat Commun
2023; 280+33 gene panel, single-cell, post-Xenium H&E with alignment file).
For a CC BY 4.0, broader-panel alternative see the Xenium Prime 5K breast dataset
(pass --base-url / --sample for that bundle once you have its file URLs).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import requests

# GigaTIME virtual channel -> candidate marker genes (first present is reported).
CHANNEL_GENES: dict[str, list[str]] = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "PD-L1": ["CD274"],
    "CK": ["KRT8", "KRT18", "KRT19", "KRT7", "EPCAM"],
    "Ki67": ["MKI67"],
    "myeloid": ["CD68", "CD163", "LYZ", "ITGAX"],
}

# Non-gene control probe prefixes in Xenium feature_name.
CONTROL_PREFIXES = ("NegControlProbe", "NegControlCodeword", "BLANK", "antisense", "UnassignedCodeword", "DeprecatedCodeword")

# Minimal probe artifacts: (suffix, required-for-probe).
PROBE_FILES: list[tuple[str, bool]] = [
    ("_gene_panel.json", True),
    ("_he_imagealignment.csv", True),
    ("_transcripts.parquet", True),
]
HE_FILE = "_he_image.ome.tif"

# Xenium morphology image pixel size (microns/pixel) for the v1 instrument.
DEFAULT_XENIUM_MPP = 0.2125


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default="Xenium_FFPE_Human_Breast_Cancer_Rep1")
    parser.add_argument(
        "--base-url",
        default="https://cf.10xgenomics.com/samples/xenium/1.0.1/{sample}/{sample}",
        help="URL prefix; '{sample}' is substituted. Each artifact is base-url + suffix.",
    )
    parser.add_argument("--out-dir", type=Path, default=None, help="Defaults to data/xenium_breast/<sample>.")
    parser.add_argument("--no-download", action="store_true", help="Probe only existing local files; do not fetch.")
    parser.add_argument("--include-he", action="store_true", help="Also download the 1.43 GB H&E OME-TIFF for the alignment cross-check.")
    parser.add_argument("--xenium-mpp", type=float, default=DEFAULT_XENIUM_MPP, help="Xenium morphology microns/pixel for the H&E extent cross-check.")
    parser.add_argument("--report-json", type=Path, default=None, help="Defaults to <out-dir>/xenium_breast_probe_report.json.")
    parser.add_argument("--out-markdown", type=Path, default=Path("docs/xenium_breast_rna_validation_probe.md"))
    parser.add_argument("--max-retries", type=int, default=4)
    parser.add_argument("--retry-sleep", type=float, default=5.0)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def require_libs():
    try:
        import numpy as np
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Run with the gigatime-tcga env, e.g.\n"
            "  conda run -n gigatime-tcga python scripts/probe_xenium_breast_rna_validation.py\n"
            "or directly: ~/miniconda3/envs/gigatime-tcga/bin/python scripts/probe_xenium_breast_rna_validation.py"
        ) from exc
    return np, pc, pq


def download_file(url: str, dest: Path, *, max_retries: int, retry_sleep: float, quiet: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, max_retries + 1):
        try:
            existing = dest.stat().st_size if dest.exists() else 0
            headers = {"Range": f"bytes={existing}-"} if existing else {}
            with requests.get(url, stream=True, timeout=(15, 300), headers=headers) as response:
                if existing and response.status_code == 200:
                    existing = 0  # server ignored Range; restart
                elif existing and response.status_code == 416:
                    return  # already complete
                response.raise_for_status()
                mode = "ab" if existing and response.status_code == 206 else "wb"
                with dest.open(mode) as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                            if not quiet:
                                gib = dest.stat().st_size / (1024 ** 3)
                                print(f"\r  {dest.name}: {gib:.2f} GiB", end="", file=sys.stderr)
            if not quiet:
                print(file=sys.stderr)
            return
        except Exception as exc:  # noqa: BLE001 - retry any transient network error
            if attempt >= max_retries:
                raise
            print(f"  retry {dest.name} after attempt {attempt}: {exc}", file=sys.stderr)
            time.sleep(retry_sleep * attempt)


def ensure_files(args: argparse.Namespace, out_dir: Path) -> dict[str, Path]:
    base = args.base_url.format(sample=args.sample)
    wanted = list(PROBE_FILES)
    if args.include_he:
        wanted.append((HE_FILE, False))
    paths: dict[str, Path] = {}
    for suffix, _required in wanted:
        dest = out_dir / f"{args.sample}{suffix}"
        paths[suffix] = dest
        if args.no_download:
            continue
        if dest.exists() and dest.stat().st_size > 0:
            if not args.quiet:
                print(f"  present: {dest.name} ({dest.stat().st_size / (1024 ** 2):.1f} MiB)")
            continue
        url = f"{base}{suffix}"
        if not args.quiet:
            print(f"  downloading: {url}")
        download_file(url, dest, max_retries=args.max_retries, retry_sleep=args.retry_sleep, quiet=args.quiet)
    return paths


def extract_panel_genes(panel: object) -> set[str]:
    """Walk the gene_panel.json structure collecting target gene symbol names."""
    genes: set[str] = set()

    def walk(node: object) -> None:
        if isinstance(node, dict):
            data = node.get("data") if isinstance(node.get("data"), dict) else None
            descriptor = node.get("descriptor") or node.get("type")
            if data and isinstance(data.get("name"), str):
                name = data["name"]
                is_gene = True
                if isinstance(descriptor, str):
                    is_gene = "gene" in descriptor.lower()
                if is_gene and not name.startswith(CONTROL_PREFIXES):
                    genes.add(name)
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(panel)
    return genes


def probe_gene_panel(path: Path) -> dict:
    if not path.exists():
        return {"available": False}
    panel = json.loads(path.read_text(encoding="utf-8"))
    genes = extract_panel_genes(panel)
    coverage = {}
    for channel, candidates in CHANNEL_GENES.items():
        present = [g for g in candidates if g in genes]
        coverage[channel] = {"present_genes": present, "covered": bool(present)}
    return {
        "available": True,
        "n_panel_genes": len(genes),
        "channel_coverage": coverage,
        "n_channels_covered": sum(1 for c in coverage.values() if c["covered"]),
        "n_channels_total": len(CHANNEL_GENES),
    }


def probe_alignment(np, path: Path) -> dict:
    if not path.exists():
        return {"available": False}
    raw = path.read_text(encoding="utf-8").strip()
    rows = []
    for line in raw.splitlines():
        parts = [p for p in line.replace(",", " ").split() if p]
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue  # skip any header
    matrix = np.array(rows, dtype=float)
    result: dict = {"available": True, "shape": list(matrix.shape), "matrix": matrix.tolist()}
    if matrix.shape == (3, 3):
        det = float(np.linalg.det(matrix))
        result["determinant"] = det
        result["invertible"] = abs(det) > 1e-9
    else:
        result["invertible"] = False
        result["note"] = "Expected a 3x3 affine; got an unexpected shape."
    return result


def probe_transcripts(pq, pc, path: Path) -> dict:
    if not path.exists():
        return {"available": False}
    schema = pq.read_schema(path)
    names = set(schema.names)
    feature_col = "feature_name" if "feature_name" in names else next((n for n in names if "feature" in n.lower()), None)
    x_col = "x_location" if "x_location" in names else next((n for n in names if n.lower().startswith("x")), None)
    y_col = "y_location" if "y_location" in names else next((n for n in names if n.lower().startswith("y")), None)
    read_cols = [c for c in [feature_col, x_col, y_col] if c]
    table = pq.read_table(path, columns=read_cols)
    n_rows = table.num_rows

    counts_by_gene: dict[str, int] = {}
    if feature_col:
        vc = pc.value_counts(table[feature_col])
        for entry in vc:
            value = entry["values"].as_py()
            count = entry["counts"].as_py()
            if value is None:
                continue
            if isinstance(value, bytes):
                value = value.decode("utf-8", "replace")
            counts_by_gene[str(value)] = int(count)
    genes_present = {g for g in counts_by_gene if not g.startswith(CONTROL_PREFIXES)}
    control_transcripts = sum(c for g, c in counts_by_gene.items() if g.startswith(CONTROL_PREFIXES))

    extent = {}
    for axis, col in (("x", x_col), ("y", y_col)):
        if col:
            mm = pc.min_max(table[col]).as_py()
            extent[axis] = {"min": float(mm["min"]), "max": float(mm["max"])}

    channel_counts = {}
    for channel, candidates in CHANNEL_GENES.items():
        per = {g: counts_by_gene.get(g, 0) for g in candidates if g in genes_present}
        channel_counts[channel] = {"transcript_counts": per, "total": int(sum(per.values()))}

    return {
        "available": True,
        "columns": list(schema.names),
        "feature_column": feature_col,
        "n_transcripts": int(n_rows),
        "n_genes_detected": len(genes_present),
        "n_control_transcripts": int(control_transcripts),
        "extent_microns": extent,
        "channel_transcript_counts": channel_counts,
    }


def probe_he_metadata(path: Path) -> dict:
    if not path.exists():
        return {"available": False, "note": "H&E not downloaded (use --include-he)."}
    try:
        import tifffile
    except ModuleNotFoundError:
        return {"available": False, "note": "tifffile not installed."}
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        levels = []
        try:
            for lvl in series.levels:
                levels.append([int(s) for s in lvl.shape])
        except Exception:  # noqa: BLE001 - single-level fallback
            levels.append([int(s) for s in series.shape])
        page = tif.pages[0]
        mpp = None
        try:
            tags = page.tags
            if "XResolution" in tags and "ResolutionUnit" in tags:
                num, den = tags["XResolution"].value
                if num:
                    px_per_unit = num / den
                    unit = tags["ResolutionUnit"].value
                    if str(unit).endswith("CENTIMETER") or unit == 3:
                        mpp = 10000.0 / px_per_unit
                    elif str(unit).endswith("INCH") or unit == 2:
                        mpp = 25400.0 / px_per_unit
        except Exception:  # noqa: BLE001
            mpp = None
    full = levels[0] if levels else None
    return {
        "available": True,
        "full_resolution_shape": full,
        "n_pyramid_levels": len(levels),
        "level_shapes": levels,
        "microns_per_pixel_tag": mpp,
    }


def cross_check_alignment(np, alignment: dict, he: dict, transcripts: dict, xenium_mpp: float) -> dict:
    if not (alignment.get("invertible") and he.get("available") and transcripts.get("available")):
        return {"performed": False, "reason": "needs invertible alignment, downloaded H&E, and transcripts."}
    shape = he.get("full_resolution_shape")
    extent = transcripts.get("extent_microns", {})
    if not shape or len(shape) < 2 or "x" not in extent or "y" not in extent:
        return {"performed": False, "reason": "missing H&E shape or transcript extent."}
    matrix = np.array(alignment["matrix"], dtype=float)
    height, width = shape[0], shape[1]
    corners = np.array([[0, 0, 1], [width, 0, 1], [0, height, 1], [width, height, 1]], dtype=float)
    mapped = corners @ matrix.T  # H&E px -> Xenium morphology px (homogeneous)
    mapped = mapped[:, :2] / mapped[:, 2:3]
    mapped_um = mapped * xenium_mpp  # morphology px -> microns
    he_x = (float(mapped_um[:, 0].min()), float(mapped_um[:, 0].max()))
    he_y = (float(mapped_um[:, 1].min()), float(mapped_um[:, 1].max()))
    tx = (extent["x"]["min"], extent["x"]["max"])
    ty = (extent["y"]["min"], extent["y"]["max"])

    def overlap_frac(a: tuple[float, float], b: tuple[float, float]) -> float:
        inter = max(0.0, min(a[1], b[1]) - max(a[0], b[0]))
        union = max(a[1], b[1]) - min(a[0], b[0])
        return inter / union if union > 0 else 0.0

    ox, oy = overlap_frac(he_x, tx), overlap_frac(he_y, ty)
    return {
        "performed": True,
        "xenium_mpp": xenium_mpp,
        "he_extent_microns": {"x": he_x, "y": he_y},
        "transcript_extent_microns": {"x": tx, "y": ty},
        "extent_overlap_fraction": {"x": ox, "y": oy},
        "applies_cleanly": ox > 0.5 and oy > 0.5,
    }


def verdict(panel: dict, alignment: dict, transcripts: dict, cross: dict) -> dict:
    gene_ok = panel.get("available") and panel.get("n_channels_covered", 0) >= panel.get("n_channels_total", 99)
    align_ok = bool(alignment.get("invertible"))
    tx_ok = transcripts.get("available") and all(
        c["total"] > 0 for ch, c in transcripts.get("channel_transcript_counts", {}).items()
    )
    cross_ok = cross.get("applies_cleanly") if cross.get("performed") else None
    feasible = bool(gene_ok and align_ok and tx_ok and (cross_ok is not False))
    return {
        "all_channels_have_panel_gene": bool(gene_ok),
        "alignment_invertible": align_ok,
        "all_channels_have_transcripts": bool(tx_ok),
        "alignment_cross_check": cross_ok,
        "feasible_for_rna_validation": feasible,
    }


def fmt_int(n) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def write_markdown(path: Path, args, report: dict) -> None:
    panel = report["gene_panel"]
    align = report["alignment"]
    tx = report["transcripts"]
    he = report["he_image"]
    cross = report["alignment_cross_check"]
    v = report["verdict"]

    lines = [
        "# Xenium Breast RNA-Validation Feasibility Probe",
        "",
        f"Status: feasibility probe for validating GigaTIME virtual channels against Xenium breast RNA. Sample `{args.sample}`.",
        "",
        "## Purpose",
        "",
        "Confirm that the 10x Xenium Human Breast dataset can serve as a within-slide RNA ground truth for",
        "GigaTIME virtual immune/tumor channels before wiring GigaTIME inference to it. Within-slide transcript",
        "co-localization sidesteps the TCGA composition/batch confound, because the correlation is computed inside",
        "one section rather than across patients.",
        "",
        "## Verdict",
        "",
        f"- Feasible for RNA validation: **{'YES' if v['feasible_for_rna_validation'] else 'NO'}**",
        f"- Every GigaTIME channel maps to a panel gene: {v['all_channels_have_panel_gene']}",
        f"- Every channel gene has transcripts: {v['all_channels_have_transcripts']}",
        f"- H&E<->Xenium alignment invertible: {v['alignment_invertible']}",
        f"- Alignment extent cross-check: {v['alignment_cross_check']}",
        "",
        "## Channel-Gene Coverage",
        "",
        f"Panel genes: {fmt_int(panel.get('n_panel_genes'))}. Channels covered: "
        f"{panel.get('n_channels_covered')} / {panel.get('n_channels_total')}.",
        "",
        "| GigaTIME channel | Panel gene(s) present | Transcript count |",
        "|---|---|---:|",
    ]
    for channel in CHANNEL_GENES:
        present = panel.get("channel_coverage", {}).get(channel, {}).get("present_genes", [])
        total = tx.get("channel_transcript_counts", {}).get(channel, {}).get("total", 0)
        lines.append(f"| {channel} | {', '.join(present) if present else '—'} | {fmt_int(total)} |")

    lines += [
        "",
        "## Transcripts",
        "",
        f"- Total transcripts: {fmt_int(tx.get('n_transcripts'))} across {fmt_int(tx.get('n_genes_detected'))} detected genes "
        f"({fmt_int(tx.get('n_control_transcripts'))} control-probe transcripts).",
        f"- Feature column: `{tx.get('feature_column')}`.",
    ]
    extent = tx.get("extent_microns", {})
    if "x" in extent and "y" in extent:
        lines.append(
            f"- Coordinate extent (microns): x [{extent['x']['min']:.1f}, {extent['x']['max']:.1f}], "
            f"y [{extent['y']['min']:.1f}, {extent['y']['max']:.1f}]."
        )

    lines += ["", "## Alignment", ""]
    if align.get("available"):
        lines.append(f"- Matrix shape: {align.get('shape')}, invertible: {align.get('invertible')}, "
                     f"determinant: {align.get('determinant', 'n/a')}.")
        lines.append("- Matrix (H&E pixel -> Xenium morphology pixel, homogeneous):")
        lines.append("")
        lines.append("```")
        for row in align.get("matrix", []):
            lines.append("  ".join(f"{val: .6g}" for val in row))
        lines.append("```")
    else:
        lines.append("- Alignment file not available.")

    lines += ["", "## H&E Image", ""]
    if he.get("available"):
        lines.append(f"- Full-resolution shape: {he.get('full_resolution_shape')}, "
                     f"pyramid levels: {he.get('n_pyramid_levels')}, "
                     f"microns/pixel tag: {he.get('microns_per_pixel_tag')}.")
    else:
        lines.append(f"- {he.get('note', 'H&E metadata not available.')}")
    if cross.get("performed"):
        lines.append(
            f"- Cross-check: H&E frame maps to transcript extent with overlap "
            f"x={cross['extent_overlap_fraction']['x']:.2f}, y={cross['extent_overlap_fraction']['y']:.2f} "
            f"(applies_cleanly={cross['applies_cleanly']}, xenium_mpp={cross['xenium_mpp']})."
        )

    lines += [
        "",
        "## Next Steps",
        "",
        "1. If feasible: tile the post-Xenium H&E at GigaTIME's expected microns/pixel and run virtual-channel inference.",
        "2. Bin transcripts to the same tile grid via the alignment transform; sum each channel gene per tile.",
        "3. Within-slide Spearman correlation of virtual-channel intensity vs transcript density per tile, per channel,",
        "   with a block-bootstrap CI over tiles. Lead on CD8A/CD3D/MKI67/keratins; report CD274/PD-L1 with the known",
        "   RNA-protein concordance caveat.",
        "4. For multi-slide breadth, repeat across Xenium breast replicates and/or HEST-1k breast Visium samples.",
        "",
        "## Sources",
        "",
        "- Janesick et al., Nat Commun 2023 (Xenium FFPE Human Breast Cancer).",
        "- 10x Xenium Human Breast Dataset Explorer: https://www.10xgenomics.com/products/xenium-in-situ/human-breast-dataset-explorer",
        "- CC BY 4.0 broader-panel alternative: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    np, pc, pq = require_libs()
    out_dir = args.out_dir or Path("data/xenium_breast") / args.sample
    out_dir.mkdir(parents=True, exist_ok=True)
    report_json = args.report_json or out_dir / "xenium_breast_probe_report.json"

    if not args.quiet:
        print(f"Sample: {args.sample}")
        print(f"Local dir: {out_dir}")
    paths = ensure_files(args, out_dir)

    panel = probe_gene_panel(paths.get("_gene_panel.json", out_dir / "missing"))
    alignment = probe_alignment(np, paths.get("_he_imagealignment.csv", out_dir / "missing"))
    transcripts = probe_transcripts(pq, pc, paths.get("_transcripts.parquet", out_dir / "missing"))
    he = probe_he_metadata(paths.get(HE_FILE, out_dir / f"{args.sample}{HE_FILE}"))
    cross = cross_check_alignment(np, alignment, he, transcripts, args.xenium_mpp)
    v = verdict(panel, alignment, transcripts, cross)

    report = {
        "sample": args.sample,
        "base_url": args.base_url.format(sample=args.sample),
        "local_dir": str(out_dir),
        "channel_gene_map": CHANNEL_GENES,
        "gene_panel": panel,
        "alignment": alignment,
        "transcripts": transcripts,
        "he_image": he,
        "alignment_cross_check": cross,
        "verdict": v,
    }
    report_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(args.out_markdown, args, report)

    if not args.quiet:
        print(f"\nChannels with a panel gene: {panel.get('n_channels_covered')}/{panel.get('n_channels_total')}")
        print(f"Transcripts: {fmt_int(transcripts.get('n_transcripts'))} over {fmt_int(transcripts.get('n_genes_detected'))} genes")
        print(f"Alignment invertible: {alignment.get('invertible')}")
        print(f"FEASIBLE for RNA validation: {v['feasible_for_rna_validation']}")
        print(f"Wrote {report_json}")
        print(f"Wrote {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
