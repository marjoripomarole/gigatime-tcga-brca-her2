#!/usr/bin/env python3
"""Audit BCNB WSI or patch files before launching image-model runs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_LABELS = Path("data/bcnb/bcnb_her2_labels.csv")
DEFAULT_WSI_DIR = Path("data/bcnb/WSIs")
DEFAULT_PATCH_ZIP = Path("data/bcnb/paper_patches.zip")
DEFAULT_PATCH_DIR = Path("data/bcnb/paper_patches")

WSI_SUFFIXES = {".svs", ".tif", ".tiff", ".ndpi", ".mrxs", ".jpg", ".jpeg"}
PATCH_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
METADATA_SUFFIXES = {".csv", ".tsv", ".xlsx", ".xls", ".json", ".txt", ".yaml", ".yml"}


@dataclass
class MatchResult:
    patient_id: str
    confidence: str
    reason: str


@dataclass
class SourceAudit:
    source: str
    path: str
    exists: bool
    total_files_seen: int = 0
    image_files_seen: int = 0
    metadata_files_seen: int = 0
    mapped_image_files: int = 0
    mapped_patients: int = 0
    ambiguous_image_files: int = 0
    unmatched_image_files: int = 0
    group_patient_counts: dict[str, int] | None = None
    suffix_counts: dict[str, int] | None = None
    metadata_examples: list[str] | None = None
    matched_examples: list[str] | None = None
    ambiguous_examples: list[str] | None = None
    unmatched_examples: list[str] | None = None
    notes: list[str] | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS, help="BCNB label CSV built from the clinical file.")
    parser.add_argument("--wsi-dir", type=Path, default=DEFAULT_WSI_DIR, help="Directory containing BCNB WSI files.")
    parser.add_argument("--patch-zip", type=Path, default=DEFAULT_PATCH_ZIP, help="BCNB paper_patches.zip path.")
    parser.add_argument("--patch-dir", type=Path, default=DEFAULT_PATCH_DIR, help="Extracted BCNB paper patches directory.")
    parser.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional path to write the machine-readable audit summary.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Maximum files to inspect per source. Use 0 for all files.",
    )
    return parser.parse_args()


def load_labels(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"BCNB label file not found: {path}")
    labels: dict[str, str] = {}
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"patient_id", "clinical_her2_group"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        for row in reader:
            patient_id = (row.get("patient_id") or "").strip()
            group = (row.get("clinical_her2_group") or "").strip()
            if patient_id:
                labels[patient_id] = group
    return labels


def normalize_patient_token(token: str) -> str:
    stripped = token.strip()
    if not stripped:
        return ""
    try:
        return str(int(stripped))
    except ValueError:
        return stripped.lstrip("0") or "0"


def named_patient_match(text: str, known_ids: set[str]) -> MatchResult | None:
    patterns = [
        r"(?:patient|case|sample|slide|bcnb|p)[_-]?0*(\d{1,4})(?:\D|$)",
        r"(?:^|/|\\)0*(\d{1,4})(?:/|\\|$)",
    ]
    for pattern in patterns:
        for match in re.finditer(pattern, text, flags=re.IGNORECASE):
            candidate = normalize_patient_token(match.group(1))
            if candidate in known_ids:
                return MatchResult(candidate, "high", f"matched pattern {pattern}")
    return None


def infer_patient_id(path_text: str, known_ids: set[str]) -> MatchResult | None:
    normalized_text = path_text.replace("\\", "/")
    named_match = named_patient_match(normalized_text, known_ids)
    if named_match:
        return named_match

    parts = [part for part in normalized_text.split("/") if part]
    for part in parts[:-1]:
        if part.isdigit():
            candidate = normalize_patient_token(part)
            if candidate in known_ids:
                return MatchResult(candidate, "high", "matched numeric directory name")

    stem = Path(parts[-1] if parts else normalized_text).stem
    start_match = re.match(r"0*(\d{1,4})(?:[_\-.]|$)", stem)
    if start_match:
        candidate = normalize_patient_token(start_match.group(1))
        if candidate in known_ids:
            stem_candidates = {
                normalize_patient_token(token)
                for token in re.findall(r"\d{1,4}", stem)
                if normalize_patient_token(token) in known_ids
            }
            if len(stem_candidates) > 1:
                return MatchResult(
                    "",
                    "ambiguous",
                    f"leading filename token plus other known ids appeared: {', '.join(sorted(stem_candidates)[:8])}",
                )
            return MatchResult(candidate, "medium", "matched leading numeric filename token")

    candidates = {
        normalize_patient_token(token)
        for token in re.findall(r"\d{1,4}", normalized_text)
        if normalize_patient_token(token) in known_ids
    }
    if len(candidates) == 1:
        return MatchResult(next(iter(candidates)), "low", "only one known patient id appeared in path")
    if len(candidates) > 1:
        return MatchResult("", "ambiguous", f"multiple known patient ids appeared: {', '.join(sorted(candidates)[:8])}")
    return None


def limited(iterable: Iterable[str], max_files: int) -> Iterable[str]:
    if max_files <= 0:
        yield from iterable
        return
    for index, value in enumerate(iterable, start=1):
        if index > max_files:
            break
        yield value


def audit_paths(
    *,
    source: str,
    path: Path,
    entries: Iterable[str],
    image_suffixes: set[str],
    labels: dict[str, str],
    max_files: int,
) -> SourceAudit:
    known_ids = set(labels)
    suffix_counts: Counter[str] = Counter()
    metadata_examples: list[str] = []
    matched_examples: list[str] = []
    ambiguous_examples: list[str] = []
    unmatched_examples: list[str] = []
    mapped_patients: defaultdict[str, Counter[str]] = defaultdict(Counter)
    mapped_image_files = 0
    ambiguous_image_files = 0
    unmatched_image_files = 0
    image_files_seen = 0
    metadata_files_seen = 0
    total_files_seen = 0

    for entry in limited(entries, max_files):
        total_files_seen += 1
        suffix = Path(entry).suffix.lower()
        suffix_counts[suffix or "<none>"] += 1

        if suffix in METADATA_SUFFIXES:
            metadata_files_seen += 1
            if len(metadata_examples) < 10:
                metadata_examples.append(entry)

        if suffix not in image_suffixes:
            continue

        image_files_seen += 1
        match = infer_patient_id(entry, known_ids)
        if match is None:
            unmatched_image_files += 1
            if len(unmatched_examples) < 10:
                unmatched_examples.append(entry)
            continue
        if match.confidence == "ambiguous":
            ambiguous_image_files += 1
            if len(ambiguous_examples) < 10:
                ambiguous_examples.append(f"{entry} ({match.reason})")
            continue

        mapped_image_files += 1
        mapped_patients[match.patient_id][labels[match.patient_id]] += 1
        if len(matched_examples) < 10:
            matched_examples.append(f"{entry} -> patient {match.patient_id} ({match.confidence})")

    group_patient_counts = Counter()
    for patient_counts in mapped_patients.values():
        group_patient_counts.update(patient_counts.keys())

    notes: list[str] = []
    if max_files > 0 and total_files_seen >= max_files:
        notes.append(f"Stopped after --max-files={max_files}; counts may be incomplete.")
    if image_files_seen and not mapped_patients:
        notes.append("Image files were present, but no patient IDs mapped unambiguously to the label table.")
    if mapped_patients:
        notes.append("Patient-ID mapping exists; inspect examples before trusting patch-level analyses.")

    return SourceAudit(
        source=source,
        path=str(path),
        exists=True,
        total_files_seen=total_files_seen,
        image_files_seen=image_files_seen,
        metadata_files_seen=metadata_files_seen,
        mapped_image_files=mapped_image_files,
        mapped_patients=len(mapped_patients),
        ambiguous_image_files=ambiguous_image_files,
        unmatched_image_files=unmatched_image_files,
        group_patient_counts=dict(sorted(group_patient_counts.items())),
        suffix_counts=dict(sorted(suffix_counts.items())),
        metadata_examples=metadata_examples,
        matched_examples=matched_examples,
        ambiguous_examples=ambiguous_examples,
        unmatched_examples=unmatched_examples,
        notes=notes,
    )


def missing_source(source: str, path: Path) -> SourceAudit:
    return SourceAudit(
        source=source,
        path=str(path),
        exists=False,
        group_patient_counts={},
        suffix_counts={},
        metadata_examples=[],
        matched_examples=[],
        ambiguous_examples=[],
        unmatched_examples=[],
        notes=["Path is not present locally."],
    )


def audit_directory(source: str, path: Path, image_suffixes: set[str], labels: dict[str, str], max_files: int) -> SourceAudit:
    if not path.exists():
        return missing_source(source, path)
    if not path.is_dir():
        audit = missing_source(source, path)
        audit.exists = True
        audit.notes = ["Path exists but is not a directory."]
        return audit
    entries = (str(file.relative_to(path)) for file in path.rglob("*") if file.is_file())
    return audit_paths(source=source, path=path, entries=entries, image_suffixes=image_suffixes, labels=labels, max_files=max_files)


def audit_zip(source: str, path: Path, image_suffixes: set[str], labels: dict[str, str], max_files: int) -> SourceAudit:
    if not path.exists():
        return missing_source(source, path)
    try:
        with zipfile.ZipFile(path) as archive:
            entries = [info.filename for info in archive.infolist() if not info.is_dir()]
    except zipfile.BadZipFile as exc:
        audit = missing_source(source, path)
        audit.exists = True
        audit.notes = [f"Zip file could not be read: {exc}"]
        return audit
    return audit_paths(source=source, path=path, entries=entries, image_suffixes=image_suffixes, labels=labels, max_files=max_files)


def print_report(labels: dict[str, str], audits: list[SourceAudit]) -> None:
    label_counts = Counter(labels.values())
    print("BCNB image input audit")
    print(f"Known patients in label table: {len(labels)}")
    print(
        "Label groups: "
        + ", ".join(f"{group}={label_counts[group]}" for group in ["HER2-zero", "HER2-low", "HER2-positive"])
    )
    print()

    for audit in audits:
        print(f"{audit.source}: {audit.path}")
        print(f"  exists: {audit.exists}")
        if not audit.exists:
            print("  note: path is not present locally")
            print()
            continue
        print(f"  total files seen: {audit.total_files_seen}")
        print(f"  image files seen: {audit.image_files_seen}")
        print(f"  metadata files seen: {audit.metadata_files_seen}")
        print(f"  mapped image files: {audit.mapped_image_files}")
        print(f"  mapped patients: {audit.mapped_patients}")
        print(f"  ambiguous image files: {audit.ambiguous_image_files}")
        print(f"  unmatched image files: {audit.unmatched_image_files}")
        if audit.group_patient_counts:
            groups = ", ".join(f"{group}={count}" for group, count in audit.group_patient_counts.items())
            print(f"  mapped patient groups: {groups}")
        if audit.metadata_examples:
            print("  metadata examples:")
            for example in audit.metadata_examples[:5]:
                print(f"    - {example}")
        if audit.matched_examples:
            print("  matched examples:")
            for example in audit.matched_examples[:5]:
                print(f"    - {example}")
        if audit.ambiguous_examples:
            print("  ambiguous examples:")
            for example in audit.ambiguous_examples[:5]:
                print(f"    - {example}")
        if audit.unmatched_examples:
            print("  unmatched examples:")
            for example in audit.unmatched_examples[:5]:
                print(f"    - {example}")
        if audit.notes:
            print("  notes:")
            for note in audit.notes:
                print(f"    - {note}")
        print()


def main() -> None:
    args = parse_args()
    labels = load_labels(args.labels)
    audits = [
        audit_directory("WSI directory", args.wsi_dir, WSI_SUFFIXES, labels, args.max_files),
        audit_zip("Patch zip", args.patch_zip, PATCH_SUFFIXES, labels, args.max_files),
        audit_directory("Extracted patch directory", args.patch_dir, PATCH_SUFFIXES, labels, args.max_files),
    ]

    print_report(labels, audits)

    if args.report_json:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {"labels": {"n": len(labels), "groups": dict(Counter(labels.values()))}, "sources": [asdict(a) for a in audits]}
        args.report_json.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
