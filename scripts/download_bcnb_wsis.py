#!/usr/bin/env python3
"""Download the gated BCNB full-WSI archive or folder into ignored local data."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import requests


DEFAULT_OUTPUT = Path("data/bcnb/WSIs")
DEFAULT_MANIFEST = Path("data/bcnb/bcnb_onedrive_wsi_manifest.csv")
DEFAULT_SHAREPOINT_SUFFIXES = ".jpg,.jpeg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--url",
        required=True,
        help=(
            "Approved BCNB WSI download URL. Google Drive URLs are handled via gdown. "
            "SharePoint/OneDrive folder links are listed through the SharePoint API and downloaded file-by-file."
        ),
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--filename", default="", help="Filename for direct file URLs. Defaults to URL basename.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="CSV manifest path for SharePoint/OneDrive listings.",
    )
    parser.add_argument(
        "--sharepoint-folder-path",
        default="",
        help="Server-relative SharePoint folder path. Defaults to the shared BCNB folder plus /WSIs.",
    )
    parser.add_argument(
        "--include-suffixes",
        default=DEFAULT_SHAREPOINT_SUFFIXES,
        help="Comma-separated suffixes to download from SharePoint/OneDrive. Defaults to WSI JPEGs only.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=0,
        help="SharePoint/OneDrive files per batch. Use 0 to download all selected files.",
    )
    parser.add_argument(
        "--batch-index",
        type=int,
        default=1,
        help="One-based SharePoint/OneDrive batch index used with --batch-size.",
    )
    parser.add_argument("--max-files", type=int, default=0, help="Maximum SharePoint/OneDrive files to download after batching.")
    parser.add_argument("--manifest-only", action="store_true", help="Only write the SharePoint/OneDrive manifest.")
    parser.add_argument(
        "--allow-nonnumeric-stems",
        action="store_true",
        help="Also download SharePoint files whose stem is not a numeric BCNB patient ID.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Redownload complete existing SharePoint/OneDrive files.")
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-sleep", type=float, default=5.0)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def looks_like_google_drive(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.netloc.endswith("google.com") or parsed.netloc.endswith("googleusercontent.com") or "drive.google.com" in parsed.netloc


def looks_like_sharepoint(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.netloc.endswith(".sharepoint.com") or "sharepoint.com" in parsed.netloc


def run_gdown(url: str, out_dir: Path, quiet: bool) -> int:
    if shutil.which("gdown") is None:
        raise SystemExit(
            "Google Drive URL detected, but gdown is not installed in this environment. "
            "Install it with: conda run -n gigatime-tcga python -m pip install gdown"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    output = str(out_dir)
    if not output.endswith("/"):
        output += "/"
    cmd = ["gdown", "--continue", "--fuzzy", url, "-O", output]
    if "/folders/" in url:
        cmd.insert(1, "--folder")
        cmd.insert(2, "--remaining-ok")
    if quiet:
        cmd.insert(1, "--quiet")
    print("Running: " + " ".join(cmd))
    return subprocess.call(cmd)


def parse_suffixes(raw: str) -> set[str]:
    return {item.strip().lower() for item in raw.split(",") if item.strip()}


def numeric_name_key(row: dict[str, str]) -> tuple[int, int | str, str]:
    stem = Path(row["name"]).stem
    try:
        return (0, int(stem), row["name"])
    except ValueError:
        return (1, stem, row["name"])


def sharepoint_site_and_folder(
    session: requests.Session,
    url: str,
    folder_override: str,
) -> tuple[str, str]:
    response = session.get(url, timeout=60)
    response.raise_for_status()
    parsed = urlparse(response.url)
    site_path = parsed.path.split("/_layouts/", 1)[0]
    if not site_path:
        raise ValueError(f"Could not derive SharePoint site path from redirected URL: {response.url}")
    site_base = f"{parsed.scheme}://{parsed.netloc}{site_path}"

    if folder_override:
        folder_path = folder_override
    else:
        query = parse_qs(parsed.query)
        folder_values = query.get("id")
        if not folder_values:
            raise ValueError(
                "Could not derive SharePoint folder path from the shared URL. "
                "Pass --sharepoint-folder-path explicitly."
            )
        folder_path = folder_values[0]
    if not folder_path.startswith("/"):
        folder_path = "/" + folder_path
    if Path(folder_path).name != "WSIs":
        folder_path = folder_path.rstrip("/") + "/WSIs"
    return site_base, folder_path


def list_sharepoint_files(
    session: requests.Session,
    site_base: str,
    folder_path: str,
) -> list[dict[str, str]]:
    folder_literal = folder_path.replace("'", "''")
    url = f"{site_base}/_api/web/GetFolderByServerRelativeUrl('{folder_literal}')/Files"
    params = {
        "$select": "Name,ServerRelativeUrl,Length,TimeLastModified,UniqueId",
        "$top": "5000",
    }
    headers = {"Accept": "application/json;odata=nometadata"}
    rows: list[dict[str, str]] = []
    while url:
        response = session.get(url, params=params, headers=headers, timeout=60)
        response.raise_for_status()
        data = response.json()
        for item in data.get("value", []):
            rows.append(
                {
                    "name": str(item.get("Name", "")),
                    "server_relative_url": str(item.get("ServerRelativeUrl", "")),
                    "length": str(item.get("Length", "")),
                    "time_last_modified": str(item.get("TimeLastModified", "")),
                    "unique_id": str(item.get("UniqueId", "")),
                }
            )
        url = data.get("odata.nextLink") or data.get("@odata.nextLink")
        params = {}
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["name", "server_relative_url", "length", "time_last_modified", "unique_id"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def local_sharepoint_path(out_dir: Path, server_relative_url: str) -> Path:
    marker = "/Documents/BCNB/"
    if marker in server_relative_url:
        relative = server_relative_url.split(marker, 1)[1]
        return out_dir / "BCNB" / relative
    return out_dir / Path(server_relative_url).name


def expected_length(row: dict[str, str]) -> int | None:
    try:
        return int(row["length"])
    except (KeyError, TypeError, ValueError):
        return None


def download_sharepoint_file(
    *,
    session: requests.Session,
    site_base: str,
    row: dict[str, str],
    destination: Path,
    overwrite: bool,
    quiet: bool,
) -> str:
    expected = expected_length(row)
    destination.parent.mkdir(parents=True, exist_ok=True)
    existing = destination.stat().st_size if destination.exists() else 0
    if overwrite and destination.exists():
        destination.unlink()
        existing = 0
    if expected is not None and existing == expected:
        return "skipped"
    if expected is not None and existing > expected:
        existing = 0

    headers = {}
    if existing:
        headers["Range"] = f"bytes={existing}-"
    url = f"{site_base}/_layouts/15/download.aspx"
    params = {"SourceUrl": row["server_relative_url"]}
    with session.get(url, params=params, headers=headers, stream=True, timeout=(15, 180)) as response:
        response.raise_for_status()
        mode = "ab" if existing and response.status_code == 206 else "wb"
        if mode == "wb":
            existing = 0
        with destination.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    final_size = destination.stat().st_size
    if expected is not None and final_size != expected:
        raise IOError(f"{destination} size mismatch: expected {expected}, got {final_size}")
    if not quiet:
        size_mb = final_size / (1024 * 1024)
        print(f"downloaded {destination} ({size_mb:.1f} MiB)")
    return "downloaded"


def run_sharepoint(args: argparse.Namespace) -> int:
    suffixes = parse_suffixes(args.include_suffixes)
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125 Safari/537.36"
            )
        }
    )
    site_base, folder_path = sharepoint_site_and_folder(session, args.url, args.sharepoint_folder_path)
    rows = list_sharepoint_files(session, site_base, folder_path)
    write_manifest(args.manifest, rows)

    selected = [row for row in rows if Path(row["name"]).suffix.lower() in suffixes]
    selected.sort(key=numeric_name_key)
    total_matching = len(selected)
    nonnumeric_selected = [row for row in selected if not Path(row["name"]).stem.isdigit()]
    if not args.allow_nonnumeric_stems:
        selected = [row for row in selected if Path(row["name"]).stem.isdigit()]
    if args.batch_size:
        if args.batch_index < 1:
            raise ValueError("--batch-index is one-based and must be >= 1")
        start = (args.batch_index - 1) * args.batch_size
        selected = selected[start : start + args.batch_size]
    if args.max_files:
        selected = selected[: args.max_files]

    print(f"SharePoint folder: {folder_path}")
    print(f"Manifest rows: {len(rows)} -> {args.manifest}")
    print(f"Selected files: {len(selected)} of {total_matching} matching {', '.join(sorted(suffixes))}")
    if nonnumeric_selected and not args.allow_nonnumeric_stems:
        examples = ", ".join(row["name"] for row in nonnumeric_selected[:5])
        print(f"Skipped nonnumeric SharePoint files by default: {len(nonnumeric_selected)} ({examples})")
    if args.manifest_only:
        return 0

    counts = {"downloaded": 0, "skipped": 0}
    for index, row in enumerate(selected, start=1):
        destination = local_sharepoint_path(args.out_dir, row["server_relative_url"])
        if not args.quiet:
            expected = expected_length(row)
            size = f"{expected / (1024 * 1024):.1f} MiB" if expected is not None else "unknown size"
            print(f"[{index}/{len(selected)}] {row['name']} -> {destination} ({size})")
        for attempt in range(1, args.max_retries + 1):
            try:
                status = download_sharepoint_file(
                    session=session,
                    site_base=site_base,
                    row=row,
                    destination=destination,
                    overwrite=args.overwrite,
                    quiet=args.quiet,
                )
                counts[status] = counts.get(status, 0) + 1
                break
            except Exception as exc:
                if attempt >= args.max_retries:
                    raise
                print(f"Retrying {row['name']} after error on attempt {attempt}: {exc}", file=sys.stderr)
                time.sleep(args.retry_sleep * attempt)
    print(f"Done. downloaded={counts.get('downloaded', 0)} skipped={counts.get('skipped', 0)}")
    return 0


def download_direct(url: str, out_dir: Path, filename: str, quiet: bool) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    if not filename:
        filename = Path(urlparse(url).path).name or "bcnb_wsi_download"
    destination = out_dir / filename
    headers = {}
    existing = destination.stat().st_size if destination.exists() else 0
    if existing:
        headers["Range"] = f"bytes={existing}-"
    with requests.get(url, stream=True, timeout=60, headers=headers) as response:
        response.raise_for_status()
        mode = "ab" if existing and response.status_code == 206 else "wb"
        if mode == "wb":
            existing = 0
        with destination.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                handle.write(chunk)
                if not quiet:
                    print(f"\r{destination.name}: {destination.stat().st_size / (1024 ** 3):.2f} GiB", end="")
    if not quiet:
        print()
    return destination


def main() -> int:
    args = parse_args()
    if looks_like_sharepoint(args.url):
        return run_sharepoint(args)
    if looks_like_google_drive(args.url):
        return run_gdown(args.url, args.out_dir, args.quiet)
    destination = download_direct(args.url, args.out_dir, args.filename, args.quiet)
    print(f"Wrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
