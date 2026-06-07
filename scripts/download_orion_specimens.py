#!/usr/bin/env python3
"""Download Orion CRC specimens (registered H&E + single-cell table + markers) from the public S3 bucket.

Deliberately SKIPS the 60-150 GB raw CyCIF image: the single-cell table (per-cell marker intensities +
X/Y centroids) is enough to build cell-level mIF targets, so each specimen is only ~2 GB.

Bucket: s3://lin-2023-orion-crc/data (anonymous / public).
Run: ~/miniconda3/envs/gigatime-tcga/bin/python scripts/download_orion_specimens.py CRC02 CRC03 CRC04 CRC05 CRC06
"""

import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

BUCKET = "https://lin-2023-orion-crc.s3.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def list_keys(prefix: str):
    xml = urllib.request.urlopen(f"{BUCKET}/?list-type=2&prefix={prefix}&max-keys=100").read()
    root = ET.fromstring(xml)
    return [(c.find(NS + "Key").text, int(c.find(NS + "Size").text)) for c in root.findall(NS + "Contents")]


def fetch(key: str, dest: str):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        print(f"  exists {dest}", flush=True)
        return
    print(f"  downloading {key} -> {dest}", flush=True)
    urllib.request.urlretrieve(f"{BUCKET}/{key}", dest)


def main(specs):
    for spec in specs:
        keys = list_keys(f"data/{spec}/")
        he = [k for k, _ in keys if k.endswith("-registered.ome.tif")]
        cells = [k for k, _ in keys if k.endswith(".csv") and not k.endswith("markers.csv")]
        markers = [k for k, _ in keys if k.endswith("markers.csv")]
        print(f"== {spec}: he={bool(he)} cells={bool(cells)} ==", flush=True)
        if he:
            fetch(he[0], f"data/orion_crc/{spec}/he.ome.tif")
        if cells:
            fetch(cells[0], f"data/orion_crc/{spec}/cells.csv")
        if markers:
            fetch(markers[0], f"data/orion_crc/{spec}/markers.csv")
        print(f"== {spec} done ==", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["CRC02", "CRC03", "CRC04", "CRC05", "CRC06"])
