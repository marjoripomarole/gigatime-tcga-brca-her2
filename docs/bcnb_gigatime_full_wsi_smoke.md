# BCNB GigaTIME Full-WSI Smoke

Date: 2026-06-05

Status: full BCNB WSI download and GigaTIME full-WSI smoke completed locally.

## Data Access

BCNB/CVSM provided approved dataset mirrors by email:

- Google Drive: `https://drive.google.com/drive/folders/1HcAgplKwbSZ7ZZl2m6PZdvVF70QJmVuR?usp=sharing`
- OneDrive/SharePoint: `https://bupteducn-my.sharepoint.com/:f:/g/personal/tangwenqi_bupt_edu_cn/EoFPuyWfkg1OpYw_QKK7aBUBKBGSazefo4qlRnaBnCtmKA?e=XPQJqh`
- Aliyun Drive and Baidu Yun mirrors were also provided.

The dataset license remains the BCNB/BALNMP non-commercial license:
`https://github.com/bupt-ai-cz/BALNMP#license`

The Google Drive folder path was brittle for a large folder download: `gdown --folder`
downloaded early files but failed on a later JSON and also hit folder page-size limits.
The robust route was the OneDrive/SharePoint mirror, using SharePoint API listing plus
resumable file-by-file downloads.

## Download

Command:

```bash
conda run -n gigatime-tcga python scripts/download_bcnb_wsis.py \
  --url 'https://bupteducn-my.sharepoint.com/:f:/g/personal/tangwenqi_bupt_edu_cn/EoFPuyWfkg1OpYw_QKK7aBUBKBGSazefo4qlRnaBnCtmKA?e=XPQJqh' \
  --manifest data/bcnb/bcnb_onedrive_wsi_manifest.csv \
  --quiet
```

Download result:

- SharePoint folder: `/personal/tangwenqi_bupt_edu_cn/Documents/BCNB/WSIs`
- Manifest rows: 2,117
- Matching JPG/JPEG files: 1,059
- Numeric patient JPG files selected by default: 1,058
- Nonnumeric JPG skipped: `47594f5f-d952-432a-aff3-e314b3d285e1.jpg`
- Local WSI directory: `data/bcnb/WSIs/BCNB/WSIs/`
- Final local image count: 1,058 `.jpg` WSIs
- Local WSI footprint: approximately 35 GB

The downloader retried transient SharePoint errors and finished with exact-size checks:

```text
Done. downloaded=1033 skipped=25
```

## Audit And Slide Table

Audit command:

```bash
conda run -n gigatime-tcga python scripts/audit_bcnb_image_inputs.py \
  --wsi-dir data/bcnb/WSIs/BCNB/WSIs \
  --report-json data/bcnb/bcnb_image_inputs_after_wsi_download.json
```

Audit result:

- Known patients in label table: 1,058
- WSI image files: 1,058
- Mapped patients: 1,058
- Mapped patient groups: HER2-low 654, HER2-positive 277, HER2-zero 127
- Ambiguous image files: 0
- Unmatched image files: 0

Low/zero slide table command:

```bash
conda run -n gigatime-tcga python scripts/build_bcnb_wsi_slide_table.py \
  --wsi-dir data/bcnb/WSIs/BCNB/WSIs \
  --groups HER2-zero,HER2-low \
  --output data/bcnb/bcnb_wsi_slide_table_low_zero.csv \
  --require-all-patients
```

Slide table result:

- Output rows: 781
- HER2-low: 654
- HER2-zero: 127
- Missing selected patients: 0

## GigaTIME Smoke Runs

The runner now supports BCNB flat JPG WSIs through `--slide-backend pil` and copies
selected clinical metadata from the input slide table into `slide_scores.csv`.

First full-WSI smoke, first 20 table rows:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/bcnb/bcnb_wsi_slide_table_low_zero.csv \
  --slide-path-column slide_local_path \
  --missing-slide-policy error \
  --out-dir results/gigatime_bcnb_wsi_low_zero_smoke20_tile64_meta \
  --tile-limit 64 \
  --tile-order random \
  --random-seed 42 \
  --batch-size 8 \
  --device auto \
  --slide-backend pil \
  --save-tile-csv \
  --max-slides 20 \
  --heatmap-channels CD3,CD8,PD-L1,CK \
  --resume
```

This run completed on MPS:

- Slide rows: 20
- Tile rows: 1,280
- Tiles per slide: 64
- Heatmaps: 80
- Backend: `pil`
- Clinical metadata present in `slide_scores.csv`
- Label distribution: all 20 were HER2-low, because this was the first 20 rows of the sorted table.

Balanced low/zero smoke table:

```text
data/bcnb/bcnb_wsi_slide_table_low_zero_balanced20.csv
```

This table contains deterministic random-seed-42 samples: 10 HER2-low and 10 HER2-zero cases.

Balanced GigaTIME smoke:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/bcnb/bcnb_wsi_slide_table_low_zero_balanced20.csv \
  --slide-path-column slide_local_path \
  --missing-slide-policy error \
  --out-dir results/gigatime_bcnb_wsi_low_zero_balanced20_tile64 \
  --tile-limit 64 \
  --tile-order random \
  --random-seed 42 \
  --batch-size 8 \
  --device auto \
  --slide-backend pil \
  --save-tile-csv \
  --heatmap-channels CD3,CD8,PD-L1,CK \
  --resume
```

Balanced run result:

- Slide rows: 20
- Tile rows: 1,280
- Tiles per slide: 64
- Heatmaps: 80
- HER2-low: 10
- HER2-zero: 10
- IHC distribution: 1+ = 5, 2+ = 5, 0 = 10
- Backend: `pil`

Small descriptive means from the balanced smoke, not inferential:

| Metric | HER2-low | HER2-zero |
|---|---:|---:|
| `mean_CD3` | 0.1101 | 0.1480 |
| `mean_CD8` | 0.0136 | 0.0208 |
| `mean_PD-L1` | 0.1183 | 0.1254 |
| `mean_CK` | 0.0660 | 0.0883 |
| `mean_Ki67` | 0.0095 | 0.0120 |

Interpretation: the full-WSI GigaTIME path works on BCNB JPG WSIs and now produces
slide-level outputs that are directly joinable to HER2 group and clinical covariates.
The smoke means are only a pipeline sanity check; the real analysis should use the
full 781-slide low/zero table with clinical and tissue/tile controls.

## Next Full-Cohort Command

When ready to launch the full low/zero WSI run:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/bcnb/bcnb_wsi_slide_table_low_zero.csv \
  --slide-path-column slide_local_path \
  --missing-slide-policy error \
  --out-dir results/gigatime_bcnb_wsi_low_zero_tile64 \
  --tile-limit 64 \
  --tile-order random \
  --random-seed 42 \
  --batch-size 8 \
  --device auto \
  --slide-backend pil \
  --save-tile-csv \
  --heatmap-channels CD3,CD8,PD-L1,CK \
  --resume
```

This would process 781 slides at 64 tiles per slide, or about 49,984 tile inferences.
Use `--resume` so interrupted runs keep completed slide rows.
