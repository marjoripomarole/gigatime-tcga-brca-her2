# Clinical HER2 Cohort Selection

This file summarizes the balanced clinical HER2-positive / HER2-low / HER2-zero pilot cohort selected for the next GigaTIME run.

## Counts

| Cohort group | Candidate cases | Selected cases | Already-downloaded slides |
|---|---:|---:|---:|
| HER2-positive | 13 | 10 | 4 |
| HER2-low | 35 | 10 | 3 |
| HER2-zero | 10 | 10 | 1 |

## Selection Priority

- Clinical HER2 group in requested groups.
- Direct clinical label before inferred label.
- Already-downloaded slide before not-yet-downloaded slide.
- Smaller slide file before larger slide file.
- Case submitter ID for deterministic tie-breaking.

## Selected Cases

| Group | Rank | Case | HER2 rule | IHC score | ISH status | ER | PR | ERBB2 TPM | Slide downloaded |
|---|---:|---|---|---|---|---|---|---:|---|
| HER2-low | 1 | TCGA-A2-A0CT | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Negative | 257.9 | yes |
| HER2-low | 2 | TCGA-A2-A04T | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 74.48 | yes |
| HER2-low | 3 | TCGA-A2-A04Q | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 35.86 | yes |
| HER2-low | 4 | TCGA-A2-A0SV | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 110.2 | no |
| HER2-low | 5 | TCGA-A1-A0SJ | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 157.4 | no |
| HER2-low | 6 | TCGA-5L-AAT0 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 156.8 | no |
| HER2-low | 7 | TCGA-A2-A0ES | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | 151.7 | no |
| HER2-low | 8 | TCGA-5T-A9QA | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Negative | 216.2 | no |
| HER2-low | 9 | TCGA-A2-A0EN | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 92.41 | no |
| HER2-low | 10 | TCGA-A2-A0T6 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 252.8 | no |
| HER2-positive | 1 | TCGA-A2-A04X | IHC score 3+ | 3+ | Positive | Positive | Positive | 877.1 | yes |
| HER2-positive | 2 | TCGA-A2-A0EY | IHC score 3+ | 3+ | Positive | Positive | Negative | 1681 | yes |
| HER2-positive | 3 | TCGA-A2-A0T1 | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | 1699 | yes |
| HER2-positive | 4 | TCGA-A1-A0SM | IHC score 3+ | 3+ | Positive | Positive | Negative | 3101 | yes |
| HER2-positive | 5 | TCGA-A2-A04U | ISH positive | 1+ | Positive | Negative | Negative | 28.21 | no |
| HER2-positive | 6 | TCGA-A2-A0D1 | IHC score 3+ | 3+ | Positive | Negative | Negative | 1123 | no |
| HER2-positive | 7 | TCGA-A2-A04W | ISH positive | [Not Available] | Positive | Negative | Negative | 1133 | no |
| HER2-positive | 8 | TCGA-A2-A0SY | ISH positive | [Not Available] | Positive | Positive | Positive | 420 | no |
| HER2-positive | 9 | TCGA-A2-A0EQ | IHC score 3+ | 3+ | Negative | Negative | Negative | 78.19 | no |
| HER2-positive | 10 | TCGA-A2-A0CX | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Negative | 1466 | no |
| HER2-zero | 1 | TCGA-A1-A0SP | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 51.02 | yes |
| HER2-zero | 2 | TCGA-A2-A0EU | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 117.1 | no |
| HER2-zero | 3 | TCGA-A2-A0T2 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 40.53 | no |
| HER2-zero | 4 | TCGA-A2-A0CM | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 42.18 | no |
| HER2-zero | 5 | TCGA-A2-A0EV | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 226.6 | no |
| HER2-zero | 6 | TCGA-A2-A0D2 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | 61.59 | no |
| HER2-zero | 7 | TCGA-A2-A0EW | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 141.1 | no |
| HER2-zero | 8 | TCGA-A2-A0D0 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | 38.39 | no |
| HER2-zero | 9 | TCGA-A1-A0SK | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 6.718 | no |
| HER2-zero | 10 | TCGA-A2-A0T0 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 54.48 | no |

## Local Outputs

- Cases CSV: `data/tcga_brca/clinical_her2_cohort_cases.csv`
- Slide table: `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- Slide manifest: `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`

These CSV/TSV files are under `data/`, so they are local reproducible outputs rather than tracked Git files.
