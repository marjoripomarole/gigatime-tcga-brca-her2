# Model Experiments Map

This folder tracks model-family experiments beyond the primary GigaTIME analysis. The detailed reports currently remain at top-level `docs/` for compatibility.

## Current Model Roles

| Model family | Role | Current status |
|---|---|---|
| GigaTIME | Primary virtual mIF/TIME feature generator | Current main result is complete for high-trust 171-slide cohort |
| H0-mini / H-Optimus-0 | Generic H&E embedding baseline | Full 171-slide H-Optimus-0 embedding control done; reproduces the HER2-low/zero separation and the source-site confound (see `clinical_her2_high_trust_tile128_hoptimus_embedding_control.md`) |
| Virchow2 | Second generic H&E embedding control | Full 171-slide Virchow2 (2560-d) embedding control done; replicates the H-Optimus result (low-vs-zero reproduced, collapses under source-site holdout) |
| Phikon | Open H&E embedding fallback | One TCGA tile embedding smoke succeeded |
| HistoPrism | Tile/spot-level virtual expression follow-up | One-vector smoke exists; true TCGA tile run waits on compatible UNI embeddings |
| DeepSpot | Virtual spatial transcriptomics follow-up | Synthetic-vector smoke and real H-Optimus one-tile TCGA smoke succeeded |

## Reports

- `clinical_her2_high_trust_tile128_hoptimus_embedding_control.md`
- `clinical_her2_high_trust_tile128_virchow2_embedding_control.md`
- `hoptimus_embedding_baseline.md`
- `histoprism_one_vector_smoke.md`
- `deepspot_one_vector_smoke.md`

## Smoke-Test Scripts

- `scripts/run_hoptimus_tcga_brca.py`
- `scripts/run_histoprism_one_vector_smoke.py`
- `scripts/run_deepspot_one_vector_smoke.py`

## Interpretation Rule

Model-family experiments should answer whether a signal generalizes across H&E representations. They should not be promoted to biological validation unless their outputs are linked to tumor-rich/pathologist review and appropriate RNA/protein evidence.
