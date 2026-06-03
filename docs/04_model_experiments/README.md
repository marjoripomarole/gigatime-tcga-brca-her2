# Model Experiments Map

This folder tracks model-family experiments beyond the primary GigaTIME analysis. The detailed reports currently remain at top-level `docs/` for compatibility.

## Current Model Roles

| Model family | Role | Current status |
|---|---|---|
| GigaTIME | Primary virtual mIF/TIME feature generator | Current main result is complete for high-trust 171-slide cohort |
| H0-mini / H-Optimus-0 | Generic H&E embedding baseline | Planned/partially set up; gated model access and large downloads matter |
| Phikon | Open H&E embedding fallback | One TCGA tile embedding smoke succeeded |
| HistoPrism | Tile/spot-level virtual expression follow-up | One-vector smoke exists; true TCGA tile run waits on compatible UNI embeddings |
| DeepSpot | Virtual spatial transcriptomics follow-up | Checkpoint smoke exists; true TCGA tile run waits on H-Optimus-0 embeddings |

## Reports

- `hoptimus_embedding_baseline.md`
- `histoprism_one_vector_smoke.md`
- `deepspot_one_vector_smoke.md`

## Smoke-Test Scripts

- `scripts/run_hoptimus_tcga_brca.py`
- `scripts/run_histoprism_one_vector_smoke.py`
- `scripts/run_deepspot_one_vector_smoke.py`

## Interpretation Rule

Model-family experiments should answer whether a signal generalizes across H&E representations. They should not be promoted to biological validation unless their outputs are linked to tumor-rich/pathologist review and appropriate RNA/protein evidence.
