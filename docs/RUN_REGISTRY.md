# Run Registry

Status: curated run-by-run map. This file is meant to answer "what is current, what is historical, and where did each result come from?"

| Run ID | Status | Cohort | Model/input | Main command or script | Main output/report | Current interpretation |
|---|---|---|---|---|---|---|
| clinical_her2_30_tile64 | Historical | 10 HER2-positive, 10 HER2-low, 10 HER2-zero | GigaTIME, 64 tiles/slide | `scripts/run_gigatime_tcga_brca.py` | `clinical_her2_gigatime_run.md` | First pilot; useful provenance only |
| clinical_her2_30_tile256 | Historical | Same 30-slide cohort | GigaTIME, 256 tiles/slide | `scripts/run_gigatime_tcga_brca.py` | `clinical_her2_tile_sampling_robustness.md` | Directional robustness check for original pilot |
| clinical_her2_60_tile256 | Historical comparator | 20/20/20 HER2 groups | GigaTIME, 256 tiles/slide | `scripts/run_gigatime_tcga_brca.py` | `clinical_her2_expanded20_results.md` | Stronger low-vs-zero signal; now a comparator |
| high_trust_171_tile128 | Current primary | 53 HER2-positive, 57 HER2-low, 61 HER2-zero | GigaTIME, 128 tiles/slide | `scripts/run_gigatime_tcga_brca.py` | `clinical_her2_high_trust_tile128_results.md` | Current top-line hypothesis-generating result |
| high_trust_tissue_composition | Current caveat | HER2-low vs HER2-zero subset | GigaTIME tile composition features | `scripts/analyze_tissue_composition_sensitivity.py` | `clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md` | Low-vs-zero signal is strongly tied to tissue composition |
| high_trust_tumor_proxy | Current caveat | HER2-low vs HER2-zero subset | Virtual CK/tumor-rich proxy views | `scripts/analyze_tumor_proxy_sensitivity.py` | `clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md` | Classifier signal persists in proxy views, but univariate channels weaken |
| high_trust_covariates | Current caveat | HER2-low vs HER2-zero subset | Clinical/source-site/slide-size covariates | `scripts/analyze_clinical_covariate_sensitivity.py` | `clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md` | Source-site and slide-size confounding are strong |
| high_trust_source_holdout | Current caveat | HER2-low vs HER2-zero subset | Leave-source-site-out validation | `scripts/analyze_source_site_generalization.py` | `clinical_her2_high_trust_tile128_source_site_generalization.md` | GigaTIME performance drops under source-site holdout |
| high_trust_erbb2_rna | Current validation | Local STAR-count cases | ERBB2 gene-level TPM | `scripts/analyze_local_erbb2_expression_validation.py` | `clinical_her2_high_trust_tile128_local_erbb2_validation.md` | ERBB2 RNA validates broad HER2-positive labels, not low-vs-zero strongly |
| isoform_feasibility | Current guardrail | Local RNA files | STAR gene-count audit | `scripts/audit_her2_isoform_validation_feasibility.py` | `her2_isoform_validation_feasibility.md` | Current files cannot directly validate isoforms or junctions |
| hoptimus_embedding_extraction | Current model baseline | 171 high-trust slides | H-Optimus-0, 1536-d, 128 tiles/slide | `scripts/run_hoptimus_tcga_brca.py` | `hoptimus_embedding_baseline.md` | Generic H&E embedding now extracted for the full high-trust cohort |
| hoptimus_embedding_control | Current control | 118 HER2-low/zero high-trust slides | H-Optimus-0 embedding classifier | `scripts/analyze_hoptimus_embedding_control.py` | `clinical_her2_high_trust_tile128_hoptimus_embedding_control.md` | Generic embedding reproduces low-vs-zero (BA 0.726, beats shuffled null p=0.005, ~matches GigaTIME 0.710) and collapses under source-site holdout (0.726->0.586) while slide-size stays portable (0.882); confound confirmed, virtual-immune framing not required |
| histoprism_one_vector | Smoke test | Official HistoPrism sample vector | HistoPrism with precomputed UNI vector | `scripts/run_histoprism_one_vector_smoke.py` | `histoprism_one_vector_smoke.md` | Model loads; true TCGA tile run waits on compatible UNI embeddings |
| deepspot_one_vector | Smoke test | Synthetic vector and one real H-Optimus-0 tile from `TCGA-A7-A26J` | DeepSpot lung Visium pretrained checkpoint | `scripts/run_deepspot_one_vector_smoke.py` | `deepspot_one_vector_smoke.md` | Model loads; real one-tile TCGA smoke succeeds; not yet a biological result |

## Current Priority

The current research priority is not another classifier score. It is reducing ambiguity around whether the HER2-low versus HER2-zero signal reflects tumor biology, tissue composition, or TCGA acquisition/source-site structure.

The 2026-06-04 generic-embedding control sharpened this: a generic H-Optimus-0 embedding with no immune interpretation reproduces the low-versus-zero separation and the same source-site collapse as GigaTIME. The remaining priority is therefore external/site-balanced validation and pathologist tumor-region review, not more GigaTIME-internal sensitivity analyses.
