# Literature Review: HER2-low vs HER2-zero Microenvironment and Confounders

*Generated via Agentic PubMed Search (June 2026).*

*Citations verified 2026-06-04: PMIDs 36967073 and 41670928 confirmed real and on-topic; 41316049 confirmed real but its specific TIL-direction claim is unconfirmed (see flag below); 42074514 removed as mis-applied. A corrected synthesis and the reframe (real grade/immune morphology entangled with batch, not HER2 protein) live in `external_validation_candidates.md`.*

This document summarizes recent scientific literature validating the hypothesis-generating results from the GigaTIME model runs, specifically focusing on morphological differences between HER2-low and HER2-zero breast cancer, as well as TCGA batch effects.

## 1. Immune & Tumor Microenvironment Differences

The GigaTIME pipeline identified differences in "virtual immune, myeloid, and checkpoint" signals between HER2-low and HER2-zero states. Recent clinical literature supports that immune infiltration is indeed a distinguishing factor:

*   **Higher TILs in HER2-null (zero):** A 2025 study (BMC Cancer) on HER2 classifications found that **HER2-null (zero) tumors were associated with significantly higher tumor-infiltrating lymphocyte (TIL) density** compared to HER2-low and HER2-ultralow tumors. 
    *   *Reference:* [PMID 41316049](https://pubmed.ncbi.nlm.nih.gov/41316049/) - "Clinical characteristics and prognostic impact of HER2-ultralow breast cancer and tumor-infiltrating lymphocytes (TILs)."
    *   ⚠️ *Verification flag (2026-06-04):* the PMID is a real BMC Cancer 2025 paper on exactly these HER2 subgroups and TILs, but the specific "HER2-null > HER2-low TIL density" direction was **not** confirmed from the abstract (which surfaced ER-positivity and prognosis). Because this is the one result that directionally matches the GigaTIME finding, confirm it in the full text before relying on it as corroboration.
*   **Proinflammatory TME in HR+ subtypes:** Looking specifically at HR+ breast cancer, researchers found that certain HER2-low/positive subgroups (with high FGFR4 co-expression) exhibit a distinctly **proinflammatory tumor microenvironment** with increased T cell, NK cell, and M1 macrophage infiltration, alongside upregulated checkpoints (like CTLA4 and LAG3). 
    *   *Reference:* [PMID 41670928](https://pubmed.ncbi.nlm.nih.gov/41670928/) - "FGFR4 and HER2 co-expression is associated with the proinflammatory tumor microenvironment in HR+ breast cancer."

## 2. Morphological Aggressiveness

The underlying visual morphology changes depending on the HER2-low vs HER2-zero state.

*   **Less Aggressive Morphology in HER2-low:** In a large cohort study of over 1,300 breast carcinomas, researchers found that **HER2-low cases had significantly less frequent grade 3 morphology** than HER2-zero cases. They concluded that HER2-low early-stage breast carcinoma may represent a biologically less aggressive morphological group.
    *   *Reference:* [PMID 36967073](https://pubmed.ncbi.nlm.nih.gov/36967073/) - "Incidence, Clinicopathologic Features, HER2 Fluorescence In Situ Hybridization Profile, and Oncotype DX Results of Human Epidermal Growth Factor Receptor 2-Low Breast Cancers..."

## 3. Confounding and Batch Effects in TCGA-BRCA

The project's internal checks identified strong source-site confounding in the pipeline.

*   **Necessity of Correction:** When integrating TCGA-BRCA data with machine learning, rigorous **batch-effect correction** is required to prevent models from learning site-specific acquisition signatures rather than true biology. This project demonstrated that risk directly: slide-size and source-site covariates out-classify the GigaTIME image model, the signal collapses under leave-source-site-out validation, and two independent generic foundation-model embeddings reproduce it.
    *   ⚠️ *Verification flag (2026-06-04):* the originally generated citation here (PMID 42074514) was **removed**. It is a real paper but a generic multi-omics ML / CHEK1 diagnostic-gene study, not about histology batch effects, so it does not support this point. Substitute a genuine WSI/histology batch-effects reference if this claim is cited.

## Impact on Current Workflow

These findings indicate that the direction of the GigaTIME signal is **directionally consistent with established biology** (HER2-zero and HER2-low do differ in immune landscape, e.g. TIL densities, and in morphological aggressiveness/grade) — but directional consistency is **not** validation, and these are tissue-context/grade differences rather than tumor-cell HER2 protein detection. 

**Critical Caveat:** Since the immune cell presence differs between these groups organically, an image-based AI model could be measuring "how much immune tissue is here" (a tissue-composition confounder) rather than directly learning a visual representation of the HER2 protein state itself. Follow-up investigations should factor this composition confounder into the validation steps.
