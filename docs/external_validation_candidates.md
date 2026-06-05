# External Validation Candidate Cohorts

Status: current reference for external (non-TCGA) validation of the HER2-low versus HER2-zero image signal. Compiled 2026-06-04 from a web/literature scout, then updated after BCNB full-clinical access confirmed IHC-score granularity. Numbers and access details outside BCNB should still be re-verified against primary sources before committing to any cohort.

## Why External Data Is Now The Bottleneck

The HER2-low versus HER2-zero signal in this project has been shown, repeatedly and from multiple angles, to be confounded by TCGA acquisition structure (slide size and source site) rather than tumor-cell HER2 biology. As of 2026-06-04 this is supported by clinical/source-site covariate baselines, matched subsets, leave-source-site-out validation, within-source-site restriction, and two independent generic foundation-model embedding controls (H-Optimus-0 and Virchow2) that reproduce the separation and the same source-site collapse. See `docs/clinical_her2_high_trust_tile128_results.md`.

TCGA-internal evidence is therefore exhausted. Pulling more TCGA-BRCA slides does not help: HER2-zero is capped at 61 cases in all of TCGA-BRCA (already fully used), and more data along a confounded axis tightens confidence intervals around a biased estimate. The only data that can change the conclusion is variation independent of HER2 status, i.e. an external cohort, ideally single-scanner / single-institution so the slide-size/source-site confound is removed by construction.

BCNB now satisfies the first external-cohort gate: full clinical data are local, HER2 IHC 0/1+/2+/3+ is recoverable, histological grade is available, and the cohort is single-scanner. The patch-input gate is also solved for a first pilot: `paper_patches.zip` is local, valid, and maps cleanly to all 1,058 patients (see `bcnb_paper_patches_audit.md`). H-Optimus-0 and Virchow2 hash-capped patient-level patch pilots found a statistically non-null but modest low-versus-zero signal (single-model BA ~0.60, AUC ~0.64; dual-model BA 0.609, AUC 0.651), comparable to clinical covariates rather than a strong standalone classifier. The two encoders are highly concordant at the patient-score level, and patch-sampling plus clinical-stratified checks do not overturn the result; the signal is uneven across grade/receptor/subtype slices. Score-driver analysis shows measured clinical covariates + patch QC explain part of the image score and reduce dual residual AUC to 0.592, so the residual signal remains weak. Visual QC shows zero-like HER2-low false positives are not obvious blank/low-tissue artifacts and often match aggressive covariate profiles. The remaining caveat is methodological rather than logistical: precomputed tumor-region patches do not preserve whole-slide acquisition and tissue-area controls, and full WSIs remain preferable for the strongest paper-grade analysis.

## The Hard Constraint: HER2-Low-vs-Zero Granularity

The exact comparison this project cares about is HER2-low (IHC 1+, or 2+/ISH-negative) versus HER2-zero (IHC 0). This is the difficult attribute to source externally, for two reasons:

1. Most public H&E + HER2 datasets label only binary HER2-positive versus HER2-negative. HER2-low is a post-2019 clinical category (driven by trastuzumab deruxtecan / DESTINY-Breast04), so older public datasets predate it and collapse 0 and 1+ into "negative."
2. Even with the IHC slide in hand, pathologist interobserver agreement on IHC 0 versus 1+ is low. The ground truth on this exact boundary is intrinsically noisy.

Practical implication: a clean external reproduction of the low-versus-zero signal would be a strong positive result; failure to reproduce would not be surprising and would be consistent with the confound interpretation. A cohort is only fully useful for the primary question if it exposes IHC score 0 separately from 1+.

## Prior Work (Closest Published Analogue)

Valieris et al., "Weakly-supervised deep learning models enable HER2-low prediction from H&E stained slides," Breast Cancer Research 2024 (PMC11331614). This is essentially the published version of this project's question. They predicted HER2-low from H&E across three cohorts: ACCCC (private, single-institution, A.C. Camargo Cancer Center, Brazil, 546 slides / 504 patients, Leica Aperio AT2), HEROHE (public), and TCGA-BRCA (535 slides). They observed external performance drop and explicitly noted that TCGA aggregates many institutions with varied protocols affecting classification, but they did not run leave-site-out or generic-embedding confound controls.

Takeaway: this project's confound analysis is more rigorous than the published state of the art on this exact task. Valieris et al. is the key citation for a cautionary-methods paper, and the ACCCC group is a natural collaborator (their cohort is the cleanest fit for the primary question).

## Literature Context: Grade And Immune Differences Are Real (And Why Grade Matters)

An internal literature scan (`docs/literature_her2_microenvironment.md`, AI-generated and citation-verified 2026-06-04) found that HER2-low and HER2-zero differ biologically in ways that are independent of TCGA and visible on H&E, in the same direction as this project's image finding. This refines the confound interpretation rather than overturning it: the low-versus-zero image signal is most plausibly real morphology (tumor grade and immune/TIL context) entangled with TCGA acquisition batch, but neither is tumor-cell HER2 protein status. "Image AI predicts HER2-low versus zero" is better read as "it predicts grade and immune context, partly confounded by site and scanner." This also explains why two generic foundation-model embeddings reproduce the signal: grade and immune infiltrate are generic visual properties.

Verified supporting literature:

- Tumor grade differs: HER2-low cases have significantly less frequent Nottingham grade 3 morphology than HER2-zero (HER2-low is the less aggressive group), among ER-positive cases. Mod Pathol 2023, PMID 36967073. (Confirmed against PubMed.)
- Immune/TIL density differs: HER2-null (zero) tumors were reported to have higher tumor-infiltrating lymphocyte (TIL) density than HER2-low. BMC Cancer 2025, PMID 41316049. (Real paper on this exact topic; the specific HER2-null > HER2-low TIL direction should be confirmed in the full text before citing, because it is the result that directionally matches the GigaTIME finding.)

Discarded citation: the source file also cited PMID 42074514 ("Integrated Multi-Omics and Machine Learning ... Breast Cancer," Genes 2026) to support the need for TCGA batch-effect correction. That paper is real but is a generic multi-omics ML / CHEK1 paper, not about histology batch effects, so it is not used here. A genuine WSI batch-effects reference should be substituted if this point is cited.

### Grade Is The Key Confounder That TCGA Cannot Give You

Tumor grade is the most biologically relevant confounder this literature surfaces: it differs between HER2-low and HER2-zero, it is visible on H&E, and it is a plausible thing an image model reads instead of HER2. The natural next analysis would be to add grade as a covariate (exactly as was done for ER/PR and slide size) and test whether the low-versus-zero signal survives grade adjustment, and whether grade alone classifies low-versus-zero the way slide-size did.

However, TCGA-BRCA does not provide histologic grade: there is no grade field in `high_trust_slides.csv` or the GDC clinical biotab (verified 2026-06-04; a known TCGA-BRCA limitation). The signal therefore cannot be adjusted for grade inside TCGA. This is an additional reason to move to external cohorts, and it raises the priority of cohorts that report grade. BCNB, ACROBAT, and ACCCC all include histologic grade.

## Candidate Shortlist

| Cohort | What it is | Single-source? | HER2 granularity | Access |
|---|---|---|---|---|
| BCNB (Early Breast Cancer Core-Needle Biopsy) | 1,058 core-biopsy H&E WSI, China, iScan Coreo 200x | Yes - one institution, one scanner | Confirmed IHC 0/1+/2+/3+; derived zero=127, low=654, positive=277; grade/ER/PR/Ki67 available | Free registration, non-commercial; full clinical file local, WSIs/patches pending |
| ACCCC (A.C. Camargo, Brazil) | 546 H&E WSI / 504 pts, Leica Aperio AT2, 0.25 um/px | Yes — one institution, one scanner | neg / low / high (the exact split) | Private; request from Valieris et al. |
| ACROBAT | 4,212 WSI / ~1,153 pts, Swedish routine diagnostics; paired H&E + IHC (ER/PR/HER2/Ki67) consecutive sections | Yes — one source | HER2 as stained IHC slide; score likely needs deriving | Public (grand-challenge, CC) |
| HEROHE | 509 cases, single scanner (3DHistech Pannoramic 1000, Ipatimup) | Yes — one scanner | Binary positive/negative only | Public |
| Yale "HER2-TUMOR-ROIS" (TCIA) | H&E + HER2 status + trastuzumab response + tumor ROI annotations | Yes — single institution | HER2-positive focus, not low/zero | Public (TCIA) |
| IMPRESS | 126 WSI (62 HER2+, 64 TNBC), neoadjuvant chemo response, multiplex IHC (PD-L1/CD8/CD163) | Yes — single-source | HER2+ vs TNBC, not low/zero | Public |

## Recommended Path

1. Keep BCNB as the immediate external validation path. The clinical gate is solved, and the first two patch pilots plus paired H-Optimus/Virchow comparison are complete: BCNB gives a single-scanner low-versus-zero cohort of 654 vs 127 with grade, ER, PR, Ki67, molecular subtype, and nodal status. Next, run multi-seed patch-sampling/PCA sensitivity or decide whether the modest replicated signal warrants full WSI download.
2. Contact the Valieris / ACCCC group in parallel. Their cohort remains the cleanest published analogue (single-institution, single-scanner, with neg/low/high labels) and they have hit the same wall. A collaboration or data request is still valuable, especially if BCNB's core-biopsy design introduces tissue-sampling concerns.
3. Use ACROBAT as the strongest independent stress test if BCNB succeeds or is limited by WSI logistics: 4,000+ WSIs from one source, with paired HER2-IHC slides to derive status.
4. Use IMPRESS multiplex IHC (real PD-L1/CD8/CD163) to validate GigaTIME's virtual immune channels against measured immune markers, closing the RNA-validation gap that never closed.

## Open Items To Verify Before Running

- BCNB: run multi-seed patch-sampling/PCA sensitivity, or decide whether the modest replicated and highly concordant H-Optimus-0/Virchow2 signal warrants full WSI download.
- BCNB: decide whether to download full WSIs after the patch pilot; full WSIs remain strongest for slide-size/tissue-area controls.
- ACROBAT: confirm whether per-case HER2 IHC scores are tabulated in metadata, or only the stained IHC slides are provided (requiring score derivation).
- ACCCC: confirm data-sharing terms / whether the Valieris group will collaborate or release.
- Confirm the HER2-null > HER2-low TIL-density direction in the full text of PMID 41316049 before citing it as biological corroboration.
- Confirm histologic grade availability and encoding in ACROBAT and ACCCC; BCNB grade is already confirmed, and TCGA-BRCA does not provide grade.

## Sources

- Valieris et al., Breast Cancer Research 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11331614/
- ACROBAT (Scientific Data 2023): https://www.nature.com/articles/s41597-023-02422-6 ; challenge: https://acrobat.grand-challenge.org/
- BCNB: https://bupt-ai-cz.github.io/BCNB/ ; challenge: https://bcnb.grand-challenge.org/
- HEROHE (J. Imaging 2022): https://www.mdpi.com/2313-433X/8/8/213 ; challenge: https://ecdp2020.grand-challenge.org/
- Yale HER2-TUMOR-ROIS (TCIA): https://www.cancerimagingarchive.net/collection/her2-tumor-rois/
- IMPRESS (npj Precision Oncology 2023): https://www.nature.com/articles/s41698-023-00352-5
- Breast H&E WSI dataset scoping review: https://arxiv.org/html/2306.01546v2
- HER2-low grade / morphology (Mod Pathol 2023): https://pubmed.ncbi.nlm.nih.gov/36967073/
- HER2 subgroups and TILs (BMC Cancer 2025): https://pubmed.ncbi.nlm.nih.gov/41316049/
