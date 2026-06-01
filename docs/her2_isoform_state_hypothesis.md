# HER2 Isoform and State Hypothesis

Status: Current hypothesis/framing document. This is not a results report; it defines careful proposal language and future validation needs.

## Why This Is a Stronger Scientific Angle

The project should not be framed only as "can H&E predict HER2 status?" That is an important diagnostic question, but it is not the most interesting biology question.

A stronger framing is:

Can image-derived features from H&E, especially GigaTIME virtual immune and tissue channels, predict or associate with HER2-related biological states that are not captured by a simple HER2-positive, HER2-low, or HER2-zero label?

This connects the image-analysis work to a deeper idea: HER2 biology is not only about how much HER2 is present. It may also depend on which HER2 molecular form, signaling state, tissue context, and antibody-accessible target state is present.

## What Would Be Especially Interesting

The most exciting future direction would be to test whether image features distinguish:

- HER2-low versus HER2-zero tumors with hidden or alternate ERBB2 transcript/isoform expression.
- HER2-positive tumors with image-derived features associated with resistance to trastuzumab or antibody-drug conjugate therapy.
- Tumors with preserved HER2 pathway signaling but reduced antibody targetability.

These are not claims from the current pilot. They are hypotheses that would make the project more biologically meaningful than a generic HER2 classifier.

## Language Guardrails

Preferred language:

- Image AI predicts or associates with HER2 isoform/state hypotheses.
- Image-derived features stratify tumors by HER2-related biological state.
- GigaTIME features may capture image-level correlates of HER2 signaling, immune context, or targetability.
- The model identifies cases that may warrant transcript-level, protein-level, or therapy-response validation.

Language to avoid:

- Image AI detects HER2 isoforms.
- H&E directly measures HER2 isoforms.
- GigaTIME diagnoses isoform status.
- The model proves antibody targetability or therapy resistance.

The core distinction is important: H&E and virtual mIF predictions can generate or prioritize biological hypotheses, but isoforms require molecular validation.

## What The Current Pilot Supports

The current pilot supports a cautious starting point:

- There is a reproducible image-derived signal separating HER2-low and HER2-zero cases, first seen in the 30-slide pilot and strengthened in the expanded 60-slide 20/20/20 run.
- In the expanded run, several HER2-low versus HER2-zero virtual immune/myeloid/checkpoint channel differences pass within-view BH correction, especially `CD3`, `CD4`, `CD11c`, and `CD68`; `PD-L1` passes in the QC-cellular view.
- The HER2-low versus HER2-zero classifier remained around balanced accuracy 0.800 in the expanded 40-case binary comparison.
- The three-group pattern is nuanced: HER2-low often appears lowest, while HER2-positive becomes highest for several broader virtual immune programs.
- HER2-positive versus HER2-negative classification remains weak, so the current pipeline does not reliably predict clinical HER2-positive disease.

This means the current study is best interpreted as evidence that image-derived tissue context may contain information related to the HER2-low versus HER2-zero boundary. It does not yet show that the images capture HER2 isoforms, and RNA validation remains weak.

## What Validation Is Needed

To test the isoform/state hypothesis, the project needs stronger orthogonal evidence:

- Transcript-level ERBB2 analysis, ideally isoform-aware or splice-aware quantification rather than only gene-level ERBB2 expression.
- Protein-level validation using IHC, ISH, real multiplex immunofluorescence, proteomics, or another antibody-based assay when available.
- Treatment-response validation in an external cohort with trastuzumab or antibody-drug conjugate exposure and outcome data.
- Pathologist review of the H&E tiles and virtual mIF-like outputs that drive the HER2-low versus HER2-zero classifier.
- Tumor purity, stromal, endothelial, and immune-composition adjustment to determine whether the image signal reflects tumor-intrinsic HER2 biology or tissue-context biology.

TCGA-BRCA is useful for histology, clinical HER2 annotation, and RNA context, but TCGA alone is unlikely to validate trastuzumab or ADC resistance because treatment-response information is limited for that purpose.

## Practical Next Analyses

1. Check whether available RNA data can support ERBB2 transcript-level or isoform-aware analysis.
2. Compare GigaTIME HER2-low versus HER2-zero predictions with ERBB2 gene expression, HER2 clinical label, and any transcript-level ERBB2 evidence available.
3. Identify high-confidence HER2-low and HER2-zero cases where the image classifier is very confident, then inspect whether those cases have unusual ERBB2 RNA or pathway features.
4. Separate tumor-epithelial, stromal, endothelial, and immune-context features so the paper can say what kind of image signal is being learned.
5. Search for an external cohort with H&E, HER2 IHC/ISH, molecular profiling, and preferably anti-HER2 therapy outcomes.

## Proposal-Safe Claim

The safest paper-proposal claim is:

This study evaluates whether GigaTIME-derived image features from breast cancer H&E slides predict or associate with HER2-related biological states, including the HER2-low versus HER2-zero boundary, and whether those image-derived states can be validated against molecular or protein-level evidence.
