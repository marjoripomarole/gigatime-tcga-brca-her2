# Advisor Brief: BRCA HER2 Pathology AI

Status: Current advisor-facing summary. For the complete documentation map, start with `docs/README.md`.

## Working Question

Can computational pathology features from TCGA-BRCA H&E slides, starting with GigaTIME virtual mIF predictions, generate biologically interpretable tumor immune microenvironment signals that vary across clinically defined HER2-positive, HER2-low, and HER2-zero breast cancers?

## Concrete Pilot

1. Query TCGA-BRCA diagnostic H&E slides, RNA-seq STAR-count files, and clinical HER2 supplement fields from GDC.
2. Assign clinical HER2-positive, HER2-low, HER2-zero, or unknown labels using IHC/ISH fields.
3. Select a balanced 30-case pilot: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero.
4. Expand to a balanced 60-case run: 20 HER2-positive, 20 HER2-low, and 20 HER2-zero.
5. Scale to the largest laptop-realistic balanced local cohort: 61 HER2-positive, 61 HER2-low, and 61 HER2-zero downloaded slides.
6. Clean the 183-slide cohort using HER2 label trust, slide metadata, file-integrity, OpenSlide readability, and a female-patient TCGA-BRCA filter. The female-patient filter follows the relevant TCGA sample-selection principle from Guardia et al., Genome Research 2025, PMID 40664477; the file-integrity and H&E slide QC steps are specific to this pathology-AI project.
7. Run the official GigaTIME model on the high-trust diagnostic slide tiles, then filter the current primary analysis to 171 strict trustworthy slides.
8. Aggregate virtual mIF activations per slide for the GigaTIME channels.
9. Compare virtual TIME markers across clinical HER2 groups, run tile-cleanup sensitivity analyses, and test baseline classifiers.

## Current Best Finding

The current result to present is the strict high-trust 171-slide TCGA-BRCA analysis, not the earlier 30- or 60-slide pilots.

The high-trust cohort contains:

| Group | High-trust slides |
|---|---:|
| HER2-positive | 53 |
| HER2-low | 57 |
| HER2-zero | 61 |

All 171 primary-analysis slides passed HER2 label, female-patient, local file-integrity, and slide-readability checks. The raw GigaTIME output contains the earlier 174-slide inference run, but the current summaries filter it to 171 strict trustworthy slides, giving 21,888 primary-analysis tile predictions.

The strongest finding is still HER2-low versus HER2-zero. In the all-sampled-tissue view, HER2-low had lower GigaTIME-predicted immune/myeloid/checkpoint and tissue-context signals than HER2-zero, including `CD68`, `PD-L1`, `PD-1`, `CD11c`, `CD4`, `CD3`, and `CK`. These differences passed within-view multiple-testing correction in the high-trust run.

The cleaned classifier comparison also supports the same boundary: the best HER2-low versus HER2-zero classifier reached balanced accuracy 0.727 and macro AUC 0.787. HER2-positive classification remained weak, so this is not yet a diagnostic HER2 classifier.

The result is also robust across GigaTIME run settings. Comparing the earlier 60-slide 256-tile run against the current strict high-trust 128-tile analysis on 56 overlapping slides showed very high slide-level channel agreement, with key-channel Spearman correlations around 0.97-0.99. All 8 tested key channels kept the same HER2-low versus HER2-zero direction across runs, and 7 of 8 kept HER2-low lower than HER2-zero.

ER/PR and HER2-detail sensitivity checks support the same cautious interpretation. In the all-sampled-tissue view, 7 of 8 tested key channels remain significant after ER/PR adjustment. The signal remains visible across HER2-low IHC `1+`, HER2-low IHC `2+`/ISH-negative, HER2-zero IHC `0`/ISH-negative, and HER2-zero IHC `0`/ISH-not-evaluated subgroups. The strict CK-enriched view still weakens, so the result is best described as tissue-context association rather than a purely epithelial HER2 phenotype.

The newest case-level driver check makes the result easier to discuss honestly. Among 118 HER2-low/HER2-zero slides, 71 matched the expected low-versus-zero profile in at least 3 of 4 cleanup views, while 47 showed the opposite profile in at least 2 views. This supports a real but imperfect case-level pattern and gives us a concrete manual pathology/QC shortlist rather than only group-average statistics.

The visual QC spot check is the key caution. In a small set of rendered H&E plus virtual mIF panels, low-like selected tiles often had high tissue fraction but almost no virtual CK, CD68, PD-L1, or CD11c, and some looked stromal/collagen-rich rather than clearly tumor-rich. This means the current result may partly reflect tissue composition. The next advisor/pathologist question should be: are the low-like driver tiles biologically meaningful tumor regions or mostly non-tumor tissue context?

We then quantified that caveat. HER2-low slides have a much higher fraction of low-marker tiles than HER2-zero slides: 0.349 versus 0.180, BH q = 0.000265. The case-driver score is strongly correlated with marker/tissue composition, and after adjusting for low-marker tile fraction, most low-versus-zero channel effects collapse. This shifts the honest interpretation toward a GigaTIME-derived tissue-context association, not tumor-cell HER2 biology.

We also ran a stricter virtual tumor-rich proxy analysis. Fixed-count CK-rich tile views weaken individual channel-test significance, but the multichannel HER2-low versus HER2-zero classifier remains above chance across proxy views. The strongest proxy result is the absolute CK-high QC view: balanced accuracy 0.761 and macro AUC 0.782 for low versus zero. This is encouraging, but still not pathologist-confirmed tumor annotation.

Finally, we ran a shuffled-label classifier sanity check. Across the selected low-versus-zero classifier views, observed repeated-CV balanced accuracy stayed around 0.693-0.729, while shuffled-label null means stayed around 0.482-0.488. All tested views beat the null distribution with empirical p = 0.0099 and BH q = 0.0099. This supports that the classifier signal is not obviously random, but it remains post-hoc exploratory evidence rather than clinical validation.

We then ran a stricter nested model-selection check. In this analysis, the classifier chooses the best feature set inside each training fold before evaluating the held-out fold. The low-versus-zero balanced accuracy still remained about 0.672-0.721 across proxy views and beat fully nested shuffled-label null tests with empirical p = 0.0323 and BH q = 0.0323. This reduces the feature-selection-bias concern, but it is still internal validation.

The clinical/source-site covariate sensitivity check is the strongest current caveat. HER2-low and HER2-zero slides are heavily imbalanced for slide size and TCGA source site. Slide-size-only and source-site-only covariates classify HER2-low versus HER2-zero better than GigaTIME features: in the top 8 CK proxy view, slide-size covariates reached balanced accuracy 0.879, source-site covariates 0.878, source-site plus slide-size 0.897, while GigaTIME mean channels reached 0.745. This means the current TCGA classifier signal may be substantially confounded by cohort construction or acquisition/site differences.

The matched HER2-low versus HER2-zero sensitivity check narrows the claim further. Exact source-site matching leaves only 12 low/zero pairs; slide-size matching gives 14-20 pairs depending on the caliper. In the top 8 CK proxy view, GigaTIME mean channels remain modestly above chance in matched subsets, with balanced accuracy 0.708 in exact source-site pairs, 0.679 in strict slide-size pairs, and 0.675 in wider slide-size pairs. However, source-site and slide-size baselines remain competitive or stronger in the larger matched subsets, and paired channel tests do not reach BH q < 0.05. This means the current result is still worth studying, but it should not be presented as independent HER2 biology from TCGA alone.

The leave-source-site-out validation gives the most direct classifier robustness test. In the top 8 CK proxy view, GigaTIME mean channels drop from balanced accuracy 0.745 under repeated stratified cross-validation to 0.669 when entire TCGA source sites are held out. Slide-size covariates remain very strong under the same source-site holdout, with balanced accuracy 0.882. This means the current classifier signal is not yet source-independent; it may include a real image pattern, but the acquisition/slide-size structure is still too strong for a biological claim.

The within-source-site sensitivity check is a smaller follow-up stress test. Only four TCGA source sites contain both HER2-low and HER2-zero strict high-trust cases: A1, A2, A8, and AO. This leaves only 51 cases, with 12 HER2-low and 39 HER2-zero. Site-fixed channel models retain 7 channel/view effects with BH q < 0.05, mostly in QC-cellular, CK-top-25%, and absolute CK-high views. In the top 8 CK proxy view, all-channel GigaTIME reaches 0.667 repeated-CV balanced accuracy and 0.628 leave-mixed-source-site-out balanced accuracy. This supports continued investigation, but the subset is too small and imbalanced to rescue the classifier as source-independent HER2 biology.

The generic H&E embedding control is the most direct test of whether the GigaTIME virtual-immune framing is even required. We extracted H-Optimus-0 embeddings (1536-d, mean-pooled over 128 random tissue tiles per slide) for the same 171 high-trust slides and classified HER2-low versus HER2-zero on identical folds, with PCA fit inside each training fold. A generic foundation-model embedding with no immune interpretation reproduces the separation about as well as GigaTIME: balanced accuracy 0.726 versus 0.710, and the embedding beats its own shuffled-label null at empirical p = 0.005 (stable across 10/20/30 PCA components). It then collapses under leave-source-site-out validation the same way GigaTIME does, from 0.726 to 0.586, while slide-size covariates remain portable at 0.882. This means the GigaTIME-specific virtual immune/myeloid/checkpoint story is not required to explain the low-versus-zero signal; generic morphology that tracks TCGA acquisition structure reproduces it. The honest interpretation is now firmly a confounded tissue-context association, and the next step is external/site-balanced data, not more TCGA-internal analysis.

The expanded local ERBB2 validation adds useful RNA context. We extracted ERBB2 gene-level TPM from all 110 local STAR count cases, including 56 strict high-trust GigaTIME/HER2 cases and 40 HER2-low/HER2-zero cases. ERBB2 RNA strongly validates the broad HER2-positive label as a sanity check: ERBB2-only AUC is 0.905 for HER2-positive versus non-positive. However, ERBB2 RNA weakly separates HER2-low from HER2-zero: AUC 0.605, median TPM 83.4 versus 62.7, and pairwise p/q 0.262/0.262. This means the GigaTIME low/zero signal is not simply a strong gene-level ERBB2 expression split, but it also does not validate HER2 isoform biology.

The HER2 isoform feasibility audit adds an important boundary around the Guardia et al. connection. Locally, we have 110 STAR augmented gene-count cases, including 56 strict high-trust cases and 40 HER2-low/HER2-zero high-trust cases. These files contain gene-level counts/TPM and can support ERBB2 expression or RNA-program context. They do not contain transcript-level isoform quantification, BAM/FASTQ reads, or junction-count outputs, so they cannot reproduce kallisto/SUPPA2/rMATS-style HER2 isoform labels. To test the isoform hypothesis directly, we need sample-level HER2 isoform labels from the Guardia/Galante workflow or appropriate RNA-seq read/junction data.

The careful presentation sentence is:

> In a strict high-trust female TCGA-BRCA diagnostic H&E cohort, using HER2 isoform/state biology from Guardia et al. as the motivating context, GigaTIME-derived virtual mIF features reproducibly associate with the HER2-low versus HER2-zero boundary, especially through lower virtual immune/myeloid/checkpoint and CK-associated signals in HER2-low tumors. Local gene-level ERBB2 RNA validates broad HER2-positive status but weakly separates HER2-low from HER2-zero, suggesting the image signal is not just a strong ERBB2 expression split. The low-versus-zero classifiers beat shuffled-label sanity checks even under nested feature-set selection, but clinical/source-site covariate, matched-subset, leave-source-site-out, and within-source-site sensitivity checks show major TCGA confounding risk: slide-size and source-site covariates classify HER2-low versus HER2-zero very well, matching does not fully remove this concern, GigaTIME performance drops when source sites are held out, and the mixed-site subset is small and imbalanced. This supports a hypothesis-generating tissue-context observation and a clear validation plan, not a clinical diagnostic claim and not evidence of tumor-cell HER2 isoform detection.

A 2026-06-04 generic H-Optimus-0 embedding control reinforces the confounding caveat: a foundation-model embedding with no immune interpretation reproduces the HER2-low versus HER2-zero separation and the same source-site collapse, so the virtual-immune framing is not required to explain the signal.

See `docs/clinical_her2_high_trust_tile128_results.md` and `docs/clinical_her2_high_trust_tile128_hoptimus_embedding_control.md`.

## Historical 30-Slide Pilot Finding

The first balanced clinical HER2 pilot processed all 30 selected slides using 64 random tissue tiles per slide.

The strongest pilot signal was not HER2-positive versus HER2-low. Instead, HER2-zero showed higher mean GigaTIME-predicted immune/checkpoint signals than HER2-low, especially:

- `CD68`
- `PD-L1`
- `CD11c`
- `CD4`
- `Ki67`

The top unadjusted three-group tests were CD68, PD-L1, and CD11c. Pairwise HER2-low versus HER2-zero tests were the strongest, but none remained significant after Benjamini-Hochberg correction. This should be framed as a hypothesis-generating signal.

## First Validation Check

We compared the GigaTIME virtual channels with matched TCGA RNA-seq marker signatures for the same 30 cases.

This did not strongly confirm the virtual immune-channel signal:

- `Ki67` had the strongest positive correlation with its RNA signature, but it was weak and not FDR-significant.
- `CD68`, `PD-L1`, and `CD11c`, the main virtual channels driving the HER2-zero versus HER2-low signal, did not show strong positive RNA-signature correlations.
- This means the current result should remain hypothesis-generating until visual QC, more tile sampling, and additional validation are done.

## Visual QC Update

We rendered H&E-versus-virtual-mIF QC panels for the top `CD68` + `PD-L1` + `CD11c` case in each HER2 group.

The high-scoring tiles were real tissue-containing, cellular H&E regions rather than obvious blank background. That supports continuing the analysis. However, the high-signal tiles were not visually unique to HER2-zero; the selected HER2-positive case also had strong high-signal tiles. The result remains a slide-level pilot trend, not a clean single-case visual phenotype.

## 256-Tile Robustness Update

We reran the same 30 selected clinical HER2 slides with up to 256 random tissue tiles per slide. This was done to test whether the original 64-tile signal was mainly a sparse-sampling artifact.

The main result persisted:

| Channel | 64-tile p | 256-tile p | 64 max-min | 256 max-min | 256 direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |

Pairwise HER2-low versus HER2-zero q values improved for CD68, PD-L1, and CD11c to 0.1133, but they still did not meet the usual 0.05 FDR threshold.

The 256-tile RNA validation remained weak. No virtual channel had an FDR-significant correlation with its matched RNA marker signature. Therefore, the current interpretation is stronger sampling robustness, not biological validation.

## RNA Program Validation Update

We then tested broader RNA immune and tissue programs rather than only single marker-channel signatures. This compared GigaTIME virtual composite programs with RNA programs for T-cell/cytotoxic, checkpoint/IFNG, myeloid/macrophage, dendritic/APC, B-cell, proliferation, epithelial, stromal, and endothelial biology.

The broader validation still did not positively confirm the virtual immune/checkpoint signal.

Key findings:

- The virtual myeloid/checkpoint composite retained the HER2-zero greater than HER2-low direction, but was not FDR-significant: Kruskal p 0.0176, BH q 0.0878.
- No broad RNA immune program showed an FDR-significant HER2-group difference.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal:
  - Virtual T cell/checkpoint versus endothelial RNA: Spearman rho -0.585, BH q 0.0309.
  - Virtual all immune/checkpoint versus endothelial RNA: Spearman rho -0.556, BH q 0.0320.

This is a cautionary result. It suggests the virtual signal is reproducible inside GigaTIME, but not yet validated against orthogonal RNA evidence. It also raises the possibility that tissue composition, stromal/endothelial context, or slide sampling may be influencing the predictions.

## First Classifier Baseline

We then moved from group-average comparisons to a first diagnostic-model style classifier. The classifier used slide-level GigaTIME features from the 256-tile run and leave-one-out cross-validation.

Three tasks were tested:

- HER2-positive versus HER2-negative.
- HER2-low versus HER2-zero.
- Full three-class HER2-positive versus HER2-low versus HER2-zero.

Main result:

| Task | Best GigaTIME/H&E feature set | Balanced accuracy | Macro AUC |
|---|---|---:|---:|
| HER2-low vs HER2-zero | GigaTIME mean + fraction channels | 0.800 | 0.870 |
| HER2-positive vs HER2-negative | GigaTIME mean + fraction channels | 0.475 | 0.430 |
| Three-class HER2 group | GigaTIME mean + fraction channels | 0.333 | 0.555 |

Interpretation:

- The HER2-low versus HER2-zero result is promising but very small-sample and potentially unstable.
- GigaTIME/H&E features do not currently classify HER2-positive status reliably.
- Full three-class prediction is at chance.
- ERBB2 RNA, included as a non-H&E reference, classified HER2-positive versus HER2-negative better than GigaTIME/H&E features. This means the labels contain molecular signal, but the current image-derived features are not capturing the clinical HER2-positive signal reliably.

## Expanded 20/20/20 Update

We expanded the balanced clinical HER2 run from 30 slides to 60 slides: 20 HER2-positive, 20 HER2-low, and 20 HER2-zero. The expanded GigaTIME run used up to 256 tissue tiles per slide and produced 15,225 tile predictions. STAR-count RNA-seq expression was downloaded for all 60 selected cases.

The strongest expanded result is that the HER2-low versus HER2-zero signal became more convincing:

- All sampled tissue HER2-low versus HER2-zero classifier: balanced accuracy 0.800, macro AUC 0.820.
- QC-cellular classifier: balanced accuracy 0.775, macro AUC 0.820.
- CK-enriched top 25% classifier: balanced accuracy 0.800, macro AUC 0.820.
- Pairwise HER2-low versus HER2-zero differences now pass within-view BH correction for several virtual immune/myeloid/checkpoint channels, especially `CD3`, `CD4`, `CD11c`, and `CD68`; `PD-L1` passes in the QC-cellular view.

The expanded three-group pattern is more nuanced than the first 30-slide run. HER2-low is often the lowest virtual immune/checkpoint group, but HER2-positive is highest for several broader virtual immune programs. RNA validation remains weak, so the result should still be presented as an image-derived, hypothesis-generating HER2-state association. This 60-slide result is now historical because the strict high-trust 171-slide analysis is the current primary result.

See `docs/clinical_her2_expanded20_results.md`.

## Sharper Paper Angle: HER2 State, Not Only HER2 Amount

The strongest proposal framing is not simply "use image AI to classify HER2." A more interesting biology question is whether H&E-derived and GigaTIME-derived features predict or associate with HER2-related molecular states that are not fully captured by the routine clinical HER2-positive, HER2-low, and HER2-zero labels.

Especially interesting future hypotheses include:

- HER2-low versus HER2-zero tumors with hidden or alternate ERBB2 transcript/isoform expression.
- HER2-positive tumors with image-derived states associated with trastuzumab or antibody-drug conjugate resistance.
- Tumors with preserved HER2 signaling but reduced antibody targetability.

This must be worded carefully. The current project can say that image AI predicts or associates with HER2 isoform/state hypotheses. It should not say that image AI detects HER2 isoforms. Isoform claims require transcript-level or protein-level validation beyond the current GigaTIME/H&E pilot.

See `docs/her2_isoform_state_hypothesis.md` for the working language guardrails and validation plan.

## Why This Is a Good First Step

- It is replication-first: the model is not retrained, only applied to public TCGA-BRCA data.
- TCGA-BRCA gives paired histology, transcriptomic context, and clinical HER2 IHC/ISH fields suitable for an exploratory HER2 axis.
- The output is easy to inspect with an advisor: slide-level CSVs, channel summaries, and figures.
- The all-channel virtual mIF figures in `docs/assets/virtual_mif_channels/` show group-level channel means, slide-by-channel relative activation, and representative HER2-high/HER2-low spatial channel grids.

## Caveats to Discuss

- Clinical HER2 labels are derived from TCGA clinical supplement IHC/ISH fields, which are incomplete and must be described carefully.
- The first balanced clinical HER2 run was small: 10 cases per group.
- The 64-tile-per-slide run is a practical pilot, not a final whole-slide sampling strategy.
- The 256-tile rerun supports robustness to denser sampling, but it is still not exhaustive whole-slide analysis.
- The current high-trust run is larger, but it still samples 128 tiles per slide rather than exhaustively processing every tissue tile.
- The high-trust slide list improves label and file trust, but it is not a substitute for pathologist-confirmed tumor-rich region review.
- Bulk RNA-seq is an indirect validation layer and did not strongly validate the current GigaTIME immune-channel pattern, even with broader RNA programs.
- The first classifier baseline is not clinically usable; it is a feasibility and failure-mode analysis.
- The current pilot does not demonstrate HER2 isoform detection. It only motivates future testing of whether image-derived features associate with HER2 state, targetability, or transcript-level biology.
- Visual QC supports that the signal is not just blank background, but it is not biological validation.
- TCGA slide quality, tissue sampling, and tumor purity need QC before strong biological claims.
- GigaTIME is research-only and not a clinical HER2 classifier.
- The virtual mIF channel images are GigaTIME predictions from H&E tiles, not real multiplex immunofluorescence measurements from TCGA.
