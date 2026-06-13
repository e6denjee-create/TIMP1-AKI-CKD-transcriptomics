# Reviewer Risk Points and Response Strategies v1

## 1. The study is cross-sectional and does not directly model AKI-to-CKD transition

**Risk:** A reviewer may argue that the title or framing exceeds the available
data because no cohort follows the same patients longitudinally from AKI to
CKD.

**Response strategy:** State consistently that the project studies
AKI-to-CKD-relevant injury-repair programs, not longitudinal transition or
prediction. Emphasize the inclusion of AKI and CKD/fibrosis cohorts as
complementary disease contexts. Add longitudinal cohorts only when outcome and
follow-up metadata are available.

## 2. TIMP1 may reflect general injury rather than kidney-specific biology

**Risk:** TIMP1 is induced in multiple injured tissues and may not be specific
to kidney disease.

**Response strategy:** Agree with the biological premise and clarify that
kidney specificity was not tested or claimed. Position TIMP1 as a candidate
associated with a kidney injury-repair and ECM-remodeling state. The value of
the study is the reproducibility and transcriptomic context within human
kidney cohorts.

## 3. Disease-control correlations may be driven by group separation

**Risk:** All-sample correlations can be inflated when both TIMP1 and pathway
scores differ between disease and control.

**Response strategy:** Lead with disease-only correlations. Report that 25 of
27 discovery disease-only signature correlations and all nine GSE180394
disease-only correlations were FDR-significant. Highlight the disease-only
module construction and diagnosis-adjusted external analyses.

## 4. The GSE180394 disease group is etiologically heterogeneous

**Risk:** Associations could be driven by lupus nephritis, diabetic
nephropathy, FSGS, or another dominant subgroup.

**Response strategy:** Present the leave-one-diagnosis-out analysis and broad
diagnosis adjustment. All 22 features remained positive after removing each
category, and 21 of 22 remained FDR-significant after adjustment. Acknowledge
that broad categories cannot replace detailed severity and treatment
covariates.

## 5. The external control group is small

**Risk:** Nine living donors provide an imprecise reference.

**Response strategy:** Report effect size and uncertainty, not P value alone.
Show the extended-control sensitivity analysis, which retained a positive,
significant effect (Hedges' g = 1.291). Also report honestly that the
six-sample tumor-nephrectomy-only comparison was positive but non-significant.

## 6. Tumor-nephrectomy tissue is not a healthy control

**Risk:** Adjacent unaffected tissue may have field effects or clinical
differences.

**Response strategy:** Keep living donors as the primary control and label
tumor-nephrectomy analyses as sensitivity tests only. Do not combine them in
the primary result. Explicitly discuss why the tumor-nephrectomy-only estimate
is secondary.

## 7. The discovery datasets use different platforms and tissue contexts

**Risk:** Absolute values and effect sizes may not be directly comparable.

**Response strategy:** Explain that every dataset was analyzed independently,
with no pooling of absolute expression. Use direction, standardized effect
size, rank correlation, and cross-dataset support counts. Figure 2 should state
that y-axis scales are cohort-specific.

## 8. Prespecified signatures overlap and may create redundant findings

**Risk:** ECM remodeling, fibrosis, collagen, and maladaptive repair share
genes and are not independent tests.

**Response strategy:** Acknowledge biological overlap and avoid counting the
nine signatures as nine independent mechanisms. Emphasize coherent program
replication. Note that TIMP1 was excluded from all scores and that module
sensitivity analyses removed signature and ECM/fibrosis genes.

## 9. The TIMP1-associated module may be circular

**Risk:** The module could simply recapitulate prespecified signatures or
known matrix genes.

**Response strategy:** Report module variants after excluding all signature
genes and after excluding ECM/fibrosis markers. The stringent module retained
389 and 414 genes, respectively, supporting a broader wound-response,
inflammatory, membrane, and cytoskeletal context.

## 10. Genome-wide correlation does not define a regulatory network

**Risk:** Network figures may be interpreted as direct regulation.

**Response strategy:** Use “correlation module,” “associated program,” and
“edge list of replicated correlations.” State that network edges do not imply
physical interaction, directionality, or regulation. Avoid “downstream genes”
and “targets.”

## 11. The 12 core genes have mixed biological functions

**Risk:** A reviewer may question whether they form a coherent mechanism.

**Response strategy:** Present the core as a reproducible transcriptomic
context spanning cytoskeletal organization, trafficking, inflammatory or
complement responses, adhesion, and matrix repair. Do not claim a single
mechanism. Use external correlation and diagnosis-adjusted support as the main
evidence.

## 12. PFN1 and FBLIM1 are weaker in external validation

**Risk:** Not all core genes reproduce at FDR < 0.05 in the unadjusted
GSE180394 disease-only analysis.

**Response strategy:** Report the exact pattern: all 12 were positive, 10 were
FDR-significant unadjusted, FBLIM1 became supported after diagnosis adjustment,
and PFN1 remained positive but non-significant. Avoid calling every gene fully
replicated.

## 13. Enrichment analysis is incomplete

**Risk:** Full Hallmark, GO, Reactome, and KEGG collections were unavailable,
and a curated subset can bias interpretation.

**Response strategy:** Label enrichment exploratory, disclose package and
network limitations, and make the tested curated gene sets available. Base
the main conclusions on expression, signature correlations, and cross-dataset
module replication rather than enrichment P values. Repeat full enrichment
before submission if the required resources become available.

## 14. Single-cell evidence is vulnerable to pseudoreplication

**Risk:** One donor cannot support patient-level inference, and cell-level
P values may be extremely small.

**Response strategy:** Keep single-cell results in a supportive or
supplementary section. State that the analysis provides exploratory
localization only. Do not present cell-level P values as clinical evidence.
Prioritize multi-donor pseudobulk validation for future work.

## 15. The TIMP1-high tubular definition is affected by dropout

**Risk:** The 75th percentile was zero, making high-versus-low groups similar
to detected-versus-undetected expression.

**Response strategy:** Disclose this directly and avoid using the comparison
as primary evidence. Treat pathway findings as hypothesis-generating. Future
analyses should use multi-donor pseudobulk, hurdle models, or robust
within-donor definitions when complete data become available.

## 16. Clinical covariates are insufficient

**Risk:** Associations may be confounded by eGFR, fibrosis severity, IFTA,
treatment, age, or disease duration.

**Response strategy:** State that only broad diagnosis adjustment was
possible. Do not imply independence from severity. Request or identify cohorts
with systematic clinical and histologic variables for future validation.

## 17. No protein-level or spatial validation is provided

**Risk:** RNA abundance may not reflect secreted TIMP1 protein or spatial
relationships with fibrosis.

**Response strategy:** Frame this as a transcriptomic database study. Propose
urinary, plasma, and tissue protein measurement only as future work, without
claiming diagnostic utility. Recommend spatial transcriptomics,
immunohistochemistry, or multiplex imaging to localize TIMP1-associated niches.

## 18. No mechanistic perturbation was performed

**Risk:** Reviewers may ask whether TIMP1 drives fibrosis or maladaptive repair.

**Response strategy:** State that causality was outside the scope of the
analysis. Suggested future experiments include TIMP1 knockdown or
overexpression in injured tubular cells, co-culture with fibroblasts or
myeloid cells, matrix-turnover assays, and in vivo injury models. These are
required before mechanistic claims.

## 19. Multiple testing and analytical flexibility

**Risk:** Numerous signatures, genes, scopes, and sensitivity analyses increase
false-positive risk.

**Response strategy:** Emphasize prespecified signatures, within-family BH
correction, explicit relaxed and stringent module definitions, independent
external validation, and complete reporting of weaker results. Make all
analysis scripts and source tables available.

## 20. The manuscript may overuse the term “validation”

**Risk:** Etiologic and platform differences mean GSE180394 is not a strict
replication of an identical cohort.

**Response strategy:** Use “independent external transcriptomic validation” or
“external tubular-compartment validation” for the predefined analysis, while
acknowledging that it validates association direction and program coherence,
not longitudinal outcomes or clinical performance.

## Recommended response posture

The strongest response strategy is transparent rather than defensive:

- lead with disease-only and external results;
- distinguish robust direction from uniform statistical significance;
- describe enrichment and single-cell findings as supportive;
- acknowledge absent longitudinal and clinical covariates;
- avoid diagnostic, kidney-specific, direct-regulation, and causal language;
- provide scripts, gene sets, source tables, and session information.
