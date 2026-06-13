# Candidate External Bulk Dataset Strategy

## Current project inventory

No independent fourth bulk expression cohort was found in the current project.
The available bulk matrices and GEO series files correspond to the three
datasets already analyzed:

- GSE139061
- GSE30718
- GSE66494

Older tables and figures in `AKI_CKD_ECM_LOXL4/` are derived from these same
cohorts and should not be counted as independent validation datasets.

## Priority dataset types

### 1. Acute kidney injury

Preferred cohorts include human kidney biopsy, nephrectomy-adjacent tissue, or
well-annotated urinary-cell transcriptomes collected during clinically defined
AKI. Useful contrasts include AKI versus non-injured controls, resolving versus
non-resolving AKI, and early versus persistent injury.

### 2. AKI-to-CKD transition

Highest-priority cohorts contain longitudinal samples or documented follow-up
after AKI, including later eGFR decline, persistent dysfunction, fibrosis, CKD
diagnosis, or maladaptive repair outcomes. These datasets are most relevant to
the project but are expected to be less common than cross-sectional cohorts.

### 3. Chronic kidney disease and renal fibrosis

Suitable cohorts include biopsy-based CKD studies with histological fibrosis,
interstitial fibrosis/tubular atrophy, eGFR, or progression information.
Datasets should allow TIMP1 associations with ECM and repair programs to be
tested without treating CKD as equivalent to prior AKI.

### 4. Diabetic kidney disease

Include glomerular or tubulointerstitial transcriptomes with diabetic kidney
disease and non-diabetic controls. Compartment-specific data are preferred
because TIMP1 associations may differ between glomerular, vascular, immune, and
tubulointerstitial compartments.

### 5. Kidney transplant injury

Relevant studies include acute tubular injury, delayed graft function,
rejection-negative injury, antibody-mediated rejection, T-cell-mediated
rejection, and interstitial fibrosis/tubular atrophy. Injury and rejection
phenotypes should be analyzed separately where possible.

### 6. Tubulointerstitial injury

Prioritize datasets enriched for tubular injury, failed repair, interstitial
inflammation, or fibrosis. Microdissected tubulointerstitial samples are
especially useful for testing whether TIMP1 associations persist outside
whole-kidney cell-composition effects.

## Inclusion criteria

1. Human kidney-related transcriptomic data with a publicly accessible
   expression matrix or raw data that can be reproducibly processed.
2. At least two clinically or pathologically interpretable groups.
3. Preferably at least 10 biological samples per principal comparison group;
   smaller cohorts may be retained as exploratory evidence.
4. Sample-level phenotype metadata sufficient to define disease, control,
   injury subtype, compartment, or outcome.
5. TIMP1 and a reasonable proportion of the project signature genes must be
   measurable on the platform.
6. Biological samples, not technical replicates, must define the statistical
   sample size.
7. Human gene identifiers must be mappable to current gene symbols.
8. Normalization and batch structure must be documented or recoverable.
9. For repeated or longitudinal samples, patient ID and time point must be
   available.
10. For biopsy cohorts, histological fibrosis, IFTA, tubular injury, rejection,
    or compartment annotation is strongly preferred.

## Exclusion criteria

1. Cell lines or isolated stimulation experiments without kidney tissue
   validation.
2. Animal-only cohorts for the primary human evidence analysis.
3. Datasets without usable sample-level metadata.
4. Studies in which cases and controls are completely confounded by platform,
   processing batch, or tissue compartment.
5. Duplicate publications or reprocessed versions of an already included
   cohort unless used solely for technical sensitivity analysis.
6. Cohorts with fewer than four independent biological samples in a required
   analysis stratum.
7. Data for which TIMP1 is absent or cannot be reliably mapped.

## Recommended extraction fields

For each candidate GEO dataset, record:

- GEO accession and publication;
- tissue and kidney compartment;
- disease and control definitions;
- number of independent donors per group;
- longitudinal or cross-sectional design;
- platform and expression-data format;
- availability of raw and processed matrices;
- phenotype, donor, batch, and clinical covariates;
- fibrosis, eGFR, IFTA, rejection, or AKI outcome measures;
- TIMP1 and signature-gene coverage;
- planned role: discovery, independent validation, or sensitivity analysis.

## Analysis principles

Each dataset should be analyzed independently. Absolute expression should not
be pooled across platforms. Report within-dataset effect sizes and pathway
associations, then integrate standardized effects or correlation coefficients
using meta-analysis where appropriate. TIMP1 should be described as a
candidate associated with kidney injury-repair and ECM-remodeling programs,
not as a kidney-specific marker or a causal factor.
