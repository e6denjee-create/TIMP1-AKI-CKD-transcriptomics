# Manuscript Results Skeleton

## TIMP1 is directionally elevated across kidney injury/fibrosis cohorts

TIMP1 expression was evaluated independently in GSE139061, GSE30718, and
GSE66494 without pooling absolute expression values across platforms. TIMP1
was directionally elevated in disease samples in all three cohorts, with
Hedges' g values of 0.415, 1.293, and 0.449, respectively. The within-cohort
BH-adjusted P values were 0.378, 0.000326, and 0.00905, indicating variation
in effect strength across datasets.

GSE180394 was then analyzed as an independent external microdissected tubular
kidney-disease cohort. TIMP1 expression was higher in 44 disease samples than
in nine healthy living-donor controls (difference = 1.318; Hedges' g = 1.569;
P = 8.01e-05). This external result supports reproducible elevation of TIMP1
in kidney injury/fibrosis-related transcriptomic datasets. It does not test
longitudinal AKI-to-CKD prediction or clinical diagnostic performance.

## TIMP1 correlates with injury-repair and ECM-remodeling signatures

Nine prespecified signatures were evaluated: ECM remodeling, collagen
formation, TGF-beta signaling, inflammation, tubular injury, maladaptive
repair, cellular senescence, fibrosis, and immune activation. TIMP1 was
excluded from signature score calculation.

Across the three discovery cohorts, all 27 all-sample correlations were
positive and FDR-significant, and 25 of 27 disease-only correlations were
FDR-significant. In GSE180394, all nine signatures were positively correlated
with TIMP1 in both all-sample and disease-only analyses. Disease-only rho
values ranged from 0.453 for collagen formation to 0.851 for fibrosis, and
all nine associations were FDR-significant.

These findings support TIMP1 as a candidate associated with injury-repair,
ECM-remodeling, inflammatory, senescence, and fibrosis-related programs.
They do not distinguish direct relationships from shared injury severity,
cell composition, or upstream regulation.

## A disease-only TIMP1-associated module is conserved across cohorts

Genome-wide disease-only correlations identified a relaxed 1,343-gene module
and a stringent 421-gene module supported by at least two discovery datasets.
In GSE180394, 399 stringent-module genes were detected. The module score was
strongly correlated with TIMP1 among all samples (rho = 0.879) and disease
samples alone (rho = 0.850).

This external replication supports a conserved TIMP1-associated transcriptomic
program. The module is a correlation-based construct and should not be
described as a regulatory network.

## Sensitivity analyses support robustness across diagnosis categories

GSE180394 disease samples were grouped into lupus nephritis, FSGS/FGGS,
diabetic nephropathy, IgA nephropathy, and other kidney disease. Each category
was removed in turn. All nine signatures, the stringent module score, and all
12 core genes retained positive correlations in every iteration. Twenty of
22 features remained FDR-significant in every iteration; FBLIM1 and PFN1 were
the weaker core-gene associations.

Diagnosis-adjusted analyses gave positive partial Spearman estimates for all
22 features. Twenty-one remained FDR-significant after adjustment; PFN1
remained positive but was not significant. Standardized linear models
including diagnosis category produced the same 21-of-22 FDR-supported pattern.

## Alternative control definitions retain the TIMP1 direction

The primary living-donor comparison produced Hedges' g = 1.569. Adding six
unaffected tumor-nephrectomy samples to the controls retained a positive and
significant effect (g = 1.291; P = 6.96e-05). Using tumor-nephrectomy samples
alone also retained a positive effect (g = 0.805), but the comparison was not
significant (P = 0.066), consistent with limited control sample size and the
imperfect nature of tumor-adjacent tissue as a healthy reference.

## A 12-gene core module is consistently correlated with TIMP1

The discovery core comprised `TMSB10`, `IFITM3`, `RAB31`, `MGP`, `PEA15`,
`TUBA1B`, `SERPING1`, `FBLIM1`, `TAGLN2`, `TGM2`, `TUBA1A`, and `PFN1`.
All 12 genes were positively correlated with TIMP1 in GSE180394 disease
samples; ten were FDR-significant in the unadjusted external analysis. After
diagnosis adjustment, 11 remained FDR-significant, with PFN1 retaining a weak
positive association.

The core genes provide candidate context involving cytoskeletal organization,
membrane trafficking, inflammation, complement regulation, wound repair, and
matrix-related biology. These associations do not establish direct regulation
by TIMP1.

## Single-cell analysis provides exploratory localization only

The available GSE210622 single-cell object contains one donor and provides
exploratory localization only. GSE267242 remains incomplete locally. These
data do not provide independent patient-level pseudobulk replication, and
single-donor cell-level P values cannot establish clinical relevance.

## Limitations

1. GSE180394 is cross-sectional and does not follow patients from AKI to CKD.
2. Its disease group is etiologically heterogeneous.
3. The primary external control group contains nine living donors.
4. Tumor-nephrectomy tissue is not an equivalent healthy reference.
5. Systematic eGFR, IFTA, fibrosis, and AKI outcome metadata are unavailable.
6. Bulk tubular expression cannot fully resolve cell-state or composition effects.
7. Signatures overlap biologically and are not statistically independent.
8. Complete database enrichment remains constrained by locally available resources.
9. Single-cell evidence remains exploratory and donor-limited.
10. Experimental perturbation and spatial validation are unavailable.

The current evidence supports TIMP1 as a candidate associated with
AKI-to-CKD-relevant injury-repair and ECM-remodeling programs. It does not
support kidney specificity, clinical diagnostic performance, longitudinal
prediction, direct regulation, or causal activity.
