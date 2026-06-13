# Manuscript First Draft v1

## Title

**Cross-cohort transcriptomic analysis associates TIMP1 with tubular
injury-repair, extracellular matrix remodeling, and a conserved disease-state
program in human kidney disease**

## Authors and affiliations

[Author 1], [Author 2], [Author 3], and [Corresponding Author]

[Affiliation 1]

[Affiliation 2]

Correspondence: [name and email]

## Abstract

### Background

Maladaptive repair after kidney injury is accompanied by persistent tubular
stress, inflammation, cellular senescence, and extracellular matrix (ECM)
remodeling. Tissue inhibitor of metalloproteinases 1 (TIMP1) is responsive to
injury in multiple tissues, but its reproducible transcriptomic context across
human kidney disease cohorts remains incompletely defined.

### Methods

We performed an integrative, platform-specific transcriptomic association
study using three discovery bulk cohorts (GSE139061, GSE30718, and GSE66494)
and an independent external cohort of microdissected human kidney tubules
(GSE180394). TIMP1 expression was compared between disease and control groups
within each cohort. Nine prespecified signatures representing ECM remodeling,
collagen formation, TGF-beta signaling, inflammation, tubular injury,
maladaptive repair, cellular senescence, fibrosis, and immune activation were
scored after excluding TIMP1. Disease-only Spearman correlations were used to
identify cross-cohort TIMP1-associated modules. External validation included
leave-one-diagnosis-out analysis, alternative control definitions, and
diagnosis-adjusted partial Spearman correlations and linear models.

### Results

TIMP1 was directionally elevated in all three discovery cohorts, with Hedges'
g values of 0.415, 1.293, and 0.449. In GSE180394, TIMP1 was higher in 44
kidney-disease samples than in nine healthy living-donor controls
(difference = 1.318; Hedges' g = 1.569; P = 8.01 x 10^-5). Across the
discovery cohorts, all 27 all-sample TIMP1-signature correlations were
positive and false discovery rate (FDR)-significant, and 25 of 27 disease-only
correlations were FDR-significant. In GSE180394, all nine signatures were
positively and significantly correlated with TIMP1 in both analysis scopes.
A 421-gene stringent disease-only module was identified in discovery
analyses; 399 genes were detected in GSE180394, where its score correlated
with TIMP1 in all samples (rho = 0.879) and disease samples (rho = 0.850).
All 12 core genes were positively correlated with TIMP1 in GSE180394, with
10 of 12 reaching disease-only FDR significance. All 22 evaluated features
remained positively correlated in leave-one-diagnosis-out analyses, and
21 of 22 remained FDR-significant after diagnosis adjustment.

### Conclusions

Across human kidney disease cohorts, TIMP1 was associated with tubular
injury-repair, ECM-remodeling, inflammatory activation, senescence,
fibrosis-related programs, and a conserved disease-only transcriptomic
module. These findings nominate TIMP1 as a candidate marker of a broader
maladaptive injury-response state and provide a reproducible framework for
future spatial and experimental validation. The study does not establish
longitudinal prediction, clinical diagnostic performance, kidney specificity,
or causality.

**Keywords:** TIMP1; acute kidney injury; chronic kidney disease; tubular
injury; maladaptive repair; extracellular matrix; fibrosis; transcriptomics

## Introduction

Acute kidney injury (AKI) and chronic kidney disease (CKD) are connected by a
continuum of incomplete or maladaptive repair. After an acute insult, injured
tubular epithelial cells may recover, remain in persistent stress states, or
adopt repair-associated phenotypes characterized by altered epithelial
identity, inflammatory signaling, cell-cycle dysregulation, senescence-like
features, and communication with immune and stromal compartments
[REFS]. These responses can coexist with extracellular matrix deposition and
remodeling, contributing to progressive tubulointerstitial dysfunction
[REFS].

Transcriptomic studies have identified recurring injury-response programs in
human kidney tissue, but individual cohorts differ in disease etiology,
compartment, platform, and control selection. Such heterogeneity complicates
the interpretation of isolated differentially expressed genes. A
cross-dataset strategy that prioritizes directionally reproducible expression,
pathway-level associations, and disease-only correlations can help distinguish
generalizable injury-associated programs from effects driven solely by
disease-control separation.

Tissue inhibitor of metalloproteinases 1 (TIMP1) is a secreted inhibitor of
several matrix metalloproteinases and participates in matrix turnover,
cell-survival signaling, inflammatory responses, and tissue repair in diverse
biological contexts [REFS]. Increased TIMP1 expression has been reported in
injured and fibrotic tissues, including the kidney [REFS]. However, increased
expression alone does not determine whether TIMP1 reflects ECM remodeling,
inflammatory activation, altered cell composition, persistent tubular injury,
or a combination of these processes. Nor does a cross-sectional association
establish direct regulation of fibrosis by TIMP1.

We therefore investigated the transcriptomic context of TIMP1 across human
kidney injury and fibrosis datasets. Three bulk cohorts were used for
discovery, and GSE180394, a microdissected tubular dataset, was used as an
independent external validation cohort. We evaluated TIMP1 expression,
prespecified injury-repair and ECM-related signatures, a cross-dataset
disease-only correlation module, and a 12-gene core module. Robustness was
examined by removing diagnosis categories in turn, changing the control
definition, and adjusting disease-only associations for diagnosis category.
Exploratory single-cell data were used only to assess localization and were
not treated as patient-level validation.

## Methods

### Study design and cohort selection

This study was designed as an integrative transcriptomic association analysis.
Publicly available human kidney datasets were selected when they provided a
processed expression matrix, sample-level disease annotations, an identifiable
control group, and sufficient gene coverage for TIMP1 and prespecified
signatures. Absolute expression values were not pooled across platforms.

The discovery analysis included GSE139061 (39 AKI and 9 control samples),
GSE30718 (28 AKI and 11 control samples), and GSE66494 (53 CKD/fibrosis and
8 control samples). GSE180394 was analyzed independently as an external
microdissected tubular cohort. Its primary analysis included 44 kidney-disease
samples and nine healthy living-donor controls. Six unaffected
tumor-nephrectomy samples were reserved for control sensitivity analyses.

### Expression preprocessing

For each discovery dataset, the locally prepared normalized gene-by-sample
matrix and matched metadata were read from the reproducible project data
structure. Expression matrices contained 20,139 genes for GSE139061, 21,755
genes for GSE30718, and 19,553 genes for GSE66494. Samples absent from either
the expression matrix or metadata were excluded and recorded in the
missing-data log.

For GSE180394, the processed GEO expression object and GPL19983 annotation
were used. Entrez gene identifiers were mapped to current human gene symbols
using the NCBI human `gene_info` table. Probe rows mapping to the same gene
symbol were averaged, yielding 24,845 mapped gene symbols. Analyses used the
provided within-dataset normalized expression values. No cross-platform
batch correction or pooling of absolute expression was performed.

### TIMP1 disease-control comparisons

TIMP1 expression was analyzed separately in each cohort. For the three
discovery cohorts, disease-control differences were evaluated using
two-sided Welch t-tests. GSE180394 comparisons used two-sided Wilcoxon
rank-sum tests because of the small and unequal control groups. Group
differences were calculated as the disease mean minus the control mean.
Standardized effects were summarized as Hedges' g using a small-sample
correction applied to the pooled-standard-deviation effect estimate.
Discovery-cohort expression P values were adjusted across the three cohort
comparisons using the Benjamini-Hochberg procedure.

### Prespecified signature scoring

Nine biologically motivated gene sets were defined before external validation:
ECM remodeling, collagen formation, TGF-beta signaling, inflammation, tubular
injury, maladaptive repair, cellular senescence, fibrosis, and immune
activation. The complete gene membership is available in
`signature_gene_sets_used.csv`.

Within each dataset, expression for every available signature gene was
standardized across samples to a gene-level z-score. Signature scores were
calculated as the mean of the available standardized genes for each sample.
At least two genes were required to calculate a score. TIMP1 was explicitly
excluded from every signature to prevent mathematical self-correlation.
Gene-set coverage was recorded for each dataset.

### TIMP1-signature correlation analysis

Spearman rank correlations were calculated between TIMP1 expression and each
signature score in two scopes: all samples and disease-only samples. The
disease-only analysis was prioritized for biological interpretation because
it reduces correlations caused only by separation of controls from disease
samples. P values were adjusted by the Benjamini-Hochberg method within each
dataset and analysis scope. FDR < 0.05 was considered statistically supported.

### Disease-only TIMP1-associated modules

Genome-wide Spearman correlations between TIMP1 and all variable genes were
calculated separately among disease samples in each discovery cohort. TIMP1
itself was removed from candidate module genes. For each gene, the correlation
direction, nominal P value, and BH-adjusted P value were retained by dataset.

The relaxed module required a positive correlation direction, nominal
P < 0.05 in at least two discovery datasets, and a positive median rho across
available cohorts. The stringent module required positive correlations with
BH-adjusted P < 0.05 in at least two datasets and a positive median rho.
These definitions produced relaxed and stringent modules of 1,343 and 421
genes, respectively. Module sensitivity analyses repeated interpretation
after removing all prespecified signature genes and after removing
ECM/fibrosis signature markers.

Over-representation analyses were attempted for Hallmark, Gene Ontology,
Reactome, and KEGG collections. Because complete online collections and
annotation packages were not available in the project environment, available
results were restricted to a documented curated offline subset and were
treated as exploratory.

### Core module definition

A 12-gene core was defined by positive TIMP1 correlations with FDR < 0.05 in
all three discovery cohorts. The genes were `TMSB10`, `IFITM3`, `RAB31`,
`MGP`, `PEA15`, `TUBA1B`, `SERPING1`, `FBLIM1`, `TAGLN2`, `TGM2`, `TUBA1A`,
and `PFN1`. Their functional descriptions were used as contextual annotations
only and not as evidence of direct regulation by TIMP1.

### GSE180394 external validation

The external workflow repeated the TIMP1 disease-control comparison, nine
signature scores, all-sample and disease-only correlations, stringent-module
scoring, and 12-core-gene correlations. Of 421 stringent genes, 399 were
available for the principal input coverage assessment. All 12 core genes were
detected. P values for correlations were BH-adjusted within each analysis
family and scope.

### Leave-one-diagnosis-out sensitivity analysis

The 44 GSE180394 disease samples were harmonized into five diagnosis
categories: lupus nephritis, FSGS/FGGS, diabetic nephropathy, IgA nephropathy,
and other kidney disease. Each category was excluded in turn. Disease-only
Spearman correlations were then recalculated for the nine signatures, the
stringent module score, and the 12 core genes. BH adjustment was applied
across the 22 features within each iteration.

### Alternative control sensitivity analysis

The primary GSE180394 reference consisted of nine healthy living-donor
samples. Two alternative definitions were assessed: an extended control group
combining the living donors with six unaffected tumor-nephrectomy samples,
and a tumor-nephrectomy-only group. For each definition, the TIMP1 group
difference, Hedges' g, approximate 95% confidence interval, Wilcoxon P value,
and BH-adjusted P value were calculated. The living-donor comparison remained
the prespecified primary analysis because tumor-adjacent tissue may not
represent a healthy reference.

### Diagnosis-adjusted disease-only associations

For each of the 22 external validation features, three disease-only estimates
were calculated: an unadjusted Spearman correlation, a partial Spearman
correlation between residuals obtained after regressing TIMP1 and the feature
separately on diagnosis category, and a standardized linear model of the form
`scale(feature) ~ scale(TIMP1) + diagnosis category`. BH correction was
applied separately to each set of 22 tests.

### Exploratory single-cell localization

The locally available GSE210622 object contained one AKI donor and
computational cell-type labels. It was used only for exploratory localization
of TIMP1 and an exploratory TIMP1-high versus TIMP1-low tubular-cell
comparison. GSE267242 lacked a complete matching expression matrix and
cell-level metadata. Because independent donor-level replication was not
available, cell-level P values were not interpreted as patient-level or
clinical evidence.

### Statistical software and reproducibility

Analyses were performed in R 4.6.0 using scripted workflows. Random seeds were
set before stochastic operations. Key packages included GEOquery, Biobase,
ggplot2, pheatmap, igraph, and msigdbr; package and session versions were
recorded in project outputs. Statistical tables and intermediate results were
saved as CSV files, and figures were saved as both PDF and 600-dpi PNG files.
Missing optional data or resources were recorded in
`missing_data_log.csv` without silently omitting the affected analysis.

## Results

### TIMP1 is directionally elevated across kidney injury/fibrosis cohorts

The three discovery cohorts represented two AKI comparisons and one
CKD/fibrosis comparison (Figure 1). TIMP1 expression was higher in disease
samples in all three datasets when each cohort was analyzed independently
(Figure 2). The disease-control differences were 1.086 in GSE139061, 1.062 in
GSE30718, and 0.417 in GSE66494. Corresponding Hedges' g values were 0.415,
1.293, and 0.449. The BH-adjusted P values were 0.378, 0.000326, and 0.00905,
respectively. Thus, TIMP1 showed a consistent direction across cohorts,
although the magnitude and statistical support varied.

### TIMP1 correlates with injury-repair and ECM-remodeling signatures

We next evaluated nine prespecified signatures spanning matrix remodeling,
collagen formation, TGF-beta signaling, inflammation, tubular injury,
maladaptive repair, senescence, fibrosis, and immune activation. TIMP1 was
excluded from all signature scores.

Across the three discovery cohorts, all 27 all-sample correlations were
positive and FDR-significant. Among disease samples alone, 25 of 27
correlations remained FDR-significant (Figure 3). The two associations not
reaching FDR < 0.05 were collagen formation and cellular senescence in
GSE139061, but both retained positive directions. The persistence of most
associations in disease-only analyses indicates that the observed pattern was
not explained solely by disease-control separation.

### A disease-only TIMP1-associated module is conserved across cohorts

Genome-wide disease-only correlations identified 1,343 relaxed module genes
with replicated nominal support and 421 stringent module genes with replicated
FDR support in at least two discovery datasets. The associated genes included
injury-response, cytoskeletal, membrane, inflammatory, adhesion, wound-response,
and matrix-related components.

Removing all prespecified signature genes retained 1,291 relaxed and 389
stringent genes. Removing ECM/fibrosis markers retained 1,323 relaxed and 414
stringent genes. In the curated enrichment subset, the residual modules
remained associated with response to wounding, inflammatory programs,
membrane and cytoskeletal biology, focal adhesion, senescence, and selected
matrix-related terms. These sensitivity results support a broader
TIMP1-associated injury-response program rather than a module defined only by
the prespecified ECM genes.

### GSE180394 validates the tubular TIMP1-associated program

GSE180394 provided independent microdissected tubular tissue from 44
kidney-disease samples and nine healthy living-donor controls. TIMP1
expression was higher in disease samples (difference = 1.318; Hedges'
g = 1.569; P = 8.01 x 10^-5; Figure 2).

All nine signatures were positively and FDR-significantly correlated with
TIMP1 in both all-sample and disease-only analyses. Disease-only rho values
ranged from 0.453 for collagen formation to 0.851 for fibrosis. The
stringent-module score was strongly correlated with TIMP1 in all samples
(rho = 0.879) and among disease samples alone (rho = 0.850; Figure 4).
These results reproduce the pathway-level and module-level associations in a
renal tubular compartment without pooling expression values with the
discovery cohorts.

### Sensitivity analyses support robustness across diagnosis and control definitions

All 22 evaluated features retained positive TIMP1 correlations in every
leave-one-diagnosis-out iteration. The nine signatures and stringent module
remained FDR-significant in all five iterations. Ten of 12 core genes were
FDR-significant in every iteration, whereas FBLIM1 and PFN1 showed weaker but
consistently positive associations. The stringent module had a median
leave-one-diagnosis-out rho of 0.851, with a range of 0.834 to 0.863.

The primary living-donor contrast yielded Hedges' g = 1.569. Adding six
unaffected tumor-nephrectomy samples to the controls retained a positive,
significant effect (difference = 1.070; g = 1.291; P = 6.96 x 10^-5).
The tumor-nephrectomy-only comparison was also positive but less precise
(difference = 0.698; g = 0.805; P = 0.0664).

After adjustment for diagnosis category, all 22 partial Spearman estimates
remained positive and 21 were FDR-significant. All nine signatures remained
supported, with partial rho values from 0.438 to 0.823. The adjusted partial
rho for the stringent module was 0.821. Standardized linear models including
diagnosis category produced the same 21-of-22 FDR-supported pattern. PFN1
remained positive but non-significant.

### A 12-gene core module supports a cytoskeletal-inflammatory-matrix repair context

Twelve genes were positively correlated with TIMP1 at FDR < 0.05 in each of
the three discovery cohorts: `TMSB10`, `IFITM3`, `RAB31`, `MGP`, `PEA15`,
`TUBA1B`, `SERPING1`, `FBLIM1`, `TAGLN2`, `TGM2`, `TUBA1A`, and `PFN1`
(Figure 5). Their discovery median rho values ranged from approximately 0.56
to 0.77.

All 12 genes retained positive disease-only correlations in GSE180394, and
10 of 12 were FDR-significant in the unadjusted external analysis. After
diagnosis adjustment, 11 remained FDR-significant; PFN1 retained a weak
positive association. Collectively, the core genes provide a reproducible
context involving actin and microtubule organization, membrane trafficking,
interferon and complement-associated responses, adhesion, wound repair, and
matrix-related biology. These correlations do not imply direct regulation by
TIMP1.

### Single-cell analysis provides exploratory localization only

The available GSE210622 single-cell object contained one AKI donor. It
provided exploratory localization of TIMP1 across computationally annotated
cell types and an exploratory tubular TIMP1-high program that included
matrix, wound-healing, maladaptive-repair, TGF-beta-related, and
senescence-related signatures. However, the analysis lacked independent
donor-level replication, and the tubular 75th percentile of TIMP1 expression
was zero, making the high-low comparison sensitive to detection and dropout.
GSE267242 was incomplete locally. These findings were therefore treated as
supportive localization only.

## Discussion

In this integrative analysis, TIMP1 was reproducibly associated with a broad
kidney injury-response context across independent human datasets. The evidence
extended beyond directional expression: TIMP1 correlated with prespecified
tubular injury, maladaptive repair, ECM remodeling, collagen, TGF-beta,
inflammation, immune activation, senescence, and fibrosis programs. A
disease-only correlation strategy identified a conserved module that was
externally reproduced in microdissected tubules. Together, these findings
position TIMP1 as a candidate indicator of an injury-repair and matrix
remodeling state rather than as an isolated disease-control marker.

Persistent tubular stress is increasingly recognized as an important feature
of maladaptive repair. Injured tubular cells may retain altered epithelial
identity, inflammatory output, and incomplete cell-cycle recovery while
communicating with infiltrating immune cells and matrix-producing stromal
cells [REFS]. The strong correlations between TIMP1 and tubular injury,
maladaptive repair, inflammation, and senescence signatures are compatible
with this framework. Nevertheless, bulk expression cannot determine whether
TIMP1 arises from a specific tubular state, other resident cells, infiltrating
cells, or coordinated changes across compartments.

The ECM-related associations were similarly coherent across analysis levels.
TIMP1 correlated with ECM remodeling, collagen formation, fibrosis, and
TGF-beta-related scores, and the cross-dataset module included matrix,
adhesion, and wound-response context. TIMP1 can inhibit metalloproteinase
activity, but transcriptomic co-expression cannot determine the balance
between matrix deposition and degradation or establish a direct effect on
fibrosis. A more conservative interpretation is that TIMP1 marks or
participates in a tissue environment in which matrix turnover and repair
programs are active.

The 12-gene core module adds biological breadth. TMSB10, PFN1, TAGLN2, and
tubulin genes are linked to cytoskeletal organization; RAB31 is associated
with membrane trafficking; IFITM3 and SERPING1 provide inflammatory,
interferon, and complement-related context; FBLIM1 and PEA15 relate to
adhesion or cellular signaling; and MGP and TGM2 provide matrix and
wound-response context. The cross-cohort consistency of these genes suggests
that TIMP1 is embedded in a coordinated cytoskeletal-inflammatory-matrix
repair program. These annotations are contextual and do not demonstrate that
TIMP1 regulates the core genes.

GSE180394 materially strengthens the study because it is independent of the
three discovery cohorts and focuses on microdissected human tubules. It
replicated TIMP1 elevation, all nine signature associations, the stringent
module score, and positive correlations for all 12 core genes. The
disease-only module correlation remained strong, and sensitivity analyses
showed that the pattern was not dependent on retaining any one major diagnosis
category. The diagnosis-adjusted analyses further reduced, but did not remove,
concern that etiologic composition explained the findings.

Disease-only analysis was central to the design. Correlations calculated
across controls and diseased tissue can be inflated when both variables
separate the same two groups. Restricting correlations to disease samples asks
whether TIMP1 tracks variation in injury-associated programs within the
disease state. The persistence of most signature, module, and core-gene
associations in this setting provides stronger association evidence than an
all-sample analysis alone, although residual confounding by severity, therapy,
cell composition, and unmeasured clinical variables remains possible.

Future work should focus on orthogonal validation. Multi-donor single-cell or
single-nucleus datasets could determine whether TIMP1-high tubular states
replicate at the patient level and how they relate to fibroblast and myeloid
programs. Spatial transcriptomics or multiplex tissue imaging could test
whether TIMP1 expression co-localizes with injured tubules, inflammatory
niches, and areas of matrix remodeling. Longitudinal cohorts with AKI
outcomes, eGFR trajectories, and histologic fibrosis measures would be needed
to assess temporal or prognostic relevance. Finally, perturbation experiments
in tubular, stromal, or multicellular systems would be required to test
whether changing TIMP1 alters matrix turnover, inflammatory signaling,
senescence-associated phenotypes, or repair outcomes.

## Limitations

First, all bulk cohorts were cross-sectional and did not longitudinally follow
the same patients from AKI to CKD. Second, the datasets differed in platform,
tissue context, disease etiology, and control definition. Platform-specific
analysis reduced inappropriate scale pooling but did not eliminate biological
heterogeneity. Third, GSE180394 contained only nine living-donor controls and
an etiologically heterogeneous disease group. The six tumor-nephrectomy
samples provided a useful sensitivity reference but are not equivalent to
healthy donor tissue.

Fourth, systematic sample-level eGFR, IFTA, fibrosis, treatment, and AKI
outcome metadata were unavailable for external adjustment. Fifth, bulk
transcriptomic associations may reflect cell composition as well as
within-cell-state regulation. Sixth, the prespecified signatures overlap
biologically and should not be considered independent pathways. Seventh,
complete Hallmark, GO, Reactome, and KEGG enrichment could not be obtained in
the available environment; curated enrichment results were exploratory.
Eighth, the single-cell analysis contained one donor and cannot provide
patient-level validation. Finally, no experimental perturbation, protein-level
validation, spatial validation, or clinical prediction analysis was performed.

## Conclusions

TIMP1 is reproducibly elevated in the evaluated human kidney injury and
fibrosis datasets and is associated with tubular injury-repair,
ECM-remodeling, inflammatory activation, senescence, fibrosis-related
programs, and a conserved disease-only transcriptomic module. Independent
validation in microdissected human kidney tubules and robustness across
diagnosis and control definitions strengthen this association-based evidence.
TIMP1 should be considered a candidate for further investigation within a
broader maladaptive repair and wound-response program. Longitudinal,
spatial, protein-level, and experimental studies are required before making
claims about prediction, clinical utility, tissue specificity, or causality.

## Figure legends

### Figure 1. TIMP1 transcriptomic validation workflow and cohort overview

Overview of the three discovery bulk cohorts, the independent GSE180394
microdissected tubular cohort, robustness analyses, and the expression,
signature, module, and core-gene evidence used for manuscript interpretation.

### Figure 2. TIMP1 expression across discovery and external cohorts

Violin plots, boxplots, and individual samples show within-cohort normalized
TIMP1 expression in GSE139061, GSE30718, GSE66494, and GSE180394. Disease and
control labels were harmonized for display. Expression scales are
platform-specific and should not be compared directly across cohorts.

### Figure 3. Disease-only correlations between TIMP1 and pathway signatures

Heatmap of Spearman correlations between TIMP1 and nine prespecified pathway
scores among disease samples in each cohort. TIMP1 was excluded from all
signature scores. Positive disease-only correlations support associations
with injury-repair, ECM-remodeling, inflammatory, senescence, and
fibrosis-related programs.

### Figure 4. TIMP1 is correlated with the stringent module score in GSE180394

Scatterplot of TIMP1 expression and the stringent module score, based on genes
detected in GSE180394. Colors indicate disease and healthy living-donor
samples. The all-sample and disease-only Spearman correlations were 0.879 and
0.850, respectively. The fitted line is descriptive and does not imply direct
regulation.

### Figure 5. Cross-cohort correlations for the 12-gene core module

Heatmap of disease-only Spearman correlations between TIMP1 and each core
module gene across the three discovery cohorts and GSE180394. All 12 genes
retained positive correlations in the external cohort, supporting a conserved
candidate program without establishing causality.

## Data availability

The public transcriptomic datasets analyzed in this study are available from
the NCBI Gene Expression Omnibus under accession numbers GSE139061,
GSE30718, GSE66494, GSE180394, GSE210622, and GSE267242. The latter two were
used or audited only for exploratory single-cell localization status. Derived
tables, feature definitions, and processed analysis outputs are stored in the
project directories `results/timp1_validation/` and
`figures/timp1_validation/`. Any repository or archival DOI for the final
processed outputs will be added before publication.

## Code availability

Reproducible R scripts are available in the project `scripts/` directory.
Principal workflows include `run_timp1_bulk_validation.R`,
`run_timp1_module_network.R`, `run_timp1_module_sensitivity_v2.R`,
`prepare_external_GSE180394.R`,
`run_external_bulk_validation_template.R`, and
`run_external_GSE180394_sensitivity_v4.R`. Package and session versions are
recorded with the results. A public code repository URL and release identifier
will be added before publication.

## Author contributions

Conceptualization: [AUTHOR INITIALS]. Data curation: [AUTHOR INITIALS].
Methodology: [AUTHOR INITIALS]. Software: [AUTHOR INITIALS]. Formal analysis:
[AUTHOR INITIALS]. Visualization: [AUTHOR INITIALS]. Writing - original draft:
[AUTHOR INITIALS]. Writing - review and editing: [AUTHOR INITIALS].
Supervision: [AUTHOR INITIALS]. Funding acquisition: [AUTHOR INITIALS].

## Funding

[Insert funding sources and grant numbers, or state that no specific funding
was received.]

## Conflict of interest

The authors declare [no competing interests / the following competing
interests: INSERT DETAILS].

## References

[Insert formatted references after target-journal selection. Literature
citations marked `[REFS]` in this first draft require final reference
curation.]
