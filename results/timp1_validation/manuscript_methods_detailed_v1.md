# Detailed Methods v1

## Study design

The analysis was designed as a reproducible, integrative transcriptomic
association study of TIMP1 in human kidney injury and fibrosis. The primary
objectives were to evaluate whether TIMP1 was directionally elevated across
independent cohorts, whether its expression tracked coherent injury-repair and
ECM-related programs within disease samples, and whether a cross-dataset
TIMP1-associated module could be reproduced in an independent tubular cohort.
The study did not evaluate longitudinal prediction, diagnostic classification,
or causal effects.

## Cohort selection and analysis roles

Datasets were eligible when they contained human kidney tissue, a processed
gene-expression matrix, sample-level disease annotations, an identifiable
control group, and sufficient gene coverage for TIMP1 and pathway scoring.
Datasets were assigned distinct analysis roles before final interpretation:

| Cohort | Analysis role | Tissue context | Disease | Disease n | Control n |
|---|---|---|---|---:|---:|
| GSE139061 | Discovery | Human kidney expression dataset | AKI | 39 | 9 |
| GSE30718 | Discovery | Renal allograft biopsy | AKI | 28 | 11 |
| GSE66494 | Discovery | Renal biopsy | CKD/fibrosis | 53 | 8 |
| GSE180394 | External validation | Microdissected human kidney tubules | Heterogeneous kidney disease | 44 | 9 |

Six unaffected tumor-nephrectomy samples in GSE180394 were excluded from the
primary contrast and retained for sensitivity analyses. GSE210622 and
GSE267242 were audited for single-cell validation. Only one GSE210622 AKI
donor was complete locally, and GSE267242 lacked a matching expression matrix
and cell-level metadata. These resources were therefore not used as
patient-level validation cohorts.

## Discovery expression matrices and metadata

The discovery workflows read normalized gene-by-sample matrices from:

- `TIMP1_AKI_CKD_project/data/processed/GSE139061_normalized_expression.csv.gz`
- `TIMP1_AKI_CKD_project/data/processed/GSE30718_normalized_expression.csv.gz`
- `TIMP1_AKI_CKD_project/data/processed/GSE66494_normalized_expression.csv.gz`

Matched metadata were read from
`TIMP1_AKI_CKD_project/data/metadata/`. The matrices contained 20,139,
21,755, and 19,553 genes, respectively. Required metadata columns were
`sample` and `group`. Samples not shared between the expression matrix and
metadata were excluded and recorded. Gene identifiers were represented as
gene symbols in the validation matrices.

No absolute expression values were merged across cohorts. All transformations,
scores, tests, and effect estimates were calculated within each dataset. This
approach avoided treating platform-dependent normalized expression scales as
directly comparable measurements.

## GSE180394 preprocessing

The GSE180394 processed GEO expression object, GPL19983 annotation, and NCBI
human `gene_info` file were used to construct a gene-symbol matrix. Probe
Entrez identifiers from GPL19983 were matched to current NCBI gene symbols.
Rows without a valid human symbol were removed. Symbols were converted to
uppercase, and probe rows mapping to the same symbol were averaged.

This procedure yielded 24,845 unique gene symbols across 59 samples. The
primary matrix retained 44 disease samples and nine healthy living-donor
controls. Six samples annotated as unaffected tissue from tumor nephrectomy
were retained in the all-sample matrix but excluded from the primary
comparison.

## TIMP1 expression analysis

TIMP1 had to be present in the gene-symbol matrix for a cohort to enter the
expression and correlation analyses. Group difference was defined as:

`mean(TIMP1 in disease) - mean(TIMP1 in control)`.

The three discovery comparisons used two-sided Welch t-tests. The GSE180394
primary and alternative-control comparisons used two-sided Wilcoxon rank-sum
tests with `exact = FALSE`.

Hedges' g was calculated from the pooled standard deviation:

`SDpooled = sqrt(((n1 - 1)SD1^2 + (n0 - 1)SD0^2) / (n1 + n0 - 2))`

`d = (mean1 - mean0) / SDpooled`

`J = 1 - 3 / (4(n1 + n0) - 9)`

`g = J x d`

For the GSE180394 control sensitivity analysis, an approximate standard error
and 95% confidence interval were also calculated. Discovery expression P
values were BH-adjusted across the three cohorts. Alternative-control P values
were adjusted across the three GSE180394 control definitions.

## Prespecified signature definitions

Nine signatures were defined before external validation:

1. ECM remodeling.
2. Collagen formation.
3. TGF-beta signaling.
4. Inflammation.
5. Tubular injury.
6. Maladaptive repair.
7. Cellular senescence.
8. Fibrosis.
9. Immune activation.

The exact gene membership is stored in
`results/timp1_validation/signature_gene_sets_used.csv`. TIMP1 was explicitly
excluded from every signature before scoring, even if introduced through a
future gene-set revision. This prevented a signature score from containing the
same expression value used as the correlation target.

## Signature score calculation

For each cohort and each signature:

1. Requested genes were intersected with the expression matrix.
2. At least two available genes were required.
3. Each available gene was standardized across samples:
   `z(gene, sample) = (expression - gene mean) / gene SD`.
4. Non-finite standardized values were set to missing.
5. The sample score was the mean standardized expression across available
   signature genes.
6. Requested, available, and missing genes were saved as a coverage table.

This scoring produced within-cohort relative scores and did not create a
shared absolute scale across platforms.

## All-sample and disease-only correlations

Spearman rank correlation was used because the objective was to test
monotonic association without assuming a linear expression relationship.
Each signature was correlated with TIMP1 in:

- all samples, including controls;
- disease samples only.

The disease-only analysis was prioritized because correlations across all
samples can arise when two variables independently separate disease from
control. Disease-only correlations test whether TIMP1 tracks program
variation within the diseased tissue set.

At least four complete, non-constant observations were required. P values were
BH-adjusted across the nine signatures within each dataset and scope.
FDR < 0.05 was used as the reporting threshold.

## Genome-wide TIMP1 correlation analysis

For each discovery cohort, metadata were restricted to disease samples.
Genes with zero variance were removed. Spearman rho was computed between
TIMP1 and every remaining gene. Ranks were generated with average handling of
ties. Approximate two-sided P values were obtained from the correlation
t-statistic, and BH correction was applied across all genes within each
dataset. TIMP1 itself was removed from the candidate gene table.

For every gene, the workflow recorded:

- rho, nominal P value, and BH-adjusted P value in each dataset;
- number of datasets with positive rho;
- number with positive rho and nominal P < 0.05;
- number with positive rho and FDR < 0.05;
- mean and median rho across available datasets.

## Relaxed and stringent module definitions

The relaxed TIMP1-associated module required:

- positive rho in at least two discovery datasets;
- positive rho with nominal P < 0.05 in at least two datasets;
- positive median rho.

The stringent module required:

- positive rho with BH-adjusted P < 0.05 in at least two datasets;
- positive median rho.

The relaxed and stringent definitions yielded 1,343 and 421 genes. These were
described as correlation modules or associated programs, not regulatory
networks.

## Module sensitivity analyses

Two exclusion-based sensitivity analyses tested whether module interpretation
was circular:

1. All genes included in any of the nine prespecified signatures were removed.
2. Genes included in the ECM remodeling, collagen formation, or fibrosis
   signatures were removed.

The resulting sizes were:

| Module variant | Genes |
|---|---:|
| Relaxed | 1,343 |
| Stringent | 421 |
| Relaxed without signature genes | 1,291 |
| Stringent without signature genes | 389 |
| Relaxed without ECM/fibrosis markers | 1,323 |
| Stringent without ECM/fibrosis markers | 414 |

## Enrichment analysis

The project attempted enrichment against Hallmark, GO biological process, GO
cellular component, Reactome, and KEGG using `msigdbr`,
`clusterProfiler`, `ReactomePA`, and `org.Hs.eg.db`. In the recorded
environment, `msigdbr` was installed but complete gene-set downloads were not
available, and the remaining annotation packages were unavailable. Failures
were recorded in `missing_data_log.csv`.

Exploratory over-representation analyses were therefore performed using a
documented curated offline subset. Hypergeometric P values were calculated
against the tested gene universe and corrected for multiple testing. These
results were used for broad contextual interpretation, not as a claim of
complete database enrichment.

## Twelve-gene core module

The core module required a positive TIMP1 correlation with FDR < 0.05 in all
three discovery cohorts. Twelve genes met this definition:

`TMSB10`, `IFITM3`, `RAB31`, `MGP`, `PEA15`, `TUBA1B`, `SERPING1`,
`FBLIM1`, `TAGLN2`, `TGM2`, `TUBA1A`, and `PFN1`.

Median rho, dataset-specific rho, support counts, and contextual functional
annotations were saved in `core_module_genes_summary.csv`. Functional
annotations were not used as statistical selection criteria.

## External validation scoring

The GSE180394 external workflow repeated the discovery scoring rules. Coverage
was assessed before analysis. TIMP1 and all 12 core genes were available.
Signature coverage ranged from 90.9% to 100%. The stringent module had high
coverage; 401 of 421 requested genes were available in the principal input
coverage table.

The stringent module score was calculated by standardizing each available
module gene across GSE180394 samples and taking the sample-level mean.
All-sample and disease-only Spearman correlations were then calculated between
TIMP1 and the score. The same approach was applied to the nine signatures.
Each core gene was correlated individually with TIMP1.

## Leave-one-diagnosis-out analysis

GSE180394 diagnoses were harmonized into:

- lupus nephritis;
- FSGS/FGGS;
- diabetic nephropathy;
- IgA nephropathy;
- other kidney disease.

For each category, all samples in that category were removed. Disease-only
Spearman correlations were recalculated for 22 features: nine signatures, one
stringent module score, and 12 core genes. P values were BH-adjusted across the
22 features within each iteration. The output retained excluded category,
excluded and retained sample counts, rho, nominal P value, and adjusted P
value.

## Alternative control definitions

Three GSE180394 contrasts were evaluated:

1. 44 disease samples versus 9 healthy living donors.
2. 44 disease samples versus 15 extended controls consisting of living donors
   plus unaffected tumor-nephrectomy samples.
3. 44 disease samples versus 6 unaffected tumor-nephrectomy samples.

The first comparison was primary. The latter two were sensitivity analyses.
The tumor-nephrectomy-only result was interpreted cautiously because adjacent
unaffected tissue may not represent a healthy baseline and the control sample
size was small.

## Diagnosis-adjusted correlations and linear models

The 44 disease samples provided sufficient representation to include the five
harmonized diagnosis categories as a covariate. For every feature, a complete
analysis frame was constructed with TIMP1, the feature, and diagnosis
category.

Partial Spearman correlation was implemented by:

1. fitting `TIMP1 ~ diagnosis category`;
2. fitting `feature ~ diagnosis category`;
3. extracting both residual vectors;
4. calculating Spearman correlation between residuals.

A complementary standardized linear model was fitted:

`scale(feature) ~ scale(TIMP1) + diagnosis category`.

The coefficient for standardized TIMP1, its standard error, and P value were
recorded. BH correction was applied across 22 features separately for raw
Spearman, partial Spearman, and linear-model tests.

## Exploratory single-cell analysis

The locally available GSE210622 resources included one expression matrix,
matching barcodes and features, a processed Seurat object, and computational
cell-type labels for one AKI donor. TIMP1 localization was summarized by
annotated cell type. The original MVP also compared tubular cells above the
within-population 75th percentile with those below the 25th percentile and
evaluated pathway-level expression differences.

This analysis was explicitly exploratory because:

- only one donor was available;
- cell-level observations were not independent patient replicates;
- the tubular TIMP1 75th percentile was zero;
- detection dropout affected the high-low definition;
- GSE267242 lacked a complete matching matrix and metadata.

No single-cell result was used as patient-level validation.

## Missing-data handling

A missing dataset, annotation, gene, optional package, or database did not
terminate independent analyses that could still be completed. The affected
unit was skipped, and the analysis step, dataset, missing item, reason, action,
and timestamp were appended to
`results/timp1_validation/missing_data_log.csv`.

## Statistical software

Analyses were performed in R 4.6.0. Recorded package versions included
GEOquery 2.80.0, limma 3.68.4, Seurat 5.5.0, fgsea 1.38.0, msigdbr 26.1.0,
pheatmap 1.0.13, ggplot2 4.0.3, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, and
Matrix 1.7.5. Each major workflow set a fixed random seed and saved
`sessionInfo()`.

## Reproducibility and output policy

All workflows were executed from the project root without hard-coded project
paths. Key intermediate and statistical results were written as CSV. Figures
were saved as PDF and 600-dpi PNG. Existing MVP outputs were preserved.
Validation outputs were written under `results/timp1_validation/` and
`figures/timp1_validation/`.

Principal commands:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_timp1_bulk_validation.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_timp1_module_network.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_timp1_module_sensitivity_v2.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/prepare_external_GSE180394.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_external_GSE180394_sensitivity_v4.R
```
