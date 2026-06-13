# Updated TIMP1 Validation Report v4

## Executive summary

GSE180394 provides an independent external microdissected tubular validation
cohort with 44 kidney-disease samples and nine healthy living-donor controls.
TIMP1 was elevated in disease, and the nine prespecified signatures, stringent
module score, and 12-gene core module showed broadly reproducible positive
associations. Diagnosis-removal and diagnosis-adjusted analyses indicate that
the main pattern is not attributable to one major disease category.

The evidence is suitable for beginning a manuscript draft centered on an
integrative transcriptomic association study. It does not support
longitudinal AKI-to-CKD prediction, clinical diagnostic performance, kidney
specificity, direct regulation, or causality.

## GSE180394 main validation

- Platform: GPL19983, Affymetrix Human Gene 2.1 ST Array.
- Tissue: microdissected human kidney tubules.
- Primary comparison: 44 disease samples versus nine living-donor controls.
- TIMP1 difference: 1.318.
- Hedges' g: 1.569.
- P value: 8.01e-05.
- All nine signatures were positively and significantly correlated with TIMP1
  in all-sample and disease-only analyses.
- Stringent module correlations were rho = 0.879 in all samples and rho =
  0.850 among disease samples.
- All 12 core genes were positively correlated with TIMP1 among disease
  samples; ten were FDR-significant before diagnosis adjustment.

## Leave-one-diagnosis-out analysis

The 44 disease samples were harmonized into five categories: lupus nephritis,
FSGS/FGGS, diabetic nephropathy, IgA nephropathy, and other kidney disease.
Each category was removed once.

- All 22 evaluated features remained positively correlated with TIMP1 in all
  five iterations.
- All nine signatures and the stringent module score remained FDR-significant
  in every iteration.
- Ten of 12 core genes remained FDR-significant in every iteration.
- FBLIM1 remained positive but was FDR-significant in one iteration.
- PFN1 remained positive but was not FDR-significant in any iteration.
- The stringent module had a median leave-one-category-out rho of 0.851
  (range 0.834-0.863).

These results reduce concern that a single major diagnosis category explains
the overall association pattern.

## Alternative-control analysis

| Control definition | Control n | Difference | Hedges' g | P value |
|---|---:|---:|---:|---:|
| Healthy living donors | 9 | 1.318 | 1.569 | 8.01e-05 |
| Living donors plus tumor-nephrectomy | 15 | 1.070 | 1.291 | 6.96e-05 |
| Tumor-nephrectomy only | 6 | 0.698 | 0.805 | 0.0664 |

The disease-control direction remained positive under all definitions. The
extended-control comparison remained significant, whereas the six-sample
tumor-nephrectomy-only comparison was imprecise and not significant. Living
donors remain the primary reference because tumor-adjacent tissue may not
represent a healthy kidney baseline.

## Robust disease-only correlations

Diagnosis category could be included as a covariate. Raw Spearman,
diagnosis-residualized partial Spearman, and standardized linear models were
calculated for the nine signatures, stringent module score, and 12 core genes.

- All 22 diagnosis-adjusted partial Spearman estimates were positive.
- Twenty-one of 22 partial Spearman associations were FDR-significant.
- Twenty-one of 22 diagnosis-adjusted linear-model coefficients were positive
  and FDR-significant.
- PFN1 remained positive but was not significant after adjustment.
- The adjusted partial rho for the stringent module was 0.821.
- Adjusted partial rho values for the nine signatures ranged from 0.438 to
  0.823 and all were FDR-significant.

Adjustment is limited to broad diagnosis categories. Disease severity,
treatment, eGFR, IFTA, fibrosis, and AKI outcomes could not be included
because systematic sample-level metadata were unavailable.

## Manuscript-ready figures

1. `v4_manuscript_figure_1_workflow_cohort_overview`
2. `v4_manuscript_figure_2_TIMP1_expression_cohorts`
3. `v4_manuscript_figure_3_signature_correlation_heatmap`
4. `v4_manuscript_figure_4_GSE180394_stringent_module_scatter`
5. `v4_manuscript_figure_5_core_gene_correlation_heatmap`

Each figure is saved as PDF and PNG. Figure legends are provided in
`manuscript_figure_legends.md`.

## Evidence boundaries

The main evidence supports TIMP1 as a candidate associated with kidney injury,
fibrosis, ECM remodeling, maladaptive repair, senescence, inflammation, and
TGF-beta-related programs across multiple transcriptomic cohorts. GSE180394
adds independent tubular-compartment validation and robustness to broad
diagnosis adjustment.

The project does not yet support:

- longitudinal AKI-to-CKD prediction;
- clinical diagnostic performance;
- kidney specificity;
- direct regulation of fibrosis or module genes;
- causal driver claims;
- patient-level single-cell replication.

## Manuscript readiness

The project is ready to begin a manuscript first draft as an integrative
transcriptomic database study. The strongest structure is: cross-cohort TIMP1
expression, pathway-level correlations, a conserved disease-only module,
external tubular validation, and sensitivity analyses. Single-cell findings
should remain exploratory localization. Experimental validation and richer
clinical outcomes would strengthen later revisions but are not prerequisites
for drafting the current database article.

## Reproducibility

Run from the project root:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_external_GSE180394_sensitivity_v4.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  tests/test_timp1_validation_v4_contracts.R
```

Key tables and reports are under `results/timp1_validation/`; PDF and PNG
figures are under `figures/timp1_validation/`. R session information is saved
in `sessionInfo_external_GSE180394_sensitivity_v4.txt`.
