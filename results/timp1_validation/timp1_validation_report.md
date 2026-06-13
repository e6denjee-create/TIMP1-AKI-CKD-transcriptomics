# TIMP1 Validation Report

## 1. Executive summary

The current evidence supports continued investigation of TIMP1 as a candidate
gene associated with kidney injury repair and extracellular matrix (ECM)
remodeling in the AKI-to-CKD context. TIMP1 expression was higher in disease
samples in all three evaluated bulk cohorts, although only two comparisons
remained significant after multiple-testing correction. Bulk correlation and
single-cell pathway analyses connect TIMP1 with injury, fibrosis, maladaptive
repair, senescence, and TGF-beta-related programs.

These results do not establish TIMP1 as a kidney-specific diagnostic biomarker
or a causal driver of fibrosis. The current single-cell analysis contains one
AKI donor, so cell-level P values are exploratory evidence and cannot support
patient-level or clinical inference.

The present report summarizes completed MVP results. Dedicated bulk signature
validation, multi-donor single-cell validation, and a cross-dataset
TIMP1-associated correlation module have not yet been completed.

## 2. What is supported

- TIMP1 was directionally elevated in disease samples in GSE139061, GSE30718,
  and GSE66494.
- The elevation was BH-significant in GSE30718 and GSE66494.
- TIMP1 showed positive disease-only correlations with multiple injury,
  fibrosis/ECM, inflammatory, and cell-cycle markers.
- In one exploratory AKI single-nucleus sample, TIMP1 expression varied among
  computationally annotated renal cell compartments.
- TIMP1-positive tubular cells were associated with ECM organization, matrix
  degradation, collagen formation, wound healing, tubular maladaptive repair,
  TGF-beta signaling, and cellular senescence gene sets.
- Taken together, TIMP1 can be described as a candidate AKI-to-CKD-related
  injury-repair and ECM-remodeling gene that warrants further validation.

## 3. What is not supported

- TIMP1 is not demonstrated to be kidney-specific.
- TIMP1 is not established as a diagnostic biomarker for AKI, CKD, fibrosis, or
  AKI-to-CKD transition.
- The current analyses do not establish diagnostic sensitivity, specificity,
  or clinically validated discrimination.
- TIMP1 is not proven to be a causal driver of fibrosis or maladaptive repair.
- Correlation does not demonstrate direct regulation of ECM or inflammatory
  genes by TIMP1.
- Single-donor cell-level P values do not prove clinical relevance.
- The available datasets do not constitute a longitudinal study of the same
  patients progressing from AKI to CKD.

## 4. Dataset-level results

| Dataset | Context | Disease | Control | TIMP1 difference | P value | BH-adjusted P value | Interpretation |
|---|---|---:|---:|---:|---:|---:|---|
| GSE139061 | AKI | 39 | 9 | 1.086 | 0.3781 | 0.3781 | Higher mean, not statistically significant |
| GSE30718 | AKI | 28 | 11 | 1.062 | 0.000109 | 0.000326 | Significant elevation |
| GSE66494 | CKD/fibrosis | 53 | 8 | 0.417 | 0.00603 | 0.00905 | Significant elevation |
| GSE210622, GSM6433706 | AKI single-nucleus RNA-seq | 1 donor | NA | Descriptive only | NA | NA | Exploratory localization |
| GSE267242 | Candidate single-cell validation dataset | NA | NA | Not analyzed | NA | NA | Expression matrix and metadata are incomplete locally |

The bulk expression differences are within-dataset contrasts. Absolute
normalized expression values were not pooled across platforms.

## 5. Bulk validation

The completed MVP includes disease-versus-control TIMP1 comparisons and
prespecified marker correlations in three independent bulk datasets.
TIMP1 was elevated in the same direction in all three cohorts, with statistical
support in two cohorts.

In disease-only sensitivity analyses, 29 positive TIMP1-marker correlations
were FDR-significant. The recurrent associations included injury markers
HAVCR1 and LCN2; ECM/fibrosis genes COL1A1, COL3A1, FN1, and ACTA2;
inflammatory genes including CCL2; and cell-state genes such as CDKN1A and
CCNB1. These findings support a shared injury/remodeling program but may still
reflect cell composition, disease severity, or other unmeasured covariates.

A dedicated validation workflow calculating ECM remodeling, collagen
formation, TGF-beta signaling, inflammation, tubular injury, maladaptive
repair, senescence, immune activation, and fibrosis signature scores has not
yet produced results. Therefore, no additional signature-score claims are made
in this report.

## 6. Single-cell validation status

Local data inspection identified:

- A complete GSE210622 Seurat object for GSM6433706.
- Raw MTX, barcode, and feature files for the same GSE210622 sample.
- Only barcode and feature files for GSE267242; the corresponding expression
  matrix and donor/cell-type metadata are not available locally.

GSE210622 currently represents one AKI donor. It permits descriptive
localization but not valid patient-level replication or between-patient
inference. The highest mean TIMP1 expression in the computational annotation
was observed in podocytes, followed by endothelial cells. These labels were
derived from marker scores and were not independently validated author
annotations.

Multi-donor pseudobulk analysis, donor-level TIMP1-positive fractions, and
patient-level comparisons remain required. UMAP and cell-level plots should be
treated as visualization rather than primary statistical evidence.

## 7. TIMP1-high tubular cell program

The current tubular-cell analysis compared cells classified as TIMP1-high and
TIMP1-low within one donor. Positive FDR-significant enrichment was observed
for:

- ECM organization
- Tubular maladaptive repair
- Matrix degradation
- Wound healing
- Collagen formation
- TGF-beta signaling
- Cellular senescence

This result is consistent with the hypothesis that TIMP1-high tubular cells may
represent a maladaptive repair-associated state.

Important qualification: the median and 75th percentile of tubular-cell TIMP1
expression were both zero. The comparison therefore approximates
TIMP1-detected versus TIMP1-undetected cells rather than a robust continuous
high-versus-low expression contrast. Dropout, sequencing depth, and tubular
subtype composition may contribute to the observed differences. TIMP1 was also
included in several curated gene sets, creating potential circularity.

## 8. TIMP1-associated module

A genome-wide, cross-dataset TIMP1-associated correlation module has not yet
been generated. The existing evidence is limited to correlations with a
prespecified marker panel.

The planned module analysis should:

1. Calculate disease-only genome-wide Spearman correlations separately in each
   bulk dataset.
2. Retain positively correlated genes with concordant direction in at least two
   datasets.
3. Perform GO, Reactome, and KEGG enrichment on the stable gene set.
4. Report the result as a TIMP1-associated correlation module or associated
   program, not as a causal regulatory network.

Until those outputs exist, no stable module size, enrichment result, or network
topology should be reported.

## 9. Limitations

1. The AKI and CKD cohorts are cross-sectional and platform-heterogeneous, not
   longitudinal AKI-to-CKD patient trajectories.
2. GSE139061 showed a non-significant TIMP1 group difference despite the same
   direction of effect.
3. Bulk correlations can be driven by tissue composition, immune infiltration,
   disease severity, and shared group effects.
4. Available metadata have not been used for comprehensive covariate
   adjustment.
5. The single-cell evidence is based on one donor and is vulnerable to
   pseudoreplication.
6. Cell identities were assigned computationally and require reference-based or
   author-annotation validation.
7. The tubular high/low definition is affected by zero inflation and dropout.
8. Including TIMP1 in TIMP1-associated pathway gene sets may inflate
   enrichment.
9. GSE267242 is incomplete locally and cannot currently provide independent
   single-cell replication.
10. Dedicated bulk signature and genome-wide module validation outputs are not
    yet available.

## 10. Recommended manuscript claim boundaries

### Appropriate claims

- "TIMP1 is reproducibly elevated in kidney injury/fibrosis datasets."
- "TIMP1 is associated with ECM remodeling, maladaptive repair, senescence,
  and TGF-beta-related programs."
- "TIMP1-high tubular cells may represent a maladaptive repair-associated
  state."
- "TIMP1 is a candidate AKI-to-CKD-related injury-repair and ECM-remodeling
  gene that warrants experimental and multi-cohort validation."

### Claims to avoid

- "TIMP1 is kidney-specific."
- "TIMP1 is a validated diagnostic biomarker."
- "TIMP1 is a proven causal driver of AKI-to-CKD transition."
- "TIMP1 directly regulates fibrosis."
- "Single-donor cell-level significance demonstrates clinical relevance."

Use "associated with," "correlated with," "candidate," "exploratory," and
"may represent" where appropriate. Reserve causal language for perturbation or
mechanistic experimental evidence.

## 11. Figure list

Existing MVP figures are stored in
`TIMP1_AKI_CKD_project/results/figures/`, each as PDF and PNG:

1. `TIMP1_bulk_expression_boxplot`
2. `TIMP1_marker_correlation_heatmap`
3. `single_cell_UMAP_cell_types`
4. `single_cell_TIMP1_FeaturePlot`
5. `single_cell_TIMP1_DotPlot`
6. `single_cell_TIMP1_ViolinPlot`
7. `TIMP1_high_tubular_enrichment`

No dedicated figures currently exist under `figures/timp1_validation/`.

## 12. Table list

Key existing MVP tables are stored in
`TIMP1_AKI_CKD_project/results/tables/`:

1. `bulk_dataset_manifest.csv`
2. `bulk_marker_gene_coverage.csv`
3. `TIMP1_bulk_expression_source_data.csv`
4. `TIMP1_bulk_expression_statistics.csv`
5. `TIMP1_marker_correlations.csv`
6. `TIMP1_marker_correlation_heatmap_source.csv`
7. `single_cell_qc_summary.csv`
8. `single_cell_cluster_marker_scores.csv`
9. `single_cell_TIMP1_cell_type_summary.csv`
10. `single_cell_TIMP1_source_data.csv`
11. `TIMP1_high_low_threshold.csv`
12. `TIMP1_high_vs_low_tubular_DE.csv`
13. `TIMP1_high_tubular_pathway_enrichment.csv`
14. `package_versions.csv`
15. `required_output_manifest.csv`

No dedicated bulk-signature, patient-level single-cell, or
TIMP1-associated-module validation tables currently exist under
`results/timp1_validation/`.

## 13. Reproducibility instructions

Run the existing MVP from the TIMP1 project directory:

```powershell
Set-Location TIMP1_AKI_CKD_project
$env:R_LIBS_USER = (Resolve-Path "renv/library/R-4.6").Path
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla scripts/run_mvp.R
```

The workflow uses a fixed random seed defined in
`config/analysis_config.R`. Standardized expression matrices are stored under
`data/processed/`, metadata under `data/metadata/`, and software versions in
`results/tables/package_versions.csv` and `sessionInfo.txt`.

Future validation scripts should be run from the repository root and must write
new tables to `results/timp1_validation/` and new PDF/PNG figures to
`figures/timp1_validation/`. Existing MVP outputs must not be overwritten.
Missing datasets or metadata should be recorded in
`results/timp1_validation/missing_data_log.csv` while allowing unaffected
analysis units to continue.
