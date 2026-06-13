# Updated TIMP1 Validation Report

## 1. Executive summary

The expanded validation supports TIMP1 as a candidate gene associated with
kidney injury repair, extracellular matrix (ECM) remodeling, fibrosis,
inflammatory activation, cellular senescence, and tubular maladaptive repair in
the AKI-to-CKD context.

Across GSE139061, GSE30718, and GSE66494, TIMP1 was directionally elevated in
disease samples. The elevation remained BH-significant in GSE30718 and
GSE66494, but not in GSE139061. TIMP1 was positively correlated with all nine
prespecified pathway signature scores in all three cohorts. In disease-only
samples, 25 of 27 TIMP1-signature associations remained FDR-significant.

A disease-only genome-wide analysis identified 1,343 genes with replicated
positive TIMP1 correlations and nominal support in at least two datasets. This
correlation module included injury, cytoskeletal, immune, wound-response, and
matrix-associated genes. Enrichment in a curated offline GO/Reactome/KEGG
subset supported ECM organization, response to wounding, collagen formation,
inflammatory signaling, focal adhesion, senescence, and TGF-beta-related
programs.

Patient-level single-cell validation cannot yet be completed. The available
GSE210622 object contains one donor, and the local GSE267242 files lack the
expression matrix and cell metadata. Existing single-donor cell-level P values
remain exploratory evidence only.

TIMP1 should not be described as kidney-specific, as a validated diagnostic
biomarker, or as a proven causal driver.

## 2. What is supported

- TIMP1 is reproducibly directionally elevated across the three evaluated
  kidney injury/fibrosis bulk cohorts.
- TIMP1 elevation is statistically supported in GSE30718 and GSE66494.
- TIMP1 is positively associated with ECM remodeling, collagen formation,
  TGF-beta signaling, inflammation, tubular injury, maladaptive repair,
  cellular senescence, fibrosis, and immune activation scores.
- These associations persist in disease-only analyses and therefore are not
  explained solely by disease-versus-control separation.
- A replicated TIMP1-associated correlation module is present across the bulk
  datasets.
- The module is compatible with ECM, wound-response, inflammatory, adhesion,
  senescence, and TGF-beta-related biological programs.
- TIMP1 can be advanced as an AKI-to-CKD-related injury-repair and
  ECM-remodeling candidate gene for further experimental and patient-level
  validation.

## 3. What is not supported

- TIMP1 is not demonstrated to be kidney-specific.
- TIMP1 is not established as a diagnostic biomarker for AKI, CKD, fibrosis, or
  AKI-to-CKD transition.
- No validated sensitivity, specificity, clinical threshold, or independent
  diagnostic model has been established.
- TIMP1 is not proven to be a causal driver of fibrosis or maladaptive repair.
- The correlation module is not a causal regulatory network.
- Single-donor cell-level P values do not establish patient-level or clinical
  relevance.
- The available cohorts do not follow the same patients longitudinally from
  AKI to CKD.

## 4. Dataset-level results

| Dataset | Context | Disease | Control | Group difference | Hedges' g | P value | BH-adjusted P value |
|---|---|---:|---:|---:|---:|---:|---:|
| GSE139061 | AKI | 39 | 9 | 1.086 | 0.415 | 0.3781 | 0.3781 |
| GSE30718 | AKI | 28 | 11 | 1.062 | 1.293 | 0.000109 | 0.000326 |
| GSE66494 | CKD/fibrosis | 53 | 8 | 0.417 | 0.449 | 0.00603 | 0.00905 |

All expression contrasts were performed within datasets. Absolute normalized
expression values were not combined across platforms.

## 5. Bulk validation

Nine signature scores were calculated separately within each dataset by
standardizing each available gene across samples and averaging the resulting
gene-level Z-scores. TIMP1 itself was excluded from the signatures to avoid
mechanical self-correlation.

All 27 all-sample TIMP1-signature correlations were positive and
FDR-significant. In disease-only samples:

- GSE139061 had 9 positive associations, of which 7 were FDR-significant.
- GSE30718 had 9 positive and FDR-significant associations.
- GSE66494 had 9 positive and FDR-significant associations.

Disease-only correlations were particularly consistent for ECM remodeling,
fibrosis, maladaptive repair, tubular injury, inflammation, and immune
activation. The two GSE139061 associations that did not reach FDR < 0.05 were
cellular senescence and collagen formation, although both retained positive
rho values.

These findings support pathway-level coherence around TIMP1. They do not show
that TIMP1 directly controls the genes in these signatures.

## 6. Single-cell validation status

### GSE210622

Available locally:

- complete MTX, barcode, and feature files for GSM6433706;
- a Seurat RDS object;
- metadata for 23,095 cells/nuclei;
- one sample ID and computational cell-type annotations.

Limitation: only one donor/sample is represented. Donor-level replication,
between-patient tests, and patient-level pseudobulk comparisons are therefore
not possible.

### GSE267242

Available locally:

- barcode file;
- feature file.

Missing locally:

- matching expression matrix;
- cell-level metadata;
- donor/sample IDs;
- disease or time-point annotations;
- cell-type annotations.

The required files and proposed pseudobulk workflow are documented in
`singlecell_validation_needed.md`. No patient-level single-cell result has been
claimed.

## 7. TIMP1-high tubular cell program

The existing one-donor analysis associated TIMP1-positive/high tubular cells
with ECM organization, tubular maladaptive repair, matrix degradation, wound
healing, collagen formation, TGF-beta signaling, and cellular senescence.

This is compatible with the statement that TIMP1-high tubular cells may
represent a maladaptive repair-associated state. However, the tubular TIMP1
median and 75th percentile were both zero, so the contrast largely reflected
TIMP1-detected versus TIMP1-undetected cells. Dropout, sequencing depth, and
tubular subtype composition may contribute. The result remains exploratory
until replicated using donor-level aggregation.

## 8. TIMP1-associated module

Genome-wide Spearman correlations were calculated using disease-only samples
in each bulk cohort. A stable gene required:

1. positive correlation direction in at least two datasets;
2. nominal positive association in at least two datasets;
3. positive median rho across datasets.

This produced 1,343 stable TIMP1-correlated genes. Among them:

- 421 had positive FDR support in at least two datasets;
- 12 had positive FDR support in all three datasets.

Highly replicated genes included TMSB10, IFITM3, RAB31, MGP, PEA15, TUBA1B,
SERPING1, FBLIM1, TAGLN2, TGM2, PFN1, CD44, MMP7, NRP1, VCAN, and VIM.
These genes are consistent with injury response, cytoskeletal remodeling,
immune activation, epithelial stress, and ECM-related processes.

The module was enriched in an explicitly labeled curated offline
GO/Reactome/KEGG subset for:

- response to wounding;
- extracellular matrix organization;
- ECM-receptor interaction;
- inflammatory response;
- collagen fibril organization and collagen formation;
- focal adhesion;
- cellular senescence;
- NF-kappa B signaling;
- TGF-beta receptor signaling;
- extracellular matrix degradation.

The full MSigDB database collections could not be downloaded because the local
environment could not resolve the Zenodo host. Therefore, this enrichment is a
targeted offline subset analysis, not a comprehensive unbiased database-wide
screen.

## 9. Limitations

1. The bulk datasets are cross-sectional, clinically heterogeneous, and
   platform-specific.
2. GSE139061 did not show a statistically significant TIMP1 group difference.
3. Bulk signature and module correlations may reflect cell composition,
   immune infiltration, disease severity, or unmeasured covariates.
4. Signature gene sets overlap biologically and are not statistically
   independent.
5. Stable-module thresholds are pragmatic and should be tested in additional
   cohorts.
6. The offline enrichment subset was prespecified around relevant biological
   themes and is not equivalent to full GO/Reactome/KEGG testing.
7. Single-cell localization is based on one donor and computational cell labels.
8. Cell-level P values are vulnerable to pseudoreplication.
9. GSE267242 is incomplete locally.
10. No perturbation experiment demonstrates direct TIMP1-mediated regulation.

## 10. Recommended manuscript claim boundaries

### Appropriate

- "TIMP1 is reproducibly elevated in kidney injury/fibrosis datasets."
- "TIMP1 is associated with ECM remodeling, maladaptive repair, senescence,
  inflammation, immune activation, and TGF-beta-related programs."
- "A replicated disease-only TIMP1-associated correlation module links TIMP1
  with injury-repair and matrix-remodeling programs."
- "TIMP1-high tubular cells may represent a maladaptive repair-associated
  state."
- "TIMP1 is a candidate AKI-to-CKD-related injury-repair and ECM-remodeling
  gene."

### Avoid

- "TIMP1 is kidney-specific."
- "TIMP1 is a validated diagnostic biomarker."
- "TIMP1 is a proven causal driver."
- "TIMP1 directly regulates fibrosis."
- "The TIMP1 correlation network is a regulatory network."
- "Single-donor cell-level significance proves clinical relevance."

Use association language such as "associated with," "correlated with,"
"candidate," "exploratory," and "may represent."

## 11. Figure list

New validation figures under `figures/timp1_validation/`, each in PDF and PNG:

1. `TIMP1_bulk_validation_expression`
2. `TIMP1_signature_correlation_heatmap`
3. `TIMP1_centered_correlation_network`
4. `TIMP1_module_enrichment_dotplot`
5. `TIMP1_correlated_gene_overlap`

Existing MVP figures remain unchanged under
`TIMP1_AKI_CKD_project/results/figures/`.

## 12. Table list

New validation tables under `results/timp1_validation/`:

1. `signature_gene_sets_used.csv`
2. `TIMP1_signature_gene_coverage.csv`
3. `TIMP1_signature_scores.csv`
4. `TIMP1_signature_correlations.csv`
5. `TIMP1_bulk_validation_source_data.csv`
6. `TIMP1_bulk_validation_statistics.csv`
7. `TIMP1_genomewide_disease_only_correlations.csv`
8. `stable_TIMP1_correlated_genes.csv`
9. `network_edge_list.csv`
10. `enrichment_TIMP1_correlated_module.csv`
11. `singlecell_data_inventory.csv`
12. `singlecell_validation_status.csv`
13. `missing_data_log.csv`

Additional documentation:

- `singlecell_validation_needed.md`
- `sessionInfo_bulk_validation.txt`
- `sessionInfo_module_network.txt`
- `sessionInfo_singlecell_audit.txt`

## 13. Reproducibility instructions

Run from the repository root:

```powershell
$env:R_LIBS_USER = (Resolve-Path `
  "TIMP1_AKI_CKD_project/renv/library/R-4.6").Path

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_timp1_bulk_validation.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_timp1_module_network.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_singlecell_validation_audit.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  tests/test_timp1_validation_contracts.R
```

All scripts use the fixed seed `20260612`. New results are written only to
`results/timp1_validation/` and new figures only to
`figures/timp1_validation/`. Existing MVP results are not overwritten.
