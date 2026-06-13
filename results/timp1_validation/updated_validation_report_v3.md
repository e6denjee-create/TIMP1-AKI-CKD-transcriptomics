# TIMP1 Validation Report v3

## 1. Executive summary

The project now contains a reproducible main-evidence layer based on bulk TIMP1
expression, nine pathway signatures, and a cross-dataset disease-only
TIMP1-associated module. It also contains a reusable workflow for screening
and analyzing an independent external bulk cohort.

Across GSE139061, GSE30718, and GSE66494, TIMP1 was directionally elevated in
disease. Two cohorts had BH-adjusted support for the expression difference.
All 27 all-sample signature correlations and 25 of 27 disease-only
correlations were FDR-significant. The disease-only module contained 1,343
relaxed and 421 stringent genes, including 12 genes positively correlated with
TIMP1 at FDR < 0.05 in all three cohorts.

No independent fourth cohort has yet been added. Complete Hallmark, GO BP, GO
CC, Reactome, and KEGG enrichment was attempted but remains unavailable.
Package installation partially succeeded, but `org.Hs.eg.db` could not be
downloaded completely and the `msigdbr` gene-set host could not be resolved.
No complete enrichment result is claimed.

Single-cell findings remain exploratory localization from one donor.

## 2. Current evidence hierarchy

### Main evidence

- Within-dataset TIMP1 disease-control expression estimates.
- All-sample and disease-only correlations with nine prespecified signatures.
- A conserved disease-only TIMP1-associated correlation module.
- Relaxed, stringent, signature-removal, and ECM/fibrosis-marker-removal
  sensitivity analyses.
- A 12-gene core module correlated with TIMP1 in all three cohorts.

### Supportive evidence

- Exploratory single-cell localization in GSE210622.
- Exploratory TIMP1-high tubular-cell-associated programs from one donor.
- Curated offline enrichment used as a reproducible sensitivity analysis.

### Evidence not yet available

- Independent replication in a fourth bulk cohort.
- Complete Hallmark/GO/Reactome/KEGG enrichment.
- Multi-donor single-cell pseudobulk validation.
- Longitudinal AKI-to-CKD patient-level validation.
- Clinical diagnostic validation, tissue-specificity testing, or experimental
  evidence of direct regulation.

## 3. External GEO screening framework

`candidate_geo_screening_table.csv` provides a standardized schema for
recording:

- GEO accession and disease context;
- kidney tissue or compartment;
- disease and control definitions;
- biological sample size;
- processed matrix, raw data, and metadata availability;
- donor/sample identifiers;
- fibrosis, eGFR, IFTA, or AKI outcome variables;
- TIMP1 detectability and signature coverage;
- inclusion recommendation, proposed use, and exclusion reason.

The table is intentionally unpopulated because no additional GEO accession was
verified in the current local project. Candidate datasets should be entered
only after checking the processed files and sample-level metadata.

Priority contexts remain AKI, longitudinal AKI-to-CKD transition, CKD/renal
fibrosis, diabetic kidney disease, kidney transplant injury, and
tubulointerstitial injury.

## 4. Reusable external bulk validation

`scripts/run_external_bulk_validation_template.R` accepts a gene-by-sample
expression matrix and sample-level metadata. It supports gene symbols,
Ensembl IDs, Entrez IDs, or probe IDs with an annotation table.

For each new cohort, the script can produce:

1. TIMP1 disease-versus-control group difference, Hedges' g, P value, and
   BH-adjusted P value.
2. Scores and coverage for all nine signatures.
3. All-sample and disease-only Spearman correlations.
4. A stringent TIMP1-module score and its TIMP1 correlations.
5. Expression and TIMP1 correlations for the 12 core genes.
6. CSV tables, PDF/PNG figures, session information, and a short markdown
   interpretation.

The script analyzes each cohort independently and does not pool absolute
expression values across platforms. Missing inputs, inadequate group sizes,
unmappable gene IDs, or absent TIMP1 are recorded in `missing_data_log.csv`.

The template was smoke-tested end to end using the already included GSE66494
files under the explicit prefix
`external_TEMPLATE_SMOKE_GSE66494_NOT_INDEPENDENT_`. These outputs verify the
workflow only and must not be counted as an independent fourth cohort.

Example:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_external_bulk_validation_template.R `
  --dataset GSEXXXX `
  --expression data/external/GSEXXXX_expression.csv.gz `
  --metadata data/external/GSEXXXX_metadata.csv `
  --gene-id-type symbol `
  --sample-col sample `
  --group-col group `
  --disease-label Disease `
  --control-label Control
```

## 5. Comprehensive enrichment attempt

`scripts/run_comprehensive_enrichment.R` was run for six module definitions:

- relaxed;
- stringent;
- relaxed without signature genes;
- stringent without signature genes;
- relaxed without ECM/fibrosis markers;
- stringent without ECM/fibrosis markers.

The requested enrichment resources were Hallmark, GO BP, GO CC, Reactome, and
KEGG. Installation and retrieval outcomes were:

| Component | Final status |
|---|---|
| msigdbr | Available, version 26.1.0 |
| AnnotationDbi | Available, version 1.74.0 |
| ggplot2 | Available, version 4.0.3 |
| clusterProfiler | Not loadable after installation attempt |
| ReactomePA | Not loadable after installation attempt |
| org.Hs.eg.db | Download incomplete after extended timeout |
| Hallmark collection | Unavailable from the msigdbr host |
| GO BP collection | Unavailable from the msigdbr host |
| GO CC collection | Unavailable from the msigdbr host |
| Reactome collection | Unavailable from the msigdbr host |
| KEGG collection | Unavailable from the msigdbr host |

The `org.Hs.eg.db` source archive was approximately 99 MB; the transfer stopped
after approximately 25 MB despite a 900-second timeout. The `msigdbr`
collection requests continued to fail because the Zenodo host could not be
resolved.

`comprehensive_enrichment_v3_results.csv` therefore contains no claimed
enrichment rows. `comprehensive_enrichment_v3_status.csv`, the paired status
figures, and `missing_data_log.csv` document the incomplete execution.

The earlier curated offline sensitivity enrichment remains the only available
pathway-level module result and must continue to be labeled as a targeted
subset rather than a complete database-wide analysis.

## 6. Bulk and module findings retained from v2

TIMP1 was directionally elevated in all three bulk cohorts. GSE30718 and
GSE66494 had BH-adjusted support, while GSE139061 did not.

All 27 all-sample TIMP1-signature correlations were positive and
FDR-significant. Twenty-five of 27 disease-only correlations were
FDR-significant. This supports TIMP1 as a candidate correlated with
injury-repair, ECM remodeling, collagen formation, TGF-beta-related,
inflammatory, tubular injury, maladaptive repair, senescence, fibrosis, and
immune activation programs.

The 421-gene stringent module and its sensitivity variants support a
cross-dataset program correlated with TIMP1. These associations do not
establish direct regulation or causal direction.

## 7. Core module

The 12 genes with positive FDR support in all three cohorts are:

`TMSB10`, `IFITM3`, `RAB31`, `MGP`, `PEA15`, `TUBA1B`, `SERPING1`, `FBLIM1`,
`TAGLN2`, `TGM2`, `TUBA1A`, and `PFN1`.

Their median rho values range from approximately 0.56 to 0.77. The core module
is associated with cytoskeletal organization, membrane trafficking,
inflammatory response, complement regulation, adhesion, wound repair, and
matrix-related biology. The annotations in `core_module_genes_summary.csv`
provide candidate context and do not imply direct TIMP1 regulation.

## 8. Manuscript results skeleton

`manuscript_results_skeleton.md` provides a results-section framework covering:

1. directional TIMP1 elevation;
2. injury-repair and ECM-remodeling signature correlations;
3. the conserved disease-only module;
4. module sensitivity analyses;
5. the 12-gene core module;
6. exploratory single-cell localization;
7. limitations.

The skeleton uses association-oriented language and separates observed results
from unresolved clinical or mechanistic questions.

## 9. Single-cell status

GSE210622 contains one donor and therefore supports exploratory localization
only. GSE267242 remains incomplete locally. Multi-donor TIMP1 expression,
positive-cell fraction, cell-type pseudobulk, and donor-level signature
comparisons remain pending.

Cell-level significance from the current one-donor object should not be used
as patient-level evidence.

## 10. Limitations and next evidence priorities

1. Screen and acquire an independent fourth human bulk cohort with processed
   expression and sample-level metadata.
2. Prioritize longitudinal AKI outcomes, biopsy fibrosis/IFTA, or
   tubulointerstitial compartments.
3. Run the external template and assess TIMP1, nine signatures, the stringent
   module score, and the 12 core genes.
4. Complete public gene-set enrichment after obtaining a stable
   `org.Hs.eg.db` installation or a locally cached MSigDB GMT resource.
5. Acquire multi-donor single-cell data for patient-level pseudobulk
   validation.
6. Treat all current module and single-cell findings as associated,
   correlated, candidate, or exploratory evidence.

## 11. New v3 files

Results and documentation:

- `candidate_geo_screening_table.csv`
- `manuscript_results_skeleton.md`
- `comprehensive_enrichment_v3_results.csv`
- `comprehensive_enrichment_v3_status.csv`
- `sessionInfo_comprehensive_enrichment_v3.txt`
- `updated_validation_report_v3.md`

Scripts:

- `scripts/run_external_bulk_validation_template.R`
- `scripts/run_comprehensive_enrichment.R`

Figures, each in PDF and PNG:

- `TIMP1_comprehensive_enrichment_v3`
- `TIMP1_comprehensive_enrichment_v3_status`

## 12. Reproducibility

Run from the project root:

```powershell
$env:R_LIBS_USER = (Resolve-Path `
  "TIMP1_AKI_CKD_project/renv/library/R-4.6").Path

$env:TIMP1_INSTALL_ENRICHMENT_PACKAGES = "true"

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_comprehensive_enrichment.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  tests/test_timp1_validation_v3_contracts.R
```

Package installation and database downloads require network access. The
script records failures and continues to produce explicit status outputs
rather than representing unavailable analyses as completed.
