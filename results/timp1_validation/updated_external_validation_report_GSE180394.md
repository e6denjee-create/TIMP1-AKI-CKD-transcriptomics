# Updated External Validation Report: GSE180394

## 1. Inclusion decision

GSE180394 is recommended for inclusion as an independent external bulk
validation cohort. It provides microdissected human kidney tubular expression
data and is directly aligned with the project's tubular injury-repair and ECM
remodeling context.

It should be presented as an external tubulointerstitial kidney-disease
validation cohort, not as a longitudinal AKI-to-CKD cohort. The disease group
contains multiple kidney disease etiologies.

## 2. Samples and platform

- GEO accession: GSE180394
- Platform: GPL19983, Affymetrix Human Gene 2.1 ST Array with custom Entrez
  gene annotation
- Tissue: microdissected tubules from human kidney biopsy
- Series size: 59 samples
- Primary disease group: 44 kidney-disease samples
- Primary control group: 9 healthy living-donor samples
- Excluded from the primary contrast: 6 unaffected tumor-nephrectomy samples

The disease samples include diabetic nephropathy, FSGS/FGGS, IgA nephropathy,
lupus nephritis, interstitial nephritis, hypertensive nephrosclerosis,
membranous nephropathy, thin basement membrane disease, and other kidney
diagnoses.

The GEO series matrix contained 25,582 Entrez-based probe rows. GPL19983 and
the NCBI human `gene_info` table were used to map and aggregate expression to
24,845 current gene symbols.

## 3. Feature coverage

| Feature set | Available | Requested | Coverage |
|---|---:|---:|---:|
| TIMP1 | 1 | 1 | 100% |
| Stringent TIMP1 module | 401 | 421 | 95.2% |
| Core module | 12 | 12 | 100% |
| ECM remodeling | 21 | 21 | 100% |
| Collagen formation | 13 | 13 | 100% |
| TGF-beta signaling | 10 | 11 | 90.9% |
| Inflammation | 13 | 13 | 100% |
| Tubular injury | 8 | 8 | 100% |
| Maladaptive repair | 13 | 13 | 100% |
| Cellular senescence | 10 | 10 | 100% |
| Fibrosis | 11 | 12 | 91.7% |
| Immune activation | 12 | 12 | 100% |

Coverage was sufficient for all planned external validation analyses.

## 4. TIMP1 expression replication

TIMP1 expression was higher in the disease group:

| Disease n | Control n | Disease mean | Control mean | Difference | Hedges' g | P value |
|---:|---:|---:|---:|---:|---:|---:|
| 44 | 9 | 5.069 | 3.752 | 1.318 | 1.569 | 8.01e-05 |

The direction and magnitude support replication of elevated TIMP1 expression
in this external tubular kidney-disease cohort. The estimate should not be
interpreted as diagnostic performance because no classifier, threshold,
sensitivity, specificity, or external clinical prediction model was tested.

## 5. Signature correlation replication

All nine TIMP1-signature correlations were positive and FDR-significant in
both analysis scopes.

| Signature | All-sample rho | Disease-only rho |
|---|---:|---:|
| ECM remodeling | 0.848 | 0.797 |
| Collagen formation | 0.573 | 0.453 |
| TGF-beta signaling | 0.437 | 0.590 |
| Inflammation | 0.828 | 0.819 |
| Tubular injury | 0.851 | 0.825 |
| Maladaptive repair | 0.806 | 0.800 |
| Cellular senescence | 0.455 | 0.703 |
| Fibrosis | 0.840 | 0.851 |
| Immune activation | 0.837 | 0.788 |

The disease-only results are particularly informative because the associations
persist after removing the living-donor versus disease separation. They
support TIMP1 as a candidate correlated with tubular injury-repair, ECM
remodeling, inflammatory activation, senescence, and fibrosis-related
programs.

## 6. Stringent module score replication

The stringent module score used 401 of the 421 module genes:

| Scope | Spearman rho | P value |
|---|---:|---:|
| All samples | 0.879 | Below numerical reporting precision |
| Disease only | 0.850 | Below numerical reporting precision |

The strong disease-only correlation supports replication of the
TIMP1-associated program within diseased tubular samples. This is a
correlation score and should not be described as a regulatory network.

## 7. Core module replication

All 12 core genes were detected and positively correlated with TIMP1 in both
analysis scopes.

- All samples: 12 of 12 positive and 12 of 12 FDR-significant.
- Disease only: 12 of 12 positive and 10 of 12 FDR-significant.
- Disease-only rho values ranged from 0.206 to 0.791.
- FBLIM1 remained positive (rho = 0.301, FDR = 0.0517).
- PFN1 remained positive (rho = 0.206, FDR = 0.179).

The strongest disease-only correlations included SERPING1, TMSB10, IFITM3,
PEA15, MGP, and TGM2. The pattern supports a conserved candidate module
associated with inflammatory, cytoskeletal, wound-repair, and matrix-related
states.

## 8. Main limitations

1. GSE180394 is cross-sectional and does not follow patients from AKI to CKD.
2. The disease group is etiologically heterogeneous.
3. The control group contains only nine living donors.
4. Six unaffected tumor-nephrectomy samples were excluded from the primary
   contrast; alternative control definitions were not used as primary
   evidence.
5. Systematic sample-level eGFR, fibrosis score, IFTA, and AKI outcome fields
   were not available in the reviewed GEO metadata.
6. Bulk tubular expression cannot resolve tubular subtype or immune/stromal
   composition.
7. Signature and module correlations may reflect shared injury severity,
   cellular composition, or upstream programs.
8. The analysis does not establish clinical diagnostic performance or direct
   regulation.

## 9. Manuscript recommendation

GSE180394 is recommended for the manuscript main results as an independent
external tubulointerstitial validation cohort.

The most defensible manuscript statement is:

> In an independent microdissected tubular kidney-disease cohort, TIMP1 was
> elevated in disease and remained correlated with injury-repair,
> ECM-remodeling, inflammatory, senescence, stringent-module, and core-module
> programs among disease samples.

The cohort materially strengthens the bulk evidence because it reproduces:

1. directional TIMP1 elevation;
2. all nine signature associations;
3. the stringent module score;
4. positive correlations for all 12 core genes.

It should not be used to claim longitudinal AKI-to-CKD prediction, kidney
specificity, clinical diagnostic performance, or direct causal activity.

## 10. Reproducibility

Downloaded source data and prepared inputs:

- `data/external/GSE180394/GSE180394_series_matrix.txt.gz`
- `data/external/GSE180394/GPL19983.soft.gz`
- `data/external/GSE180394/Homo_sapiens.gene_info.gz`
- `data/external/GSE180394/external_GSE180394_expression_gene_symbol.csv.gz`
- `data/external/GSE180394/external_GSE180394_metadata.csv`
- `data/external/GSE180394/external_GSE180394_metadata_all_samples.csv`

Run from the project root:

```powershell
$env:R_LIBS_USER = (Resolve-Path `
  "TIMP1_AKI_CKD_project/renv/library/R-4.6").Path

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/prepare_external_GSE180394.R

& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla `
  scripts/run_external_bulk_validation_template.R `
  --dataset GSE180394 `
  --expression data/external/GSE180394/external_GSE180394_expression_gene_symbol.csv.gz `
  --metadata data/external/GSE180394/external_GSE180394_metadata.csv `
  --gene-id-type symbol `
  --sample-col sample `
  --group-col group `
  --disease-label Disease `
  --control-label Control
```

All statistical outputs use the `external_GSE180394` prefix. Figures are saved
as both PDF and PNG, and session information is recorded.
