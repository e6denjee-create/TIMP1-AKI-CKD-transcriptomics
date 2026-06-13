# Single-cell Validation Data Needed

## Current status

Patient-level or sample-level TIMP1 validation cannot currently be completed.
The local GSE210622 object contains one AKI donor/sample (GSM6433706).
GSE267242 has barcode and feature files, but its expression matrix and
cell-level metadata are absent locally.

The existing single-donor cell-level results remain exploratory and must not
be presented as clinical or patient-level evidence.

## GSE210622

Available:

- `data/raw/GSE210622/GSM6433706/matrix.mtx.gz`
- matching barcode and feature files
- `data/processed/GSE210622_GSM6433706_seurat.rds`
- cell metadata with a sample column and computational cell-type labels

Still needed:

- expression matrices for additional GSE210622 donors/samples
- donor ID and clinical group for every cell
- consistent cell-type annotations across donors
- preferably author-provided annotations or a reproducible reference mapping

## GSE267242

Available:

- `data/raw/GSE267242/barcodes.tsv.gz`
- `data/raw/GSE267242/features.tsv.gz`

Still needed:

- the matching expression matrix, such as `matrix.mtx.gz`, H5, RDS, or H5AD
- cell-level metadata linking barcodes to donor/sample IDs
- disease/control or time-point annotation
- cell-type annotation, including tubular, podocyte, fibroblast, and
  myeloid/macrophage compartments

## Required validation workflow after data completion

1. Construct or import one consistently processed object containing all donors.
2. Verify unique donor/sample IDs and harmonized cell-type labels.
3. Aggregate counts by donor and cell type to create pseudobulk profiles.
4. Report donor-level TIMP1 mean expression and TIMP1-positive cell fraction.
5. Analyze tubular, podocyte, fibroblast, and myeloid/macrophage compartments.
6. Define TIMP1-high and TIMP1-low tubular cells within donor and cell type
   using the 75th and 25th percentiles.
7. Aggregate signature scores to the donor level before statistical testing.
8. Use UMAP only for visualization; use donors as the independent units.

## Re-run

After adding the missing matrices and metadata, run:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" --vanilla \
  scripts/run_singlecell_validation_audit.R
```

Proceed to a patient-level validation script only when at least two
independent donors and usable cell-type annotations are available.
