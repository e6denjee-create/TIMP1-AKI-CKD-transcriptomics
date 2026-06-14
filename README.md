# TIMP1 kidney injury transcriptomic analysis

This repository contains the reproducible code, derived tables, intermediate
outputs, session information, source data, and figures for:

**Cross-Cohort Transcriptomic Analysis Associates TIMP1 with Tubular Injury-Repair, Extracellular Matrix Remodeling, and a Reproducible Disease-State Program in Human Kidney Disease**

Authors: Yanzhao Ji and Zhihong Gao.

The analysis positions TIMP1 as a candidate kidney injury and extracellular
matrix remodeling-associated gene. Results are association-based and do not
establish diagnostic utility, kidney specificity, longitudinal prediction, or
causality.

## Contents

- `scripts/`: R and Python workflows.
- `results/timp1_validation/`: statistical tables, intermediate results,
  session information, interpretations, and missing-data records.
- `figures/timp1_validation/`: PDF and PNG figure files.
- `data/processed/`: processed expression matrices small enough for GitHub.
- `data/metadata/`: sample metadata.
- `data/external/GSE180394/`: prepared external expression matrices and metadata.
- `INPUT_MANIFEST.csv`: SHA-256 checksums for every packaged input.

The complete versioned archive is deposited on Zenodo.

## Archived release

- GitHub release: https://github.com/e6denjee-create/TIMP1-AKI-CKD-transcriptomics/releases/tag/v1.2.1
- Zenodo DOI: https://doi.org/10.5281/zenodo.20689824

## Public datasets

NCBI GEO: GSE139061, GSE30718, GSE66494, and GSE180394.

## Clean-room reproduction

R 4.6.0 and Python 3.11 were used. Required R packages are checked by each
workflow. Optional enrichment is skipped and logged when `msigdbr` or its
database cache is unavailable. Python figure dependencies are listed in
`requirements-reproduction.txt`.

From a fresh checkout, run:

```powershell
python -m pip install -r requirements-reproduction.txt
python scripts/run_cleanroom_reproduction.py
```

The runner deletes packaged result and figure copies, rebuilds the core
discovery, external-validation, robustness, bootstrap, and manuscript-figure
outputs, then checks prespecified sample counts and numerical invariants. A
successful run ends with `CLEANROOM_REPRODUCTION_OK`. Package versions are
recorded in `sessionInfo*.txt` files.

## Missing data

Unavailable optional datasets or resources are recorded in
`results/timp1_validation/missing_data_log.csv`; affected units are skipped
without silently terminating the remaining workflow.
