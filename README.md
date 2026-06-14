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

The complete versioned archive is deposited on Zenodo.

## Archived release

- GitHub release: https://github.com/e6denjee-create/TIMP1-AKI-CKD-transcriptomics/releases/tag/v1.2.0
- Zenodo DOI: https://doi.org/10.5281/zenodo.20688976

## Public datasets

NCBI GEO: GSE139061, GSE30718, GSE66494, and GSE180394.

## Reproduction

R 4.6.0 was used for the principal analyses. Run the scripts in numeric or
workflow order from the repository root. Package versions are recorded in the
`sessionInfo*.txt` files under `results/timp1_validation/`.

## Missing data

Unavailable optional datasets or resources are recorded in
`results/timp1_validation/missing_data_log.csv`; affected units are skipped
without silently terminating the remaining workflow.
