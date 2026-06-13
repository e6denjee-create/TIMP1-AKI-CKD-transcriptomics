# TIMP1 kidney injury transcriptomic analysis

This repository contains the reproducible code, derived tables, intermediate
outputs, session information, source data, and figures for:

**Cross-Cohort Transcriptomic Analysis Associates TIMP1 with Tubular Injury-Repair, Extracellular Matrix Remodeling, and a Conserved Disease-State Program in Human Kidney Disease**

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

The complete archive, including the >100 MB exploratory single-cell RDS object,
is intended for Zenodo because GitHub rejects individual files above 100 MB.

## Public datasets

NCBI GEO: GSE139061, GSE30718, GSE66494, GSE180394, GSE210622, and GSE267242.

## Reproduction

R 4.6.0 was used for the principal analyses. Run the scripts in numeric or
workflow order from the repository root. Package versions are recorded in the
`sessionInfo*.txt` files under `results/timp1_validation/`.

## Missing data

Unavailable optional datasets or resources are recorded in
`results/timp1_validation/missing_data_log.csv`; affected units are skipped
without silently terminating the remaining workflow.
