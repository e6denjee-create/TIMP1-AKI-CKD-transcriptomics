#!/usr/bin/env Rscript

# Audit local GSE210622 and GSE267242 resources before any patient-level
# single-cell validation. This script does not substitute cell-level tests for
# missing donor-level replication.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260612)

root <- validation_root()
project <- timp1_project_dir(root)

file_specs <- data.frame(
  dataset = c(
    rep("GSE210622", 4),
    rep("GSE267242", 4)
  ),
  resource = rep(
    c("expression_matrix", "barcodes", "features", "cell_metadata"), 2
  ),
  path = c(
    file.path(
      project, "data", "raw", "GSE210622", "GSM6433706", "matrix.mtx.gz"
    ),
    file.path(
      project, "data", "raw", "GSE210622", "GSM6433706", "barcodes.tsv.gz"
    ),
    file.path(
      project, "data", "raw", "GSE210622", "GSM6433706", "features.tsv.gz"
    ),
    file.path(
      project, "data", "metadata",
      "GSE210622_GSM6433706_cell_metadata.csv"
    ),
    file.path(project, "data", "raw", "GSE267242", "matrix.mtx.gz"),
    file.path(project, "data", "raw", "GSE267242", "barcodes.tsv.gz"),
    file.path(project, "data", "raw", "GSE267242", "features.tsv.gz"),
    file.path(project, "data", "metadata", "GSE267242_cell_metadata.csv")
  ),
  stringsAsFactors = FALSE
)
file_specs$exists <- file.exists(file_specs$path)
file_specs$bytes <- ifelse(
  file_specs$exists, file.info(file_specs$path)$size, NA_real_
)

rds_path <- file.path(
  project, "data", "processed", "GSE210622_GSM6433706_seurat.rds"
)
rds_row <- data.frame(
  dataset = "GSE210622",
  resource = "Seurat_RDS",
  path = rds_path,
  exists = file.exists(rds_path),
  bytes = ifelse(file.exists(rds_path), file.info(rds_path)$size, NA_real_),
  stringsAsFactors = FALSE
)
inventory <- rbind(file_specs, rds_row)
write_validation_csv(inventory, "singlecell_data_inventory.csv")

status_rows <- list()

audit_metadata <- function(dataset, metadata_path) {
  if (!file.exists(metadata_path)) {
    append_missing_data(
      "singlecell_validation", dataset, metadata_path,
      "Cell metadata file is absent.",
      "Patient-level analysis was not attempted."
    )
    return(data.frame(
      dataset = dataset,
      metadata_available = FALSE,
      cells = NA_integer_,
      donor_column = NA_character_,
      donor_count = NA_integer_,
      sample_column = NA_character_,
      sample_count = NA_integer_,
      cell_type_column = NA_character_,
      cell_type_count = NA_integer_,
      patient_level_ready = FALSE,
      reason = "Cell metadata is absent.",
      stringsAsFactors = FALSE
    ))
  }
  metadata <- tryCatch(
    read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(error) {
      append_missing_data(
        "singlecell_validation", dataset, metadata_path,
        conditionMessage(error),
        "Patient-level analysis was not attempted."
      )
      NULL
    }
  )
  if (is.null(metadata)) {
    return(data.frame(
      dataset = dataset,
      metadata_available = FALSE,
      cells = NA_integer_,
      donor_column = NA_character_,
      donor_count = NA_integer_,
      sample_column = NA_character_,
      sample_count = NA_integer_,
      cell_type_column = NA_character_,
      cell_type_count = NA_integer_,
      patient_level_ready = FALSE,
      reason = "Cell metadata could not be read.",
      stringsAsFactors = FALSE
    ))
  }
  donor_candidates <- intersect(
    c("donor", "donor_id", "patient", "patient_id", "individual"),
    colnames(metadata)
  )
  sample_candidates <- intersect(
    c("sample", "sample_id", "orig.ident"), colnames(metadata)
  )
  type_candidates <- intersect(
    c("cell_type", "celltype", "annotation", "cell_type_annotation"),
    colnames(metadata)
  )
  donor_column <- if (length(donor_candidates)) donor_candidates[[1]] else NA
  sample_column <- if (length(sample_candidates)) sample_candidates[[1]] else NA
  type_column <- if (length(type_candidates)) type_candidates[[1]] else NA
  donor_count <- if (!is.na(donor_column)) {
    length(unique(metadata[[donor_column]]))
  } else {
    NA_integer_
  }
  sample_count <- if (!is.na(sample_column)) {
    length(unique(metadata[[sample_column]]))
  } else {
    NA_integer_
  }
  type_count <- if (!is.na(type_column)) {
    length(unique(metadata[[type_column]]))
  } else {
    NA_integer_
  }
  effective_replicates <- if (!is.na(donor_count)) donor_count else sample_count
  ready <- !is.na(effective_replicates) && effective_replicates >= 2 &&
    !is.na(type_column)
  reasons <- character()
  if (is.na(donor_column)) {
    reasons <- c(reasons, "No explicit donor ID column")
  }
  if (is.na(sample_column)) {
    reasons <- c(reasons, "No sample ID column")
  }
  if (is.na(type_column)) {
    reasons <- c(reasons, "No cell-type annotation column")
  }
  if (!is.na(effective_replicates) && effective_replicates < 2) {
    reasons <- c(reasons, "Only one donor/sample is represented")
  }
  data.frame(
    dataset = dataset,
    metadata_available = TRUE,
    cells = nrow(metadata),
    donor_column = donor_column,
    donor_count = donor_count,
    sample_column = sample_column,
    sample_count = sample_count,
    cell_type_column = type_column,
    cell_type_count = type_count,
    patient_level_ready = ready,
    reason = paste(reasons, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

status_rows[["GSE210622"]] <- audit_metadata(
  "GSE210622",
  file.path(
    project, "data", "metadata",
    "GSE210622_GSM6433706_cell_metadata.csv"
  )
)
status_rows[["GSE267242"]] <- audit_metadata(
  "GSE267242",
  file.path(project, "data", "metadata", "GSE267242_cell_metadata.csv")
)
status <- do.call(rbind, status_rows)

for (dataset in unique(inventory$dataset)) {
  absent <- inventory[
    inventory$dataset == dataset & !inventory$exists, , drop = FALSE
  ]
  for (index in seq_len(nrow(absent))) {
    append_missing_data(
      "singlecell_validation", dataset, absent$resource[[index]],
      paste("Required local file is absent:", absent$path[[index]]),
      "Patient-level analysis was not attempted."
    )
  }
}
if (!status$patient_level_ready[status$dataset == "GSE210622"]) {
  append_missing_data(
    "singlecell_validation", "GSE210622", "independent donors",
    "Only GSM6433706 is represented in the current object.",
    "Retained existing cell-level results as exploratory evidence only."
  )
}
write_validation_csv(status, "singlecell_validation_status.csv")

needed_lines <- c(
  "# Single-cell Validation Data Needed",
  "",
  "## Current status",
  "",
  "Patient-level or sample-level TIMP1 validation cannot currently be completed.",
  "The local GSE210622 object contains one AKI donor/sample (GSM6433706).",
  "GSE267242 has barcode and feature files, but its expression matrix and",
  "cell-level metadata are absent locally.",
  "",
  "The existing single-donor cell-level results remain exploratory and must not",
  "be presented as clinical or patient-level evidence.",
  "",
  "## GSE210622",
  "",
  "Available:",
  "",
  "- `data/raw/GSE210622/GSM6433706/matrix.mtx.gz`",
  "- matching barcode and feature files",
  "- `data/processed/GSE210622_GSM6433706_seurat.rds`",
  "- cell metadata with a sample column and computational cell-type labels",
  "",
  "Still needed:",
  "",
  "- expression matrices for additional GSE210622 donors/samples",
  "- donor ID and clinical group for every cell",
  "- consistent cell-type annotations across donors",
  "- preferably author-provided annotations or a reproducible reference mapping",
  "",
  "## GSE267242",
  "",
  "Available:",
  "",
  "- `data/raw/GSE267242/barcodes.tsv.gz`",
  "- `data/raw/GSE267242/features.tsv.gz`",
  "",
  "Still needed:",
  "",
  "- the matching expression matrix, such as `matrix.mtx.gz`, H5, RDS, or H5AD",
  "- cell-level metadata linking barcodes to donor/sample IDs",
  "- disease/control or time-point annotation",
  "- cell-type annotation, including tubular, podocyte, fibroblast, and",
  "  myeloid/macrophage compartments",
  "",
  "## Required validation workflow after data completion",
  "",
  "1. Construct or import one consistently processed object containing all donors.",
  "2. Verify unique donor/sample IDs and harmonized cell-type labels.",
  "3. Aggregate counts by donor and cell type to create pseudobulk profiles.",
  "4. Report donor-level TIMP1 mean expression and TIMP1-positive cell fraction.",
  "5. Analyze tubular, podocyte, fibroblast, and myeloid/macrophage compartments.",
  "6. Define TIMP1-high and TIMP1-low tubular cells within donor and cell type",
  "   using the 75th and 25th percentiles.",
  "7. Aggregate signature scores to the donor level before statistical testing.",
  "8. Use UMAP only for visualization; use donors as the independent units.",
  "",
  "## Re-run",
  "",
  "After adding the missing matrices and metadata, run:",
  "",
  "```powershell",
  "& \"C:\\Program Files\\R\\R-4.6.0\\bin\\Rscript.exe\" --vanilla \\",
  "  scripts/run_singlecell_validation_audit.R",
  "```",
  "",
  "Proceed to a patient-level validation script only when at least two",
  "independent donors and usable cell-type annotations are available."
)
writeLines(
  needed_lines,
  file.path(validation_result_dir(root), "singlecell_validation_needed.md")
)

write_validation_session_info("sessionInfo_singlecell_audit.txt")
message("Single-cell validation audit completed; patient-level data are incomplete.")
