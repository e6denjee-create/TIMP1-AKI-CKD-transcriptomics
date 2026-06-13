#!/usr/bin/env Rscript

# Prepare GSE180394 for independent TIMP1 external bulk validation.
# Main contrast: heterogeneous kidney disease versus living-donor controls.
# Tumor-nephrectomy unaffected tissue is excluded from the primary contrast.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260613)

required_packages <- c("Biobase", "GEOquery")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

input_dir <- file.path("data", "external", "GSE180394")
eset_path <- file.path(input_dir, "GSE180394_eset_1.rds")
gpl_path <- file.path(input_dir, "GPL19983.soft.gz")
gene_info_path <- file.path(input_dir, "Homo_sapiens.gene_info.gz")
required_inputs <- c(eset_path, gpl_path, gene_info_path)
if (!all(file.exists(required_inputs))) {
  missing <- required_inputs[!file.exists(required_inputs)]
  append_missing_data(
    "prepare_external_bulk", "GSE180394", paste(missing, collapse = "; "),
    "Required GEO or NCBI annotation input is absent.",
    "Preparation stopped without claiming external validation."
  )
  stop("Missing GSE180394 preparation inputs.")
}

eset <- readRDS(eset_path)
expression_probe <- Biobase::exprs(eset)
phenotype <- Biobase::pData(eset)
platform <- GEOquery::getGEO(filename = gpl_path)
platform_table <- GEOquery::Table(platform)

gene_info <- read.delim(
  gzfile(gene_info_path), check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
tax_id_column <- if ("#tax_id" %in% colnames(gene_info)) {
  "#tax_id"
} else {
  "tax_id"
}
if (!tax_id_column %in% colnames(gene_info)) {
  stop("NCBI gene_info does not contain a tax_id column.")
}
gene_info <- gene_info[gene_info[[tax_id_column]] == 9606, , drop = FALSE]
gene_info <- gene_info[!duplicated(gene_info$GeneID), , drop = FALSE]

probe_entrez <- platform_table$ENTREZ_GENE_ID[
  match(rownames(expression_probe), platform_table$ID)
]
gene_symbol <- gene_info$Symbol[
  match(as.character(probe_entrez), as.character(gene_info$GeneID))
]
keep_genes <- !is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != "-"
expression_probe <- expression_probe[keep_genes, , drop = FALSE]
gene_symbol <- toupper(gene_symbol[keep_genes])

# Average probes mapping to the same current human gene symbol.
expression_symbol <- rowsum(
  expression_probe, group = gene_symbol, reorder = FALSE
)
probe_counts <- as.numeric(table(factor(
  gene_symbol, levels = unique(gene_symbol)
)))
expression_symbol <- expression_symbol / probe_counts
if (!nrow(expression_symbol) || !"TIMP1" %in% rownames(expression_symbol)) {
  stop("Entrez-to-symbol mapping failed or TIMP1 is absent after mapping.")
}

sample_group <- as.character(phenotype[["sample group:ch1"]])
analysis_group <- ifelse(
  sample_group == "Living donor", "Control",
  ifelse(
    grepl("Tumor Nephrectomy", sample_group, ignore.case = TRUE),
    "Excluded_tumor_nephrectomy_unaffected", "Disease"
  )
)
metadata_all <- data.frame(
  sample = as.character(phenotype$geo_accession),
  group = analysis_group,
  original_diagnosis = sample_group,
  title = as.character(phenotype$title),
  tissue = sub(
    "^tissue:\\s*", "", as.character(phenotype[["tissue:ch1"]]),
    ignore.case = TRUE
  ),
  platform = as.character(phenotype$platform_id),
  organism = as.character(phenotype$organism_ch1),
  primary_analysis_included = analysis_group %in% c("Disease", "Control"),
  exclusion_reason = ifelse(
    analysis_group == "Excluded_tumor_nephrectomy_unaffected",
    paste(
      "Unaffected tumor-nephrectomy tissue excluded to keep living donors",
      "as the primary control definition."
    ),
    ""
  ),
  stringsAsFactors = FALSE
)

analysis_metadata <- metadata_all[
  metadata_all$primary_analysis_included, , drop = FALSE
]
analysis_expression <- expression_symbol[
  , analysis_metadata$sample, drop = FALSE
]

expression_path <- file.path(
  input_dir, "external_GSE180394_expression_gene_symbol.csv.gz"
)
expression_all_path <- file.path(
  input_dir, "external_GSE180394_expression_all_samples_gene_symbol.csv.gz"
)
metadata_path <- file.path(
  input_dir, "external_GSE180394_metadata.csv"
)
metadata_all_path <- file.path(
  input_dir, "external_GSE180394_metadata_all_samples.csv"
)
write.csv(
  analysis_expression, gzfile(expression_path),
  row.names = TRUE, na = ""
)
write.csv(
  expression_symbol, gzfile(expression_all_path),
  row.names = TRUE, na = ""
)
write.csv(analysis_metadata, metadata_path, row.names = FALSE, na = "")
write.csv(metadata_all, metadata_all_path, row.names = FALSE, na = "")

signature_table <- read.csv(
  file.path(validation_result_dir(), "signature_gene_sets_used.csv"),
  stringsAsFactors = FALSE
)
stringent <- read.csv(
  file.path(validation_result_dir(), "stringent_TIMP1_correlated_module.csv"),
  stringsAsFactors = FALSE
)
core <- read.csv(
  file.path(validation_result_dir(), "core_module_genes_summary.csv"),
  stringsAsFactors = FALSE
)
available <- rownames(analysis_expression)

coverage_rows <- list(
  data.frame(
    feature_set = "TIMP1",
    requested_genes = 1L,
    available_genes = as.integer("TIMP1" %in% available),
    coverage_fraction = as.numeric("TIMP1" %in% available),
    available_gene_symbols = if ("TIMP1" %in% available) "TIMP1" else "",
    stringsAsFactors = FALSE
  ),
  data.frame(
    feature_set = "stringent_TIMP1_module",
    requested_genes = length(unique(stringent$gene)),
    available_genes = length(intersect(stringent$gene, available)),
    coverage_fraction = length(intersect(stringent$gene, available)) /
      length(unique(stringent$gene)),
    available_gene_symbols = paste(
      intersect(stringent$gene, available), collapse = ";"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    feature_set = "core_12_gene_module",
    requested_genes = length(unique(core$gene)),
    available_genes = length(intersect(core$gene, available)),
    coverage_fraction = length(intersect(core$gene, available)) /
      length(unique(core$gene)),
    available_gene_symbols = paste(intersect(core$gene, available), collapse = ";"),
    stringsAsFactors = FALSE
  )
)

for (signature in unique(signature_table$signature)) {
  genes <- unique(signature_table$gene[
    signature_table$signature == signature
  ])
  present <- intersect(genes, available)
  coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
    feature_set = paste0("signature_", signature),
    requested_genes = length(genes),
    available_genes = length(present),
    coverage_fraction = length(present) / length(genes),
    available_gene_symbols = paste(present, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
coverage <- do.call(rbind, coverage_rows)
write_validation_csv(
  coverage, "external_GSE180394_input_feature_coverage.csv"
)

preparation_summary <- data.frame(
  GEO_accession = "GSE180394",
  platform = unique(analysis_metadata$platform),
  tissue_compartment = "Microdissected tubules from human kidney biopsy",
  total_series_samples = nrow(metadata_all),
  disease_samples = sum(analysis_metadata$group == "Disease"),
  living_donor_controls = sum(analysis_metadata$group == "Control"),
  excluded_tumor_nephrectomy_unaffected = sum(
    metadata_all$group == "Excluded_tumor_nephrectomy_unaffected"
  ),
  probe_rows = nrow(expression_probe),
  mapped_unique_gene_symbols = nrow(analysis_expression),
  expression_min = min(analysis_expression, na.rm = TRUE),
  expression_max = max(analysis_expression, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_validation_csv(
  preparation_summary, "external_GSE180394_preparation_summary.csv"
)
write_validation_session_info(
  "external_GSE180394_preparation_sessionInfo.txt"
)

message(
  "Prepared GSE180394: ", nrow(analysis_expression), " genes, ",
  sum(analysis_metadata$group == "Disease"), " disease samples, ",
  sum(analysis_metadata$group == "Control"), " living-donor controls."
)
