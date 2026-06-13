#!/usr/bin/env Rscript

# Reusable external bulk-cohort validation template for TIMP1.
# Expected expression layout: genes in rows, samples in columns, first column
# containing gene IDs. Expected metadata: one row per biological sample.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260612)

parse_args <- function(args) {
  defaults <- list(
    dataset = "external_cohort",
    expression = "",
    metadata = "",
    annotation = "",
    gene_id_type = "auto",
    sample_col = "sample",
    group_col = "group",
    disease_label = "Disease",
    control_label = "Control",
    help = FALSE
  )
  index <- 1
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("--help", "-h")) {
      defaults$help <- TRUE
      index <- index + 1
      next
    }
    if (!startsWith(argument, "--") || index == length(args)) {
      stop("Arguments must use --name value syntax.")
    }
    key <- gsub("-", "_", substring(argument, 3), fixed = TRUE)
    if (!key %in% names(defaults)) stop("Unknown argument: ", argument)
    defaults[[key]] <- args[[index + 1]]
    index <- index + 2
  }
  defaults
}

print_help <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/run_external_bulk_validation_template.R",
      "--dataset GSEXXXX --expression expression.csv[.gz]",
      "--metadata metadata.csv --gene-id-type auto|symbol|ensembl|entrez|probe",
      "[--annotation gene_annotation.csv]",
      "[--sample-col sample --group-col group]",
      "[--disease-label Disease --control-label Control]\n"
    )
  )
}

read_delimited <- function(path, row_names = FALSE) {
  separator <- if (grepl("\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$", path, TRUE)) {
    "\t"
  } else {
    ","
  }
  connection <- if (grepl("\\.gz$", path, TRUE)) gzfile(path) else path
  read.table(
    connection, header = TRUE, sep = separator, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "\"", comment.char = "",
    row.names = if (row_names) 1 else NULL
  )
}

stop_unit <- function(dataset, item, reason) {
  append_missing_data(
    "external_bulk_validation", dataset, item, reason,
    "Recorded missing or invalid input; external cohort analysis was skipped."
  )
  message("External cohort skipped: ", reason)
  quit(save = "no", status = 0)
}

collapse_duplicate_genes <- function(expression, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[keep, , drop = FALSE]
  symbols <- toupper(symbols[keep])
  rowsum(expression, group = symbols, reorder = FALSE) /
    as.numeric(table(factor(symbols, levels = unique(symbols))))
}

map_gene_ids <- function(expression, gene_ids, gene_id_type, annotation_path,
                         dataset) {
  gene_ids <- as.character(gene_ids)
  inferred_type <- gene_id_type
  if (gene_id_type == "auto") {
    inferred_type <- if (mean(grepl("^ENSG", gene_ids)) > 0.5) {
      "ensembl"
    } else if (mean(grepl("^[0-9]+$", gene_ids)) > 0.8) {
      "entrez"
    } else {
      "symbol"
    }
  }

  if (inferred_type == "symbol") {
    symbols <- gene_ids
  } else if (nzchar(annotation_path)) {
    if (!file.exists(annotation_path)) {
      stop_unit(dataset, annotation_path, "Gene annotation file is absent.")
    }
    annotation <- read_delimited(annotation_path)
    required <- c("gene_id", "gene_symbol")
    if (!all(required %in% colnames(annotation))) {
      stop_unit(
        dataset, annotation_path,
        "Annotation must contain gene_id and gene_symbol columns."
      )
    }
    symbols <- annotation$gene_symbol[
      match(gene_ids, as.character(annotation$gene_id))
    ]
  } else if (
    inferred_type %in% c("ensembl", "entrez") &&
      requireNamespace("AnnotationDbi", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)
  ) {
    clean_ids <- if (inferred_type == "ensembl") {
      sub("\\..*$", "", gene_ids)
    } else {
      gene_ids
    }
    keytype <- if (inferred_type == "ensembl") "ENSEMBL" else "ENTREZID"
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db, keys = unique(clean_ids),
      columns = "SYMBOL", keytype = keytype
    )
    mapping <- mapping[!duplicated(mapping[[keytype]]), , drop = FALSE]
    symbols <- mapping$SYMBOL[match(clean_ids, mapping[[keytype]])]
  } else {
    stop_unit(
      dataset, "gene ID conversion",
      paste(
        "Cannot map", inferred_type,
        "IDs without annotation or AnnotationDbi/org.Hs.eg.db."
      )
    )
  }

  mapped <- collapse_duplicate_genes(expression, symbols)
  if (!nrow(mapped)) {
    stop_unit(dataset, "gene ID conversion", "No genes mapped to symbols.")
  }
  mapped
}

z_score_rows <- function(matrix) {
  scaled <- t(scale(t(matrix)))
  scaled[!is.finite(scaled)] <- NA_real_
  scaled
}

score_gene_set <- function(expression, genes) {
  available <- intersect(toupper(genes), rownames(expression))
  if (!length(available)) return(rep(NA_real_, ncol(expression)))
  colMeans(z_score_rows(expression[available, , drop = FALSE]), na.rm = TRUE)
}

spearman_row <- function(x, y, dataset, scope, feature) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4 || length(unique(x[keep])) < 2 ||
      length(unique(y[keep])) < 2) {
    return(data.frame(
      dataset = dataset, scope = scope, feature = feature, rho = NA_real_,
      p_value = NA_real_, n = sum(keep), stringsAsFactors = FALSE
    ))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
  data.frame(
    dataset = dataset, scope = scope, feature = feature,
    rho = unname(test$estimate), p_value = test$p.value, n = sum(keep),
    stringsAsFactors = FALSE
  )
}

hedges_g <- function(disease, control) {
  n1 <- length(disease)
  n0 <- length(control)
  pooled_sd <- sqrt(
    ((n1 - 1) * var(disease) + (n0 - 1) * var(control)) /
      (n1 + n0 - 2)
  )
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  correction * (mean(disease) - mean(control)) / pooled_sd
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help) || !nzchar(args$expression) || !nzchar(args$metadata)) {
  print_help()
  if (!isTRUE(args$help)) {
    append_missing_data(
      "external_bulk_validation", args$dataset, "expression/metadata paths",
      "Required command-line inputs were not supplied.",
      "Printed template usage; no cohort result was claimed."
    )
  }
  quit(save = "no", status = 0)
}

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop_unit(
    args$dataset, paste(missing_packages, collapse = "; "),
    "Required plotting package is not installed."
  )
}
if (!file.exists(args$expression)) {
  stop_unit(args$dataset, args$expression, "Expression matrix is absent.")
}
if (!file.exists(args$metadata)) {
  stop_unit(args$dataset, args$metadata, "Metadata file is absent.")
}

raw_expression <- read_delimited(args$expression, row_names = TRUE)
expression <- as.matrix(raw_expression)
storage.mode(expression) <- "numeric"
expression <- map_gene_ids(
  expression, rownames(expression), tolower(args$gene_id_type),
  args$annotation, args$dataset
)
metadata <- read_delimited(args$metadata)

required_metadata <- c(args$sample_col, args$group_col)
if (!all(required_metadata %in% colnames(metadata))) {
  stop_unit(
    args$dataset, paste(setdiff(required_metadata, colnames(metadata)),
                        collapse = "; "),
    "Required metadata columns are absent."
  )
}
metadata$sample <- as.character(metadata[[args$sample_col]])
metadata$analysis_group <- as.character(metadata[[args$group_col]])
samples <- intersect(metadata$sample, colnames(expression))
if (length(samples) < 8) {
  stop_unit(
    args$dataset, "matched samples",
    "Fewer than eight metadata-expression matched samples are available."
  )
}
metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]
expression <- expression[, samples, drop = FALSE]
keep_groups <- metadata$analysis_group %in%
  c(args$disease_label, args$control_label)
metadata <- metadata[keep_groups, , drop = FALSE]
expression <- expression[, metadata$sample, drop = FALSE]

if (!"TIMP1" %in% rownames(expression)) {
  stop_unit(args$dataset, "TIMP1", "TIMP1 is not detectable after gene mapping.")
}
group_counts <- table(factor(
  metadata$analysis_group,
  levels = c(args$disease_label, args$control_label)
))
if (any(group_counts < 4)) {
  stop_unit(
    args$dataset, "group sample size",
    "At least one required group contains fewer than four biological samples."
  )
}

safe_dataset <- gsub("[^A-Za-z0-9_.-]", "_", args$dataset)
prefix <- paste0("external_", safe_dataset, "_")
reference_files <- c(
  "signature_gene_sets_used.csv",
  "stringent_TIMP1_correlated_module.csv",
  "core_module_genes_summary.csv"
)
missing_references <- reference_files[
  !file.exists(file.path(validation_result_dir(), reference_files))
]
if (length(missing_references)) {
  stop_unit(
    args$dataset, paste(missing_references, collapse = "; "),
    "Required project reference tables are absent."
  )
}
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

timp1 <- as.numeric(expression["TIMP1", ])
disease_index <- metadata$analysis_group == args$disease_label
control_index <- metadata$analysis_group == args$control_label
group_test <- suppressWarnings(wilcox.test(
  timp1[disease_index], timp1[control_index], exact = FALSE
))
group_statistics <- data.frame(
  dataset = args$dataset,
  contrast = paste(args$disease_label, "vs", args$control_label),
  n_disease = sum(disease_index), n_control = sum(control_index),
  disease_mean = mean(timp1[disease_index]),
  control_mean = mean(timp1[control_index]),
  group_difference = mean(timp1[disease_index]) - mean(timp1[control_index]),
  hedges_g = hedges_g(timp1[disease_index], timp1[control_index]),
  p_value = group_test$p.value,
  bh_adjusted_p_value = group_test$p.value,
  stringsAsFactors = FALSE
)
write_validation_csv(
  group_statistics, paste0(prefix, "TIMP1_group_statistics.csv")
)

signature_names <- unique(signature_table$signature)
signature_scores <- data.frame(
  dataset = args$dataset, sample = metadata$sample,
  group = metadata$analysis_group, stringsAsFactors = FALSE
)
coverage <- list()
for (signature in signature_names) {
  genes <- setdiff(
    signature_table$gene[signature_table$signature == signature], "TIMP1"
  )
  available <- intersect(toupper(genes), rownames(expression))
  signature_scores[[signature]] <- score_gene_set(expression, genes)
  coverage[[signature]] <- data.frame(
    dataset = args$dataset, signature = signature,
    requested_genes = length(unique(genes)),
    available_genes = length(available),
    coverage_fraction = length(available) / length(unique(genes)),
    available_gene_symbols = paste(available, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
coverage <- do.call(rbind, coverage)
write_validation_csv(
  coverage, paste0(prefix, "signature_gene_coverage.csv")
)
write_validation_csv(
  signature_scores, paste0(prefix, "signature_scores.csv")
)

correlation_rows <- list()
row_index <- 1
for (scope in c("all_samples", "disease_only")) {
  selected <- if (scope == "all_samples") {
    rep(TRUE, nrow(metadata))
  } else {
    disease_index
  }
  for (signature in signature_names) {
    correlation_rows[[row_index]] <- spearman_row(
      timp1[selected], signature_scores[[signature]][selected],
      args$dataset, scope, signature
    )
    row_index <- row_index + 1
  }
}
signature_correlations <- do.call(rbind, correlation_rows)
signature_correlations$bh_adjusted_p_value <- ave(
  signature_correlations$p_value, signature_correlations$scope,
  FUN = function(x) p.adjust(x, "BH")
)
write_validation_csv(
  signature_correlations, paste0(prefix, "signature_correlations.csv")
)

stringent_genes <- intersect(toupper(stringent$gene), rownames(expression))
module_score <- score_gene_set(expression, stringent_genes)
module_scores <- data.frame(
  dataset = args$dataset, sample = metadata$sample,
  group = metadata$analysis_group, stringent_TIMP1_module_score = module_score,
  genes_available = length(stringent_genes), stringsAsFactors = FALSE
)
write_validation_csv(
  module_scores, paste0(prefix, "stringent_module_scores.csv")
)
module_correlations <- do.call(rbind, lapply(
  c("all_samples", "disease_only"),
  function(scope) {
    selected <- if (scope == "all_samples") rep(TRUE, nrow(metadata)) else {
      disease_index
    }
    spearman_row(
      timp1[selected], module_score[selected], args$dataset, scope,
      "stringent_TIMP1_module_score"
    )
  }
))
module_correlations$bh_adjusted_p_value <- p.adjust(
  module_correlations$p_value, "BH"
)
write_validation_csv(
  module_correlations, paste0(prefix, "stringent_module_correlations.csv")
)

core_genes <- intersect(toupper(core$gene), rownames(expression))
core_expression <- if (length(core_genes)) {
  data.frame(
    dataset = args$dataset,
    gene = rep(core_genes, ncol(expression)),
    sample = rep(colnames(expression), each = length(core_genes)),
    group = rep(metadata$analysis_group, each = length(core_genes)),
    expression = as.vector(expression[core_genes, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    dataset = character(), gene = character(), sample = character(),
    group = character(), expression = numeric(), stringsAsFactors = FALSE
  )
}
write_validation_csv(
  core_expression, paste0(prefix, "core_module_expression.csv")
)
core_correlations <- if (length(core_genes)) {
  output <- do.call(rbind, lapply(core_genes, function(gene) {
    do.call(rbind, lapply(c("all_samples", "disease_only"), function(scope) {
      selected <- if (scope == "all_samples") rep(TRUE, nrow(metadata)) else {
        disease_index
      }
      spearman_row(
        timp1[selected], as.numeric(expression[gene, selected]),
        args$dataset, scope, gene
      )
    }))
  }))
  output$bh_adjusted_p_value <- ave(
    output$p_value, output$scope,
    FUN = function(x) p.adjust(x, "BH")
  )
  output
} else {
  append_missing_data(
    "external_bulk_validation", args$dataset, "12 core module genes",
    "No core module gene was detected after gene mapping.",
    "Wrote empty core-gene tables and continued."
  )
  data.frame(
    dataset = character(), scope = character(), feature = character(),
    rho = numeric(), p_value = numeric(), n = integer(),
    bh_adjusted_p_value = numeric(), stringsAsFactors = FALSE
  )
}
write_validation_csv(
  core_correlations, paste0(prefix, "core_module_correlations.csv")
)

plot_data <- data.frame(
  group = metadata$analysis_group, TIMP1_expression = timp1
)
expression_plot <- ggplot2::ggplot(
  plot_data, ggplot2::aes(group, TIMP1_expression, fill = group)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.55) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.75) +
  ggplot2::labs(
    title = paste("TIMP1 expression in", args$dataset),
    x = NULL, y = "Within-cohort normalized expression"
  ) +
  validation_theme() +
  ggplot2::theme(legend.position = "none")
save_validation_plot(
  expression_plot, paste0(prefix, "TIMP1_expression"), width = 6.5, height = 5
)

correlation_plot <- ggplot2::ggplot(
  signature_correlations,
  ggplot2::aes(scope, feature, fill = rho)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = paste("TIMP1-signature correlations in", args$dataset),
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  validation_theme()
save_validation_plot(
  correlation_plot, paste0(prefix, "signature_correlation_heatmap"),
  width = 7.5, height = 6
)

core_plot <- if (nrow(core_correlations)) {
  ggplot2::ggplot(
    core_correlations,
    ggplot2::aes(scope, feature, fill = rho)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", rho)), size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4C78A8", mid = "white", high = "#C44E52",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = paste("Core-module correlations with TIMP1 in", args$dataset),
      x = NULL, y = NULL, fill = "Spearman rho"
    ) +
    validation_theme()
} else {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = "No core module genes were detectable."
    ) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(
      title = paste("Core-module correlations in", args$dataset)
    ) +
    validation_theme()
}
save_validation_plot(
  core_plot, paste0(prefix, "core_module_correlation_heatmap"),
  width = 7, height = 7
)

interpretation <- c(
  paste0("# External bulk validation: ", args$dataset),
  "",
  paste(
    "This cohort was analyzed independently using within-cohort normalized",
    "expression values."
  ),
  paste(
    "TIMP1 was evaluated as a candidate associated with injury-repair and",
    "ECM-remodeling programs; correlation does not establish direct regulation."
  ),
  paste(
    "Available stringent-module genes:", length(stringent_genes), "of",
    nrow(stringent), "."
  ),
  paste("Available core genes:", length(core_genes), "of", nrow(core), "."),
  "",
  "Review the CSV effect estimates, FDR values, gene coverage, and figures",
  "before deciding whether this cohort provides independent support."
)
writeLines(
  interpretation,
  file.path(validation_result_dir(), paste0(prefix, "interpretation.md"))
)
write_validation_session_info(paste0(prefix, "sessionInfo.txt"))
message("External bulk validation completed for ", args$dataset, ".")
