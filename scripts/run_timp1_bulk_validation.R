#!/usr/bin/env Rscript

# Validate TIMP1 expression and pathway-level associations in independent bulk
# kidney cohorts without pooling absolute expression across platforms.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260612)

required_packages <- c("ggplot2", "pheatmap")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

signature_gene_sets <- list(
  ECM_remodeling = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
    "COL6A1", "COL6A2", "FN1", "SPARC", "POSTN", "VCAN", "LUM",
    "DCN", "MMP2", "MMP7", "MMP9", "MMP14", "TIMP2", "LOX", "LOXL2"
  ),
  Collagen_formation = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
    "COL5A2", "COL6A1", "P4HA1", "P4HA2", "PLOD1", "PLOD2", "SERPINH1"
  ),
  TGF_beta_signaling = c(
    "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
    "SMAD7", "CTGF", "SERPINE1", "ACTA2"
  ),
  Inflammation = c(
    "IL1B", "IL6", "TNF", "NFKB1", "NFKBIA", "CCL2", "CCL5",
    "CXCL8", "CXCL10", "ICAM1", "VCAM1", "STAT1", "STAT3"
  ),
  Tubular_injury = c(
    "HAVCR1", "LCN2", "KRT8", "KRT18", "KRT19", "VIM", "SPP1", "MMP7"
  ),
  Maladaptive_repair = c(
    "HAVCR1", "LCN2", "VCAM1", "VIM", "KRT8", "KRT18", "KRT19",
    "SOX9", "MMP7", "SPP1", "CCL2", "IL6", "CDKN1A"
  ),
  Cellular_senescence = c(
    "CDKN1A", "CDKN2A", "TP53", "SERPINE1", "GDF15", "IL6", "CCL2",
    "CXCL8", "MMP3", "MMP9"
  ),
  Fibrosis = c(
    "COL1A1", "COL1A2", "COL3A1", "FN1", "ACTA2", "CTGF", "TGFB1",
    "POSTN", "DCN", "LUM", "SPARC", "SERPINE1"
  ),
  Immune_activation = c(
    "PTPRC", "TYROBP", "LST1", "FCER1G", "CTSS", "CD68", "HLA-DRA",
    "HLA-DPA1", "HLA-DPB1", "C1QA", "C1QB", "C1QC"
  )
)

gene_set_table <- do.call(rbind, lapply(names(signature_gene_sets), function(x) {
  data.frame(
    signature = x,
    gene = signature_gene_sets[[x]],
    stringsAsFactors = FALSE
  )
}))
write_validation_csv(gene_set_table, "signature_gene_sets_used.csv")

zscore_signature <- function(expression, genes) {
  available <- intersect(genes, rownames(expression))
  if (length(available) < 2) return(NULL)
  matrix <- expression[available, , drop = FALSE]
  scaled <- t(scale(t(matrix)))
  scaled[!is.finite(scaled)] <- NA_real_
  list(
    score = colMeans(scaled, na.rm = TRUE),
    available = available
  )
}

datasets <- c("GSE139061", "GSE30718", "GSE66494")
score_rows <- list()
coverage_rows <- list()
correlation_rows <- list()
expression_rows <- list()
expression_stats <- list()

for (dataset_name in datasets) {
  dataset <- read_validation_bulk(dataset_name)
  if (is.null(dataset)) next
  if (!"TIMP1" %in% rownames(dataset$expression)) {
    append_missing_data(
      "bulk_signature", dataset_name, "TIMP1",
      "TIMP1 is absent from the expression matrix."
    )
    next
  }

  metadata <- dataset$metadata
  timp1 <- as.numeric(dataset$expression["TIMP1", metadata$sample])
  disease_groups <- setdiff(unique(metadata$group), "Control")
  if (length(disease_groups) != 1) {
    append_missing_data(
      "bulk_expression", dataset_name, "group",
      "Expected exactly one disease group and one Control group."
    )
  } else {
    disease <- disease_groups[[1]]
    disease_values <- timp1[metadata$group == disease]
    control_values <- timp1[metadata$group == "Control"]
    if (length(disease_values) >= 2 && length(control_values) >= 2) {
      test <- stats::t.test(disease_values, control_values)
      pooled_sd <- sqrt(
        ((length(disease_values) - 1) * stats::var(disease_values) +
           (length(control_values) - 1) * stats::var(control_values)) /
          (length(disease_values) + length(control_values) - 2)
      )
      cohen_d <- (mean(disease_values) - mean(control_values)) / pooled_sd
      correction <- 1 - 3 / (
        4 * (length(disease_values) + length(control_values)) - 9
      )
      expression_stats[[dataset_name]] <- data.frame(
        dataset = dataset_name,
        contrast = paste0(disease, "_vs_Control"),
        n_disease = length(disease_values),
        n_control = length(control_values),
        disease_mean = mean(disease_values),
        control_mean = mean(control_values),
        group_difference = mean(disease_values) - mean(control_values),
        hedges_g = cohen_d * correction,
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  expression_rows[[dataset_name]] <- data.frame(
    dataset = dataset_name,
    sample = metadata$sample,
    group = metadata$group,
    TIMP1_expression = timp1,
    stringsAsFactors = FALSE
  )

  dataset_scores <- data.frame(
    dataset = dataset_name,
    sample = metadata$sample,
    group = metadata$group,
    TIMP1_expression = timp1,
    stringsAsFactors = FALSE
  )
  for (signature in names(signature_gene_sets)) {
    scored <- zscore_signature(
      dataset$expression, signature_gene_sets[[signature]]
    )
    available <- intersect(
      signature_gene_sets[[signature]], rownames(dataset$expression)
    )
    missing <- setdiff(signature_gene_sets[[signature]], available)
    coverage_rows[[paste(dataset_name, signature)]] <- data.frame(
      dataset = dataset_name,
      signature = signature,
      genes_requested = length(signature_gene_sets[[signature]]),
      genes_available = length(available),
      available_genes = paste(available, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
    if (is.null(scored)) {
      append_missing_data(
        "bulk_signature", dataset_name, signature,
        "Fewer than two signature genes were available."
      )
      next
    }
    dataset_scores[[signature]] <- unname(scored$score[metadata$sample])
  }
  score_rows[[dataset_name]] <- dataset_scores

  for (scope in c("all_samples", "disease_only")) {
    keep <- if (scope == "all_samples") {
      rep(TRUE, nrow(dataset_scores))
    } else {
      dataset_scores$group != "Control"
    }
    scoped <- dataset_scores[keep, , drop = FALSE]
    signature_columns <- intersect(
      names(signature_gene_sets), colnames(scoped)
    )
    for (signature in signature_columns) {
      valid <- is.finite(scoped$TIMP1_expression) &
        is.finite(scoped[[signature]])
      if (sum(valid) < 4 ||
          stats::sd(scoped$TIMP1_expression[valid]) == 0 ||
          stats::sd(scoped[[signature]][valid]) == 0) {
        append_missing_data(
          "bulk_signature_correlation", dataset_name,
          paste(scope, signature, sep = ":"),
          "Insufficient non-constant observations for Spearman correlation."
        )
        next
      }
      test <- suppressWarnings(stats::cor.test(
        scoped$TIMP1_expression[valid],
        scoped[[signature]][valid],
        method = "spearman", exact = FALSE
      ))
      correlation_rows[[paste(dataset_name, scope, signature)]] <- data.frame(
        dataset = dataset_name,
        scope = scope,
        signature = signature,
        rho = unname(test$estimate),
        p_value = test$p.value,
        n = sum(valid),
        stringsAsFactors = FALSE
      )
    }
  }
}

scores <- do.call(rbind, score_rows)
coverage <- do.call(rbind, coverage_rows)
correlations <- do.call(rbind, correlation_rows)
statistics <- do.call(rbind, expression_stats)
expression_source <- do.call(rbind, expression_rows)

if (!is.null(statistics) && nrow(statistics)) {
  statistics$bh_adjusted_p_value <- p.adjust(statistics$p_value, "BH")
  write_validation_csv(statistics, "TIMP1_bulk_validation_statistics.csv")
}
write_validation_csv(expression_source, "TIMP1_bulk_validation_source_data.csv")
write_validation_csv(scores, "TIMP1_signature_scores.csv")
write_validation_csv(coverage, "TIMP1_signature_gene_coverage.csv")

correlations$bh_adjusted_p_value <- ave(
  correlations$p_value, correlations$dataset, correlations$scope,
  FUN = function(x) p.adjust(x, "BH")
)
write_validation_csv(correlations, "TIMP1_signature_correlations.csv")

expression_source$group <- factor(
  expression_source$group,
  levels = c("Control", "AKI", "CKD")
)
expression_plot <- ggplot2::ggplot(
  expression_source,
  ggplot2::aes(group, TIMP1_expression, fill = group)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.35, colour = NA) +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75) +
  ggplot2::geom_jitter(width = 0.10, size = 1.4, alpha = 0.65) +
  ggplot2::facet_wrap(~dataset, scales = "free", nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(Control = "#4C78A8", AKI = "#E45756", CKD = "#72B7B2"),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "TIMP1 expression in independent kidney cohorts",
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  validation_theme() +
  ggplot2::theme(legend.position = "none")
save_validation_plot(
  expression_plot, "TIMP1_bulk_validation_expression",
  width = 10, height = 4.6
)

heatmap_source <- correlations[, c(
  "signature", "dataset", "scope", "rho"
)]
heatmap_source$dataset_scope <- paste(
  heatmap_source$dataset, heatmap_source$scope, sep = " | "
)
heatmap_wide <- reshape(
  heatmap_source[, c("signature", "dataset_scope", "rho")],
  idvar = "signature", timevar = "dataset_scope", direction = "wide"
)
rownames(heatmap_wide) <- heatmap_wide$signature
heatmap_wide$signature <- NULL
heatmap_matrix <- as.matrix(heatmap_wide)
colnames(heatmap_matrix) <- sub("^rho\\.", "", colnames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[names(signature_gene_sets), , drop = FALSE]

draw_heatmap <- function(path, type) {
  if (type == "png") {
    grDevices::png(path, width = 10, height = 7, units = "in", res = 600)
  } else {
    grDevices::cairo_pdf(path, width = 10, height = 7)
  }
  pheatmap::pheatmap(
    heatmap_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      c("#3B6FB6", "white", "#C44E52")
    )(101),
    breaks = seq(-1, 1, length.out = 102),
    display_numbers = TRUE,
    number_format = "%.2f",
    border_color = "white",
    main = "TIMP1 correlation with pathway signature scores"
  )
  grDevices::dev.off()
}
figure_base <- file.path(
  validation_figure_dir(), "TIMP1_signature_correlation_heatmap"
)
draw_heatmap(paste0(figure_base, ".png"), "png")
draw_heatmap(paste0(figure_base, ".pdf"), "pdf")

write_validation_session_info("sessionInfo_bulk_validation.txt")
message("Bulk TIMP1 signature validation completed.")
