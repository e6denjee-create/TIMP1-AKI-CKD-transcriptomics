#!/usr/bin/env Rscript

# Reviewer-driven robustness analyses for the TIMP1 manuscript.
# New outputs supplement, rather than replace, the immutable MVP results.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260614)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

result_dir <- validation_result_dir()
figure_dir <- validation_figure_dir()
datasets <- c("GSE139061", "GSE30718", "GSE66494")

signature_table <- read.csv(
  file.path(result_dir, "signature_gene_sets_used.csv"),
  stringsAsFactors = FALSE
)
signature_names <- unique(signature_table$signature)
signature_list <- split(signature_table$gene, signature_table$signature)
signature_list <- lapply(signature_list, unique)

signature_provenance <- data.frame(
  signature = signature_names,
  construction = c(
    "Investigator-curated composite of canonical matrix structural, cross-linking, protease, and matricellular genes",
    "Investigator-curated collagen-chain and collagen-maturation program",
    "Investigator-curated canonical ligand-receptor-SMAD and profibrotic target program",
    "Investigator-curated cytokine, chemokine, NF-kappaB, adhesion, and STAT response program",
    "Kidney tubular injury marker panel",
    "Kidney failed-repair marker panel integrating injury, epithelial-state, inflammatory, and cell-cycle genes",
    "Cell-cycle arrest and senescence-associated secretory phenotype panel",
    "Investigator-curated renal fibrosis and matrix-deposition panel",
    "Leukocyte, myeloid, antigen-presentation, and complement-associated activation panel"
  ),
  source_basis = c(
    "Kidney fibrosis reviews and TIMP/MMP literature",
    "Reactome-like collagen biosynthesis concepts and kidney fibrosis literature",
    "Canonical TGF-beta/SMAD signaling and renal fibrosis literature",
    "AKI inflammation reviews and canonical inflammatory markers",
    "Human and experimental kidney injury literature",
    "AKI-to-CKD failed-repair and single-cell kidney literature",
    "Kidney senescence and G2/M arrest literature",
    "Renal fibrosis reviews and biopsy transcriptomic literature",
    "Kidney immune-injury and single-cell atlas literature"
  ),
  manuscript_reference_numbers = c(
    "5;23;32-40", "5;23;34", "6;7;34-40", "8-21",
    "3;8;15;22", "1;2;6;7;15;16;22;26", "6;7;17-19",
    "5;23-27;33-40", "8;15;22;24;26;28-31"
  ),
  exact_database_export = "No",
  prespecification_boundary = paste(
    "Gene membership was frozen before analysis of GSE180394;",
    "the sets were not selected using external-cohort TIMP1 correlations."
  ),
  stringsAsFactors = FALSE
)
write_validation_csv(signature_provenance, "signature_provenance_v5.csv")

overlap_rows <- list()
index <- 1L
for (i in seq_along(signature_names)) {
  for (j in i:length(signature_names)) {
    a <- signature_names[[i]]
    b <- signature_names[[j]]
    genes_a <- signature_list[[a]]
    genes_b <- signature_list[[b]]
    intersection <- intersect(genes_a, genes_b)
    union_genes <- union(genes_a, genes_b)
    overlap_rows[[index]] <- data.frame(
      signature_a = a,
      signature_b = b,
      genes_a = length(genes_a),
      genes_b = length(genes_b),
      overlap_count = length(intersection),
      jaccard_index = length(intersection) / length(union_genes),
      overlapping_genes = paste(intersection, collapse = ";"),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
signature_overlap <- do.call(rbind, overlap_rows)
write_validation_csv(signature_overlap, "signature_overlap_audit_v5.csv")

score_discovery <- read.csv(
  file.path(result_dir, "TIMP1_signature_scores.csv"),
  stringsAsFactors = FALSE
)
score_external <- read.csv(
  file.path(result_dir, "external_GSE180394_signature_scores.csv"),
  stringsAsFactors = FALSE
)
external_expression_path <- file.path(
  "data", "external", "GSE180394",
  "external_GSE180394_expression_gene_symbol.csv.gz"
)
external_metadata_path <- file.path(
  "data", "external", "GSE180394", "external_GSE180394_metadata.csv"
)
external_expression <- as.matrix(read.csv(
  gzfile(external_expression_path), row.names = 1, check.names = FALSE
))
storage.mode(external_expression) <- "numeric"
external_metadata <- read.csv(
  external_metadata_path, stringsAsFactors = FALSE, check.names = FALSE
)
external_metadata <- external_metadata[
  match(colnames(external_expression), external_metadata$sample), ,
  drop = FALSE
]
external_scores <- merge(
  data.frame(
    dataset = "GSE180394",
    sample = external_metadata$sample,
    group = external_metadata$group,
    TIMP1_expression = as.numeric(external_expression["TIMP1", ]),
    stringsAsFactors = FALSE
  ),
  score_external,
  by = c("dataset", "sample", "group"),
  all.x = TRUE,
  sort = FALSE
)
all_scores <- rbind(
  score_discovery[, colnames(external_scores), drop = FALSE],
  external_scores
)

inter_signature_rows <- list()
index <- 1L
for (dataset_name in unique(all_scores$dataset)) {
  data <- all_scores[
    all_scores$dataset == dataset_name & all_scores$group != "Control", ,
    drop = FALSE
  ]
  for (i in seq_along(signature_names)) {
    for (j in i:length(signature_names)) {
      a <- signature_names[[i]]
      b <- signature_names[[j]]
      valid <- is.finite(data[[a]]) & is.finite(data[[b]])
      rho <- if (sum(valid) >= 4) {
        suppressWarnings(cor(data[[a]][valid], data[[b]][valid],
          method = "spearman"
        ))
      } else {
        NA_real_
      }
      inter_signature_rows[[index]] <- data.frame(
        dataset = dataset_name,
        scope = "disease_only",
        signature_a = a,
        signature_b = b,
        rho = rho,
        n = sum(valid),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
}
inter_signature <- do.call(rbind, inter_signature_rows)
write_validation_csv(
  inter_signature, "signature_intercorrelations_disease_only_v5.csv"
)

bootstrap_cor <- function(x, y, reps = 2000L) {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  boot <- replicate(reps, {
    selected <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[selected], y[selected], method = "spearman"))
  })
  boot <- boot[is.finite(boot)]
  c(
    estimate = observed,
    ci_lower = unname(quantile(boot, 0.025, na.rm = TRUE)),
    ci_upper = unname(quantile(boot, 0.975, na.rm = TRUE))
  )
}

bootstrap_g <- function(disease, control, reps = 2000L) {
  hedges_g <- function(x, y) {
    pooled_sd <- sqrt(
      ((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
        (length(x) + length(y) - 2)
    )
    correction <- 1 - 3 / (4 * (length(x) + length(y)) - 9)
    ((mean(x) - mean(y)) / pooled_sd) * correction
  }
  observed <- hedges_g(disease, control)
  boot <- replicate(reps, hedges_g(
    sample(disease, replace = TRUE),
    sample(control, replace = TRUE)
  ))
  c(
    hedges_g = observed,
    ci_lower = unname(quantile(boot, 0.025, na.rm = TRUE)),
    ci_upper = unname(quantile(boot, 0.975, na.rm = TRUE))
  )
}

effect_ci_rows <- list()
index <- 1L
for (dataset_name in datasets) {
  dataset <- read_validation_bulk(dataset_name)
  if (is.null(dataset) || !"TIMP1" %in% rownames(dataset$expression)) next
  values <- as.numeric(dataset$expression["TIMP1", dataset$metadata$sample])
  disease <- values[dataset$metadata$group != "Control"]
  control <- values[dataset$metadata$group == "Control"]
  ci <- bootstrap_g(disease, control)
  effect_ci_rows[[index]] <- data.frame(
    dataset = dataset_name,
    contrast = "Disease_vs_Control",
    n_disease = length(disease),
    n_control = length(control),
    hedges_g = ci[["hedges_g"]],
    hedges_g_ci_lower = ci[["ci_lower"]],
    hedges_g_ci_upper = ci[["ci_upper"]],
    bootstrap_repetitions = 2000L,
    stringsAsFactors = FALSE
  )
  index <- index + 1L
}
external_disease <- as.numeric(
  external_expression["TIMP1", external_metadata$group == "Disease"]
)
external_control <- as.numeric(
  external_expression["TIMP1", external_metadata$group == "Control"]
)
ci <- bootstrap_g(external_disease, external_control)
effect_ci_rows[[index]] <- data.frame(
  dataset = "GSE180394",
  contrast = "Disease_vs_Living_Donor",
  n_disease = length(external_disease),
  n_control = length(external_control),
  hedges_g = ci[["hedges_g"]],
  hedges_g_ci_lower = ci[["ci_lower"]],
  hedges_g_ci_upper = ci[["ci_upper"]],
  bootstrap_repetitions = 2000L,
  stringsAsFactors = FALSE
)
effect_ci <- do.call(rbind, effect_ci_rows)
write_validation_csv(effect_ci, "TIMP1_expression_effect_size_ci_v5.csv")

correlation_ci_rows <- list()
index <- 1L
for (dataset_name in unique(all_scores$dataset)) {
  data <- all_scores[
    all_scores$dataset == dataset_name & all_scores$group != "Control", ,
    drop = FALSE
  ]
  for (signature in signature_names) {
    ci <- bootstrap_cor(data$TIMP1_expression, data[[signature]])
    correlation_ci_rows[[index]] <- data.frame(
      dataset = dataset_name,
      scope = "disease_only",
      feature_type = "signature",
      feature = signature,
      rho = ci[["estimate"]],
      rho_ci_lower = ci[["ci_lower"]],
      rho_ci_upper = ci[["ci_upper"]],
      n = sum(is.finite(data$TIMP1_expression) & is.finite(data[[signature]])),
      bootstrap_repetitions = 2000L,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
write_validation_csv(
  do.call(rbind, correlation_ci_rows),
  "TIMP1_signature_correlation_ci_v5.csv"
)

score_gene_set <- function(expression, genes) {
  available <- intersect(toupper(genes), rownames(expression))
  if (length(available) < 2) return(rep(NA_real_, ncol(expression)))
  matrix <- expression[available, , drop = FALSE]
  scaled <- t(scale(t(matrix)))
  scaled[!is.finite(scaled)] <- NA_real_
  colMeans(scaled, na.rm = TRUE)
}

broad_factor_rows <- list()
index <- 1L
for (dataset_name in c(datasets, "GSE180394")) {
  if (dataset_name == "GSE180394") {
    expression <- external_expression
    metadata <- data.frame(
      sample = external_metadata$sample,
      group = external_metadata$group,
      stringsAsFactors = FALSE
    )
  } else {
    dataset <- read_validation_bulk(dataset_name)
    if (is.null(dataset)) next
    expression <- dataset$expression
    metadata <- dataset$metadata
  }
  disease_samples <- metadata$sample[metadata$group != "Control"]
  expression <- expression[, disease_samples, drop = FALSE]
  timp1 <- as.numeric(expression["TIMP1", ])
  for (signature in signature_names) {
    focal_genes <- signature_list[[signature]]
    background_genes <- setdiff(
      unique(unlist(signature_list, use.names = FALSE)),
      c("TIMP1", focal_genes)
    )
    background_genes <- intersect(background_genes, rownames(expression))
    focal_score <- score_gene_set(expression, focal_genes)
    if (length(background_genes) < 5) next
    background_matrix <- t(scale(t(
      expression[background_genes, , drop = FALSE]
    )))
    background_matrix[!is.finite(background_matrix)] <- 0
    pc1 <- prcomp(t(background_matrix), center = FALSE, scale. = FALSE)$x[, 1]
    timp1_residual <- residuals(lm(timp1 ~ pc1))
    focal_residual <- residuals(lm(focal_score ~ pc1))
    test <- suppressWarnings(cor.test(
      timp1_residual, focal_residual, method = "spearman", exact = FALSE
    ))
    broad_factor_rows[[index]] <- data.frame(
      dataset = dataset_name,
      signature = signature,
      n = length(timp1),
      background_gene_count = length(background_genes),
      unadjusted_rho = suppressWarnings(cor(
        timp1, focal_score, method = "spearman"
      )),
      broad_injury_pc1_adjusted_rho = unname(test$estimate),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
broad_factor <- do.call(rbind, broad_factor_rows)
broad_factor$bh_adjusted_p_value <- ave(
  broad_factor$p_value, broad_factor$dataset,
  FUN = function(x) p.adjust(x, method = "BH")
)
write_validation_csv(
  broad_factor, "signature_broad_injury_factor_adjustment_v5.csv"
)

genomewide <- read.csv(
  file.path(result_dir, "TIMP1_genomewide_disease_only_correlations.csv"),
  stringsAsFactors = FALSE
)
stringent <- read.csv(
  file.path(result_dir, "stringent_TIMP1_correlated_module.csv"),
  stringsAsFactors = FALSE
)
external_disease_expression <- external_expression[
  , external_metadata$group == "Disease", drop = FALSE
]
external_timp1 <- as.numeric(external_disease_expression["TIMP1", ])

loco_rows <- list()
index <- 1L
for (left_out in datasets) {
  retained <- setdiff(datasets, left_out)
  wide <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), lapply(
    retained,
    function(dataset_name) {
      x <- genomewide[genomewide$dataset == dataset_name, c(
        "gene", "rho", "bh_adjusted_p_value"
      )]
      colnames(x)[-1] <- paste0(
        c("rho_", "fdr_"), dataset_name
      )
      x
    }
  ))
  selected <- wide[
    is.finite(wide[[paste0("rho_", retained[[1]])]]) &
      is.finite(wide[[paste0("rho_", retained[[2]])]]) &
      wide[[paste0("rho_", retained[[1]])]] > 0 &
      wide[[paste0("rho_", retained[[2]])]] > 0 &
      wide[[paste0("fdr_", retained[[1]])]] < 0.05 &
      wide[[paste0("fdr_", retained[[2]])]] < 0.05,
    "gene"
  ]
  available <- intersect(toupper(selected), rownames(external_disease_expression))
  score <- score_gene_set(external_disease_expression, available)
  test <- suppressWarnings(cor.test(
    external_timp1, score, method = "spearman", exact = FALSE
  ))
  loco_rows[[index]] <- data.frame(
    left_out_discovery_cohort = left_out,
    retained_discovery_cohorts = paste(retained, collapse = ";"),
    selected_gene_count = length(unique(selected)),
    external_available_gene_count = length(available),
    overlap_with_full_module = length(intersect(selected, stringent$gene)),
    jaccard_with_full_module = length(intersect(selected, stringent$gene)) /
      length(union(selected, stringent$gene)),
    external_disease_only_rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
  index <- index + 1L
}
loco <- do.call(rbind, loco_rows)
loco$bh_adjusted_p_value <- p.adjust(loco$p_value, method = "BH")
write_validation_csv(loco, "module_leave_one_discovery_out_v5.csv")

module_genes <- intersect(
  toupper(stringent$gene), rownames(external_disease_expression)
)
module_score <- score_gene_set(external_disease_expression, module_genes)
observed_rho <- suppressWarnings(cor(
  external_timp1, module_score, method = "spearman"
))
candidate_genes <- setdiff(
  rownames(external_disease_expression), c("TIMP1", module_genes)
)
gene_mean <- rowMeans(external_disease_expression, na.rm = TRUE)
gene_sd <- apply(external_disease_expression, 1, sd, na.rm = TRUE)
valid_genes <- names(gene_mean)[
  is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
]
candidate_genes <- intersect(candidate_genes, valid_genes)
module_genes <- intersect(module_genes, valid_genes)
mean_breaks <- unique(quantile(gene_mean[valid_genes], seq(0, 1, 0.1),
  na.rm = TRUE
))
sd_breaks <- unique(quantile(gene_sd[valid_genes], seq(0, 1, 0.1),
  na.rm = TRUE
))
mean_bin <- cut(gene_mean, breaks = mean_breaks, include.lowest = TRUE)
sd_bin <- cut(gene_sd, breaks = sd_breaks, include.lowest = TRUE)
stratum <- interaction(mean_bin, sd_bin, drop = TRUE)
names(stratum) <- names(gene_mean)

random_rho <- numeric(1000L)
external_scaled <- t(scale(t(external_disease_expression[valid_genes, , drop = FALSE])))
external_scaled[!is.finite(external_scaled)] <- 0
module_stratum_counts <- table(stratum[module_genes])
candidate_by_stratum <- split(candidate_genes, stratum[candidate_genes])
for (iteration in seq_along(random_rho)) {
  sampled <- character()
  for (stratum_name in names(module_stratum_counts)) {
    needed <- unname(module_stratum_counts[[stratum_name]])
    pool <- setdiff(candidate_by_stratum[[stratum_name]], sampled)
    if (length(pool) >= needed) {
      sampled <- c(sampled, sample(pool, needed, replace = FALSE))
    } else {
      sampled <- c(sampled, pool)
      remaining <- needed - length(pool)
      fallback <- setdiff(candidate_genes, sampled)
      sampled <- c(sampled, sample(fallback, remaining, replace = FALSE))
    }
  }
  random_score <- colMeans(external_scaled[sampled, , drop = FALSE])
  random_rho[[iteration]] <- suppressWarnings(cor(
    external_timp1, random_score, method = "spearman"
  ))
}
random_module_summary <- data.frame(
  dataset = "GSE180394",
  scope = "disease_only",
  observed_module_gene_count = length(module_genes),
  observed_rho = observed_rho,
  random_repetitions = length(random_rho),
  random_rho_mean = mean(random_rho, na.rm = TRUE),
  random_rho_sd = sd(random_rho, na.rm = TRUE),
  random_rho_95th_percentile = unname(quantile(random_rho, 0.95, na.rm = TRUE)),
  empirical_one_sided_p = (
    1 + sum(random_rho >= observed_rho, na.rm = TRUE)
  ) / (1 + sum(is.finite(random_rho))),
  matching_variables = "External-cohort mean-expression decile and variability decile",
  stringsAsFactors = FALSE
)
write_validation_csv(
  random_module_summary, "module_random_matched_benchmark_v5.csv"
)
write_validation_csv(
  data.frame(iteration = seq_along(random_rho), random_rho = random_rho),
  "module_random_matched_null_distribution_v5.csv"
)

null_plot <- ggplot2::ggplot(
  data.frame(random_rho = random_rho),
  ggplot2::aes(random_rho)
) +
  ggplot2::geom_histogram(
    bins = 35, fill = "#8FB3D9", colour = "white"
  ) +
  ggplot2::geom_vline(
    xintercept = observed_rho, colour = "#B22222", linewidth = 1.1
  ) +
  ggplot2::annotate(
    "text", x = observed_rho, y = Inf,
    label = sprintf("Observed rho = %.2f", observed_rho),
    hjust = 1.05, vjust = 1.5, colour = "#B22222", fontface = "bold"
  ) +
  ggplot2::labs(
    title = "Matched random-module benchmark in GSE180394",
    x = "Disease-only Spearman rho with TIMP1",
    y = "Random modules"
  ) +
  validation_theme()
save_validation_plot(
  null_plot, "v5_module_random_matched_benchmark",
  width = 7.2, height = 5.2
)

signature_corr_discovery <- read.csv(
  file.path(result_dir, "TIMP1_signature_correlations.csv"),
  stringsAsFactors = FALSE
)
signature_corr_external <- read.csv(
  file.path(result_dir, "external_GSE180394_signature_correlations.csv"),
  stringsAsFactors = FALSE
)
colnames(signature_corr_external)[
  colnames(signature_corr_external) == "feature"
] <- "signature"
signature_corr <- rbind(signature_corr_discovery, signature_corr_external)
signature_corr <- signature_corr[signature_corr$scope == "disease_only", ]
signature_corr$dataset <- factor(
  signature_corr$dataset,
  levels = c(datasets, "GSE180394")
)
signature_corr$signature <- factor(
  signature_corr$signature, levels = rev(signature_names)
)
signature_corr$label <- sprintf(
  "%.2f%s", signature_corr$rho,
  ifelse(signature_corr$bh_adjusted_p_value < 0.05, "*", "")
)
signature_plot <- ggplot2::ggplot(
  signature_corr,
  ggplot2::aes(dataset, signature, fill = rho)
) +
  ggplot2::geom_tile(colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C94F57",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Disease-only TIMP1-signature correlations",
    subtitle = "* BH-adjusted P < 0.05 within dataset and scope",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  validation_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )
save_validation_plot(
  signature_plot, "v5_manuscript_figure_3_signature_correlation_heatmap",
  width = 9.5, height = 6.8
)

core_summary <- read.csv(
  file.path(result_dir, "core_module_genes_summary.csv"),
  stringsAsFactors = FALSE
)
core_external <- read.csv(
  file.path(result_dir, "external_GSE180394_core_module_correlations.csv"),
  stringsAsFactors = FALSE
)
core_long <- rbind(
  data.frame(
    gene = core_summary$gene, dataset = "GSE139061",
    rho = core_summary$rho.GSE139061, significant = TRUE
  ),
  data.frame(
    gene = core_summary$gene, dataset = "GSE30718",
    rho = core_summary$rho.GSE30718, significant = TRUE
  ),
  data.frame(
    gene = core_summary$gene, dataset = "GSE66494",
    rho = core_summary$rho.GSE66494, significant = TRUE
  ),
  data.frame(
    gene = core_external$feature,
    dataset = "GSE180394",
    rho = core_external$rho,
    significant = core_external$bh_adjusted_p_value < 0.05
  )
)
core_long$dataset <- factor(
  core_long$dataset, levels = c(datasets, "GSE180394")
)
core_long$gene <- factor(core_long$gene, levels = rev(core_summary$gene))
core_long$label <- sprintf(
  "%.2f%s", core_long$rho, ifelse(core_long$significant, "*", "")
)
core_plot <- ggplot2::ggplot(
  core_long, ggplot2::aes(dataset, gene, fill = rho)
) +
  ggplot2::geom_tile(colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.1) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C94F57",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Disease-only TIMP1 correlations across the 12-gene core",
    subtitle = "* BH-adjusted P < 0.05",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  validation_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )
save_validation_plot(
  core_plot, "v5_manuscript_figure_5_core_gene_correlation_heatmap",
  width = 9.0, height = 7.2
)

expression_source <- read.csv(
  file.path(result_dir, "TIMP1_bulk_validation_source_data.csv"),
  stringsAsFactors = FALSE
)
external_source <- data.frame(
  dataset = "GSE180394",
  sample = external_metadata$sample,
  group = external_metadata$group,
  TIMP1_expression = as.numeric(external_expression["TIMP1", ]),
  stringsAsFactors = FALSE
)
expression_source$group[expression_source$group != "Control"] <- "Disease"
expression_source <- rbind(expression_source, external_source)
expression_source$dataset <- factor(
  expression_source$dataset, levels = c(datasets, "GSE180394")
)
expression_source$group <- factor(
  expression_source$group, levels = c("Control", "Disease")
)
expression_plot <- ggplot2::ggplot(
  expression_source,
  ggplot2::aes(group, TIMP1_expression, colour = group)
) +
  ggplot2::geom_boxplot(
    width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.6
  ) +
  ggplot2::geom_jitter(
    width = 0.12, height = 0, alpha = 0.72, size = 1.7
  ) +
  ggplot2::facet_wrap(~dataset, scales = "free_y", nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(Control = "#356D9E", Disease = "#C94F57")
  ) +
  ggplot2::labs(
    title = "TIMP1 expression across kidney injury and fibrosis cohorts",
    subtitle = "Points represent biological samples; expression scales are cohort-specific",
    x = NULL, y = "Within-cohort normalized expression", colour = NULL
  ) +
  validation_theme() +
  ggplot2::theme(legend.position = "top")
save_validation_plot(
  expression_plot, "v5_manuscript_figure_2_TIMP1_expression_cohorts",
  width = 12.5, height = 5.8
)

interpretation <- c(
  "# Reviewer-driven robustness analyses (v5)",
  "",
  paste0(
    "- The nine signatures were investigator-curated and frozen before ",
    "GSE180394 analysis; they are not exact exports from a single database."
  ),
  paste0(
    "- Signature overlap and inter-signature correlations confirm that the ",
    "scores represent related, non-independent injury-response programs."
  ),
  paste0(
    "- Broad-injury PC1 adjustment is a deliberately stringent sensitivity ",
    "analysis. Attenuation is interpreted as evidence of shared injury-state ",
    "variance, not as evidence against the unadjusted association."
  ),
  sprintf(
    paste0(
      "- The external stringent-module correlation was %.3f. Against 1000 ",
      "mean/variability-matched random modules, the empirical one-sided P ",
      "value was %.4g."
    ),
    as.numeric(observed_rho),
    as.numeric(random_module_summary$empirical_one_sided_p)
  ),
  paste0(
    "- Leave-one-discovery-cohort-out reconstruction tests whether external ",
    "module replication depends on any single discovery cohort."
  ),
  "",
  "Limitations: these analyses do not remove unmeasured clinical severity,",
  "cell-composition, treatment, or technical confounding. They strengthen an",
  "association claim but do not establish TIMP1 specificity or causality."
)
writeLines(
  interpretation,
  file.path(result_dir, "revision_robustness_v5_interpretation.md")
)
write_validation_session_info("sessionInfo_revision_robustness_v5.txt")

message("Revision robustness v5 analyses completed.")
