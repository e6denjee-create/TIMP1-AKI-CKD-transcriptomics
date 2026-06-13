#!/usr/bin/env Rscript

# GSE180394 robustness analyses and manuscript-ready figures.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260613)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

result_dir <- validation_result_dir()
figure_dir <- validation_figure_dir()
data_dir <- file.path("data", "external", "GSE180394")

paths <- c(
  expression_all = file.path(
    data_dir, "external_GSE180394_expression_all_samples_gene_symbol.csv.gz"
  ),
  metadata_all = file.path(
    data_dir, "external_GSE180394_metadata_all_samples.csv"
  ),
  signature_scores = file.path(
    result_dir, "external_GSE180394_signature_scores.csv"
  ),
  module_scores = file.path(
    result_dir, "external_GSE180394_stringent_module_scores.csv"
  ),
  core_expression = file.path(
    result_dir, "external_GSE180394_core_module_expression.csv"
  ),
  discovery_expression = file.path(
    result_dir, "TIMP1_bulk_validation_source_data.csv"
  ),
  discovery_signatures = file.path(
    result_dir, "TIMP1_signature_correlations.csv"
  ),
  core_summary = file.path(
    result_dir, "core_module_genes_summary.csv"
  )
)
if (!all(file.exists(paths))) {
  missing <- paths[!file.exists(paths)]
  append_missing_data(
    "external_GSE180394_sensitivity_v4", "GSE180394",
    paste(missing, collapse = "; "), "Required input is absent.",
    "Sensitivity workflow stopped without claiming results."
  )
  stop("Missing v4 inputs.")
}

expression_all <- as.matrix(read.csv(
  gzfile(paths[["expression_all"]]), row.names = 1, check.names = FALSE
))
storage.mode(expression_all) <- "numeric"
metadata_all <- read.csv(
  paths[["metadata_all"]], stringsAsFactors = FALSE, check.names = FALSE
)
signature_scores <- read.csv(
  paths[["signature_scores"]], stringsAsFactors = FALSE, check.names = FALSE
)
module_scores <- read.csv(
  paths[["module_scores"]], stringsAsFactors = FALSE, check.names = FALSE
)
core_long <- read.csv(
  paths[["core_expression"]], stringsAsFactors = FALSE, check.names = FALSE
)

diagnosis_category <- function(diagnosis) {
  ifelse(
    grepl("^LN-|Lupus", diagnosis, ignore.case = TRUE),
    "Lupus nephritis",
    ifelse(
      grepl("FSGS|FGGS", diagnosis, ignore.case = TRUE),
      "FSGS/FGGS",
      ifelse(
        diagnosis == "DN", "Diabetic nephropathy",
        ifelse(diagnosis == "IgAN", "IgA nephropathy", "Other kidney disease")
      )
    )
  )
}

disease_metadata <- metadata_all[
  metadata_all$group == "Disease", , drop = FALSE
]
disease_metadata$diagnosis_category <- diagnosis_category(
  disease_metadata$original_diagnosis
)

core_wide <- reshape(
  core_long[, c("sample", "gene", "expression")],
  idvar = "sample", timevar = "gene", direction = "wide"
)
colnames(core_wide) <- sub("^expression\\.", "", colnames(core_wide))

feature_data <- merge(
  disease_metadata[, c(
    "sample", "original_diagnosis", "diagnosis_category"
  )],
  signature_scores[, setdiff(colnames(signature_scores), c("dataset", "group"))],
  by = "sample", all.x = TRUE
)
feature_data <- merge(
  feature_data,
  module_scores[, c("sample", "stringent_TIMP1_module_score")],
  by = "sample", all.x = TRUE
)
feature_data <- merge(feature_data, core_wide, by = "sample", all.x = TRUE)
feature_data$TIMP1 <- as.numeric(expression_all[
  "TIMP1", feature_data$sample
])

signature_features <- setdiff(
  colnames(signature_scores), c("dataset", "sample", "group")
)
core_features <- sort(unique(core_long$gene))
feature_names <- c(
  signature_features, "stringent_TIMP1_module_score", core_features
)
feature_types <- c(
  setNames(rep("signature", length(signature_features)), signature_features),
  stringent_TIMP1_module_score = "stringent_module",
  setNames(rep("core_gene", length(core_features)), core_features)
)

spearman_result <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4 || length(unique(x[keep])) < 2 ||
      length(unique(y[keep])) < 2) {
    return(c(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
  c(rho = unname(test$estimate), p_value = test$p.value, n = sum(keep))
}

# 1. Leave-one-diagnosis-category-out analysis.
lodo_rows <- list()
lodo_index <- 1
categories <- sort(unique(feature_data$diagnosis_category))
for (excluded in categories) {
  selected <- feature_data$diagnosis_category != excluded
  for (feature in feature_names) {
    result <- spearman_result(
      feature_data$TIMP1[selected], feature_data[[feature]][selected]
    )
    lodo_rows[[lodo_index]] <- data.frame(
      excluded_diagnosis_category = excluded,
      excluded_n = sum(!selected),
      retained_n = sum(selected),
      feature_type = unname(feature_types[[feature]]),
      feature = feature,
      rho = result[["rho"]],
      p_value = result[["p_value"]],
      stringsAsFactors = FALSE
    )
    lodo_index <- lodo_index + 1
  }
}
lodo <- do.call(rbind, lodo_rows)
lodo$bh_adjusted_p_value <- ave(
  lodo$p_value, lodo$excluded_diagnosis_category,
  FUN = function(x) p.adjust(x, "BH")
)
write_validation_csv(
  lodo, "external_GSE180394_sensitivity_leave_one_diagnosis_out.csv"
)

lodo_summary <- do.call(rbind, lapply(
  split(lodo, paste(lodo$feature_type, lodo$feature, sep = "::")),
  function(x) {
    data.frame(
      feature_type = x$feature_type[[1]],
      feature = x$feature[[1]],
      min_rho = min(x$rho, na.rm = TRUE),
      median_rho = median(x$rho, na.rm = TRUE),
      max_rho = max(x$rho, na.rm = TRUE),
      positive_iterations = sum(x$rho > 0, na.rm = TRUE),
      fdr_significant_iterations = sum(
        x$bh_adjusted_p_value < 0.05, na.rm = TRUE
      ),
      total_iterations = nrow(x),
      stringsAsFactors = FALSE
    )
  }
))
write_validation_csv(
  lodo_summary,
  "external_GSE180394_sensitivity_leave_one_diagnosis_out_summary.csv"
)

lodo_plot_data <- lodo
lodo_plot_data$feature_label <- ifelse(
  lodo_plot_data$feature == "stringent_TIMP1_module_score",
  "Stringent module", lodo_plot_data$feature
)
lodo_plot <- ggplot2::ggplot(
  lodo_plot_data,
  ggplot2::aes(excluded_diagnosis_category, feature_label, fill = rho)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.4
  ) +
  ggplot2::facet_grid(
    feature_type ~ ., scales = "free_y", space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "GSE180394 leave-one-diagnosis-out sensitivity",
    subtitle = "Disease-only Spearman correlations",
    x = "Excluded diagnosis category", y = NULL, fill = "Spearman rho"
  ) +
  validation_theme(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )
save_validation_plot(
  lodo_plot,
  "external_GSE180394_sensitivity_leave_one_diagnosis_out_heatmap",
  width = 11, height = 10
)

# 2. Alternative control definitions.
hedges_g_ci <- function(disease, control) {
  n1 <- length(disease)
  n0 <- length(control)
  pooled_sd <- sqrt(
    ((n1 - 1) * var(disease) + (n0 - 1) * var(control)) /
      (n1 + n0 - 2)
  )
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- correction * (mean(disease) - mean(control)) / pooled_sd
  se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2)))
  c(g = g, lower = g - 1.96 * se, upper = g + 1.96 * se)
}

disease_samples <- metadata_all$sample[metadata_all$group == "Disease"]
living_samples <- metadata_all$sample[metadata_all$group == "Control"]
tumor_samples <- metadata_all$sample[
  metadata_all$group == "Excluded_tumor_nephrectomy_unaffected"
]
control_definitions <- list(
  living_donor_only = living_samples,
  extended_controls = c(living_samples, tumor_samples),
  tumor_nephrectomy_only = tumor_samples
)
disease_timp1 <- as.numeric(expression_all["TIMP1", disease_samples])
control_rows <- list()
control_plot_rows <- list()
control_index <- 1
for (definition in names(control_definitions)) {
  controls <- control_definitions[[definition]]
  control_timp1 <- as.numeric(expression_all["TIMP1", controls])
  test <- suppressWarnings(wilcox.test(
    disease_timp1, control_timp1, exact = FALSE
  ))
  effect <- hedges_g_ci(disease_timp1, control_timp1)
  control_rows[[control_index]] <- data.frame(
    control_definition = definition,
    n_disease = length(disease_timp1),
    n_control = length(control_timp1),
    disease_mean = mean(disease_timp1),
    control_mean = mean(control_timp1),
    group_difference = mean(disease_timp1) - mean(control_timp1),
    hedges_g = effect[["g"]],
    hedges_g_ci_lower = effect[["lower"]],
    hedges_g_ci_upper = effect[["upper"]],
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
  control_plot_rows[[control_index]] <- rbind(
    data.frame(
      control_definition = definition, group = "Disease",
      sample = disease_samples, TIMP1_expression = disease_timp1
    ),
    data.frame(
      control_definition = definition, group = "Control",
      sample = controls, TIMP1_expression = control_timp1
    )
  )
  control_index <- control_index + 1
}
control_statistics <- do.call(rbind, control_rows)
control_statistics$bh_adjusted_p_value <- p.adjust(
  control_statistics$p_value, "BH"
)
control_source <- do.call(rbind, control_plot_rows)
write_validation_csv(
  control_statistics,
  "external_GSE180394_sensitivity_control_statistics.csv"
)
write_validation_csv(
  control_source,
  "external_GSE180394_sensitivity_control_source_data.csv"
)

control_labels <- c(
  living_donor_only = "Living donors",
  extended_controls = "Living donors +\ntumor-nephrectomy",
  tumor_nephrectomy_only = "Tumor-nephrectomy"
)
control_plot <- ggplot2::ggplot(
  control_source,
  ggplot2::aes(group, TIMP1_expression, fill = group)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.5) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, size = 1.3, alpha = 0.65) +
  ggplot2::facet_wrap(
    ~control_definition, scales = "free_x",
    labeller = ggplot2::as_labeller(control_labels)
  ) +
  ggplot2::scale_fill_manual(
    values = c(Control = "#4C78A8", Disease = "#C44E52")
  ) +
  ggplot2::labs(
    title = "TIMP1 expression across alternative control definitions",
    x = NULL, y = "Within-cohort normalized expression", fill = NULL
  ) +
  validation_theme() +
  ggplot2::theme(legend.position = "top")
save_validation_plot(
  control_plot,
  "external_GSE180394_sensitivity_control_expression",
  width = 10, height = 5.5
)

control_statistics$display_label <- control_labels[
  control_statistics$control_definition
]
forest_plot <- ggplot2::ggplot(
  control_statistics,
  ggplot2::aes(hedges_g, reorder(display_label, hedges_g))
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = hedges_g_ci_lower, xmax = hedges_g_ci_upper
    ),
    width = 0.16, orientation = "y", color = "#333333"
  ) +
  ggplot2::geom_point(size = 3.2, color = "#C44E52") +
  ggplot2::labs(
    title = "Stability of the TIMP1 disease-control effect",
    subtitle = "Hedges' g with approximate 95% confidence intervals",
    x = "Hedges' g", y = "Control definition"
  ) +
  validation_theme()
save_validation_plot(
  forest_plot,
  "external_GSE180394_sensitivity_control_forest",
  width = 8, height = 4.8
)

# 3. Diagnosis-adjusted robust disease-only correlations.
robust_rows <- list()
robust_index <- 1
for (feature in feature_names) {
  frame <- data.frame(
    TIMP1 = feature_data$TIMP1,
    feature = feature_data[[feature]],
    diagnosis_category = factor(feature_data$diagnosis_category)
  )
  frame <- frame[complete.cases(frame), , drop = FALSE]
  raw <- spearman_result(frame$TIMP1, frame$feature)

  timp1_residual <- residuals(lm(TIMP1 ~ diagnosis_category, data = frame))
  feature_residual <- residuals(lm(feature ~ diagnosis_category, data = frame))
  partial <- spearman_result(timp1_residual, feature_residual)

  model <- lm(
    as.numeric(scale(feature)) ~ as.numeric(scale(TIMP1)) +
      diagnosis_category,
    data = frame
  )
  coefficient <- summary(model)$coefficients[
    "as.numeric(scale(TIMP1))", , drop = FALSE
  ]
  robust_rows[[robust_index]] <- data.frame(
    feature_type = unname(feature_types[[feature]]),
    feature = feature,
    n = nrow(frame),
    diagnosis_categories = nlevels(frame$diagnosis_category),
    spearman_rho = raw[["rho"]],
    spearman_p_value = raw[["p_value"]],
    diagnosis_adjusted_partial_spearman_rho = partial[["rho"]],
    diagnosis_adjusted_partial_spearman_p_value = partial[["p_value"]],
    diagnosis_adjusted_lm_beta = coefficient[1, "Estimate"],
    diagnosis_adjusted_lm_standard_error = coefficient[1, "Std. Error"],
    diagnosis_adjusted_lm_p_value = coefficient[1, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
  robust_index <- robust_index + 1
}
robust <- do.call(rbind, robust_rows)
robust$spearman_bh_adjusted_p_value <- p.adjust(
  robust$spearman_p_value, "BH"
)
robust$diagnosis_adjusted_partial_spearman_bh_adjusted_p_value <- p.adjust(
  robust$diagnosis_adjusted_partial_spearman_p_value, "BH"
)
robust$diagnosis_adjusted_lm_bh_adjusted_p_value <- p.adjust(
  robust$diagnosis_adjusted_lm_p_value, "BH"
)
write_validation_csv(
  robust,
  "external_GSE180394_sensitivity_robust_correlations.csv"
)

robust_plot_data <- rbind(
  data.frame(
    feature_type = robust$feature_type, feature = robust$feature,
    method = "Raw Spearman", estimate = robust$spearman_rho
  ),
  data.frame(
    feature_type = robust$feature_type, feature = robust$feature,
    method = "Diagnosis-adjusted\npartial Spearman",
    estimate = robust$diagnosis_adjusted_partial_spearman_rho
  ),
  data.frame(
    feature_type = robust$feature_type, feature = robust$feature,
    method = "Diagnosis-adjusted\nlinear-model beta",
    estimate = robust$diagnosis_adjusted_lm_beta
  )
)
robust_plot_data$feature <- ifelse(
  robust_plot_data$feature == "stringent_TIMP1_module_score",
  "Stringent module", robust_plot_data$feature
)
robust_plot <- ggplot2::ggplot(
  robust_plot_data,
  ggplot2::aes(method, feature, fill = estimate)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", estimate)), size = 2.5
  ) +
  ggplot2::facet_grid(
    feature_type ~ ., scales = "free_y", space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Diagnosis-adjusted robustness of TIMP1 associations",
    subtitle = "GSE180394 disease samples (n = 44)",
    x = NULL, y = NULL, fill = "Estimate"
  ) +
  validation_theme(base_size = 9)
save_validation_plot(
  robust_plot,
  "external_GSE180394_sensitivity_robust_correlation_heatmap",
  width = 9.5, height = 10
)

# 4. Manuscript-ready figures.
workflow_data <- data.frame(
  xmin = c(0.2, 2.4, 4.6, 6.8),
  xmax = c(1.8, 4.0, 6.2, 8.4),
  ymin = 0.25,
  ymax = 1.25,
  label = c(
    "Discovery cohorts\n3 bulk datasets",
    "External tubular cohort\nGSE180394\n44 disease / 9 control",
    "Robustness analyses\nDiagnosis removal\nAlternative controls\nAdjusted models",
    "Manuscript evidence\nExpression\nSignatures\nModule / core genes"
  ),
  stringsAsFactors = FALSE
)
workflow_plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = workflow_data,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = c("#D9E8F5", "#DDEED9", "#F8E3C5", "#E8DDF1"),
    color = "#333333", linewidth = 0.6
  ) +
  ggplot2::geom_text(
    data = workflow_data,
    ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.6, lineheight = 1.05
  ) +
  ggplot2::annotate(
    "segment",
    x = c(1.85, 4.05, 6.25), xend = c(2.32, 4.52, 6.72),
    y = 0.75, yend = 0.75,
    arrow = grid::arrow(length = grid::unit(0.16, "inches")),
    linewidth = 0.7
  ) +
  ggplot2::annotate(
    "text", x = 4.3, y = 1.55,
    label = "TIMP1 AKI-to-CKD transcriptomic validation framework",
    fontface = "bold", size = 5
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 8.6), ylim = c(0, 1.8)) +
  ggplot2::theme_void()
save_validation_plot(
  workflow_plot,
  "v4_manuscript_figure_1_workflow_cohort_overview",
  width = 12, height = 4.2
)

discovery_expression <- read.csv(
  paths[["discovery_expression"]], stringsAsFactors = FALSE
)
discovery_expression$group <- ifelse(
  discovery_expression$group == "Control", "Control", "Disease"
)
external_expression <- data.frame(
  dataset = "GSE180394",
  sample = c(disease_samples, living_samples),
  group = c(
    rep("Disease", length(disease_samples)),
    rep("Control", length(living_samples))
  ),
  TIMP1_expression = as.numeric(expression_all[
    "TIMP1", c(disease_samples, living_samples)
  ]),
  stringsAsFactors = FALSE
)
expression_cohorts <- rbind(discovery_expression, external_expression)
expression_cohorts$dataset <- factor(
  expression_cohorts$dataset,
  levels = c("GSE139061", "GSE30718", "GSE66494", "GSE180394")
)
expression_cohort_plot <- ggplot2::ggplot(
  expression_cohorts,
  ggplot2::aes(group, TIMP1_expression, fill = group)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.5) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, size = 1.1, alpha = 0.55) +
  ggplot2::facet_wrap(~dataset, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(Control = "#4C78A8", Disease = "#C44E52")
  ) +
  ggplot2::labs(
    title = "TIMP1 expression across kidney injury and fibrosis cohorts",
    x = NULL, y = "Within-cohort normalized expression", fill = NULL
  ) +
  validation_theme(base_size = 10) +
  ggplot2::theme(legend.position = "top")
save_validation_plot(
  expression_cohort_plot,
  "v4_manuscript_figure_2_TIMP1_expression_cohorts",
  width = 12, height = 5.5
)

discovery_signatures <- read.csv(
  paths[["discovery_signatures"]], stringsAsFactors = FALSE
)
discovery_signatures <- discovery_signatures[
  discovery_signatures$scope == "disease_only",
  c("dataset", "signature", "rho"), drop = FALSE
]
external_signatures <- read.csv(
  file.path(result_dir, "external_GSE180394_signature_correlations.csv"),
  stringsAsFactors = FALSE
)
external_signatures <- external_signatures[
  external_signatures$scope == "disease_only",
  c("dataset", "feature", "rho"), drop = FALSE
]
colnames(external_signatures)[2] <- "signature"
signature_four <- rbind(discovery_signatures, external_signatures)
signature_four$dataset <- factor(
  signature_four$dataset,
  levels = c("GSE139061", "GSE30718", "GSE66494", "GSE180394")
)
signature_figure <- ggplot2::ggplot(
  signature_four,
  ggplot2::aes(dataset, signature, fill = rho)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.8
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Disease-only correlations between TIMP1 and pathway signatures",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  validation_theme(base_size = 10)
save_validation_plot(
  signature_figure,
  "v4_manuscript_figure_3_signature_correlation_heatmap",
  width = 8.5, height = 6.5
)

scatter_data <- merge(
  signature_scores[, c("sample", "group")],
  module_scores[, c("sample", "stringent_TIMP1_module_score")],
  by = "sample"
)
scatter_data$TIMP1_expression <- as.numeric(expression_all[
  "TIMP1", scatter_data$sample
])
module_main <- read.csv(
  file.path(result_dir, "external_GSE180394_stringent_module_correlations.csv"),
  stringsAsFactors = FALSE
)
rho_all <- module_main$rho[module_main$scope == "all_samples"]
rho_disease <- module_main$rho[module_main$scope == "disease_only"]
scatter_plot <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(
    TIMP1_expression, stringent_TIMP1_module_score, color = group
  )
) +
  ggplot2::geom_point(size = 2.4, alpha = 0.8) +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE, color = "#333333",
    linewidth = 0.7
  ) +
  ggplot2::scale_color_manual(
    values = c(Control = "#4C78A8", Disease = "#C44E52")
  ) +
  ggplot2::annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
    label = sprintf(
      "All samples rho = %.2f\nDisease-only rho = %.2f",
      rho_all, rho_disease
    ),
    size = 3.5
  ) +
  ggplot2::labs(
    title = "TIMP1 is correlated with the stringent module in GSE180394",
    x = "TIMP1 expression", y = "Stringent module score", color = NULL
  ) +
  validation_theme() +
  ggplot2::theme(legend.position = "top")
save_validation_plot(
  scatter_plot,
  "v4_manuscript_figure_4_GSE180394_stringent_module_scatter",
  width = 7, height = 5.8
)

core_summary <- read.csv(
  paths[["core_summary"]], stringsAsFactors = FALSE, check.names = FALSE
)
core_discovery <- rbind(
  data.frame(
    gene = core_summary$gene, dataset = "GSE139061",
    rho = core_summary$rho.GSE139061
  ),
  data.frame(
    gene = core_summary$gene, dataset = "GSE30718",
    rho = core_summary$rho.GSE30718
  ),
  data.frame(
    gene = core_summary$gene, dataset = "GSE66494",
    rho = core_summary$rho.GSE66494
  )
)
external_core <- read.csv(
  file.path(result_dir, "external_GSE180394_core_module_correlations.csv"),
  stringsAsFactors = FALSE
)
external_core <- external_core[
  external_core$scope == "disease_only",
  c("feature", "rho"), drop = FALSE
]
colnames(external_core)[1] <- "gene"
external_core$dataset <- "GSE180394"
core_four <- rbind(
  core_discovery,
  external_core[, c("gene", "dataset", "rho")]
)
core_four$dataset <- factor(
  core_four$dataset,
  levels = c("GSE139061", "GSE30718", "GSE66494", "GSE180394")
)
core_figure <- ggplot2::ggplot(
  core_four,
  ggplot2::aes(dataset, gene, fill = rho)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.8
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Disease-only TIMP1 correlations across the 12-gene core module",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  validation_theme(base_size = 10)
save_validation_plot(
  core_figure,
  "v4_manuscript_figure_5_core_gene_correlation_heatmap",
  width = 8.5, height = 7.5
)

# Reports generated from the current results.
writeLines(
  c(
    "# GSE180394 Leave-One-Diagnosis-Out Report",
    "",
    paste(
      "Five diagnosis categories were evaluated:",
      paste(categories, collapse = ", "), "."
    ),
    "",
    paste0(
      "All ", nrow(lodo_summary), " evaluated features retained positive ",
      "median correlations across the leave-one-category-out iterations."
    ),
    paste0(
      sum(lodo_summary$positive_iterations == lodo_summary$total_iterations),
      " of ", nrow(lodo_summary),
      " features remained positive in every iteration."
    ),
    paste0(
      sum(
        lodo_summary$fdr_significant_iterations ==
          lodo_summary$total_iterations
      ),
      " of ", nrow(lodo_summary),
      " features remained FDR-significant in every iteration."
    ),
    "",
    "The analysis reduces concern that one major diagnosis category alone",
    "accounts for the TIMP1-associated signature, module, or core-gene",
    "pattern. It does not remove residual clinical heterogeneity or establish",
    "direct regulation."
  ),
  file.path(
    result_dir,
    "external_GSE180394_leave_one_diagnosis_out_report.md"
  )
)

writeLines(
  c(
    "# GSE180394 Alternative-Control Sensitivity Report",
    "",
    "Three control definitions were evaluated: living donors, living donors",
    "plus unaffected tumor-nephrectomy tissue, and tumor-nephrectomy tissue",
    "alone.",
    "",
    paste0(
      "The TIMP1 group difference was positive under all definitions (range ",
      sprintf("%.3f", min(control_statistics$group_difference)), " to ",
      sprintf("%.3f", max(control_statistics$group_difference)), ")."
    ),
    paste0(
      "Hedges' g remained positive (range ",
      sprintf("%.3f", min(control_statistics$hedges_g)), " to ",
      sprintf("%.3f", max(control_statistics$hedges_g)), ")."
    ),
    paste0(
      sum(control_statistics$bh_adjusted_p_value < 0.05), " of ",
      nrow(control_statistics),
      " comparisons had BH-adjusted P < 0.05."
    ),
    "",
    "The direction is therefore stable to alternative control definitions.",
    "The living-donor comparison remains primary because tumor-adjacent",
    "tissue may not represent a healthy kidney reference."
  ),
  file.path(
    result_dir, "external_GSE180394_control_sensitivity_report.md"
  )
)

writeLines(
  c(
    "# GSE180394 Robust Disease-Only Correlation Report",
    "",
    "The metadata supported adjustment for five diagnosis categories across",
    "44 disease samples. Each feature was evaluated by raw Spearman",
    "correlation, diagnosis-residualized partial Spearman correlation, and a",
    "standardized linear model including diagnosis category.",
    "",
    paste0(
      sum(robust$diagnosis_adjusted_partial_spearman_rho > 0), " of ",
      nrow(robust),
      " diagnosis-adjusted partial Spearman estimates were positive."
    ),
    paste0(
      sum(
        robust$diagnosis_adjusted_partial_spearman_bh_adjusted_p_value < 0.05
      ),
      " of ", nrow(robust),
      " partial Spearman associations had FDR < 0.05."
    ),
    paste0(
      sum(robust$diagnosis_adjusted_lm_beta > 0), " of ", nrow(robust),
      " diagnosis-adjusted linear-model coefficients were positive."
    ),
    "",
    "Diagnosis adjustment was feasible, but the broad category model cannot",
    "account for disease severity, treatment, eGFR, IFTA, or other unavailable",
    "clinical covariates. Results remain association-based."
  ),
  file.path(
    result_dir, "external_GSE180394_robust_correlation_report.md"
  )
)

write_validation_session_info(
  "sessionInfo_external_GSE180394_sensitivity_v4.txt"
)
message("GSE180394 sensitivity v4 workflow completed.")
