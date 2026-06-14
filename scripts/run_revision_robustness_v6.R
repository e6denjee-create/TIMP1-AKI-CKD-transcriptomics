#!/usr/bin/env Rscript

# Second-round reviewer analyses. These outputs preserve the v5 files and add
# scoring-consistent module benchmarks and leave-one-signature-out injury PCs.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260614)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}

result_dir <- validation_result_dir()
external_expression <- as.matrix(read.csv(
  gzfile(file.path(
    "data", "external", "GSE180394",
    "external_GSE180394_expression_gene_symbol.csv.gz"
  )),
  row.names = 1, check.names = FALSE
))
storage.mode(external_expression) <- "numeric"
external_metadata <- read.csv(
  file.path(
    "data", "external", "GSE180394",
    "external_GSE180394_metadata.csv"
  ),
  stringsAsFactors = FALSE
)
external_metadata <- external_metadata[
  match(colnames(external_expression), external_metadata$sample), ,
  drop = FALSE
]
disease_index <- external_metadata$group == "Disease"

score_gene_set <- function(expression, genes) {
  available <- intersect(toupper(genes), rownames(expression))
  if (length(available) < 2) return(rep(NA_real_, ncol(expression)))
  scaled <- t(scale(t(expression[available, , drop = FALSE])))
  scaled[!is.finite(scaled)] <- NA_real_
  colMeans(scaled, na.rm = TRUE)
}

stringent <- read.csv(
  file.path(result_dir, "stringent_TIMP1_correlated_module.csv"),
  stringsAsFactors = FALSE
)
module_genes <- intersect(toupper(stringent$gene), rownames(external_expression))
module_score <- score_gene_set(external_expression, module_genes)
timp1 <- as.numeric(external_expression["TIMP1", ])
observed_test <- suppressWarnings(cor.test(
  timp1[disease_index], module_score[disease_index],
  method = "spearman", exact = FALSE
))

genomewide <- read.csv(
  file.path(result_dir, "TIMP1_genomewide_disease_only_correlations.csv"),
  stringsAsFactors = FALSE
)
datasets <- c("GSE139061", "GSE30718", "GSE66494")
loco_rows <- list()
for (left_out in datasets) {
  retained <- setdiff(datasets, left_out)
  wide <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), lapply(
    retained,
    function(dataset_name) {
      x <- genomewide[
        genomewide$dataset == dataset_name,
        c("gene", "rho", "bh_adjusted_p_value")
      ]
      colnames(x)[-1] <- paste0(c("rho_", "fdr_"), dataset_name)
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
  available <- intersect(toupper(selected), rownames(external_expression))
  score <- score_gene_set(external_expression, available)
  test <- suppressWarnings(cor.test(
    timp1[disease_index], score[disease_index],
    method = "spearman", exact = FALSE
  ))
  loco_rows[[left_out]] <- data.frame(
    left_out_discovery_cohort = left_out,
    retained_discovery_cohorts = paste(retained, collapse = ";"),
    selected_gene_count = length(unique(selected)),
    external_available_gene_count = length(available),
    overlap_with_full_module = length(intersect(selected, stringent$gene)),
    jaccard_with_full_module = length(intersect(selected, stringent$gene)) /
      length(union(selected, stringent$gene)),
    external_disease_only_rho = unname(test$estimate),
    p_value = test$p.value,
    scoring_standardization_scope = "All 53 GSE180394 samples",
    stringsAsFactors = FALSE
  )
}
loco <- do.call(rbind, loco_rows)
loco$bh_adjusted_p_value <- p.adjust(loco$p_value, "BH")
write_validation_csv(loco, "module_leave_one_discovery_out_v6.csv")

gene_mean <- rowMeans(external_expression, na.rm = TRUE)
gene_sd <- apply(external_expression, 1, sd, na.rm = TRUE)
valid_genes <- names(gene_mean)[
  is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
]
module_genes <- intersect(module_genes, valid_genes)
candidate_genes <- setdiff(valid_genes, c("TIMP1", module_genes))
mean_breaks <- unique(quantile(
  gene_mean[valid_genes], seq(0, 1, 0.1), na.rm = TRUE
))
sd_breaks <- unique(quantile(
  gene_sd[valid_genes], seq(0, 1, 0.1), na.rm = TRUE
))
stratum <- interaction(
  cut(gene_mean, mean_breaks, include.lowest = TRUE),
  cut(gene_sd, sd_breaks, include.lowest = TRUE),
  drop = TRUE
)
names(stratum) <- names(gene_mean)
scaled_all <- t(scale(t(external_expression[valid_genes, , drop = FALSE])))
scaled_all[!is.finite(scaled_all)] <- 0
module_stratum_counts <- table(stratum[module_genes])
candidate_by_stratum <- split(candidate_genes, stratum[candidate_genes])

random_rho <- numeric(1000L)
for (iteration in seq_along(random_rho)) {
  sampled <- character()
  for (stratum_name in names(module_stratum_counts)) {
    needed <- unname(module_stratum_counts[[stratum_name]])
    pool <- setdiff(candidate_by_stratum[[stratum_name]], sampled)
    if (length(pool) >= needed) {
      sampled <- c(sampled, sample(pool, needed, replace = FALSE))
    } else {
      sampled <- c(sampled, pool)
      fallback <- setdiff(candidate_genes, sampled)
      sampled <- c(sampled, sample(fallback, needed - length(pool), FALSE))
    }
  }
  random_score <- colMeans(scaled_all[sampled, , drop = FALSE])
  random_rho[[iteration]] <- suppressWarnings(cor(
    timp1[disease_index], random_score[disease_index], method = "spearman"
  ))
}
benchmark <- data.frame(
  dataset = "GSE180394",
  scope = "disease_only",
  observed_module_gene_count = length(module_genes),
  observed_rho = unname(observed_test$estimate),
  observed_p_value = observed_test$p.value,
  random_repetitions = length(random_rho),
  random_rho_mean = mean(random_rho),
  random_rho_sd = sd(random_rho),
  random_rho_95th_percentile = unname(quantile(random_rho, 0.95)),
  empirical_one_sided_p = (1 + sum(random_rho >= observed_test$estimate)) /
    (1 + length(random_rho)),
  matching_variables = paste(
    "Mean-expression and variability deciles calculated across all",
    "53 GSE180394 samples"
  ),
  scoring_standardization_scope = "All 53 GSE180394 samples",
  stringsAsFactors = FALSE
)
write_validation_csv(benchmark, "module_random_matched_benchmark_v6.csv")
write_validation_csv(
  data.frame(iteration = seq_along(random_rho), random_rho = random_rho),
  "module_random_matched_null_distribution_v6.csv"
)

plot <- ggplot2::ggplot(
  data.frame(random_rho = random_rho),
  ggplot2::aes(random_rho)
) +
  ggplot2::geom_histogram(
    bins = 35, fill = "#8FB3D9", colour = "white"
  ) +
  ggplot2::geom_vline(
    xintercept = unname(observed_test$estimate),
    colour = "#B22222", linewidth = 1
  ) +
  ggplot2::annotate(
    "text", x = unname(observed_test$estimate), y = Inf,
    label = sprintf(
      "Observed rho = %.3f\nEmpirical P = %.3f",
      observed_test$estimate, benchmark$empirical_one_sided_p
    ),
    hjust = 1.05, vjust = 1.3, colour = "#B22222"
  ) +
  ggplot2::labs(
    title = "Matched random-program benchmark in GSE180394",
    subtitle = "All scores standardized across the same 53-sample cohort",
    x = "Disease-only Spearman rho", y = "Random programs"
  ) +
  validation_theme()
save_validation_plot(
  plot, "v6_module_random_matched_benchmark",
  width = 7.2, height = 5.4
)

signature_table <- read.csv(
  file.path(result_dir, "signature_gene_sets_used.csv"),
  stringsAsFactors = FALSE
)
signature_names <- unique(signature_table$signature)
signature_list <- split(signature_table$gene, signature_table$signature)
score_discovery <- read.csv(
  file.path(result_dir, "TIMP1_signature_scores.csv"),
  stringsAsFactors = FALSE
)
score_external <- read.csv(
  file.path(result_dir, "external_GSE180394_signature_scores.csv"),
  stringsAsFactors = FALSE
)
score_external$TIMP1_expression <- timp1[
  match(score_external$sample, external_metadata$sample)
]
all_scores <- rbind(
  score_discovery[, c(
    "dataset", "sample", "group", "TIMP1_expression", signature_names
  )],
  score_external[, c(
    "dataset", "sample", "group", "TIMP1_expression", signature_names
  )]
)

loso_rows <- list()
loading_rows <- list()
index <- 1L
loading_index <- 1L
for (dataset_name in unique(all_scores$dataset)) {
  data <- all_scores[
    all_scores$dataset == dataset_name & all_scores$group != "Control", ,
    drop = FALSE
  ]
  for (signature in signature_names) {
    background_names <- setdiff(signature_names, signature)
    background <- scale(data[, background_names, drop = FALSE])
    background[!is.finite(background)] <- 0
    pc <- prcomp(background, center = FALSE, scale. = FALSE)
    pc1 <- pc$x[, 1]
    if (cor(pc1, rowMeans(background), use = "pairwise.complete.obs") < 0) {
      pc1 <- -pc1
      pc$rotation[, 1] <- -pc$rotation[, 1]
    }
    timp1_residual <- residuals(lm(data$TIMP1_expression ~ pc1))
    signature_residual <- residuals(lm(data[[signature]] ~ pc1))
    test <- suppressWarnings(cor.test(
      timp1_residual, signature_residual,
      method = "spearman", exact = FALSE
    ))
    loso_rows[[index]] <- data.frame(
      dataset = dataset_name,
      signature = signature,
      n = nrow(data),
      background_signature_count = length(background_names),
      pc1_variance_explained = pc$sdev[[1]]^2 / sum(pc$sdev^2),
      unadjusted_rho = cor(
        data$TIMP1_expression, data[[signature]], method = "spearman"
      ),
      leave_one_signature_out_pc1_adjusted_rho = unname(test$estimate),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
    for (background_signature in background_names) {
      loading_rows[[loading_index]] <- data.frame(
        dataset = dataset_name,
        focal_signature_excluded = signature,
        background_signature = background_signature,
        pc1_loading = pc$rotation[background_signature, 1],
        stringsAsFactors = FALSE
      )
      loading_index <- loading_index + 1L
    }
  }
}
loso <- do.call(rbind, loso_rows)
loso$bh_adjusted_p_value <- ave(
  loso$p_value, loso$dataset,
  FUN = function(x) p.adjust(x, "BH")
)
write_validation_csv(
  loso, "signature_leave_one_out_injury_pc_adjustment_v6.csv"
)
write_validation_csv(
  do.call(rbind, loading_rows),
  "signature_leave_one_out_injury_pc_loadings_v6.csv"
)

writeLines(
  c(
    "# Second-round robustness interpretation (v6)",
    "",
    sprintf(
      "- The scoring-consistent external program correlation was %.3f.",
      benchmark$observed_rho
    ),
    sprintf(
      "- Against 1,000 matched random programs, empirical P = %.4f.",
      benchmark$empirical_one_sided_p
    ),
    "- Leave-one-discovery-cohort scores remain correlated externally, but low Jaccard overlap in two iterations shows that gene membership is not compositionally conserved.",
    "- Leave-one-signature-out injury PCs avoid using the focal signature to define its own adjustment factor; variance explained and loadings are reported for every analysis.",
    "- These analyses support a reproducible injury-state score, not a fixed conserved gene module or a TIMP1-specific causal program."
  ),
  file.path(result_dir, "revision_robustness_v6_interpretation.md")
)
write_validation_session_info("sessionInfo_revision_robustness_v6.txt")
message("Revision robustness v6 analyses completed.")
