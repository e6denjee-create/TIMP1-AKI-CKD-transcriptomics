#!/usr/bin/env Rscript

# Authoritative effect-size confidence intervals for manuscript version 7.
# Each contrast uses stratified, within-group nonparametric bootstrap sampling.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()

result_dir <- validation_result_dir()
repetitions <- 2000L
base_seed <- 20260614L

hedges_g <- function(disease, control) {
  pooled_sd <- sqrt(
    ((length(disease) - 1) * var(disease) +
      (length(control) - 1) * var(control)) /
      (length(disease) + length(control) - 2)
  )
  correction <- 1 - 3 / (4 * (length(disease) + length(control)) - 9)
  ((mean(disease) - mean(control)) / pooled_sd) * correction
}

bootstrap_contrast <- function(
    dataset, contrast, disease, control, seed, p_value, p_adjustment_family) {
  set.seed(seed)
  estimates <- replicate(
    repetitions,
    hedges_g(
      sample(disease, length(disease), replace = TRUE),
      sample(control, length(control), replace = TRUE)
    )
  )
  data.frame(
    dataset = dataset,
    contrast = contrast,
    n_disease = length(disease),
    n_control = length(control),
    disease_mean = mean(disease),
    control_mean = mean(control),
    group_difference = mean(disease) - mean(control),
    hedges_g = hedges_g(disease, control),
    hedges_g_ci_lower = unname(quantile(estimates, 0.025, na.rm = TRUE)),
    hedges_g_ci_upper = unname(quantile(estimates, 0.975, na.rm = TRUE)),
    p_value = p_value,
    p_adjustment_family = p_adjustment_family,
    bootstrap_method = paste(
      "Stratified nonparametric percentile bootstrap;",
      "disease and control groups resampled independently"
    ),
    bootstrap_repetitions = repetitions,
    random_seed = seed,
    stringsAsFactors = FALSE
  )
}

bulk_source <- read.csv(
  file.path(result_dir, "TIMP1_bulk_validation_source_data.csv"),
  stringsAsFactors = FALSE
)
bulk_stats <- read.csv(
  file.path(result_dir, "TIMP1_bulk_validation_statistics.csv"),
  stringsAsFactors = FALSE
)

rows <- list()
index <- 1L
for (dataset_name in c("GSE139061", "GSE30718", "GSE66494")) {
  data <- bulk_source[bulk_source$dataset == dataset_name, , drop = FALSE]
  stats <- bulk_stats[bulk_stats$dataset == dataset_name, , drop = FALSE]
  disease_values <- data$TIMP1_expression[data$group != "Control"]
  rows[[index]] <- bootstrap_contrast(
    dataset_name,
    stats$contrast[[1]],
    disease_values,
    data$TIMP1_expression[data$group == "Control"],
    base_seed + index,
    stats$p_value[[1]],
    "BH across the three discovery-cohort expression comparisons"
  )
  index <- index + 1L
}

external_source <- read.csv(
  file.path(
    result_dir,
    "external_GSE180394_sensitivity_control_source_data.csv"
  ),
  stringsAsFactors = FALSE
)
external_stats <- read.csv(
  file.path(
    result_dir,
    "external_GSE180394_sensitivity_control_statistics.csv"
  ),
  stringsAsFactors = FALSE
)
contrast_labels <- c(
  living_donor_only = "Disease_vs_Living_Donor",
  extended_controls = "Disease_vs_Living_Donor_plus_Tumor_Nephrectomy",
  tumor_nephrectomy_only = "Disease_vs_Tumor_Nephrectomy"
)
for (definition in names(contrast_labels)) {
  data <- external_source[
    external_source$control_definition == definition, ,
    drop = FALSE
  ]
  stats <- external_stats[
    external_stats$control_definition == definition, ,
    drop = FALSE
  ]
  rows[[index]] <- bootstrap_contrast(
    "GSE180394",
    contrast_labels[[definition]],
    data$TIMP1_expression[data$group == "Disease"],
    data$TIMP1_expression[data$group == "Control"],
    base_seed + index,
    stats$p_value[[1]],
    if (definition == "living_donor_only") {
      "Prespecified single primary external expression comparison; unadjusted"
    } else {
      "Exploratory alternative-control sensitivity family"
    }
  )
  index <- index + 1L
}

unified <- do.call(rbind, rows)
discovery <- unified$dataset != "GSE180394"
unified$bh_adjusted_p_value <- NA_real_
unified$bh_adjusted_p_value[discovery] <- p.adjust(
  unified$p_value[discovery], method = "BH"
)
sensitivity <- unified$dataset == "GSE180394" &
  unified$contrast != "Disease_vs_Living_Donor"
unified$bh_adjusted_p_value[sensitivity] <- p.adjust(
  unified$p_value[sensitivity], method = "BH"
)
write_validation_csv(
  unified,
  "TIMP1_expression_statistics_unified_bootstrap_v7.csv"
)

policy <- data.frame(
  dataset = "GSE210622",
  current_manuscript_status = "Excluded from all inferential and exploratory results",
  reason = paste(
    "Only one donor was locally reconstructed; selective use would not reproduce",
    "the original multi-donor study."
  ),
  treatment_of_local_single_donor_files = paste(
    "Retained only as historical audit material outside the submission evidence;",
    "not summarized, plotted, tested, or cited as a study result."
  ),
  manuscript_scope = paste(
    "The dataset is mentioned only to document exclusion and as a future",
    "multi-donor validation opportunity."
  ),
  supersedes = "Earlier exploratory single-donor audit wording dated 12 June 2026",
  effective_date = "2026-06-14",
  stringsAsFactors = FALSE
)
write_validation_csv(policy, "singlecell_exclusion_policy_v7.csv")

writeLines(
  c(
    "# Unified expression bootstrap and single-cell exclusion policy (v7)",
    "",
    "- All Hedges' g intervals use stratified within-group nonparametric percentile bootstrapping.",
    "- Each contrast uses 2,000 resamples and a recorded deterministic seed.",
    "- The unified CSV is the sole source for expression effect sizes and confidence intervals in the manuscript and Supplementary Table S2.",
    "- GSE210622 is excluded from all inferential and exploratory manuscript results because only one donor was locally reconstructed.",
    "- Historical single-donor audit files are not submission evidence and are superseded by the v7 exclusion policy."
  ),
  file.path(result_dir, "unified_bootstrap_and_singlecell_policy_v7.md")
)
write_validation_session_info("sessionInfo_unified_bootstrap_v7.txt")
