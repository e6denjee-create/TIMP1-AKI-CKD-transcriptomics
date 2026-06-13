#!/usr/bin/env Rscript

# Comprehensive enrichment attempt for TIMP1-associated module variants.
# Package installation, collection retrieval, and failed analyses are recorded
# explicitly so unavailable resources are never represented as completed work.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260612)
options(timeout = 900)

package_names <- c(
  "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "msigdbr",
  "AnnotationDbi", "ggplot2"
)
install_requested <- tolower(Sys.getenv(
  "TIMP1_INSTALL_ENRICHMENT_PACKAGES", "true"
)) %in% c("1", "true", "yes")

record_failure <- function(component, reason, action) {
  append_missing_data(
    "comprehensive_enrichment_v3", "TIMP1_modules", component,
    reason, action
  )
}

install_missing_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (!length(missing) || !install_requested) return(invisible(NULL))

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    tryCatch(
      install.packages("BiocManager", repos = "https://cloud.r-project.org"),
      error = function(error) {
        record_failure(
          "BiocManager installation", conditionMessage(error),
          "Package installation failed; dependent enrichments were skipped."
        )
      }
    )
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) return(invisible(NULL))

  tryCatch(
    BiocManager::install(missing, ask = FALSE, update = FALSE),
    error = function(error) {
      record_failure(
        paste("package installation", paste(missing, collapse = "; ")),
        conditionMessage(error),
        "Installation failure recorded; available methods continued."
      )
    }
  )
  invisible(NULL)
}

install_missing_packages(package_names)

status_rows <- list()
status_index <- 1
add_status <- function(component, type, available, status) {
  status_rows[[status_index]] <<- data.frame(
    component = component, type = type, available = available,
    status = status, stringsAsFactors = FALSE
  )
  status_index <<- status_index + 1
}

for (package in package_names) {
  available <- requireNamespace(package, quietly = TRUE)
  add_status(
    package, "R_package", available,
    if (available) {
      paste0("installed version ", utils::packageVersion(package))
    } else {
      "not installed after installation attempt"
    }
  )
  if (!available) {
    record_failure(
      package, "Package is unavailable after the installation attempt.",
      "Analyses requiring this package were skipped."
    )
  }
}

membership_path <- file.path(
  validation_result_dir(), "module_sensitivity_membership.csv"
)
correlation_path <- file.path(
  validation_result_dir(),
  "TIMP1_genomewide_disease_only_correlations.csv"
)
if (!file.exists(membership_path) || !file.exists(correlation_path)) {
  record_failure(
    "module inputs",
    "module_sensitivity_membership.csv or genome-wide correlations are absent.",
    "Comprehensive enrichment was skipped."
  )
  quit(save = "no", status = 0)
}

membership <- read.csv(
  membership_path, check.names = FALSE, stringsAsFactors = FALSE
)
correlations <- read.csv(
  correlation_path, check.names = FALSE, stringsAsFactors = FALSE
)
module_names <- c(
  "relaxed", "stringent",
  "relaxed_without_signature_genes",
  "stringent_without_signature_genes",
  "relaxed_without_ECM_fibrosis_markers",
  "stringent_without_ECM_fibrosis_markers"
)
modules <- split(membership$gene, membership$analysis)
modules <- modules[intersect(module_names, names(modules))]
symbol_universe <- Reduce(
  intersect, split(correlations$gene, correlations$dataset)
)

empty_results <- data.frame(
  module = character(), database = character(), ID = character(),
  Description = character(), GeneRatio = character(), BgRatio = character(),
  pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
  geneID = character(), Count = integer(), stringsAsFactors = FALSE
)
result_rows <- list()
result_index <- 1

standardize_result <- function(result, module, database) {
  if (is.null(result)) return(NULL)
  result <- as.data.frame(result)
  if (!nrow(result)) return(NULL)
  required <- colnames(empty_results)[-c(1, 2)]
  for (column in setdiff(required, colnames(result))) result[[column]] <- NA
  output <- result[, required, drop = FALSE]
  output$module <- module
  output$database <- database
  output[, colnames(empty_results), drop = FALSE]
}

append_result <- function(result, module, database) {
  standardized <- standardize_result(result, module, database)
  if (!is.null(standardized)) {
    result_rows[[result_index]] <<- standardized
    result_index <<- result_index + 1
  }
}

symbol_to_entrez <- NULL
if (
  requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  symbol_to_entrez <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db, keys = unique(symbol_universe),
      columns = "ENTREZID", keytype = "SYMBOL"
    ),
    error = function(error) {
      record_failure(
        "SYMBOL-to-ENTREZ mapping", conditionMessage(error),
        "GO, Reactome, and KEGG analyses were skipped."
      )
      NULL
    }
  )
  if (!is.null(symbol_to_entrez)) {
    symbol_to_entrez <- symbol_to_entrez[
      !is.na(symbol_to_entrez$ENTREZID) &
        !duplicated(symbol_to_entrez$SYMBOL),
      , drop = FALSE
    ]
  }
}

msig_specs <- list(
  Hallmark = c("H", NA_character_),
  GO_BP = c("C5", "GO:BP"),
  GO_CC = c("C5", "GO:CC"),
  Reactome = c("C2", "CP:REACTOME"),
  KEGG = c("C2", "CP:KEGG_MEDICUS")
)
msig_collections <- list()
if (requireNamespace("msigdbr", quietly = TRUE)) {
  cache_dir <- file.path(timp1_project_dir(), ".cache", "R")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(R_USER_CACHE_DIR = cache_dir)
  for (database in names(msig_specs)) {
    spec <- msig_specs[[database]]
    collection <- tryCatch(
      {
        if (is.na(spec[[2]])) {
          msigdbr::msigdbr(
            species = "Homo sapiens", collection = spec[[1]]
          )
        } else {
          msigdbr::msigdbr(
            species = "Homo sapiens", collection = spec[[1]],
            subcollection = spec[[2]]
          )
        }
      },
      error = function(error) {
        record_failure(
          paste("msigdbr", database), conditionMessage(error),
          "Collection retrieval failed; no result was claimed."
        )
        NULL
      }
    )
    available <- is.data.frame(collection) && nrow(collection) > 0
    add_status(
      paste0("msigdbr_", database), "gene_set_collection", available,
      if (available) {
        paste("retrieved", nrow(collection), "membership rows")
      } else {
        "collection unavailable"
      }
    )
    if (available) msig_collections[[database]] <- collection
  }
}

run_attempt <- function(module, database, expression) {
  tryCatch(
    expression,
    error = function(error) {
      record_failure(
        paste(module, database, sep = "::"), conditionMessage(error),
        "This module-database analysis was skipped."
      )
      add_status(
        paste(module, database, sep = "::"), "analysis", FALSE,
        conditionMessage(error)
      )
      NULL
    }
  )
}

for (module in names(modules)) {
  symbols <- intersect(unique(modules[[module]]), symbol_universe)

  if (
    requireNamespace("clusterProfiler", quietly = TRUE) &&
      "Hallmark" %in% names(msig_collections)
  ) {
    hallmark <- msig_collections$Hallmark
    hallmark_result <- run_attempt(
      module, "Hallmark",
      clusterProfiler::enricher(
        symbols, universe = symbol_universe,
        TERM2GENE = hallmark[, c("gs_name", "gene_symbol")],
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
      )
    )
    append_result(hallmark_result, module, "Hallmark")
    add_status(
      paste(module, "Hallmark", sep = "::"), "analysis",
      !is.null(hallmark_result),
      if (is.null(hallmark_result)) {
        "analysis failed; see missing_data_log.csv"
      } else {
        paste("completed;", nrow(as.data.frame(hallmark_result)), "terms")
      }
    )
  } else {
    add_status(
      paste(module, "Hallmark", sep = "::"), "analysis", FALSE,
      "clusterProfiler or Hallmark collection unavailable"
    )
  }

  entrez <- if (!is.null(symbol_to_entrez)) {
    unique(symbol_to_entrez$ENTREZID[
      symbol_to_entrez$SYMBOL %in% symbols
    ])
  } else {
    character()
  }
  universe_entrez <- if (!is.null(symbol_to_entrez)) {
    unique(symbol_to_entrez$ENTREZID)
  } else {
    character()
  }

  for (ontology in c(BP = "GO_BP", CC = "GO_CC")) {
    database <- unname(ontology)
    if (
      length(entrez) &&
        requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
    ) {
      go_result <- run_attempt(
        module, database,
        clusterProfiler::enrichGO(
          gene = entrez, universe = universe_entrez,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = "ENTREZID", ont = names(ontology),
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
          readable = TRUE
        )
      )
      append_result(go_result, module, database)
      add_status(
        paste(module, database, sep = "::"), "analysis",
        !is.null(go_result),
        if (is.null(go_result)) {
          "analysis failed; see missing_data_log.csv"
        } else {
          paste("completed;", nrow(as.data.frame(go_result)), "terms")
        }
      )
    } else {
      add_status(
        paste(module, database, sep = "::"), "analysis", FALSE,
        "clusterProfiler, org.Hs.eg.db, or ENTREZ mapping unavailable"
      )
    }
  }

  if (
    length(entrez) && requireNamespace("ReactomePA", quietly = TRUE)
  ) {
    reactome_result <- run_attempt(
      module, "Reactome",
      ReactomePA::enrichPathway(
        gene = entrez, universe = universe_entrez, organism = "human",
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        readable = TRUE
      )
    )
    append_result(reactome_result, module, "Reactome")
    add_status(
      paste(module, "Reactome", sep = "::"), "analysis",
      !is.null(reactome_result),
      if (is.null(reactome_result)) {
        "analysis failed; see missing_data_log.csv"
      } else {
        paste("completed;", nrow(as.data.frame(reactome_result)), "terms")
      }
    )
  } else {
    add_status(
      paste(module, "Reactome", sep = "::"), "analysis", FALSE,
      "ReactomePA or ENTREZ mapping unavailable"
    )
  }

  if (
    length(entrez) && requireNamespace("clusterProfiler", quietly = TRUE)
  ) {
    kegg_result <- run_attempt(
      module, "KEGG",
      clusterProfiler::enrichKEGG(
        gene = entrez, universe = universe_entrez, organism = "hsa",
        keyType = "ncbi-geneid", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, use_internal_data = FALSE
      )
    )
    append_result(kegg_result, module, "KEGG")
    kegg_available <- !is.null(kegg_result)
    add_status(
      paste(module, "KEGG", sep = "::"), "analysis", kegg_available,
      if (kegg_available) {
        paste("completed;", nrow(as.data.frame(kegg_result)), "terms")
      } else {
        "KEGG analysis unavailable"
      }
    )
  } else {
    add_status(
      paste(module, "KEGG", sep = "::"), "analysis", FALSE,
      "clusterProfiler or ENTREZ mapping unavailable"
    )
  }
}

results <- if (length(result_rows)) {
  do.call(rbind, result_rows)
} else {
  empty_results
}
status <- do.call(rbind, status_rows)
write_validation_csv(results, "comprehensive_enrichment_v3_results.csv")
write_validation_csv(status, "comprehensive_enrichment_v3_status.csv")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  record_failure(
    "ggplot2", "ggplot2 is unavailable.",
    "Enrichment figures could not be generated."
  )
} else {
  significant <- results[
    is.finite(results$p.adjust) & results$p.adjust < 0.05,
    , drop = FALSE
  ]
  if (nrow(significant)) {
    significant <- do.call(rbind, lapply(
      split(significant, paste(significant$module, significant$database)),
      function(x) head(x[order(x$p.adjust), , drop = FALSE], 4)
    ))
    significant$Description <- factor(
      significant$Description,
      levels = rev(unique(significant$Description))
    )
    enrichment_plot <- ggplot2::ggplot(
      significant,
      ggplot2::aes(
        x = Count, y = Description, size = Count,
        color = -log10(pmax(p.adjust, 1e-300))
      )
    ) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::facet_grid(database ~ module, scales = "free_y", space = "free") +
      ggplot2::scale_color_gradient(low = "#4C78A8", high = "#C44E52") +
      ggplot2::labs(
        title = "Comprehensive enrichment of TIMP1-associated modules",
        x = "Overlap genes", y = NULL, color = "-log10(FDR)"
      ) +
      validation_theme(base_size = 9)
  } else {
    enrichment_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = paste(
          "No comprehensive enrichment result passed FDR < 0.05.",
          "See status CSV and missing_data_log.csv."
        )
      ) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::labs(
        title = "Comprehensive enrichment status"
      ) +
      validation_theme()
  }
  save_validation_plot(
    enrichment_plot, "TIMP1_comprehensive_enrichment_v3",
    width = 14, height = 10
  )

  status_summary <- aggregate(
    available ~ type, status, function(x) sum(as.logical(x), na.rm = TRUE)
  )
  totals <- aggregate(component ~ type, status, length)
  status_summary$total <- totals$component[
    match(status_summary$type, totals$type)
  ]
  status_summary$unavailable <- status_summary$total - status_summary$available
  status_long <- rbind(
    data.frame(
      type = status_summary$type, outcome = "Available",
      count = status_summary$available
    ),
    data.frame(
      type = status_summary$type, outcome = "Unavailable",
      count = status_summary$unavailable
    )
  )
  status_plot <- ggplot2::ggplot(
    status_long, ggplot2::aes(type, count, fill = outcome)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c(Available = "#4C78A8", Unavailable = "#C44E52")
    ) +
    ggplot2::labs(
      title = "Comprehensive enrichment execution status",
      x = NULL, y = "Components", fill = NULL
    ) +
    validation_theme()
  save_validation_plot(
    status_plot, "TIMP1_comprehensive_enrichment_v3_status",
    width = 8, height = 5.5
  )
}

write_validation_session_info(
  "sessionInfo_comprehensive_enrichment_v3.txt"
)
message(
  "Comprehensive enrichment v3 finished with ", nrow(results),
  " result rows. Review the status table for unavailable components."
)
