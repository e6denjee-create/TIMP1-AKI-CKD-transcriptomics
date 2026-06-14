#!/usr/bin/env Rscript

# Build a cross-dataset TIMP1-associated correlation module from disease-only
# bulk samples. The resulting network represents association, not causality.

source(file.path("scripts", "_timp1_validation_utils.R"))
activate_timp1_library()
ensure_validation_dirs()
ensure_missing_log()
set.seed(20260612)

required_packages <- c("ggplot2", "igraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

datasets <- c("GSE139061", "GSE30718", "GSE66494")
cache_dir <- file.path(validation_root(), ".cache", "R")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = cache_dir)

vectorized_spearman <- function(expression, target) {
  target_rank <- rank(target, ties.method = "average")
  gene_ranks <- t(apply(
    expression, 1, rank, ties.method = "average"
  ))
  rho <- as.numeric(stats::cor(t(gene_ranks), target_rank))
  names(rho) <- rownames(expression)
  n <- length(target)
  rho_for_p <- pmin(pmax(rho, -0.999999), 0.999999)
  statistic <- rho_for_p * sqrt((n - 2) / (1 - rho_for_p^2))
  p_value <- 2 * stats::pt(-abs(statistic), df = n - 2)
  p_value[!is.finite(rho)] <- NA_real_
  data.frame(
    gene = rownames(expression),
    rho = rho,
    p_value = p_value,
    n = n,
    stringsAsFactors = FALSE
  )
}

correlation_tables <- list()
positive_sets <- list()

correlation_cache <- file.path(
  validation_result_dir(),
  "TIMP1_genomewide_disease_only_correlations.csv"
)
if (file.exists(correlation_cache)) {
  all_correlations <- read.csv(
    correlation_cache, check.names = FALSE, stringsAsFactors = FALSE
  )
  correlation_tables <- split(all_correlations, all_correlations$dataset)
  positive_sets <- lapply(correlation_tables, function(result) {
    result$gene[is.finite(result$rho) & result$rho > 0]
  })
} else {
  for (dataset_name in datasets) {
    dataset <- read_validation_bulk(dataset_name)
    if (is.null(dataset)) next
    disease_metadata <- dataset$metadata[
      dataset$metadata$group != "Control", , drop = FALSE
    ]
    if (nrow(disease_metadata) < 4) {
      append_missing_data(
        "module_correlation", dataset_name, "disease-only samples",
        "Fewer than four disease samples were available."
      )
      next
    }
    if (!"TIMP1" %in% rownames(dataset$expression)) {
      append_missing_data(
        "module_correlation", dataset_name, "TIMP1",
        "TIMP1 is absent from the expression matrix."
      )
      next
    }
    expression <- dataset$expression[, disease_metadata$sample, drop = FALSE]
    variable <- apply(expression, 1, stats::sd, na.rm = TRUE) > 0
    expression <- expression[variable, , drop = FALSE]
    result <- vectorized_spearman(
      expression, as.numeric(expression["TIMP1", ])
    )
    result$bh_adjusted_p_value <- p.adjust(result$p_value, "BH")
    result$dataset <- dataset_name
    result <- result[result$gene != "TIMP1", , drop = FALSE]
    correlation_tables[[dataset_name]] <- result
    positive_sets[[dataset_name]] <- result$gene[
      is.finite(result$rho) & result$rho > 0
    ]
  }
  all_correlations <- do.call(rbind, correlation_tables)
  write_validation_csv(
    all_correlations, "TIMP1_genomewide_disease_only_correlations.csv"
  )
}

correlation_wide <- reshape(
  all_correlations[, c(
    "gene", "dataset", "rho", "p_value", "bh_adjusted_p_value"
  )],
  idvar = "gene", timevar = "dataset", direction = "wide"
)

rho_columns <- paste0("rho.", datasets)
p_columns <- paste0("p_value.", datasets)
fdr_columns <- paste0("bh_adjusted_p_value.", datasets)
for (column in c(rho_columns, p_columns, fdr_columns)) {
  if (!column %in% colnames(correlation_wide)) {
    correlation_wide[[column]] <- NA_real_
  }
}
correlation_wide$positive_dataset_count <- rowSums(
  correlation_wide[, rho_columns, drop = FALSE] > 0, na.rm = TRUE
)
correlation_wide$nominal_positive_dataset_count <- rowSums(
  correlation_wide[, rho_columns, drop = FALSE] > 0 &
    correlation_wide[, p_columns, drop = FALSE] < 0.05,
  na.rm = TRUE
)
correlation_wide$fdr_positive_dataset_count <- rowSums(
  correlation_wide[, rho_columns, drop = FALSE] > 0 &
    correlation_wide[, fdr_columns, drop = FALSE] < 0.05,
  na.rm = TRUE
)
correlation_wide$mean_rho <- rowMeans(
  correlation_wide[, rho_columns, drop = FALSE], na.rm = TRUE
)
correlation_wide$median_rho <- apply(
  correlation_wide[, rho_columns, drop = FALSE], 1,
  stats::median, na.rm = TRUE
)

# Directional replication in at least two datasets is mandatory. Requiring
# nominal support in at least two datasets prevents near-zero positive
# correlations from defining the associated module.
stable <- correlation_wide[
  correlation_wide$positive_dataset_count >= 2 &
    correlation_wide$nominal_positive_dataset_count >= 2 &
    correlation_wide$median_rho > 0,
  , drop = FALSE
]
stable <- stable[order(
  -stable$fdr_positive_dataset_count,
  -stable$nominal_positive_dataset_count,
  -stable$mean_rho
), ]
write_validation_csv(stable, "stable_TIMP1_correlated_genes.csv")

edges <- data.frame(
  from = "TIMP1",
  to = stable$gene,
  weight = stable$mean_rho,
  median_rho = stable$median_rho,
  positive_dataset_count = stable$positive_dataset_count,
  nominal_positive_dataset_count = stable$nominal_positive_dataset_count,
  fdr_positive_dataset_count = stable$fdr_positive_dataset_count,
  relationship = "positive disease-only Spearman association",
  stringsAsFactors = FALSE
)
write_validation_csv(edges, "network_edge_list.csv")

get_msigdb_collection <- function(collection, subcollection = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    append_missing_data(
      "module_enrichment", "MSigDB",
      paste(collection, subcollection, sep = ":"),
      "Optional package msigdbr is unavailable.",
      "Skipped optional enrichment database and continued."
    )
    return(NULL)
  }
  tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(
          species = "Homo sapiens", collection = collection
        )
      } else {
        msigdbr::msigdbr(
          species = "Homo sapiens", collection = collection,
          subcollection = subcollection
        )
      }
    },
    error = function(error) {
      append_missing_data(
        "module_enrichment", "MSigDB",
        paste(collection, subcollection, sep = ":"),
        conditionMessage(error),
        "Skipped unavailable enrichment database and continued."
      )
      NULL
    }
  )
}

enrich_gene_sets <- function(module_genes, universe, msig, database) {
  if (is.null(msig) || !nrow(msig)) return(NULL)
  gene_column <- if ("gene_symbol" %in% colnames(msig)) {
    "gene_symbol"
  } else {
    "gene_symbol"
  }
  gene_sets <- split(msig[[gene_column]], msig$gs_name)
  rows <- lapply(names(gene_sets), function(pathway) {
    pathway_genes <- intersect(unique(gene_sets[[pathway]]), universe)
    overlap <- intersect(module_genes, pathway_genes)
    if (length(pathway_genes) < 5 || length(overlap) < 2) return(NULL)
    p_value <- stats::phyper(
      length(overlap) - 1,
      length(pathway_genes),
      length(universe) - length(pathway_genes),
      length(module_genes),
      lower.tail = FALSE
    )
    data.frame(
      database = database,
      pathway = pathway,
      enrichment_source = "MSigDB_database_collection",
      overlap_count = length(overlap),
      pathway_size = length(pathway_genes),
      module_size = length(module_genes),
      universe_size = length(universe),
      gene_ratio = length(overlap) / length(module_genes),
      fold_enrichment = (
        length(overlap) / length(module_genes)
      ) / (
        length(pathway_genes) / length(universe)
      ),
      p_value = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (is.null(result)) return(NULL)
  result$bh_adjusted_p_value <- p.adjust(result$p_value, "BH")
  result
}

universe <- Reduce(
  intersect, lapply(correlation_tables, function(x) x$gene)
)
module_genes <- intersect(stable$gene, universe)
collections <- list(
  GO_BP = get_msigdb_collection("C5", "GO:BP"),
  Reactome = get_msigdb_collection("C2", "CP:REACTOME"),
  KEGG = get_msigdb_collection("C2", "CP:KEGG_MEDICUS")
)
enrichment_list <- lapply(names(collections), function(database) {
  enrich_gene_sets(
    module_genes, universe, collections[[database]], database
  )
})
enrichment_list <- Filter(Negate(is.null), enrichment_list)
offline_pathways <- list(
  GO_BP = list(
    `GO:0030198 extracellular matrix organization` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL5A2", "COL6A1", "COL6A2", "FN1", "SPARC", "POSTN", "VCAN",
      "LUM", "DCN", "MMP2", "MMP7", "MMP9", "MMP14", "LOX", "LOXL2"
    ),
    `GO:0030199 collagen fibril organization` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "DCN",
      "LUM", "P4HA1", "P4HA2", "PLOD1", "PLOD2", "SERPINH1"
    ),
    `GO:0009611 response to wounding` = c(
      "FN1", "SPP1", "VIM", "KRT8", "KRT18", "KRT19", "VCAM1",
      "ICAM1", "CCL2", "IL6", "STAT3", "JUN", "FOS", "TGM2"
    ),
    `GO:0006954 inflammatory response` = c(
      "IL1B", "IL6", "TNF", "NFKB1", "NFKBIA", "CCL2", "CCL5",
      "CXCL8", "CXCL10", "ICAM1", "VCAM1", "STAT1", "STAT3",
      "TYROBP", "FCER1G", "CTSS"
    ),
    `GO:0090398 cellular senescence` = c(
      "CDKN1A", "CDKN2A", "TP53", "SERPINE1", "GDF15", "IL6",
      "CCL2", "CXCL8", "MMP3", "MMP9"
    )
  ),
  Reactome = list(
    `R-HSA-1474244 Extracellular matrix organization` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL6A1", "COL6A2", "FN1", "SPARC", "VCAN", "LUM", "DCN",
      "MMP2", "MMP7", "MMP9", "MMP14", "TIMP2", "LOX"
    ),
    `R-HSA-1474290 Collagen formation` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL5A2", "COL6A1", "P4HA1", "P4HA2", "PLOD1", "PLOD2",
      "SERPINH1"
    ),
    `R-HSA-1474228 Degradation of the extracellular matrix` = c(
      "MMP2", "MMP3", "MMP7", "MMP9", "MMP14", "TIMP2", "PLAU",
      "PLAUR", "CTSB", "CTSL"
    ),
    `R-HSA-170834 Signaling by TGF-beta receptor complex` = c(
      "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
      "SMAD4", "SMAD7", "CTGF", "SERPINE1", "ACTA2"
    ),
    `R-HSA-2559583 Cellular senescence` = c(
      "CDKN1A", "CDKN2A", "TP53", "SERPINE1", "GDF15", "IL6",
      "CCL2", "MMP3", "MMP9"
    )
  ),
  KEGG = list(
    `hsa04512 ECM-receptor interaction` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL6A1", "COL6A2", "FN1", "SPP1", "LAMA1", "LAMB1",
      "ITGA1", "ITGA2", "ITGA5", "ITGB1", "CD44"
    ),
    `hsa04510 Focal adhesion` = c(
      "COL1A1", "COL1A2", "COL3A1", "FN1", "ITGA1", "ITGA2",
      "ITGA5", "ITGB1", "PTK2", "PXN", "VCL", "ACTN1", "TLN1",
      "PIK3CA", "AKT1", "MAPK1"
    ),
    `hsa04350 TGF-beta signaling pathway` = c(
      "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
      "SMAD4", "SMAD7", "BMP2", "BMP4", "ID1", "ID2", "SERPINE1"
    ),
    `hsa04064 NF-kappa B signaling pathway` = c(
      "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA", "IKBKB", "TNF",
      "TNFRSF1A", "IL1B", "IL1R1", "CCL2", "CXCL8", "ICAM1"
    ),
    `hsa04218 Cellular senescence` = c(
      "CDKN1A", "CDKN2A", "TP53", "RB1", "E2F1", "CCNB1", "IL6",
      "CXCL8", "SERPINE1", "MAPK1", "MAPK14"
    )
  )
)
if (!length(enrichment_list)) {
  enrichment_list <- lapply(names(offline_pathways), function(database) {
    pseudo_msig <- do.call(rbind, lapply(
      names(offline_pathways[[database]]),
      function(pathway) {
        data.frame(
          gs_name = pathway,
          gene_symbol = offline_pathways[[database]][[pathway]],
          stringsAsFactors = FALSE
        )
      }
    ))
    result <- enrich_gene_sets(
      module_genes, universe, pseudo_msig, database
    )
    if (!is.null(result)) {
      result$enrichment_source <- "curated_offline_subset"
    }
    result
  })
  enrichment_list <- Filter(Negate(is.null), enrichment_list)
}
if (length(enrichment_list)) {
  enrichment <- do.call(rbind, enrichment_list)
  enrichment$global_bh_adjusted_p_value <- p.adjust(enrichment$p_value, "BH")
  enrichment <- enrichment[order(
    enrichment$global_bh_adjusted_p_value,
    -enrichment$overlap_count
  ), ]
} else {
  enrichment <- data.frame(
    database = character(), pathway = character(),
    enrichment_source = character(),
    overlap_count = integer(), pathway_size = integer(),
    module_size = integer(), universe_size = integer(),
    gene_ratio = numeric(), fold_enrichment = numeric(), p_value = numeric(),
    overlap_genes = character(), bh_adjusted_p_value = numeric(),
    global_bh_adjusted_p_value = numeric()
  )
}
write_validation_csv(
  enrichment, "enrichment_TIMP1_correlated_module.csv"
)

draw_network <- function(path, type) {
  if (type == "png") {
    grDevices::png(path, width = 9, height = 8, units = "in", res = 600)
  } else {
    grDevices::cairo_pdf(path, width = 9, height = 8)
  }
  top_edges <- head(edges, 40)
  if (nrow(top_edges)) {
    graph <- igraph::graph_from_data_frame(top_edges, directed = FALSE)
    vertex_size <- ifelse(igraph::V(graph)$name == "TIMP1", 24, 9)
    vertex_color <- ifelse(
      igraph::V(graph)$name == "TIMP1", "#C44E52", "#4C78A8"
    )
    edge_width <- 1 + 6 * igraph::E(graph)$weight
    set.seed(20260612)
    plot(
      graph,
      layout = igraph::layout_with_fr(graph),
      vertex.size = vertex_size,
      vertex.color = vertex_color,
      vertex.frame.color = "white",
      vertex.label.cex = 0.72,
      vertex.label.color = "black",
      edge.width = edge_width,
      edge.color = grDevices::adjustcolor("#777777", alpha.f = 0.65),
      main = "TIMP1-associated disease-only correlation module"
    )
    graphics::mtext(
      "Edges represent replicated positive Spearman associations",
      side = 1, line = -1.2, cex = 0.8
    )
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No stable correlated genes met criteria.")
  }
  grDevices::dev.off()
}
network_base <- file.path(
  validation_figure_dir(), "TIMP1_centered_correlation_network"
)
draw_network(paste0(network_base, ".png"), "png")
draw_network(paste0(network_base, ".pdf"), "pdf")

plot_enrichment <- enrichment[
  is.finite(enrichment$global_bh_adjusted_p_value) &
    enrichment$global_bh_adjusted_p_value < 0.25,
  , drop = FALSE
]
if (!nrow(plot_enrichment)) plot_enrichment <- head(enrichment, 15)
plot_enrichment <- head(plot_enrichment, 20)
if (nrow(plot_enrichment)) {
  plot_enrichment$pathway_label <- gsub(
    "^(GOBP_|REACTOME_|KEGG_MEDICUS_)", "", plot_enrichment$pathway
  )
  plot_enrichment$pathway_label <- gsub("_", " ", plot_enrichment$pathway_label)
  plot_enrichment$pathway_label <- paste0(
    "[", plot_enrichment$database, "] ", plot_enrichment$pathway_label
  )
  plot_enrichment$pathway_label <- factor(
    plot_enrichment$pathway_label,
    levels = rev(unique(plot_enrichment$pathway_label))
  )
  enrichment_plot <- ggplot2::ggplot(
    plot_enrichment,
    ggplot2::aes(
      x = fold_enrichment, y = pathway_label,
      size = overlap_count,
      color = -log10(pmax(global_bh_adjusted_p_value, 1e-300))
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(low = "#4C78A8", high = "#C44E52") +
    ggplot2::labs(
      title = "Enrichment of the TIMP1-associated correlation module",
      subtitle = "Curated offline GO/Reactome/KEGG subset",
      x = "Fold enrichment", y = NULL,
      size = "Overlap genes", color = "-log10(global FDR)"
    ) +
    validation_theme()
} else {
  enrichment_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = "No enrichment result was available."
    ) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(
      title = "Enrichment of the TIMP1-associated correlation module"
    ) +
    validation_theme()
}
save_validation_plot(
  enrichment_plot, "TIMP1_module_enrichment_dotplot",
  width = 10, height = 7
)

draw_overlap <- function(path, type) {
  if (type == "png") {
    grDevices::png(path, width = 8, height = 7, units = "in", res = 600)
  } else {
    grDevices::cairo_pdf(path, width = 8, height = 7)
  }
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 10), ylim = c(0, 9))
  colors <- grDevices::adjustcolor(
    c("#4C78A8", "#E45756", "#72B7B2"), alpha.f = 0.35
  )
  symbols(
    c(4.1, 5.9, 5.0), c(5.2, 5.2, 3.8),
    circles = c(2.45, 2.45, 2.45), inches = FALSE,
    add = TRUE, bg = colors, fg = c("#4C78A8", "#E45756", "#3A8D88")
  )
  a <- positive_sets[[datasets[1]]]
  b <- positive_sets[[datasets[2]]]
  cset <- positive_sets[[datasets[3]]]
  only_a <- length(setdiff(a, union(b, cset)))
  only_b <- length(setdiff(b, union(a, cset)))
  only_c <- length(setdiff(cset, union(a, b)))
  ab <- length(setdiff(intersect(a, b), cset))
  ac <- length(setdiff(intersect(a, cset), b))
  bc <- length(setdiff(intersect(b, cset), a))
  abc <- length(Reduce(intersect, list(a, b, cset)))
  graphics::text(2.7, 5.6, only_a, cex = 1.1)
  graphics::text(7.3, 5.6, only_b, cex = 1.1)
  graphics::text(5.0, 2.3, only_c, cex = 1.1)
  graphics::text(5.0, 6.0, ab, cex = 1.1)
  graphics::text(4.0, 4.1, ac, cex = 1.1)
  graphics::text(6.0, 4.1, bc, cex = 1.1)
  graphics::text(5.0, 4.8, abc, cex = 1.2, font = 2)
  graphics::text(2.0, 8.2, datasets[1], font = 2)
  graphics::text(8.0, 8.2, datasets[2], font = 2)
  graphics::text(5.0, 0.6, datasets[3], font = 2)
  graphics::title(
    "Overlap of positively TIMP1-correlated genes",
    sub = "Disease-only Spearman direction; stable module applies stricter criteria"
  )
  grDevices::dev.off()
}
overlap_base <- file.path(
  validation_figure_dir(), "TIMP1_correlated_gene_overlap"
)
draw_overlap(paste0(overlap_base, ".png"), "png")
draw_overlap(paste0(overlap_base, ".pdf"), "pdf")

write_validation_session_info("sessionInfo_module_network.txt")
message(
  "TIMP1-associated module completed with ", nrow(stable),
  " stable genes."
)
