#!/usr/bin/env Rscript

# Publication-oriented sensitivity analyses for the cross-dataset
# TIMP1-associated bulk correlation module.

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

stable_path <- file.path(
  validation_result_dir(), "stable_TIMP1_correlated_genes.csv"
)
correlation_path <- file.path(
  validation_result_dir(),
  "TIMP1_genomewide_disease_only_correlations.csv"
)
signature_path <- file.path(
  validation_result_dir(), "signature_gene_sets_used.csv"
)
if (!all(file.exists(c(stable_path, correlation_path, signature_path)))) {
  stop("Required v1 module or signature outputs are missing.")
}

stable <- read.csv(stable_path, check.names = FALSE, stringsAsFactors = FALSE)
correlations <- read.csv(
  correlation_path, check.names = FALSE, stringsAsFactors = FALSE
)
signature_table <- read.csv(
  signature_path, check.names = FALSE, stringsAsFactors = FALSE
)

relaxed <- stable[
  stable$positive_dataset_count >= 2 &
    stable$nominal_positive_dataset_count >= 2 &
    stable$median_rho > 0,
  , drop = FALSE
]
relaxed$module_definition <- paste(
  "positive direction and nominal P < 0.05 in at least two datasets"
)

stringent <- stable[
  stable$fdr_positive_dataset_count >= 2 &
    stable$median_rho > 0,
  , drop = FALSE
]
stringent$module_definition <- paste(
  "positive direction and BH-adjusted P < 0.05 in at least two datasets"
)

write_validation_csv(relaxed, "relaxed_TIMP1_correlated_module.csv")
write_validation_csv(stringent, "stringent_TIMP1_correlated_module.csv")

signature_genes <- unique(signature_table$gene)
ecm_fibrosis_signatures <- c(
  "ECM_remodeling", "Collagen_formation", "Fibrosis"
)
ecm_fibrosis_genes <- unique(signature_table$gene[
  signature_table$signature %in% ecm_fibrosis_signatures
])

module_variants <- list(
  relaxed = relaxed$gene,
  stringent = stringent$gene,
  relaxed_without_signature_genes = setdiff(relaxed$gene, signature_genes),
  stringent_without_signature_genes = setdiff(stringent$gene, signature_genes),
  relaxed_without_ECM_fibrosis_markers = setdiff(
    relaxed$gene, ecm_fibrosis_genes
  ),
  stringent_without_ECM_fibrosis_markers = setdiff(
    stringent$gene, ecm_fibrosis_genes
  )
)

membership <- do.call(rbind, lapply(names(module_variants), function(name) {
  data.frame(
    analysis = name,
    gene = module_variants[[name]],
    stringsAsFactors = FALSE
  )
}))
write_validation_csv(membership, "module_sensitivity_membership.csv")

module_sizes <- data.frame(
  analysis = names(module_variants),
  gene_count = lengths(module_variants),
  stringsAsFactors = FALSE
)
write_validation_csv(module_sizes, "module_sensitivity_sizes.csv")

package_status <- data.frame(
  component = c(
    "msigdbr", "clusterProfiler", "ReactomePA", "org.Hs.eg.db",
    "Hallmark", "GO_BP", "GO_CC", "Reactome", "KEGG"
  ),
  type = c(
    rep("R_package", 4), rep("gene_set_collection", 5)
  ),
  available = FALSE,
  status = "",
  stringsAsFactors = FALSE
)
for (package in package_status$component[package_status$type == "R_package"]) {
  installed <- requireNamespace(package, quietly = TRUE)
  package_status$available[package_status$component == package] <- installed
  package_status$status[package_status$component == package] <- if (installed) {
    paste0("installed version ", as.character(utils::packageVersion(package)))
  } else {
    "not installed in the project library"
  }
}

cache_dir <- file.path(validation_root(), ".cache", "R")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = cache_dir)

collection_specs <- list(
  Hallmark = c("H", NA_character_),
  GO_BP = c("C5", "GO:BP"),
  GO_CC = c("C5", "GO:CC"),
  Reactome = c("C2", "CP:REACTOME"),
  KEGG = c("C2", "CP:KEGG_MEDICUS")
)

full_collections <- list()
if (requireNamespace("msigdbr", quietly = TRUE)) {
  for (collection_name in names(collection_specs)) {
    spec <- collection_specs[[collection_name]]
    result <- tryCatch(
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
        append_missing_data(
          "comprehensive_module_enrichment", "MSigDB",
          collection_name, conditionMessage(error),
          "Recorded network/database failure; used curated offline subset."
        )
        attr(error, "message")
      }
    )
    if (is.data.frame(result) && nrow(result)) {
      full_collections[[collection_name]] <- result
      package_status$available[
        package_status$component == collection_name
      ] <- TRUE
      package_status$status[
        package_status$component == collection_name
      ] <- paste("downloaded", nrow(result), "gene-set membership rows")
    } else {
      package_status$status[
        package_status$component == collection_name
      ] <- "unavailable: MSigDB download/cache could not be completed"
    }
  }
} else {
  package_status$status[package_status$type == "gene_set_collection"] <-
    "unavailable because msigdbr is not installed"
}
write_validation_csv(
  package_status, "comprehensive_enrichment_status.csv"
)

curated_sets <- list(
  Hallmark = list(
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = c(
      "COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "VCAN", "SPARC",
      "TGM2", "MMP2", "MMP14", "TAGLN", "ACTA2", "SERPINE1"
    ),
    HALLMARK_TGF_BETA_SIGNALING = c(
      "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
      "SMAD7", "SERPINE1", "CTGF"
    ),
    HALLMARK_INFLAMMATORY_RESPONSE = c(
      "IL1B", "IL6", "CCL2", "CCL5", "CXCL8", "ICAM1", "VCAM1",
      "NFKB1", "STAT1", "STAT3", "TYROBP", "FCER1G"
    ),
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = c(
      "TNF", "NFKB1", "RELA", "NFKBIA", "CCL2", "CXCL8", "ICAM1",
      "JUN", "FOS"
    ),
    HALLMARK_P53_PATHWAY = c(
      "TP53", "CDKN1A", "GDF15", "SERPINE1", "BAX", "BBC3", "MDM2"
    )
  ),
  GO_BP = list(
    `GO:0030198 extracellular matrix organization` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL6A1", "COL6A2", "FN1", "SPARC", "POSTN", "VCAN", "LUM",
      "DCN", "MMP2", "MMP7", "MMP9", "MMP14", "LOX", "LOXL2"
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
  GO_CC = list(
    `GO:0031012 extracellular matrix` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "FN1",
      "SPARC", "VCAN", "LUM", "DCN", "MGP"
    ),
    `GO:0005925 focal adhesion` = c(
      "ACTN1", "VCL", "PXN", "TLN1", "ITGA5", "ITGB1", "FBLIM1",
      "PFN1", "FLNA"
    ),
    `GO:0005886 plasma membrane` = c(
      "CD44", "VCAM1", "ICAM1", "TGFBR1", "TGFBR2", "ITGA1",
      "ITGA5", "ITGB1", "PTPRE"
    ),
    `GO:0016020 membrane` = c(
      "IFITM1", "IFITM2", "IFITM3", "RAB31", "NRP1", "MPEG1",
      "MS4A7", "C1R"
    )
  ),
  Reactome = list(
    `R-HSA-1474244 Extracellular matrix organization` = c(
      "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
      "COL6A1", "COL6A2", "FN1", "SPARC", "VCAN", "LUM", "DCN",
      "MMP2", "MMP7", "MMP9", "MMP14", "LOX"
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

to_gene_sets <- function(collection) {
  if (!is.data.frame(collection) || !nrow(collection)) return(NULL)
  split(collection$gene_symbol, collection$gs_name)
}

gene_sets <- if (length(full_collections) == length(collection_specs)) {
  lapply(full_collections, to_gene_sets)
} else {
  curated_sets
}
enrichment_source <- if (
  length(full_collections) == length(collection_specs)
) {
  "MSigDB_full_collections"
} else {
  "curated_offline_subset"
}

universe <- Reduce(
  intersect, split(correlations$gene, correlations$dataset)
)

enrich_module <- function(module_genes, database, pathway, pathway_genes) {
  module_genes <- intersect(module_genes, universe)
  pathway_genes <- intersect(unique(pathway_genes), universe)
  overlap <- intersect(module_genes, pathway_genes)
  if (length(module_genes) < 2 || length(pathway_genes) < 5 ||
      length(overlap) < 2) {
    return(NULL)
  }
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
    enrichment_source = enrichment_source,
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
}

enrichment_rows <- list()
index <- 1
for (analysis in names(module_variants)) {
  for (database in names(gene_sets)) {
    for (pathway in names(gene_sets[[database]])) {
      row <- enrich_module(
        module_variants[[analysis]], database, pathway,
        gene_sets[[database]][[pathway]]
      )
      if (!is.null(row)) {
        row$analysis <- analysis
        enrichment_rows[[index]] <- row
        index <- index + 1
      }
    }
  }
}
enrichment <- if (length(enrichment_rows)) {
  do.call(rbind, enrichment_rows)
} else {
  data.frame(
    database = character(), pathway = character(),
    enrichment_source = character(), overlap_count = integer(),
    pathway_size = integer(), module_size = integer(),
    universe_size = integer(), gene_ratio = numeric(),
    fold_enrichment = numeric(), p_value = numeric(),
    overlap_genes = character(), analysis = character()
  )
}
if (nrow(enrichment)) {
  enrichment$bh_adjusted_p_value <- ave(
    enrichment$p_value, enrichment$analysis,
    FUN = function(x) p.adjust(x, "BH")
  )
  enrichment$global_bh_adjusted_p_value <- p.adjust(
    enrichment$p_value, "BH"
  )
  enrichment <- enrichment[order(
    enrichment$analysis,
    enrichment$bh_adjusted_p_value,
    -enrichment$fold_enrichment
  ), ]
}
write_validation_csv(
  enrichment, "module_sensitivity_enrichment.csv"
)

core <- stable[
  stable$fdr_positive_dataset_count == 3, , drop = FALSE
]
core_annotations <- data.frame(
  gene = c(
    "TMSB10", "IFITM3", "RAB31", "MGP", "PEA15", "TUBA1B",
    "SERPING1", "FBLIM1", "TAGLN2", "TGM2", "TUBA1A", "PFN1"
  ),
  possible_biological_function = c(
    "Actin dynamics and cytoskeletal organization",
    "Interferon-associated membrane response and cellular stress",
    "Endosomal trafficking and membrane transport",
    "Extracellular matrix mineralization and vascular matrix regulation",
    "Apoptosis and MAPK-associated signaling scaffold",
    "Microtubule organization and proliferative or repair-associated state",
    "Complement regulation and inflammatory protease inhibition",
    "Focal adhesion and actin-cytoskeleton organization",
    "Actin-associated cytoskeletal remodeling in activated cells",
    "Protein cross-linking, wound repair, and matrix remodeling",
    "Microtubule organization and cellular repair state",
    "Actin polymerization, migration, and cytoskeletal remodeling"
  ),
  ECM_related = c(
    "no", "no", "no", "yes", "no", "no",
    "indirect", "adhesion-associated", "indirect", "yes", "no", "indirect"
  ),
  inflammation_related = c(
    "indirect", "yes", "indirect", "indirect", "indirect", "no",
    "yes", "no", "immune-cell-associated", "yes", "no", "indirect"
  ),
  injury_repair_related = c(
    "yes", "yes", "yes", "yes", "yes", "yes",
    "yes", "yes", "yes", "yes", "yes", "yes"
  ),
  interpretation_boundary = paste(
    "Functional category is contextual annotation; the gene is correlated",
    "with TIMP1 and is not inferred to be directly regulated by TIMP1."
  ),
  stringsAsFactors = FALSE
)
core_summary <- merge(
  core[, c(
    "gene", "median_rho", "mean_rho", "positive_dataset_count",
    "nominal_positive_dataset_count", "fdr_positive_dataset_count",
    "rho.GSE139061", "rho.GSE30718", "rho.GSE66494"
  )],
  core_annotations, by = "gene", all.x = TRUE, sort = FALSE
)
core_summary <- core_summary[
  match(core$gene, core_summary$gene), , drop = FALSE
]
write_validation_csv(core_summary, "core_module_genes_summary.csv")

size_plot <- ggplot2::ggplot(
  module_sizes,
  ggplot2::aes(
    x = reorder(analysis, gene_count), y = gene_count,
    fill = grepl("^stringent", analysis)
  )
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = gene_count), hjust = -0.1, size = 3.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c(`FALSE` = "#4C78A8", `TRUE` = "#C44E52"),
    labels = c(`FALSE` = "Relaxed", `TRUE` = "Stringent")
  ) +
  ggplot2::expand_limits(y = max(module_sizes$gene_count) * 1.12) +
  ggplot2::labs(
    title = "TIMP1-associated module sensitivity definitions",
    x = NULL, y = "Gene count", fill = "Threshold"
  ) +
  validation_theme()
save_validation_plot(
  size_plot, "TIMP1_module_sensitivity_sizes", width = 9, height = 5.8
)

plot_enrichment <- enrichment[
  enrichment$bh_adjusted_p_value < 0.05 &
    enrichment$analysis %in% c(
      "relaxed", "stringent",
      "relaxed_without_signature_genes",
      "relaxed_without_ECM_fibrosis_markers"
    ),
  , drop = FALSE
]
plot_enrichment <- do.call(rbind, lapply(
  split(plot_enrichment, plot_enrichment$analysis),
  function(x) head(x[order(x$bh_adjusted_p_value), ], 8)
))
if (nrow(plot_enrichment)) {
  plot_enrichment$pathway_label <- paste0(
    "[", plot_enrichment$database, "] ", plot_enrichment$pathway
  )
  enrichment_plot <- ggplot2::ggplot(
    plot_enrichment,
    ggplot2::aes(
      x = fold_enrichment, y = pathway_label,
      size = overlap_count,
      color = -log10(pmax(bh_adjusted_p_value, 1e-300))
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_wrap(~analysis, scales = "free_y", ncol = 2) +
    ggplot2::scale_color_gradient(low = "#4C78A8", high = "#C44E52") +
    ggplot2::labs(
      title = "Sensitivity of TIMP1-module pathway associations",
      subtitle = enrichment_source,
      x = "Fold enrichment", y = NULL,
      size = "Overlap genes", color = "-log10(FDR)"
    ) +
    validation_theme(base_size = 10)
} else {
  enrichment_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = "No sensitivity enrichment passed FDR < 0.05."
    ) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    validation_theme()
}
save_validation_plot(
  enrichment_plot, "TIMP1_module_sensitivity_enrichment",
  width = 13, height = 10
)

rho_matrix <- as.matrix(core_summary[, c(
  "rho.GSE139061", "rho.GSE30718", "rho.GSE66494"
)])
rownames(rho_matrix) <- core_summary$gene
colnames(rho_matrix) <- c("GSE139061", "GSE30718", "GSE66494")
draw_core_heatmap <- function(path, type) {
  if (type == "png") {
    grDevices::png(path, width = 7.5, height = 7.5, units = "in", res = 600)
  } else {
    grDevices::cairo_pdf(path, width = 7.5, height = 7.5)
  }
  pheatmap::pheatmap(
    rho_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      c("white", "#F0B5B5", "#C44E52")
    )(100),
    breaks = seq(0, 1, length.out = 101),
    display_numbers = TRUE,
    number_format = "%.2f",
    border_color = "white",
    main = "Core TIMP1-module correlations"
  )
  grDevices::dev.off()
}
core_base <- file.path(
  validation_figure_dir(), "TIMP1_core_module_gene_correlations"
)
draw_core_heatmap(paste0(core_base, ".png"), "png")
draw_core_heatmap(paste0(core_base, ".pdf"), "pdf")

write_validation_session_info("sessionInfo_module_sensitivity_v2.txt")
message(
  "TIMP1 module sensitivity v2 completed: ",
  nrow(relaxed), " relaxed genes and ", nrow(stringent),
  " stringent genes."
)
