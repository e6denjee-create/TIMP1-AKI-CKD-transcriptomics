## LOXL4 in AKI-to-CKD transition: GEO analysis workflow
## Datasets: GSE139061, GSE30718, GSE66494
## Outputs: Figure2-6 publication-ready plots and CSV result files

set.seed(20260607)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

cran_pkgs <- c("ggplot2", "pROC", "pheatmap", "dplyr", "tidyr", "readr", "scales")
bioc_pkgs <- c("GEOquery", "limma", "fgsea")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(NULL))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing)
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(fgsea)
  library(ggplot2)
  library(pROC)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(readr)
})

target_gene <- "LOXL4"
fibrosis_genes <- c("COL1A1", "COL3A1", "FN1", "TGFB1")

gene_sets <- list(
  "Fibrosis core" = c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","ACTA2","TGFB1","TGFBR1","TGFBR2","CTGF","POSTN","SPARC","MMP2","TIMP1","LOX","LOXL1","LOXL2","LOXL4","VIM","SNAI1","SNAI2"),
  "ECM organization" = c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL6A1","COL6A2","COL6A3","FN1","LAMB1","LAMC1","ITGA5","ITGB1","MMP2","MMP9","TIMP1","SPARC","POSTN","THBS1","VCAN"),
  "TGF beta signaling" = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7","CTGF","SERPINE1","ID1","ID2","JUN","FOS","SNAI1","SNAI2","ZEB1","ZEB2"),
  "Kidney injury repair" = c("HAVCR1","LCN2","VCAM1","PROM1","SOX9","VIM","KRT8","KRT18","KRT19","MMP7","TIMP1","IL6","CCL2","CXCL1","CXCL8","STAT3","JUN","FOS","EGR1"),
  "Inflammation" = c("IL1B","IL6","TNF","NFKB1","NFKBIA","CCL2","CCL5","CXCL8","CXCL10","ICAM1","VCAM1","TLR2","TLR4","STAT1","STAT3","IRF1","HLA-DRA","CD68","LST1"),
  "Hypoxia oxidative stress" = c("HIF1A","VEGFA","SLC2A1","LDHA","ENO1","NQO1","HMOX1","SOD2","CAT","GPX1","TXN","TXNRD1","DDIT4","BNIP3","P4HA1","PLOD2")
)

download_if_missing <- function(url, path) {
  if (!file.exists(path) || file.info(path)$size < 100) {
    download.file(url, path, mode = "wb", quiet = FALSE)
  }
}

download_geo_files <- function() {
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139061/matrix/GSE139061_series_matrix.txt.gz", file.path(data_dir, "GSE139061_series_matrix.txt.gz"))
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30718/matrix/GSE30718_series_matrix.txt.gz", file.path(data_dir, "GSE30718_series_matrix.txt.gz"))
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66494/matrix/GSE66494_series_matrix.txt.gz", file.path(data_dir, "GSE66494_series_matrix.txt.gz"))
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139061/suppl/GSE139061_Eadon_processed_QN_101419.csv.gz", file.path(data_dir, "GSE139061_Eadon_processed_QN_101419.csv.gz"))
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz", file.path(data_dir, "GPL570.annot.gz"))
  download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz", file.path(data_dir, "GPL6480.annot.gz"))
}

read_platform <- function(gpl) {
  gpl_obj <- getGEO(gpl, AnnotGPL = TRUE, destdir = data_dir)
  tab <- Table(gpl_obj)
  sym_col <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol"), colnames(tab))[1]
  out <- tab[, c("ID", sym_col)]
  colnames(out) <- c("ID", "Gene")
  out$Gene <- trimws(vapply(strsplit(as.character(out$Gene), "///", fixed = TRUE), `[`, character(1), 1))
  out <- out[!is.na(out$Gene) & out$Gene != "", ]
  out[!duplicated(out$ID), ]
}

collapse_to_gene <- function(expr, anno = NULL) {
  if (!is.null(anno)) {
    expr <- expr[rownames(expr) %in% anno$ID, , drop = FALSE]
    gene <- anno$Gene[match(rownames(expr), anno$ID)]
  } else {
    gene <- rownames(expr)
  }
  gene <- trimws(vapply(strsplit(as.character(gene), "///", fixed = TRUE), `[`, character(1), 1))
  expr <- expr[!is.na(gene) & gene != "", , drop = FALSE]
  gene <- gene[!is.na(gene) & gene != ""]
  rowsum(expr, group = gene, reorder = FALSE) / as.vector(table(gene))
}

prepare_dataset <- function(gse) {
  if (gse == "GSE139061") {
    raw <- read_csv(file.path(data_dir, "GSE139061_Eadon_processed_QN_101419.csv.gz"), show_col_types = FALSE)
    expr <- as.matrix(raw[, -1])
    rownames(expr) <- raw[[1]]
    expr <- normalizeBetweenArrays(log2(expr + 1), method = "quantile")
    meta <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
    meta$group <- ifelse(grepl("^AKI", meta$sample), "AKI", "Control")
    meta$comparison <- "AKI_vs_Control"
    gene_expr <- collapse_to_gene(expr)
  } else {
    gset <- getGEO(gse, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = data_dir)[[1]]
    expr <- exprs(gset)
    meta <- pData(gset)
    meta$sample <- rownames(meta)
    if (max(expr, na.rm = TRUE) > 100) expr <- log2(expr + 1)
    expr <- normalizeBetweenArrays(expr, method = "quantile")
    if (gse == "GSE30718") {
      anno <- read_platform("GPL570")
      txt <- tolower(apply(meta[, grep("^characteristics", colnames(meta)), drop = FALSE], 1, paste, collapse = " "))
      meta$group <- ifelse(grepl("acute kidney injury", txt), "AKI", ifelse(grepl("protocol", tolower(meta$source_name_ch1)), "Control", "Exclude"))
      meta$comparison <- "AKI_vs_Control"
    } else {
      anno <- read_platform("GPL6480")
      txt <- tolower(apply(meta[, grep("^characteristics", colnames(meta)), drop = FALSE], 1, paste, collapse = " "))
      meta$group <- ifelse(grepl("chronic kidney disease", txt), "CKD", ifelse(grepl("normal kidney", txt), "Control", "Exclude"))
      meta$comparison <- "CKD_vs_Control"
    }
    keep <- meta$group != "Exclude"
    meta <- meta[keep, ]
    expr <- expr[, meta$sample, drop = FALSE]
    gene_expr <- collapse_to_gene(expr, anno)
  }
  write.csv(gene_expr, file.path(result_dir, paste0(gse, "_normalized_gene_expression.csv")))
  write.csv(meta, file.path(result_dir, paste0(gse, "_sample_metadata.csv")), row.names = FALSE)
  list(expr = gene_expr, meta = meta)
}

run_limma <- function(gse, expr, meta) {
  group <- factor(meta$group, levels = c("Control", setdiff(unique(meta$group), "Control")))
  design <- model.matrix(~ group)
  fit <- eBayes(lmFit(expr[, meta$sample, drop = FALSE], design))
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$Gene <- rownames(res)
  res <- res[, c("Gene", setdiff(colnames(res), "Gene"))]
  write.csv(res, file.path(result_dir, paste0(gse, "_differential_expression.csv")), row.names = FALSE)
  write.csv(res[res$Gene == target_gene, ], file.path(result_dir, paste0(gse, "_LOXL4_differential_expression.csv")), row.names = FALSE)
  res
}

run_roc <- function(gse, expr, meta) {
  y <- ifelse(meta$group == "Control", 0, 1)
  r <- roc(y, as.numeric(expr[target_gene, meta$sample]), quiet = TRUE)
  out <- data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, threshold = r$thresholds, AUC = as.numeric(auc(r)))
  write.csv(out, file.path(result_dir, paste0(gse, "_LOXL4_ROC.csv")), row.names = FALSE)
  out
}

run_correlation <- function(gse, expr, meta) {
  out <- lapply(fibrosis_genes, function(g) {
    ct <- suppressWarnings(cor.test(as.numeric(expr[target_gene, meta$sample]), as.numeric(expr[g, meta$sample]), method = "spearman"))
    data.frame(Dataset = gse, Gene = g, rho = unname(ct$estimate), P.Value = ct$p.value, n = length(meta$sample))
  }) |> bind_rows()
  out$adj.P.Val <- p.adjust(out$P.Value, method = "BH")
  write.csv(out, file.path(result_dir, paste0(gse, "_LOXL4_fibrosis_correlations.csv")), row.names = FALSE)
  out
}

run_gsea <- function(gse, de) {
  ranks <- sign(de$logFC) * -log10(pmax(de$P.Value, 1e-300))
  names(ranks) <- de$Gene
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  fg <- fgsea(pathways = gene_sets, stats = ranks, minSize = 5, maxSize = 500, nperm = 10000)
  out <- as.data.frame(fg) |> select(pathway, pval, padj, ES, NES, size, leadingEdge)
  out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse = ";")
  colnames(out) <- c("Pathway", "P.Value", "adj.P.Val", "ES", "NES", "Size", "LeadingEdge")
  out$Dataset <- gse
  write.csv(out, file.path(result_dir, paste0(gse, "_GSEA_preranked.csv")), row.names = FALSE)
  out
}

save_pub <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 600)
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height)
}

download_geo_files()

datasets <- c("GSE139061", "GSE30718", "GSE66494")
prepared <- lapply(datasets, prepare_dataset)
names(prepared) <- datasets
de <- Map(run_limma, datasets, lapply(prepared, `[[`, "expr"), lapply(prepared, `[[`, "meta"))
roc_res <- Map(run_roc, datasets, lapply(prepared, `[[`, "expr"), lapply(prepared, `[[`, "meta"))
corr_res <- Map(run_correlation, datasets, lapply(prepared, `[[`, "expr"), lapply(prepared, `[[`, "meta"))
gsea_res <- Map(run_gsea, datasets, de)

fig2 <- bind_rows(lapply(datasets, function(gse) {
  expr <- prepared[[gse]]$expr
  meta <- prepared[[gse]]$meta
  data.frame(Dataset = gse, Group = meta$group, LOXL4 = as.numeric(expr[target_gene, meta$sample]))
}))
write.csv(fig2, file.path(result_dir, "Figure2_LOXL4_expression.csv"), row.names = FALSE)
p2 <- ggplot(fig2, aes(Dataset, LOXL4, fill = Group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 1.4, alpha = 0.55) +
  scale_fill_manual(values = c(Control = "#4C78A8", AKI = "#D65F5F", CKD = "#B279A2")) +
  labs(title = "Figure 2. LOXL4 expression in AKI and CKD datasets", x = NULL, y = "LOXL4 expression (normalized log2)") +
  theme_bw(base_size = 12) + theme(legend.title = element_blank())
save_pub(p2, "Figure2_LOXL4_differential_expression", 8, 4.8)

fig3 <- bind_rows(roc_res, .id = "Dataset")
write.csv(bind_rows(lapply(names(roc_res), function(x) data.frame(Dataset = x, AUC = unique(roc_res[[x]]$AUC)))), file.path(result_dir, "Figure3_LOXL4_ROC_summary.csv"), row.names = FALSE)
p3 <- ggplot(fig3, aes(FPR, TPR, color = Dataset)) +
  geom_line(linewidth = 1.1) + geom_abline(linetype = 2, color = "grey45") +
  labs(title = "Figure 3. Diagnostic ROC for LOXL4", x = "False positive rate", y = "True positive rate") +
  theme_bw(base_size = 12)
save_pub(p3, "Figure3_LOXL4_ROC", 5.7, 5.2)

fig4 <- bind_rows(corr_res)
write.csv(fig4, file.path(result_dir, "Figure4_LOXL4_fibrosis_correlations.csv"), row.names = FALSE)
p4 <- ggplot(fig4, aes(Gene, Dataset, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.4) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#D65F5F", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Figure 4. LOXL4 correlation with fibrosis genes", x = NULL, y = NULL, fill = "Spearman rho") +
  theme_bw(base_size = 12)
save_pub(p4, "Figure4_LOXL4_fibrosis_correlation_heatmap", 7.6, 4.5)

fig5 <- bind_rows(gsea_res)
write.csv(fig5, file.path(result_dir, "Figure5_GSEA_summary.csv"), row.names = FALSE)
p5 <- ggplot(fig5, aes(Dataset, Pathway, color = NES, size = -log10(pmax(adj.P.Val, 1e-300)))) +
  geom_point() +
  scale_color_gradient2(low = "#4C78A8", mid = "white", high = "#D65F5F", midpoint = 0) +
  labs(title = "Figure 5. Preranked GSEA across LOXL4-associated disease contrasts", x = NULL, y = NULL, size = "-log10(FDR)") +
  theme_bw(base_size = 12)
save_pub(p5, "Figure5_GSEA_dotplot", 8.2, 5.4)

fig6 <- bind_rows(lapply(datasets, function(gse) {
  lox <- de[[gse]][de[[gse]]$Gene == target_gene, ]
  data.frame(Dataset = gse, LOXL4_logFC = lox$logFC, LOXL4_P.Value = lox$P.Value,
             LOXL4_adj.P.Val = lox$adj.P.Val, LOXL4_AUC = unique(roc_res[[gse]]$AUC),
             Mean_fibrosis_rho = mean(corr_res[[gse]]$rho, na.rm = TRUE))
}))
write.csv(fig6, file.path(result_dir, "Figure6_integrated_LOXL4_evidence.csv"), row.names = FALSE)
p6 <- fig6 |>
  pivot_longer(cols = c(LOXL4_logFC, LOXL4_AUC, Mean_fibrosis_rho), names_to = "Metric", values_to = "Value") |>
  mutate(Metric = recode(Metric, LOXL4_logFC = "logFC", LOXL4_AUC = "AUC", Mean_fibrosis_rho = "mean rho")) |>
  ggplot(aes(Value, Dataset)) +
  geom_col(fill = "#5B8C85") +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  labs(title = "Figure 6. Integrated evidence for LOXL4 in AKI-to-CKD transition", x = NULL, y = NULL) +
  theme_bw(base_size = 12)
save_pub(p6, "Figure6_integrated_LOXL4_evidence", 9, 4.3)

message("Analysis complete.")
message("Results: ", result_dir)
message("Figures: ", figure_dir)
