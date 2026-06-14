# Shared helpers for non-destructive TIMP1 validation workflows.

validation_root <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

validation_result_dir <- function(root = validation_root()) {
  file.path(root, "results", "timp1_validation")
}

validation_figure_dir <- function(root = validation_root()) {
  file.path(root, "figures", "timp1_validation")
}

activate_timp1_library <- function(root = validation_root()) {
  library_dir <- file.path(root, "renv", "library", "R-4.6")
  if (dir.exists(library_dir)) .libPaths(c(library_dir, .libPaths()))
  invisible(.libPaths())
}

ensure_validation_dirs <- function(root = validation_root()) {
  dirs <- c(validation_result_dir(root), validation_figure_dir(root))
  invisible(vapply(
    dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE
  ))
}

append_missing_data <- function(step, dataset, missing_item, reason,
                                action = "Skipped affected analysis unit",
                                root = validation_root()) {
  ensure_validation_dirs(root)
  path <- file.path(validation_result_dir(root), "missing_data_log.csv")
  row <- data.frame(
    analysis_step = step,
    dataset = dataset,
    missing_item = missing_item,
    reason = reason,
    action_taken = action,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    existing <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    all_columns <- union(colnames(existing), colnames(row))
    for (column in setdiff(all_columns, colnames(existing))) {
      existing[[column]] <- ""
    }
    for (column in setdiff(all_columns, colnames(row))) {
      row[[column]] <- ""
    }
    existing <- existing[, all_columns, drop = FALSE]
    row <- row[, all_columns, drop = FALSE]
    row <- rbind(existing, row)
  }
  write.csv(unique(row), path, row.names = FALSE, na = "")
  invisible(path)
}

ensure_missing_log <- function(root = validation_root()) {
  ensure_validation_dirs(root)
  path <- file.path(validation_result_dir(root), "missing_data_log.csv")
  if (!file.exists(path)) {
    write.csv(
      data.frame(
        analysis_step = character(),
        dataset = character(),
        missing_item = character(),
        reason = character(),
        action_taken = character(),
        timestamp = character()
      ),
      path, row.names = FALSE
    )
  }
  invisible(path)
}

read_validation_bulk <- function(dataset, root = validation_root()) {
  expression_path <- file.path(
    root, "data", "processed",
    paste0(dataset, "_normalized_expression.csv.gz")
  )
  metadata_path <- file.path(
    root, "data", "metadata", paste0(dataset, "_metadata.csv")
  )
  missing <- c(
    if (!file.exists(expression_path)) expression_path,
    if (!file.exists(metadata_path)) metadata_path
  )
  if (length(missing)) {
    append_missing_data(
      "bulk_input", dataset, paste(missing, collapse = "; "),
      "Required expression matrix or metadata file is absent.", root = root
    )
    return(NULL)
  }

  expression <- tryCatch(
    as.matrix(read.csv(
      gzfile(expression_path), row.names = 1, check.names = FALSE
    )),
    error = function(error) {
      append_missing_data(
        "bulk_input", dataset, expression_path, conditionMessage(error),
        root = root
      )
      NULL
    }
  )
  if (is.null(expression)) return(NULL)
  storage.mode(expression) <- "numeric"

  metadata <- tryCatch(
    read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(error) {
      append_missing_data(
        "bulk_input", dataset, metadata_path, conditionMessage(error),
        root = root
      )
      NULL
    }
  )
  if (is.null(metadata)) return(NULL)
  required_columns <- c("sample", "group")
  absent_columns <- setdiff(required_columns, colnames(metadata))
  if (length(absent_columns)) {
    append_missing_data(
      "bulk_input", dataset, paste(absent_columns, collapse = "; "),
      "Required metadata columns are absent.", root = root
    )
    return(NULL)
  }
  available_samples <- intersect(metadata$sample, colnames(expression))
  missing_samples <- setdiff(metadata$sample, colnames(expression))
  if (length(missing_samples)) {
    append_missing_data(
      "bulk_input", dataset, paste(missing_samples, collapse = "; "),
      "Metadata samples are absent from the expression matrix.",
      "Excluded unmatched samples and continued.", root
    )
  }
  metadata <- metadata[match(available_samples, metadata$sample), , drop = FALSE]
  expression <- expression[, available_samples, drop = FALSE]
  list(dataset = dataset, expression = expression, metadata = metadata)
}

save_validation_plot <- function(plot, stem, width = 8, height = 6,
                                 root = validation_root(), dpi = 600) {
  ensure_validation_dirs(root)
  base <- file.path(validation_figure_dir(root), stem)
  ggplot2::ggsave(
    paste0(base, ".png"), plot, width = width, height = height,
    dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot, width = width, height = height,
    device = grDevices::cairo_pdf, bg = "white"
  )
  invisible(base)
}

write_validation_csv <- function(x, filename, root = validation_root()) {
  ensure_validation_dirs(root)
  path <- file.path(validation_result_dir(root), filename)
  write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

validation_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

write_validation_session_info <- function(filename, root = validation_root()) {
  path <- file.path(validation_result_dir(root), filename)
  capture.output(sessionInfo(), file = path)
  invisible(path)
}
