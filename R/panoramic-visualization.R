#' Volcano plot for differential spatial colocalization
#'
#' Updated volcano helper based on the CRC TMA manuscript workflow.
#' It uses \code{beta_diff} for the x-axis and \code{-log10(p_diff)} for
#' the y-axis, and colors points by significance direction.
#'
#' @param se_diff A \code{SummarizedExperiment} with contrast columns in
#'   \code{rowData(se_diff)} (at least \code{ct1}, \code{ct2},
#'   \code{beta_diff}, \code{p_diff}; and \code{fdr_diff} when
#'   \code{sig_col = "fdr_diff"}).
#' @param sig_col Character, one of \code{"fdr_diff"} or \code{"p_diff"}.
#'   Controls which column determines significance coloring.
#' @param alpha Numeric threshold for significance based on \code{sig_col}.
#' @param x_scale Character, one of \code{"beta_diff"} or \code{"log2fc"}.
#'   \code{"log2fc"} uses \code{sign(beta_diff) * log2(1 + abs(beta_diff))}.
#' @param label_top Integer number of significant hits to label. Use
#'   \code{NULL} or \code{0} to disable labels.
#' @param label_text_pt Numeric font size (points) for significant-hit labels.
#' @param p_floor Numeric floor applied to \code{p_diff} before
#'   \code{-log10()}.
#' @param title Optional character title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' se_diff <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(dummy = matrix(0, nrow = 2, ncol = 1)),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = c("A", "A"),
#'     ct2 = c("B", "C"),
#'     radius_um = c(10, 10),
#'     beta_diff = c(0.4, -0.3),
#'     p_diff = c(0.01, 0.20),
#'     fdr_diff = c(0.02, 0.30)
#'   )
#' )
#' m <- S4Vectors::metadata(se_diff)
#' m$panoramic <- list(contrast = list(control = "control", case = "case"))
#' S4Vectors::metadata(se_diff) <- m
#' plot_volcano(se_diff, sig_col = "fdr_diff")
#' @export
plot_volcano <- function(
    se_diff,
    sig_col = c("fdr_diff", "p_diff"),
    alpha = 0.05,
    x_scale = c("beta_diff", "log2fc"),
    label_top = 12L,
    label_text_pt = 4,
    p_floor = 1e-300,
    title = NULL
) {
  sig_col <- match.arg(sig_col)
  x_scale <- match.arg(x_scale)

  rd <- as.data.frame(SummarizedExperiment::rowData(se_diff))
  required <- c("ct1", "ct2", "beta_diff", "p_diff")
  if (identical(sig_col, "fdr_diff")) {
    required <- c(required, "fdr_diff")
  }
  missing <- setdiff(required, colnames(rd))
  if (length(missing) > 0L) {
    stop(
      "Missing required columns in rowData(se_diff): ",
      paste(missing, collapse = ", ")
    )
  }

  contrast_meta <- S4Vectors::metadata(se_diff)$panoramic$contrast
  control_group <- contrast_meta$control %||% "control"
  case_group <- contrast_meta$case %||% "case"

  df <- rd
  df$coloc_source <- as.character(df$ct2)
  df$coloc_target <- as.character(df$ct1)
  df$feature_label <- paste0(df$coloc_source, " -> ", df$coloc_target)
  df$p_plot <- pmax(df$p_diff, p_floor)
  df$neg_log10_p <- -log10(df$p_plot)
  df$sig_value <- df[[sig_col]]
  df$is_sig <- is.finite(df$sig_value) & (df$sig_value <= alpha)
  df$direction_group <- ifelse(df$beta_diff >= 0, case_group, control_group)
  df$color_class <- ifelse(df$is_sig, df$direction_group, "Not significant")
  df$x_value <- if (identical(x_scale, "log2fc")) {
    sign(df$beta_diff) * log2(1 + abs(df$beta_diff))
  } else {
    df$beta_diff
  }
  df <- df[is.finite(df$x_value) & is.finite(df$neg_log10_p), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No finite rows available for volcano plotting.")
  }

  top_labels <- df[0, , drop = FALSE]
  if (!is.null(label_top) && is.finite(label_top) && label_top > 0L) {
    top_labels <- df[df$is_sig, , drop = FALSE]
    if (nrow(top_labels) > 0L) {
      ord <- order(top_labels$sig_value, top_labels$p_diff, na.last = TRUE)
      top_labels <- top_labels[ord, , drop = FALSE]
      top_labels <- utils::head(top_labels, as.integer(label_top))
    }
  }
  label_text_pt <- suppressWarnings(as.numeric(label_text_pt))
  if (!is.finite(label_text_pt) || label_text_pt <= 0) {
    label_text_pt <- 4
  }
  label_size <- label_text_pt / ggplot2::.pt

  col_map <- if (identical(sig_col, "fdr_diff")) {
    stats::setNames(
      c("grey80", "#2F5597", "#C00000"),
      c("Not significant", control_group, case_group)
    )
  } else {
    stats::setNames(
      c("grey80", "#6BAED6", "#FB6A4A"),
      c("Not significant", control_group, case_group)
    )
  }

  radii_vals <- sort(unique(df$radius_um[is.finite(df$radius_um)]))
  radius_text <- if (length(radii_vals) == 0L) {
    "(r=NA um)"
  } else {
    paste0(
      "(r=",
      paste(format(radii_vals, trim = TRUE, scientific = FALSE), collapse = ", "),
      " um)"
    )
  }

  if (is.null(title)) {
    title <- paste("Differential spatial colocalization", radius_text)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = x_value,
      y = neg_log10_p,
      color = color_class
    )
  ) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    ggplot2::geom_point(alpha = 0.9, size = 1.8) +
    ggplot2::scale_color_manual(values = col_map, breaks = c("Not significant", control_group, case_group)) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Left: higher in ", control_group, " | Right: higher in ", case_group),
      x = if (identical(x_scale, "log2fc")) {
        paste0("log2FC (", case_group, " - ", control_group, ")")
      } else {
        paste0("Group contrast (beta_diff = ", case_group, " - ", control_group, ")")
      },
      y = expression(-log[10](raw~p)),
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 11)

  if (nrow(top_labels) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = top_labels,
        ggplot2::aes(label = feature_label),
        size = label_size,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = top_labels,
        ggplot2::aes(label = feature_label),
        size = label_size,
        vjust = -0.2,
        check_overlap = TRUE,
        show.legend = FALSE
      )
    }
  }

  p
}

.panoramic_make_feature_key <- function(df) {
  paste(df$ct1, df$ct2, df$radius_um, sep = "|")
}

.panoramic_choose_representative_sample <- function(values, target) {
  ok <- is.finite(values)
  if (!any(ok)) {
    return(NA_integer_)
  }
  idx_ok <- which(ok)
  if (!is.finite(target)) {
    return(idx_ok[1])
  }
  idx_ok[which.min(abs(values[idx_ok] - target))]
}

#' Plot Representative Samples for Significant Differential Hits
#'
#' Select representative samples per group for significant differential
#' colocalization features, then build side-by-side spatial panels
#' highlighting source and target cell types.
#'
#' @param se_stats A \code{SummarizedExperiment} from
#'   \code{panoramic_spatialstats()} containing sample-level \code{yi}.
#' @param se_meta A \code{SummarizedExperiment} from
#'   \code{panoramic_meta_mv()} containing contrast columns.
#' @param spe_list Named list of \code{SpatialExperiment} objects.
#' @param top_hits Optional data.frame of selected features. If \code{NULL},
#'   features are selected from \code{se_meta} using \code{sig_col},
#'   \code{alpha}, and \code{top_n}.
#' @param sig_col One of \code{"fdr_diff"} or \code{"p_diff"} used when
#'   selecting \code{top_hits = NULL}.
#' @param alpha Numeric threshold for feature selection when
#'   \code{top_hits = NULL}.
#' @param top_n Integer number of selected hits when \code{top_hits = NULL}.
#' @param group_col Character group column in \code{colData(se_stats)}.
#' @param cell_type_col Character cell-type column in \code{colData(spe)}.
#' @param sample_col Character sample-id column in \code{colData(se_stats)}.
#' @param out_prefix Optional output prefix. If provided, each panel is saved as
#'   PNG/PDF and the returned index includes file paths.
#'
#' @return A list with \code{plots} (named list of ggplot objects) and
#'   \code{index} (data.frame summarizing selected hits and samples).
#'
#' @details
#' This helper currently supports exactly two groups in \code{group_col}.
#'
#' @examples
#' sample_ids <- paste0("sample_", seq_len(4L))
#' spe_list <- stats::setNames(
#'   lapply(sample_ids, function(sid) {
#'     panoramic_simulate_spe(
#'       n_cells = 40,
#'       sample_id = sid,
#'       cell_types = c("A", "B"),
#'       scenario = "random",
#'       seed = 1
#'     )
#'   }),
#'   sample_ids
#' )
#'
#' yi <- matrix(c(0.10, 0.15, 0.25, 0.30), nrow = 1)
#' vi <- matrix(rep(0.05, 4), nrow = 1)
#' colnames(yi) <- colnames(vi) <- sample_ids
#' se_stats <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(yi = yi, vi = vi),
#'   rowData = S4Vectors::DataFrame(ct1 = "A", ct2 = "B", radius_um = 10),
#'   colData = S4Vectors::DataFrame(
#'     sample = sample_ids,
#'     group = c("control", "control", "case", "case")
#'   )
#' )
#'
#' se_meta <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(dummy = matrix(0, nrow = 1, ncol = 1)),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = "A",
#'     ct2 = "B",
#'     radius_um = 10,
#'     beta_diff = 0.25,
#'     se_diff = 0.10,
#'     z_diff = 2.5,
#'     p_diff = 0.012,
#'     fdr_diff = 0.02,
#'     control_mu_hat = 0.125,
#'     case_mu_hat = 0.275
#'   )
#' )
#'
#' out <- plot_representative_samples(
#'   se_stats = se_stats,
#'   se_meta = se_meta,
#'   spe_list = spe_list,
#'   sig_col = "p_diff",
#'   alpha = 0.05,
#'   top_n = 1
#' )
#' names(out)
#' @export
plot_representative_samples <- function(
    se_stats,
    se_meta,
    spe_list,
    top_hits = NULL,
    sig_col = c("fdr_diff", "p_diff"),
    alpha = 0.05,
    top_n = 10L,
    group_col = "group",
    cell_type_col = "cell_type",
    sample_col = "sample",
    out_prefix = NULL
) {
  sig_col <- match.arg(sig_col)
  if (!inherits(se_stats, "SummarizedExperiment")) {
    stop("`se_stats` must be a SummarizedExperiment.")
  }
  if (!inherits(se_meta, "SummarizedExperiment")) {
    stop("`se_meta` must be a SummarizedExperiment.")
  }
  if (!is.list(spe_list) || is.null(names(spe_list))) {
    stop("`spe_list` must be a named list of SpatialExperiment objects.")
  }

  rd_stats <- as.data.frame(SummarizedExperiment::rowData(se_stats))
  rd_meta <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  cd_stats <- as.data.frame(SummarizedExperiment::colData(se_stats))
  yi <- SummarizedExperiment::assay(se_stats, "yi")

  needed_stats_cols <- c("ct1", "ct2", "radius_um")
  needed_meta_cols <- c("ct1", "ct2", "radius_um")
  if (!all(needed_stats_cols %in% colnames(rd_stats))) {
    stop("rowData(se_stats) must include: ct1, ct2, radius_um.")
  }
  if (!all(needed_meta_cols %in% colnames(rd_meta))) {
    stop("rowData(se_meta) must include: ct1, ct2, radius_um.")
  }
  if (!all(c(group_col, sample_col) %in% colnames(cd_stats))) {
    stop("colData(se_stats) must include: ", group_col, ", ", sample_col)
  }

  if (is.null(top_hits)) {
    top_hits <- panoramic_extract_contrast(se_meta)
    if (!all(c("ct1", "ct2", "radius_um", sig_col, "p_diff") %in% colnames(top_hits))) {
      stop("Extracted contrast table is missing required columns.")
    }
    top_hits <- top_hits[
      is.finite(top_hits[[sig_col]]) & (top_hits[[sig_col]] <= alpha),
      ,
      drop = FALSE
    ]
    if (nrow(top_hits) > 0L) {
      ord <- order(top_hits[[sig_col]], top_hits$p_diff, na.last = TRUE)
      top_hits <- top_hits[ord, , drop = FALSE]
      if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
        top_hits <- utils::head(top_hits, as.integer(top_n))
      }
    }
  }

  if (nrow(top_hits) == 0L) {
    return(list(plots = list(), index = data.frame()))
  }

  group_raw <- cd_stats[[group_col]]
  groups <- if (is.factor(group_raw)) {
    lev <- levels(group_raw)
    lev[lev %in% as.character(group_raw[!is.na(group_raw)])]
  } else {
    unique(as.character(group_raw))
  }
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) != 2L) {
    stop(
      "Representative sample plotting currently expects exactly 2 groups in `",
      group_col, "`. Found: ", paste(groups, collapse = ", "), "."
    )
  }

  key_stats <- .panoramic_make_feature_key(rd_stats)
  key_meta <- .panoramic_make_feature_key(rd_meta)

  sanitize_slug <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    ifelse(nzchar(x), x, "NA")
  }

  plots <- list()
  index_rows <- list()
  k_out <- 0L

  for (i in seq_len(nrow(top_hits))) {
    hit <- top_hits[i, , drop = FALSE]
    key <- .panoramic_make_feature_key(hit)
    idx_stats <- which(key_stats == key)
    idx_meta <- which(key_meta == key)
    if (length(idx_stats) != 1L || length(idx_meta) != 1L) {
      next
    }
    row_stat <- idx_stats[1]
    row_meta <- idx_meta[1]
    y_row <- yi[row_stat, ]

    ct1 <- as.character(hit$ct1)
    ct2 <- as.character(hit$ct2)
    source_ct <- ct2
    target_ct <- ct1

    rep_pts <- list()
    panel_levels <- character(0)
    picked_samples <- stats::setNames(rep(NA_character_, length(groups)), groups)
    j_out <- 0L

    for (g in groups) {
      idx_group <- which(as.character(cd_stats[[group_col]]) == g)
      if (length(idx_group) == 0L) next

      mu_col <- paste0(make.names(g), "_mu_hat")
      mu_target <- if (mu_col %in% colnames(rd_meta)) rd_meta[row_meta, mu_col] else NA_real_
      pick_local <- .panoramic_choose_representative_sample(y_row[idx_group], mu_target)
      if (is.na(pick_local)) next

      pick_col <- idx_group[pick_local]
      sample_id <- as.character(cd_stats[[sample_col]][pick_col])
      if (!sample_id %in% names(spe_list)) next

      spe <- spe_list[[sample_id]]
      coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
      if (ncol(coords) < 2) next
      colnames(coords)[seq_len(2L)] <- c("x", "y")
      ct <- as.character(SummarizedExperiment::colData(spe)[[cell_type_col]])
      role <- ifelse(ct == source_ct, source_ct, ifelse(ct == target_ct, target_ct, "Other"))

      panel_label <- paste0(g, " | ", sample_id)
      panel_levels <- c(panel_levels, panel_label)
      picked_samples[g] <- sample_id
      j_out <- j_out + 1L
      rep_pts[[j_out]] <- data.frame(
        x = coords$x,
        y = coords$y,
        role = role,
        group = g,
        sample = sample_id,
        panel = panel_label,
        stringsAsFactors = FALSE
      )
    }

    if (length(rep_pts) == 0L) next
    rep_df <- dplyr::bind_rows(rep_pts)
    rep_df$panel <- factor(rep_df$panel, levels = unique(panel_levels))

    hit_title <- paste0(
      "Representative samples: ", source_ct, " -> ", target_ct,
      " (r=", format(hit$radius_um, trim = TRUE, scientific = FALSE), " um)"
    )
    col_map <- stats::setNames(c("grey80", "#1F77B4", "#FF7F0E"), c("Other", source_ct, target_ct))
    p <- ggplot2::ggplot(rep_df, ggplot2::aes(x = x, y = y, color = role)) +
      ggplot2::geom_point(size = 0.55, alpha = 0.9) +
      ggplot2::scale_color_manual(values = col_map) +
      ggplot2::coord_equal() +
      ggplot2::facet_wrap(~panel, nrow = 1, scales = "fixed") +
      ggplot2::labs(
        title = hit_title,
        subtitle = "Highlighted cell types in blue/orange; all other cells in gray",
        x = "x",
        y = "y",
        color = "Cell type"
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8), legend.position = "bottom")

    file_png <- NA_character_
    file_pdf <- NA_character_
    if (!is.null(out_prefix) && nzchar(out_prefix)) {
      radius_slug <- gsub("\\.", "p", format(hit$radius_um, trim = TRUE, scientific = FALSE))
      file_stub <- paste0(
        out_prefix, "_",
        sanitize_slug(source_ct), "_to_", sanitize_slug(target_ct),
        "_r", radius_slug, "um"
      )
      file_png <- paste0(file_stub, ".png")
      file_pdf <- paste0(file_stub, ".pdf")
      ggplot2::ggsave(file_png, p, width = 10, height = 4.6, dpi = 300)
      ggplot2::ggsave(file_pdf, p, width = 10, height = 4.6)
    }

    k_out <- k_out + 1L
    plot_name <- paste0(source_ct, "_to_", target_ct, "_r", as.character(hit$radius_um))
    plots[[plot_name]] <- p
    index_rows[[k_out]] <- data.frame(
      hit_rank = i,
      ct1 = ct1,
      ct2 = ct2,
      radius_um = hit$radius_um,
      sample_group1 = picked_samples[groups[1]],
      sample_group2 = picked_samples[groups[2]],
      png = file_png,
      pdf = file_pdf,
      stringsAsFactors = FALSE
    )
  }

  index_df <- if (length(index_rows) == 0L) {
    data.frame()
  } else {
    dplyr::bind_rows(index_rows)
  }
  list(plots = plots, index = index_df)
}


#' Forest plot for PANORAMIC spatial colocalization
#'
#' Create a forest plot showing individual-sample spatial statistics and
#' pooled random-effects estimates for a single cell-type pair at a given
#' radius. Styling and text columns can be customized for figure tuning.
#'
#' @param se_meta A SummarizedExperiment returned by
#'   \code{panoramic_meta_mv()} or a downstream testing function, containing
#'   assays "yi" (effect estimates) and "vi" (variances), and rowData
#'   columns such as ct1, ct2, radius_um, and, if group-level pooling was
#'   performed, mu_control, se_control, mu_case, se_case.
#' @param ct1 Character. First cell type (must match rowData(se_meta)$ct1).
#' @param ct2 Character. Second cell type (must match rowData(se_meta)$ct2).
#' @param radius_um Numeric radius (in microns) to plot. If NULL,
#'   the first available radius for this pair is used.
#' @param group_col Character. Column name in colData(se_meta)
#'   containing group labels (e.g. treatment vs control). Default "group".
#' @param group_colors Optional named character vector mapping group labels
#'   to colors. If NULL, a default palette is used.
#' @param text_size Numeric text size for row annotations.
#' @param header_text_size Numeric text size for header annotations.
#' @param title_text_size Numeric text size for plot title.
#' @param base_size Numeric base font size passed to \code{theme_classic()}.
#' @param show_est_se Logical. If TRUE, show the Est (SE) text column.
#' @param show_ci Logical. If TRUE, show the 95% CI text column.
#'
#' @return A ggplot object representing the forest plot.
#'
#' @details
#' This helper currently supports exactly two groups in \code{group_col}.
#' 
#' @examples
#' yi <- matrix(c(0.12, 0.18, 0.35, 0.30), nrow = 1)
#' vi <- matrix(c(0.03, 0.04, 0.05, 0.05), nrow = 1)
#' colnames(yi) <- colnames(vi) <- paste0("s", seq_len(4L))
#'
#' se_meta <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(yi = yi, vi = vi),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = "A",
#'     ct2 = "B",
#'     radius_um = 10,
#'     control_mu_hat = 0.15,
#'     control_se_mu = 0.10,
#'     case_mu_hat = 0.32,
#'     case_se_mu = 0.11
#'   ),
#'   colData = S4Vectors::DataFrame(
#'     group = c("control", "control", "case", "case")
#'   )
#' )
#'
#' p_forest <- plot_forest(
#'   se_meta,
#'   ct1 = "A",
#'   ct2 = "B",
#'   radius_um = 10,
#'   group_col = "group"
#' )
#' p_forest
#' @export
plot_forest <- function(se_meta, ct1, ct2, radius_um = NULL,
                        group_col = "group",
                        group_colors = NULL,
                        text_size = 2,
                        header_text_size = 2,
                        title_text_size = 8,
                        base_size = 5,
                        show_est_se = TRUE,
                        show_ci = TRUE) {

  yi <- SummarizedExperiment::assay(se_meta, "yi")
  vi <- SummarizedExperiment::assay(se_meta, "vi")
  rd <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  cd <- as.data.frame(SummarizedExperiment::colData(se_meta))

  if (is.null(radius_um)) {
    radius_um <- rd$radius_um[1]
    message("Using radius: ", radius_um, " um")
  }

  row_idx <- which(rd$ct1 == ct1 & rd$ct2 == ct2 & rd$radius_um == radius_um)
  if (length(row_idx) == 0L) {
    stop("Cell pair ", ct1, ":", ct2, " at radius ", radius_um, " not found")
  }
  row_idx <- row_idx[1]

  y_vals <- yi[row_idx, ]
  v_vals <- vi[row_idx, ]

  sample_data <- data.frame(
    sample = colnames(yi),
    estimate = y_vals,
    variance = v_vals,
    se = sqrt(v_vals),
    group = cd[[group_col]],
    type = "Individual",
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(estimate) & !is.na(se)) |>
    dplyr::mutate(
      ci_lower = estimate - 1.96 * se,
      ci_upper = estimate + 1.96 * se,
      label = sample,
      stats_text = sprintf("%.2f (%.2f)", estimate, se)
    )

  group_raw <- cd[[group_col]]
  groups <- if (is.factor(group_raw)) {
    lev <- levels(group_raw)
    lev[lev %in% as.character(group_raw[!is.na(group_raw)])]
  } else {
    unique(as.character(group_raw))
  }
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) != 2L) {
    stop(
      "plot_forest currently expects exactly two groups in `", group_col,
      "`. Found: ", paste(groups, collapse = ", "), "."
    )
  }
  grp_prefix <- make.names(groups)
  ctrl_pref <- grp_prefix[1]
  case_pref <- grp_prefix[2]

  ctrl_mu_col <- paste0(ctrl_pref, "_mu_hat")
  case_mu_col <- paste0(case_pref, "_mu_hat")
  ctrl_se_col <- paste0(ctrl_pref, "_se_mu")
  case_se_col <- paste0(case_pref, "_se_mu")

  has_pooled <- all(c(ctrl_mu_col, case_mu_col, ctrl_se_col, case_se_col) %in% colnames(rd))

  pooled_data <- NULL
  if (has_pooled) {
    pooled_data <- data.frame(
      sample = c("Pooled (RE)", "Pooled (RE)"),
      estimate = c(rd[[ctrl_mu_col]][row_idx], rd[[case_mu_col]][row_idx]),
      se = c(rd[[ctrl_se_col]][row_idx], rd[[case_se_col]][row_idx]),
      group = c(groups[1], groups[2]),
      type = "Pooled",
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        ci_lower = estimate - 1.96 * se,
        ci_upper = estimate + 1.96 * se,
        label = "Pooled (RE)",
        stats_text = sprintf("%.2f (%.2f)", estimate, se)
      )
  }

  if (is.null(pooled_data)) {
    warning("No pooled group estimates found for this feature; showing individual samples only.")
  }

  plot_data <- dplyr::bind_rows(sample_data, pooled_data) |>
    dplyr::arrange(group, dplyr::desc(type), estimate) |>
    dplyr::group_by(group) |>
    dplyr::mutate(row_within_group = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      type = factor(type, levels = c("Individual", "Pooled")),
      group_offset = ifelse(group == unique(group)[1], 0,
                            max(row_within_group[group == unique(group)[1]]) + 2),
      y_pos = row_within_group + group_offset + ifelse(type == "Pooled", 0.2, 0)
    )

  if (is.null(group_colors)) {
    groups <- unique(plot_data$group)
    group_colors <- stats::setNames(c("#1f77b4", "#d62728")[seq_along(groups)], groups)
  }

  x_min <- min(plot_data$ci_lower, na.rm = TRUE)
  x_max <- max(plot_data$ci_upper, na.rm = TRUE)
  x_range <- x_max - x_min
  if (!is.finite(x_range) || x_range <= 0) x_range <- 1
  x_est <- x_max + x_range * 0.02
  x_ci <- if (isTRUE(show_est_se)) x_max + x_range * 0.20 else x_max + x_range * 0.02

  right_expand <- if (isTRUE(show_est_se) && isTRUE(show_ci)) {
    0.35
  } else if (isTRUE(show_est_se) || isTRUE(show_ci)) {
    0.22
  } else {
    0.10
  }
  right_margin <- if (isTRUE(show_est_se) && isTRUE(show_ci)) {
    50
  } else if (isTRUE(show_est_se) || isTRUE(show_ci)) {
    35
  } else {
    20
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = y_pos)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lower, xend = ci_upper, yend = y_pos, color = group),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, type == "Individual"),
      ggplot2::aes(x = estimate, color = group),
      size = 1.7, shape = 22, stroke = 0.5, fill = "white",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, type == "Pooled"),
      ggplot2::aes(x = estimate, color = group, fill = group),
      size = 2.4, shape = 23, stroke = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = x_min - x_range * 0.02,
        label = label,
        color = group,
        fontface = ifelse(type == "Pooled", "bold", "plain")
      ),
      hjust = 1, size = text_size, show.legend = FALSE
    ) +
    {
      if (isTRUE(show_est_se)) {
        ggplot2::geom_text(
          ggplot2::aes(
            x = x_est,
            label = stats_text,
            color = group,
            fontface = ifelse(type == "Pooled", "bold", "plain")
          ),
          hjust = 0, size = text_size, family = "mono", show.legend = FALSE
        )
      }
    } +
    {
      if (isTRUE(show_ci)) {
        ggplot2::geom_text(
          ggplot2::aes(
            x = x_ci,
            label = sprintf("[%.2f, %.2f]", ci_lower, ci_upper),
            color = group,
            fontface = ifelse(type == "Pooled", "bold", "plain")
          ),
          hjust = 0, size = text_size, family = "mono", show.legend = FALSE
        )
      }
    } +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray30",
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.15, right_expand))) +
    ggplot2::scale_y_continuous(breaks = NULL, labels = NULL) +
    ggplot2::annotate(
      "text", x = x_min - x_range * 0.02,
      y = max(plot_data$y_pos) + 1,
      label = "Sample", hjust = 1,
      fontface = "bold", size = header_text_size
    ) +
    {
      if (isTRUE(show_est_se)) {
        ggplot2::annotate(
          "text", x = x_est,
          y = max(plot_data$y_pos) + 1,
          label = "Est (SE)", hjust = 0,
          fontface = "bold", size = header_text_size
        )
      }
    } +
    {
      if (isTRUE(show_ci)) {
        ggplot2::annotate(
          "text", x = x_ci,
          y = max(plot_data$y_pos) + 1,
          label = "95% CI", hjust = 0,
          fontface = "bold", size = header_text_size
        )
      }
    } +
    ggplot2::labs(
      x = paste0("Colocalization at ", format(radius_um, trim = TRUE, scientific = FALSE), " um"),
      y = NULL,
      color = NULL,
      fill = NULL,
      title = paste0(gsub("[(){}]", "", ct2), " -> ", gsub("[(){}]", "", ct1)),
      subtitle = NULL
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = title_text_size),
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, right_margin, 4, 20)
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        override.aes = list(shape = NA, linetype = 1, linewidth = 1.2)
      )
    ) +
    ggplot2::coord_cartesian(
      clip = "off",
      ylim = c(0.5, max(plot_data$y_pos) + 1.5)
    )

  p
}



#' Construct spatial colocalization network from PANORAMIC results
#'
#' Build a cell-type network where edges represent differential
#' spatial colocalization between cell-type pairs, based on the output of
#' \code{panoramic_meta_mv(...)} (or a compatible function). Edges are
#' filtered by FDR and optional \code{z_diff} sign, and weighted by
#' \code{|z_diff|}.
#'
#' @param se_diff A SummarizedExperiment returned by
#'   \code{panoramic_meta_mv(...)}, where
#'   rowData(se_diff) contains
#'   at least ct1, ct2, z_diff, p_diff, and fdr_diff.
#' @param fdr_threshold Numeric scalar. FDR threshold used to define
#'   significance for \code{edge_class} and, when
#'   \code{include_nonsig = FALSE}, to filter edges (default 0.05).
#' @param leiden_resolution Numeric. Clustering resolution for the Leiden algorithm.
#' @param z_sign Character sign filter for \code{z_diff}, one of
#'   \code{"both"} (default), \code{"positive"}, or \code{"negative"}.
#'   This enables group-specific networks by retaining only positive or
#'   only negative differential effects.
#' @param include_nonsig Logical. If TRUE, include non-significant edges
#'   (up to \code{nonsig_max_fdr}) and mark them as dotted in plotting.
#'   If FALSE (default), keep only \code{fdr_diff < fdr_threshold}.
#' @param nonsig_max_fdr Numeric upper bound for retained edges when
#'   \code{include_nonsig = TRUE} (default 1.0).
#' @param directed Logical. If \code{TRUE}, construct a directed network
#'  (and directed centrality metrics). If \code{FALSE}, construct an
#'  undirected network.
#' @param sig_operator One of \code{"lt"} (default) or \code{"gt"} controlling
#'  significance threshold direction for \code{fdr_diff} when
#'  \code{include_nonsig = FALSE}. \code{"lt"} uses standard
#'  \code{fdr_diff < fdr_threshold}; \code{"gt"} inverts the threshold
#'  (diagnostic use only).
#'
#' @return A list with components:
#' \itemize{
#'   \item graph: an igraph object (directed or undirected) with edge attributes
#'     \code{weight}, \code{fdr}, \code{pval}, \code{edge_sig},
#'     \code{edge_class} and vertex attributes
#'     \code{cluster}, \code{cluster_id}, \code{degree},
#'     \code{betweenness}, \code{strength}.
#'   \item clusters: the community structure object from
#'     \code{igraph::cluster_leiden()}.
#'   \item n_clusters: number of detected clusters.
#'   \item modularity: modularity score of the clustering.
#'   \item z_sign: the applied sign filter.
#'   \item fdr_threshold: the applied significance threshold.
#'   \item include_nonsig: whether non-significant edges were included.
#'   \item nonsig_max_fdr: maximum retained FDR when non-significant edges
#'     are included.
#' }
#'
#' @examples
#' se_diff <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(dummy = matrix(0, nrow = 3, ncol = 1)),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = c("A", "A", "B"),
#'     ct2 = c("B", "C", "C"),
#'     z_diff = c(2.0, -1.5, 1.2),
#'     p_diff = c(0.01, 0.03, 0.20),
#'     fdr_diff = c(0.02, 0.05, 0.30)
#'   )
#' )
#'
#' net <- create_spatial_network(
#'   se_diff,
#'   fdr_threshold = 0.2,
#'   z_sign = "both",
#'   include_nonsig = TRUE
#' )
#'
#' net$graph
#' @export
create_spatial_network <- function(se_diff,
                                   fdr_threshold = 0.05,
                                   leiden_resolution = 1.0,
                                   z_sign = c("both", "positive", "negative"),
                                   include_nonsig = FALSE,
                                   nonsig_max_fdr = 1.0,
                                   directed = FALSE,
                                   sig_operator = c("lt", "gt")) {
  if (!inherits(se_diff, "SummarizedExperiment")) {
    stop("`se_diff` must be a SummarizedExperiment.")
  }
  if (!is.finite(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1) {
    stop("`fdr_threshold` must be a finite number in [0, 1].")
  }
  if (!is.finite(leiden_resolution) || leiden_resolution <= 0) {
    stop("`leiden_resolution` must be a positive finite number.")
  }
  z_sign <- match.arg(z_sign)
  sig_operator <- match.arg(sig_operator)
  if (identical(sig_operator, "gt")) {
    warning(
      "`sig_operator = \"gt\"` inverts the usual FDR rule and treats larger FDR values as significant. Use only for diagnostic/inverted-threshold views.",
      call. = FALSE
    )
  }
  if (!is.finite(nonsig_max_fdr) || nonsig_max_fdr <= 0) nonsig_max_fdr <- 1.0
  nonsig_max_fdr <- min(nonsig_max_fdr, 1.0)

  results_df <- as.data.frame(SummarizedExperiment::rowData(se_diff))
  required_cols <- c("ct1", "ct2", "z_diff", "p_diff", "fdr_diff")
  missing_cols <- setdiff(required_cols, colnames(results_df))
  if (length(missing_cols) > 0L) {
    stop(
      "rowData(se_diff) is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "."
    )
  }

  edges_df <- results_df |>
    dplyr::filter(
      !is.na(ct1),
      !is.na(ct2),
      is.finite(fdr_diff),
      is.finite(z_diff)
    )
  if (identical(z_sign, "positive")) {
    edges_df <- dplyr::filter(edges_df, z_diff > 0)
  } else if (identical(z_sign, "negative")) {
    edges_df <- dplyr::filter(edges_df, z_diff < 0)
  }
  if (isTRUE(include_nonsig)) {
    edges_df <- dplyr::filter(edges_df, fdr_diff <= nonsig_max_fdr)
  } else {
    if (identical(sig_operator, "gt")) {
      edges_df <- dplyr::filter(edges_df, fdr_diff > fdr_threshold)
    } else {
      edges_df <- dplyr::filter(edges_df, fdr_diff < fdr_threshold)
    }
  }

  sig_edge <- if (identical(sig_operator, "gt")) {
    function(x, thr) as.integer(x > thr)
  } else {
    function(x, thr) as.integer(x < thr)
  }
  fdr_comb <- if (identical(sig_operator, "gt")) "max" else "min"
  pval_comb <- if (identical(sig_operator, "gt")) "max" else "min"

  edges <- edges_df |>
    dplyr::mutate(
      from = gsub("[(){}|]", "", ct1),
      to   = gsub("[(){}|]", "", ct2),
      weight = abs(z_diff),
      fdr    = fdr_diff,
      pval   = p_diff,
      edge_sig = sig_edge(fdr_diff, fdr_threshold)
    ) |>
    dplyr::select(from, to, weight, fdr, pval, edge_sig)

  if (nrow(edges) == 0) {
    stop(
      "No edges found for network at FDR/sign filter with z_sign = '", z_sign, "'."
    )
  }

  g_net <- igraph::graph_from_data_frame(edges, directed = directed)
  g_net <- igraph::simplify(
    g_net,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "max", fdr = fdr_comb, pval = pval_comb, edge_sig = "max")
  )
  igraph::E(g_net)$edge_class <- ifelse(igraph::E(g_net)$edge_sig >= 1, "Significant", "Not significant")

  # Cluster on an undirected projection to keep community detection stable.
  g_cluster <- if (isTRUE(directed)) {
    igraph::as.undirected(
      g_net,
      mode = "collapse",
      edge.attr.comb = list(weight = "max", fdr = fdr_comb, pval = pval_comb, edge_sig = "max")
    )
  } else {
    g_net
  }
  clusters <- igraph::cluster_leiden(
    g_cluster,
    objective_function = "modularity",
    weights = igraph::E(g_cluster)$weight,
    resolution = leiden_resolution
  )

  mod <- igraph::modularity(g_cluster, membership = clusters$membership,
                            weights = igraph::E(g_cluster)$weight)

  igraph::V(g_net)$cluster <- as.factor(clusters$membership)
  igraph::V(g_net)$cluster_id <- clusters$membership
  igraph::V(g_net)$degree <- igraph::degree(g_net, mode = "all")
  igraph::V(g_net)$betweenness <- igraph::betweenness(
    g_net,
    directed = directed,
    weights = 1 / pmax(igraph::E(g_net)$weight, 1e-8)
  )
  igraph::V(g_net)$strength <- igraph::strength(g_net, mode = "all", weights = igraph::E(g_net)$weight)

  list(
    graph = g_net,
    clusters = clusters,
    n_clusters = max(clusters$membership),
    modularity = mod,
    z_sign = z_sign,
    fdr_threshold = fdr_threshold,
    include_nonsig = include_nonsig,
    nonsig_max_fdr = nonsig_max_fdr,
    directed = directed,
    sig_operator = sig_operator
  )
}


#' Plot PANORAMIC spatial colocalization network
#'
#' Wrapper to construct and visualize a PANORAMIC spatial colocalization network.
#' It builds the directed network via \code{create_spatial_network()} and
#' then renders it with ggraph/tidygraph.
#'
#' @param se_diff A SummarizedExperiment returned by
#'   \code{panoramic_meta_mv()}.
#' @param fdr_threshold Numeric threshold used for significance and
#'   edge classification.
#' @param leiden_resolution Numeric Leiden resolution passed to
#'   \code{create_spatial_network()}.
#' @param z_sign Character sign filter for \code{z_diff}
#'   (\code{"both"}, \code{"positive"}, \code{"negative"}).
#' @param include_nonsig Logical; include non-significant edges.
#' @param nonsig_max_fdr Numeric max FDR retained when
#'   \code{include_nonsig = TRUE}.
#' @param directed Logical; pass through to \code{create_spatial_network()}.
#' @param sig_operator One of \code{"lt"} or \code{"gt"} passed through to
#'  \code{create_spatial_network()} for threshold direction.
#' @param layout Character string specifying the graph layout passed to
#'   ggraph (e.g. "fr", "kk", "stress"). Default "fr".
#' @param node_size_by Character name of a vertex attribute used to scale
#'   node sizes (e.g. "degree", "strength", "betweenness").
#'   Default "degree".
#' @param label_repel Logical; if TRUE use repelled node labels.
#' @param label_box_padding Numeric box padding for repelled labels.
#' @param label_point_padding Numeric point padding for repelled labels.
#' @param label_force Numeric repulsion force for labels.
#' @param return_net Logical; if TRUE return a list with
#'   \code{plot} and \code{net}.
#'
#' @return A ggplot object by default, or a list with \code{plot} and
#'   \code{net} when \code{return_net = TRUE}.
#' 
#' @examples
#' se_diff <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(dummy = matrix(0, nrow = 3, ncol = 1)),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = c("A", "A", "B"),
#'     ct2 = c("B", "C", "C"),
#'     z_diff = c(2.0, -1.5, 1.2),
#'     p_diff = c(0.01, 0.03, 0.20),
#'     fdr_diff = c(0.02, 0.05, 0.30)
#'   )
#' )
#'
#' if (requireNamespace("ggraph", quietly = TRUE) &&
#'     requireNamespace("tidygraph", quietly = TRUE)) {
#'   p_net <- plot_spatial_network(
#'     se_diff,
#'     fdr_threshold = 0.2,
#'     z_sign = "both",
#'     include_nonsig = TRUE,
#'     layout = "fr",
#'     node_size_by = "degree"
#'   )
#'   print(p_net)
#' }
#' @export
plot_spatial_network <- function(se_diff,
                                 fdr_threshold = 0.05,
                                 leiden_resolution = 1.0,
                                 z_sign = c("both", "positive", "negative"),
                                 include_nonsig = FALSE,
                                 nonsig_max_fdr = 1.0,
                                 directed = FALSE,
                                 sig_operator = c("lt", "gt"),
                                 layout = "fr",
                                 node_size_by = "degree",
                                 label_repel = TRUE,
                                 label_box_padding = 0.2,
                                 label_point_padding = 0.1,
                                 label_force = 0.8,
                                 return_net = FALSE) {
  net_result <- create_spatial_network(
    se_diff = se_diff,
    fdr_threshold = fdr_threshold,
    leiden_resolution = leiden_resolution,
    z_sign = z_sign,
    include_nonsig = include_nonsig,
    nonsig_max_fdr = nonsig_max_fdr,
    directed = directed,
    sig_operator = sig_operator
  )
  g <- net_result$graph

  edge_list <- igraph::as_data_frame(g, what = "edges")
  node_list <- igraph::as_data_frame(g, what = "vertices")

  tg <- tidygraph::tbl_graph(nodes = node_list, edges = edge_list, directed = directed)

  p <- ggraph::ggraph(tg, layout = layout) +
    {
      if (isTRUE(directed)) {
        ggraph::geom_edge_link(
          ggplot2::aes(edge_width = weight, edge_linetype = edge_class),
          color = "gray50",
          arrow = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed"),
          end_cap = ggraph::circle(2.2, "mm")
        )
      } else {
        ggraph::geom_edge_link(
          ggplot2::aes(edge_width = weight, edge_linetype = edge_class),
          color = "gray50"
        )
      }
    } +
    ggraph::scale_edge_width_continuous(
      name = "Edge weight (|Z-score|)",
      range = c(0.5, 3)
    ) +
    ggraph::scale_edge_linetype_manual(
      name = "Edge",
      values = c("Significant" = "solid", "Not significant" = "dotted")
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = !!rlang::sym(node_size_by), color = cluster),
      alpha = 0.9
    ) +
    {
      if (isTRUE(label_repel)) {
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          size = 3.5,
          fontface = "bold",
          repel = TRUE,
          box.padding = label_box_padding,
          point.padding = label_point_padding,
          force = label_force
        )
      } else {
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          size = 3.5,
          fontface = "bold",
          repel = FALSE
        )
      }
    } +
    ggplot2::scale_color_brewer(palette = "Set2", name = "Cluster") +
    ggplot2::scale_size_continuous(range = c(6, 20), guide = "none") +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::guides(
      edge_width = ggplot2::guide_legend(order = 1),
      edge_linetype = ggplot2::guide_legend(order = 2),
      color = ggplot2::guide_legend(order = 3)
    )

  if (isTRUE(return_net)) {
    return(list(plot = p, net = net_result))
  }
  p
}
