#' Volcano plot for differential spatial colocalization
#'
#' Create a volcano plot summarizing differential spatial colocalization
#' results from PANORAMIC. Points represent cell-type pairs, with the
#' x-axis showing effect sizes and the y-axis showing -log10(p)
#' from \code{p_diff}. Significant pairs can be highlighted and labeled.
#'
#' @param se_diff A SummarizedExperiment returned by
#'   \code{pano_compare_groups()}, where rowData(se_diff) contains at least
#'   ct1, ct2, beta_diff, p_diff, and fdr_diff.
#' @param fdr_threshold Numeric scalar. FDR threshold used to define
#'   significance for highlighting and labeling. Default is 0.05.
#' @param effect_threshold Numeric scalar. Absolute effect size threshold
#'   |beta| used to highlight strong effects. Default is 5.
#' @param label_top Optional integer. If non-NULL, only the top
#'   label_top pairs (smallest p_diff) among FDR-significant
#'   features are labeled. If NULL (default), all FDR-significant
#'   pairs are labeled.
#' @param title Optional character string for a custom plot title. If
#'   NULL, a title is constructed from metadata(se_diff)$panoramic
#'   if available.
#'
#' @return A ggplot object representing the volcano plot.
#' @export
plot_volcano <- function(se_diff,
                         fdr_threshold = 0.05,
                         effect_threshold = 5,
                         label_top = NULL,
                         title = NULL) {

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(se_diff))

  required_cols <- c("ct1", "ct2", "beta_diff", "p_diff", "fdr_diff")
  if (!all(required_cols %in% colnames(rd))) {
    stop("se_diff must have columns: ", paste(required_cols, collapse = ", "))
  }

  plot_data <- dplyr::mutate(
    rd,
    pair_label = paste0(
      gsub("[(){}]", "", ct1),
      ":",
      gsub("[(){}]", "", ct2)
    ),
    neg_log10_p = -log10(p_diff),
    significance = dplyr::case_when(
      fdr_diff < 0.01 & abs(beta_diff) > effect_threshold ~
        paste0("FDR < 0.01, |beta| > ", effect_threshold),
      fdr_diff < fdr_threshold & abs(beta_diff) > effect_threshold ~
        paste0("FDR < ", fdr_threshold, ", |beta| > ", effect_threshold),
      fdr_diff < fdr_threshold ~ paste0("FDR < ", fdr_threshold),
      p_diff < 0.05 ~ "p < 0.05",
      TRUE ~ "Not significant"
    ),
    direction = ifelse(beta_diff > 0, "Increased", "Decreased")
  )
  plot_data <- dplyr::filter(plot_data, !is.na(p_diff) & !is.na(beta_diff))

  if (!is.null(label_top)) {
    label_data <- plot_data |>
      dplyr::filter(fdr_diff < fdr_threshold) |>
      dplyr::arrange(p_diff) |>
      dplyr::slice_head(n = label_top)
  } else {
    label_data <- dplyr::filter(plot_data, fdr_diff < fdr_threshold)
  }

  meta <- S4Vectors::metadata(se_diff)$panoramic
  radius <- if (!is.null(rd$radius_um)) unique(rd$radius_um)[1] else "unknown"

  if (is.null(title)) {
    if (!is.null(meta$contrast)) {
      title <- paste0(
        "Differential spatial colocalization: ",
        meta$contrast$case, " vs. ", meta$contrast$control
      )
    } else {
      title <- "Differential spatial colocalization"
    }
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = beta_diff, y = neg_log10_p)) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray50",
      alpha = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = c(-effect_threshold, effect_threshold),
      linetype = "dashed",
      color = "gray50",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = significance, size = abs(beta_diff)),
      alpha = 0.7
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = pair_label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FDR < 0.01, |beta| > 5" = "#d62728",
        "FDR < 0.05, |beta| > 5" = "#ff7f0e",
        "FDR < 0.05" = "#2ca02c",
        "p < 0.05" = "#1f77b4",
        "Not significant" = "gray70"
      ),
      breaks = c(
        paste0("FDR < 0.01, |beta| > ", effect_threshold),
        paste0("FDR < ", fdr_threshold, ", |beta| > ", effect_threshold),
        paste0("FDR < ", fdr_threshold),
        "p < 0.05",
        "Not significant"
      )
    ) +
    ggplot2::scale_size_continuous(range = c(1, 5)) +
    ggplot2::labs(
      x = "Effect size (beta: case - control)",
      y = expression(-log[10](italic(p) * "-value")),
      color = "Significance",
      size = "|Effect size|",
      title = title,
      subtitle = paste0("Radius: ", radius, " um")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold")
    )

  p
}


#' Forest plot for PANORAMIC spatial colocalization
#'
#' Create a forest plot showing individual-sample spatial statistics and
#' pooled random-effects estimates for a single cell-type pair at a given
#' radius. Optionally annotates heterogeneity statistics (tau2, I2) for
#' the pooled effects.
#'
#' @param se_meta A SummarizedExperiment returned by
#'   \code{panoramic_meta()} or a downstream testing function, containing
#'   assays "yi" (effect estimates) and "vi" (variances), and rowData
#'   columns such as ct1, ct2, radius_um, and, if group-level pooling was
#'   performed, mu_control, se_control, mu_case, se_case, and optionally
#'   tau2_control, tau2_case, I2.
#' @param ct1 Character. First cell type (must match rowData(se_meta)$ct1).
#' @param ct2 Character. Second cell type (must match rowData(se_meta)$ct2).
#' @param radius_um Numeric radius (in microns) to plot. If NULL,
#'   the first available radius for this pair is used.
#' @param group_col Character. Column name in colData(se_meta)
#'   containing group labels (e.g. treatment vs control). Default "group".
#' @param group_colors Optional named character vector mapping group labels
#'   to colors. If NULL, a default palette is used.
#' @param show_heterogeneity Logical. If TRUE, show tau2 annotations
#'   for pooled effects when available. Default TRUE.
#'
#' @return A ggplot object representing the forest plot.
#' @export
plot_forest <- function(se_meta, ct1, ct2, radius_um = NULL,
                        group_col = "group",
                        group_colors = NULL,
                        show_heterogeneity = TRUE) {

  yi <- SummarizedExperiment::assay(se_meta, "yi")
  vi <- SummarizedExperiment::assay(se_meta, "vi")
  rd <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  cd <- as.data.frame(SummarizedExperiment::colData(se_meta))

  if (is.null(radius_um)) {
    radius_um <- rd$radius_um[1]
    message("Using radius: ", radius_um, " um")
  }

  row_idx <- which(rd$ct1 == ct1 & rd$ct2 == ct2 & rd$radius_um == radius_um)
  if (length(row_idx) == 0) {
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

  pooled_data <- NULL
  tau2_control <- NA
  tau2_case <- NA
  I2_control <- NA
  I2_case <- NA

  if (all(c("mu_control", "se_control", "mu_case", "se_case") %in% colnames(rd))) {
    groups <- unique(cd[[group_col]])

    if ("tau2_control" %in% colnames(rd)) {
      tau2_control <- rd$tau2_control[row_idx]
      tau2_case <- rd$tau2_case[row_idx]
    } else if ("tau2" %in% colnames(rd)) {
      tau2_control <- tau2_case <- rd$tau2[row_idx]
    }

    if ("I2" %in% colnames(rd)) {
      I2_control <- I2_case <- rd$I2[row_idx]
    }

    pooled_data <- data.frame(
      sample = c("Pooled (RE)", "Pooled (RE)"),
      estimate = c(rd$mu_control[row_idx], rd$mu_case[row_idx]),
      se = c(rd$se_control[row_idx], rd$se_case[row_idx]),
      group = c(groups[1], groups[2]),
      type = "Pooled",
      tau2 = c(tau2_control, tau2_case),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        ci_lower = estimate - 1.96 * se,
        ci_upper = estimate + 1.96 * se,
        label = "Pooled (RE)",
        stats_text = sprintf("%.2f (%.2f)", estimate, se),
        tau2_text = ifelse(!is.na(tau2), sprintf("tau2 = %.2f", tau2), "")
      )
  }

  plot_data <- dplyr::bind_rows(sample_data, pooled_data) |>
    dplyr::arrange(group, dplyr::desc(type), estimate) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      row_within_group = dplyr::row_number()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      group_offset = ifelse(group == unique(group)[1], 0,
                            max(row_within_group[group == unique(group)[1]]) + 2),
      y_pos = row_within_group + group_offset
    )

  if (is.null(group_colors)) {
    groups <- unique(plot_data$group)
    group_colors <- stats::setNames(c("#1f77b4", "#d62728"), groups)
  }

  x_min <- min(plot_data$ci_lower, na.rm = TRUE)
  x_max <- max(plot_data$ci_upper, na.rm = TRUE)
  x_range <- x_max - x_min

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = y_pos)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lower, xend = ci_upper,
                   yend = y_pos, color = group),
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, type == "Individual"),
      ggplot2::aes(x = estimate, color = group),
      size = 2.5, shape = 15
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, type == "Pooled"),
      ggplot2::aes(x = estimate, color = group),
      size = 4.5, shape = 18
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_min - x_range * 0.02,
                   label = label,
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 1, size = 3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_max + x_range * 0.02,
                   label = stats_text,
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 0, size = 3, family = "mono"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_max + x_range * 0.20,
                   label = sprintf("[%.2f, %.2f]", ci_lower, ci_upper),
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 0, size = 3, family = "mono"
    ) +
    {
      if (show_heterogeneity && !is.null(pooled_data)) {
        ggplot2::geom_text(
          data = subset(plot_data, type == "Pooled" & tau2_text != ""),
          ggplot2::aes(x = x_max + x_range * 0.38,
                       label = tau2_text,
                       color = group),
          hjust = 0, size = 2.8, fontface = "italic"
        )
      }
    } +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray30", alpha = 0.5) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.50))
    ) +
    ggplot2::annotate("text", x = x_min - x_range * 0.02,
                      y = max(plot_data$y_pos) + 1,
                      label = "Sample", hjust = 1,
                      fontface = "bold", size = 3.5) +
    ggplot2::annotate("text", x = x_max + x_range * 0.02,
                      y = max(plot_data$y_pos) + 1,
                      label = "Est (SE)", hjust = 0,
                      fontface = "bold", size = 3.5) +
    ggplot2::annotate("text", x = x_max + x_range * 0.20,
                      y = max(plot_data$y_pos) + 1,
                      label = "95% CI", hjust = 0,
                      fontface = "bold", size = 3.5) +
    {
      if (show_heterogeneity && !is.null(pooled_data)) {
        ggplot2::annotate("text", x = x_max + x_range * 0.38,
                          y = max(plot_data$y_pos) + 1,
                          label = "Heterogeneity", hjust = 0,
                          fontface = "bold", size = 3.5)
      }
    } +
    ggplot2::labs(
      x = paste0("L-cross statistic at ", radius_um, " um"),
      y = NULL,
      color = NULL,
      title = paste0(gsub("[(){}]", "", ct1), " : ",
                     gsub("[(){}]", "", ct2)),
      subtitle = "Random-effects meta-analysis with individual sample estimates"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      plot.margin = ggplot2::margin(10, 120, 10, 40)
    ) +
    ggplot2::coord_cartesian(
      clip = "off",
      ylim = c(0.5, max(plot_data$y_pos) + 1.5)
    )

  p
}

#' Construct spatial colocalization network from PANORAMIC results
#'
#' Build a cell-type network where edges represent significant differential
#' spatial colocalization between cell-type pairs, based on the output of
#' \code{pano_compare_groups()} (or a compatible function). Negative
#' \code{z_diff} values (reduced colocalization in the case group) are
#' encoded as edges, weighted by |z|.
#'
#' @param se_diff A SummarizedExperiment returned by
#'   \code{pano_compare_groups()}, where rowData(se_diff) contains
#'   at least ct1, ct2, z_diff, p_diff, and fdr_diff.
#' @param fdr_threshold Numeric scalar. FDR threshold used to define
#'   significant edges (default 0.05).
#' @param directed Logical. If TRUE, keep edges directed as
#'   (ct1 -> ct2). If FALSE (default), create a symmetric
#'   undirected network by adding reversed edges before simplification.
#' @param leiden_resolution Numeric. Clustering resolution for the Leiden algorithm 
#'
#' @return A list with components:
#' \itemize{
#'   \item graph: an igraph object with edge attributes
#'     weight, fdr, pval and vertex attributes
#'     cluster, cluster_id, degree, betweenness, strength.
#'   \item clusters: the community structure object from
#'     \code{igraph::cluster_leiden()}.
#'   \item n_clusters: number of detected clusters.
#'   \item modularity: modularity score of the clustering.
#' }
#' @export
create_spatial_network <- function(se_diff,
                                   fdr_threshold = 0.05,
                                   directed = FALSE, 
                                   leiden_resolution = 1.0) {

  results_df <- as.data.frame(SummarizedExperiment::rowData(se_diff))

  edges <- results_df |>
    dplyr::filter(
      fdr_diff < fdr_threshold
    ) |>
    dplyr::mutate(
      from = gsub("[(){}|]", "", ct1),
      to   = gsub("[(){}|]", "", ct2),
      weight = abs(z_diff),
      fdr    = fdr_diff,
      pval   = p_diff
    ) |>
    dplyr::select(from, to, weight, fdr, pval)

  if (nrow(edges) == 0) {
    stop("No significant edges found at FDR < ", fdr_threshold)
  }

  if (!directed) {
    edges_rev <- edges |>
      dplyr::rename(from = to, to = from)
    edges <- dplyr::bind_rows(edges, edges_rev)
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  g <- igraph::simplify(
    g,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "max", fdr = "min", pval = "min")
  )

  set.seed(123)
  clusters <- igraph::cluster_leiden(
    g,
    objective_function   = "modularity",
    weights              = igraph::E(g)$weight,
    resolution = leiden_resolution
  )

  mod <- igraph::modularity(g, membership = clusters$membership,
                            weights = igraph::E(g)$weight)

  igraph::V(g)$cluster    <- as.factor(clusters$membership)
  igraph::V(g)$cluster_id <- clusters$membership
  igraph::V(g)$degree     <- igraph::degree(g)
  igraph::V(g)$betweenness <- igraph::betweenness(g, weights = 1 / igraph::E(g)$weight)
  igraph::V(g)$strength    <- igraph::strength(g,  weights = igraph::E(g)$weight)

  list(
    graph     = g,
    clusters  = clusters,
    n_clusters = max(clusters$membership),
    modularity = mod
  )
}


#' Plot PANORAMIC spatial colocalization network
#'
#' Visualize the spatial colocalization network constructed by
#' \code{create_spatial_network()} using ggraph and tidygraph.
#' Nodes represent cell types and edges represent significant differential
#' colocalization, with edge width proportional to |z| and node size
#' based on a chosen centrality measure.
#'
#' @param net_result A list returned by \code{create_spatial_network()},
#'   containing at least a graph element (an igraph object)
#'   and summary statistics such as n_clusters and modularity.
#' @param layout Character string specifying the graph layout passed to
#'   ggraph (e.g. "fr", "kk", "stress"). Default "fr".
#' @param node_size_by Character name of a vertex attribute used to scale
#'   node sizes (e.g. "degree", "strength", "betweenness").
#'   Default "degree".
#'
#' @return A ggplot object representing the network plot.
#' @export
plot_spatial_network <- function(net_result,
                                 layout = "fr",
                                 node_size_by = "degree") {

  g <- net_result$graph

  edge_list <- igraph::as_data_frame(g, what = "edges")
  node_list <- igraph::as_data_frame(g, what = "vertices")

  tg <- tidygraph::tbl_graph(nodes = node_list, edges = edge_list, directed = FALSE)

  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_width = weight, edge_alpha = fdr),
      color = "gray50"
    ) +
    ggraph::scale_edge_width_continuous(
      name = "|z-score|",
      range = c(0.5, 3)
    ) +
    ggraph::scale_edge_alpha_continuous(
      name = "FDR",
      range = c(0.9, 0.3),
      trans = "reverse"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = !!rlang::sym(node_size_by), color = cluster),
      alpha = 0.9
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      size = 3.5,
      fontface = "bold",
      repel = TRUE,
      box.padding = 0.5
    ) +
    ggplot2::scale_color_brewer(palette = "Set2", name = "Cluster") +
    ggplot2::scale_size_continuous(range = c(6, 20), name = "Degree") +
    ggplot2::labs(
      title = "Spatial Colocalization Network",
      subtitle = sprintf(
        "%d clusters | Q = %.3f | %d nodes, %d edges",
        net_result$n_clusters, net_result$modularity,
        igraph::vcount(g), igraph::ecount(g)
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    )

  p
}