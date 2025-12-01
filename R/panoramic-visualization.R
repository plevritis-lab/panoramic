#' Volcano plot for differential spatial colocalization
#'
#' Create a volcano plot summarizing differential spatial colocalization
#' results from PANORAMIC. Points represent cell-type pairs, with the
#' x-axis showing effect sizes and the y-axis showing \eqn{-\log_{10}(p)}
#' from \code{p_diff}. Significant pairs can be highlighted and labeled.
#'
#' @param se_diff A \code{SummarizedExperiment} returned by
#'   \code{pano_compare_groups()} (or an analogous function), where
#'   \code{rowData(se_diff)} contains at least \code{ct1}, \code{ct2},
#'   \code{beta_diff}, \code{p_diff}, and \code{fdr_diff}.
#' @param fdr_threshold Numeric scalar. FDR threshold used to define
#'   significance for highlighting and labeling. Default is \code{0.05}.
#' @param effect_threshold Numeric scalar. Absolute effect size threshold
#'   \eqn{|β|} used to highlight strong effects. Default is \code{5}.
#' @param label_top Optional integer. If non-NULL, only the top
#'   \code{label_top} pairs (smallest \code{p_diff}) among FDR-significant
#'   features are labeled. If \code{NULL} (default), all FDR-significant
#'   pairs are labeled.
#' @param title Optional character string for a custom plot title. If
#'   \code{NULL}, a title is constructed from \code{metadata(se_diff)$panoramic}
#'   if available.
#'
#' @return A \code{ggplot} object representing the volcano plot, which can
#'   be further customized or printed.
#'
#' @details
#' The function expects that PANORAMIC differential testing has already
#' been performed and stored in \code{rowData(se_diff)}. It categorizes
#' each feature by FDR and effect size, and overlays reference lines at
#' \eqn{p = 0.05} and \eqn{|β| =} \code{effect_threshold}. Labels are
#' optionally added using \code{ggrepel} for readability.
#'
#' @examples
#' \dontrun{
#' p <- plot_volcano(se_diff, fdr_threshold = 0.05, effect_threshold = 3)
#' print(p)
#' }
#'
#' @export
plot_volcano <- function(se_diff, 
                         fdr_threshold = 0.05,
                         effect_threshold = 5,
                         label_top = NULL,
                         title = NULL) {
  
  # Check required packages
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  }
  
  # Extract data from SummarizedExperiment
  rd <- as.data.frame(SummarizedExperiment::rowData(se_diff))
  
  # Check required columns
  required_cols <- c("ct1", "ct2", "beta_diff", "p_diff", "fdr_diff")
  if (!all(required_cols %in% colnames(rd))) {
    stop("se_diff must have columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Prepare plot data
  plot_data <- rd %>%
    dplyr::mutate(
      # Create clean label for cell type pairs
      pair_label = paste0(
        gsub("[(){}]", "", ct1), 
        ":", 
        gsub("[(){}]", "", ct2)
      ),
      # -log10(p-value) for y-axis
      neg_log10_p = -log10(p_diff),
      # Significance categories
      significance = dplyr::case_when(
        fdr_diff < 0.01 & abs(beta_diff) > effect_threshold ~ 
          paste0("FDR < 0.01, |β| > ", effect_threshold),
        fdr_diff < fdr_threshold & abs(beta_diff) > effect_threshold ~ 
          paste0("FDR < ", fdr_threshold, ", |β| > ", effect_threshold),
        fdr_diff < fdr_threshold ~ paste0("FDR < ", fdr_threshold),
        p_diff < 0.05 ~ "p < 0.05",
        TRUE ~ "Not significant"
      ),
      # Direction (not currently mapped to aesthetics but available)
      direction = ifelse(beta_diff > 0, "Increased", "Decreased")
    ) %>%
    dplyr::filter(!is.na(p_diff) & !is.na(beta_diff))
  
  # Determine which pairs to label
  if (!is.null(label_top)) {
    label_data <- plot_data %>%
      dplyr::filter(fdr_diff < fdr_threshold) %>%
      dplyr::arrange(p_diff) %>%
      dplyr::slice_head(n = label_top)
  } else {
    label_data <- plot_data %>%
      dplyr::filter(fdr_diff < fdr_threshold)
  }
  
  # Get metadata for subtitle
  meta <- S4Vectors::metadata(se_diff)$panoramic
  radius <- if (!is.null(rd$radius_um)) unique(rd$radius_um)[1] else "unknown"
  
  # Default title
  if (is.null(title)) {
    if (!is.null(meta$contrast)) {
      title <- paste0("Differential spatial colocalization: ",
                      meta$contrast$case, " vs. ", meta$contrast$control)
    } else {
      title <- "Differential spatial colocalization"
    }
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = beta_diff, y = neg_log10_p)) +
    # Reference lines
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                        color = "gray50", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(-effect_threshold, effect_threshold), 
                        linetype = "dashed", color = "gray50", alpha = 0.5) +
    # Points
    ggplot2::geom_point(ggplot2::aes(color = significance, size = abs(beta_diff)), 
                        alpha = 0.7) +
    # Labels
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
    # Color scheme
    ggplot2::scale_color_manual(
      values = c(
        "FDR < 0.01, |β| > 5" = "#d62728",
        "FDR < 0.05, |β| > 5" = "#ff7f0e",
        "FDR < 0.05" = "#2ca02c",
        "p < 0.05" = "#1f77b4",
        "Not significant" = "gray70"
      ),
      breaks = c(
        paste0("FDR < 0.01, |β| > ", effect_threshold),
        paste0("FDR < ", fdr_threshold, ", |β| > ", effect_threshold),
        paste0("FDR < ", fdr_threshold),
        "p < 0.05",
        "Not significant"
      )
    ) +
    ggplot2::scale_size_continuous(range = c(1, 5)) +
    # Labels
    ggplot2::labs(
      x = "Effect size (β: case - control)",
      y = expression(-log[10](italic(p)*"-value")),
      color = "Significance",
      size = "|Effect size|",
      title = title,
      subtitle = paste0("Radius: ", radius, " μm")
    ) +
    # Theme
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}

#' Forest plot for PANORAMIC spatial colocalization
#'
#' Create a forest plot showing individual-sample spatial statistics and
#' pooled random-effects estimates for a single cell-type pair at a given
#' radius. Optionally annotates heterogeneity statistics (tau², I²) for
#' the pooled effects.
#'
#' @param se_meta A \code{SummarizedExperiment} returned by
#'   \code{panoramic_meta()} or a downstream testing function, containing
#'   assays \code{"yi"} (effect estimates) and \code{"vi"} (variances), and
#'   rowData columns such as \code{ct1}, \code{ct2}, \code{radius_um},
#'   and, if group-level pooling was performed, \code{mu_control},
#'   \code{se_control}, \code{mu_case}, \code{se_case}, and optionally
#'   \code{tau2_control}, \code{tau2_case}, \code{I2}.
#' @param ct1 Character. First cell type (must match \code{rowData(se_meta)$ct1}).
#' @param ct2 Character. Second cell type (must match \code{rowData(se_meta)$ct2}).
#' @param radius_um Numeric radius (in microns) to plot. If \code{NULL},
#'   the first available radius for this pair is used.
#' @param group_col Character. Column name in \code{colData(se_meta)}
#'   containing group labels (e.g. treatment vs control). Default \code{"group"}.
#' @param group_colors Optional named character vector mapping group labels
#'   to colors. If \code{NULL}, a default palette is used.
#' @param show_heterogeneity Logical. If \code{TRUE}, show tau² annotations
#'   for pooled effects when available. Default \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the forest plot, which can
#'   be further customized or printed.
#'
#' @details
#' For the specified cell-type pair and radius, the function extracts
#' individual-sample estimates and 95% confidence intervals from
#' \code{assay(se_meta, "yi")} and \code{assay(se_meta, "vi")}, grouped
#' by \code{group_col}. If pooled random-effects estimates are present
#' in \code{rowData(se_meta)} (e.g. \code{mu_control}, \code{mu_case}),
#' these are shown as group-level diamonds with optional tau² annotations.
#'
#' The x-axis is centered at zero, representing no deviation from the
#' theoretical spatial summary (e.g. L(r) - r). The y-axis lists samples
#' within each group, with pooled estimates shown at the bottom of each group.
#'
#' @examples
#' \dontrun{
#' p <- plot_forest(se_meta, ct1 = "T cells", ct2 = "B cells", radius_um = 50)
#' print(p)
#' }
#'
#' @export
plot_forest <- function(se_meta, ct1, ct2, radius_um = NULL, 
                        group_col = "group",
                        group_colors = NULL,
                        show_heterogeneity = TRUE) {
  
  # Get data
  yi <- SummarizedExperiment::assay(se_meta, "yi")
  vi <- SummarizedExperiment::assay(se_meta, "vi")
  rd <- as.data.frame(SummarizedExperiment::rowData(se_meta))
  cd <- as.data.frame(SummarizedExperiment::colData(se_meta))
  
  # Select radius if not specified
  if (is.null(radius_um)) {
    radius_um <- rd$radius_um[1]
    message("Using radius: ", radius_um, " μm")
  }
  
  # Find matching row
  row_idx <- which(rd$ct1 == ct1 & rd$ct2 == ct2 & rd$radius_um == radius_um)
  
  if (length(row_idx) == 0) {
    stop("Cell pair ", ct1, ":", ct2, " at radius ", radius_um, " not found")
  }
  row_idx <- row_idx[1]
  
  # Extract data for this pair
  y_vals <- yi[row_idx, ]
  v_vals <- vi[row_idx, ]
  
  # Build individual sample data
  sample_data <- data.frame(
    sample = colnames(yi),
    estimate = y_vals,
    variance = v_vals,
    se = sqrt(v_vals),
    group = cd[[group_col]],
    type = "Individual",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(estimate) & !is.na(se)) %>%
    dplyr::mutate(
      ci_lower = estimate - 1.96 * se,
      ci_upper = estimate + 1.96 * se,
      label = sample,
      stats_text = sprintf("%.2f (%.2f)", estimate, se)
    )
  
  # Get pooled estimates with heterogeneity stats (if available)
  pooled_data <- NULL
  tau2_control <- NA
  tau2_case <- NA
  I2_control <- NA
  I2_case <- NA
  
  if (all(c("mu_control", "se_control", "mu_case", "se_case") %in% colnames(rd))) {
    groups <- unique(cd[[group_col]])
    
    # Get tau2 if available
    if ("tau2_control" %in% colnames(rd)) {
      tau2_control <- rd$tau2_control[row_idx]
      tau2_case <- rd$tau2_case[row_idx]
    } else if ("tau2" %in% colnames(rd)) {
      tau2_control <- tau2_case <- rd$tau2[row_idx]
    }
    
    # Get I2 if available (not currently shown, but could be used)
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
    ) %>%
      dplyr::mutate(
        ci_lower = estimate - 1.96 * se,
        ci_upper = estimate + 1.96 * se,
        label = "Pooled (RE)",
        stats_text = sprintf("%.2f (%.2f)", estimate, se),
        tau2_text = ifelse(!is.na(tau2), sprintf("τ² = %.2f", tau2), "")
      )
  }
  
  # Combine data
  plot_data <- dplyr::bind_rows(sample_data, pooled_data) %>%
    dplyr::arrange(group, dplyr::desc(type), estimate) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      row_within_group = dplyr::row_number()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      group_offset = ifelse(group == unique(group)[1], 0, 
                            max(row_within_group[group == unique(group)[1]]) + 2),
      y_pos = row_within_group + group_offset
    )
  
  # Set colors
  if (is.null(group_colors)) {
    groups <- unique(plot_data$group)
    group_colors <- setNames(c("#1f77b4", "#d62728"), groups)
  }
  
  # Calculate x-axis limits
  x_min <- min(plot_data$ci_lower, na.rm = TRUE)
  x_max <- max(plot_data$ci_upper, na.rm = TRUE)
  x_range <- x_max - x_min
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = y_pos)) +
    # Confidence intervals
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lower, xend = ci_upper, 
                   yend = y_pos, color = group),
      linewidth = 0.5
    ) +
    # Point estimates
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(type == "Individual"),
      ggplot2::aes(x = estimate, color = group),
      size = 2.5, shape = 15  # Squares
    ) +
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(type == "Pooled"),
      ggplot2::aes(x = estimate, color = group),
      size = 4.5, shape = 18  # Diamonds
    ) +
    # Sample labels (left)
    ggplot2::geom_text(
      ggplot2::aes(x = x_min - x_range * 0.02, 
                   label = label, 
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 1, size = 3
    ) +
    # Estimate (SE) - middle right
    ggplot2::geom_text(
      ggplot2::aes(x = x_max + x_range * 0.02,
                   label = stats_text,
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 0, size = 3, family = "mono"
    ) +
    # 95% CI - far right
    ggplot2::geom_text(
      ggplot2::aes(x = x_max + x_range * 0.20,
                   label = sprintf("[%.2f, %.2f]", ci_lower, ci_upper),
                   color = group,
                   fontface = ifelse(type == "Pooled", "bold", "plain")),
      hjust = 0, size = 3, family = "mono"
    ) +
    # Tau² annotation for pooled
    {
      if (show_heterogeneity && !is.null(pooled_data)) {
        ggplot2::geom_text(
          data = plot_data %>% dplyr::filter(type == "Pooled" & tau2_text != ""),
          ggplot2::aes(x = x_max + x_range * 0.38,
                       label = tau2_text,
                       color = group),
          hjust = 0, size = 2.8, fontface = "italic"
        )
      }
    } +
    # Reference line at 0
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", 
                        color = "gray30", alpha = 0.5) +
    # Colors
    ggplot2::scale_color_manual(values = group_colors) +
    # Expand x-axis for annotations
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.50))
    ) +
    # Column headers
    ggplot2::annotate("text", x = x_min - x_range * 0.02, 
                      y = max(plot_data$y_pos) + 1,
                      label = "Sample", hjust = 1, fontface = "bold", size = 3.5) +
    ggplot2::annotate("text", x = x_max + x_range * 0.02,
                      y = max(plot_data$y_pos) + 1,
                      label = "Est (SE)", hjust = 0, fontface = "bold", size = 3.5) +
    ggplot2::annotate("text", x = x_max + x_range * 0.20,
                      y = max(plot_data$y_pos) + 1,
                      label = "95% CI", hjust = 0, fontface = "bold", size = 3.5) +
    {
      if (show_heterogeneity && !is.null(pooled_data)) {
        ggplot2::annotate("text", x = x_max + x_range * 0.38,
                          y = max(plot_data$y_pos) + 1,
                          label = "Heterogeneity", hjust = 0, 
                          fontface = "bold", size = 3.5)
      }
    } +
    # Labels
    ggplot2::labs(
      x = paste0("L-cross statistic at ", radius_um, " μm"),
      y = NULL,
      color = NULL,
      title = paste0(gsub("[(){}]", "", ct1), " : ", 
                     gsub("[(){}]", "", ct2)),
      subtitle = "Random-effects meta-analysis with individual sample estimates"
    ) +
    # Theme
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
    ggplot2::coord_cartesian(clip = "off", 
                             ylim = c(0.5, max(plot_data$y_pos) + 1.5))
  
  return(p)
}

#' Construct spatial colocalization network from PANORAMIC results
#'
#' Build a cell-type network where edges represent significant differential
#' spatial colocalization between cell-type pairs, based on the output of
#' \code{pano_compare_groups()} (or a compatible function). Negative
#' \code{z_diff} values (reduced colocalization in the case group) are
#' encoded as edges, weighted by \eqn{|z|}.
#'
#' @param se_diff A \code{SummarizedExperiment} returned by
#'   \code{pano_compare_groups()}, where \code{rowData(se_diff)} contains
#'   at least \code{ct1}, \code{ct2}, \code{z_diff}, \code{p_diff},
#'   and \code{fdr_diff}.
#' @param fdr_threshold Numeric scalar. FDR threshold used to define
#'   significant edges (default \code{0.05}).
#' @param directed Logical. If \code{TRUE}, keep edges directed as
#'   (ct1 → ct2). If \code{FALSE} (default), create a symmetric
#'   undirected network by adding reversed edges before simplification.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{graph}: an \code{igraph} object with edge attributes
#'     \code{weight}, \code{fdr}, \code{pval} and vertex attributes
#'     \code{cluster}, \code{cluster_id}, \code{degree}, \code{betweenness},
#'     \code{strength}.
#'   \item \code{clusters}: the community structure object from
#'     \code{igraph::cluster_leiden()}.
#'   \item \code{n_clusters}: number of detected clusters.
#'   \item \code{modularity}: modularity score of the clustering.
#' }
#'
#' @details
#' Only features with \code{fdr_diff < fdr_threshold} and \code{z_diff < 0}
#' (i.e. reduced colocalization in the case group) are included as edges.
#' Edge weights are set to \eqn{|z_diff|}, and Leiden clustering is applied
#' with modularity as the objective. The resulting network can be visualized
#' with \code{plot_spatial_network()}.
#'
#' @examples
#' \dontrun{
#' net <- create_spatial_network(se_diff, fdr_threshold = 0.05)
#' p   <- plot_spatial_network(net)
#' print(p)
#' }
#'
#' @export
create_spatial_network <- function(se_diff, 
                                   fdr_threshold = 0.05,
                                   directed = FALSE) {
  
  results_df <- as.data.frame(SummarizedExperiment::rowData(se_diff))
  
  # Filter for significant negative associations
  edges <- results_df %>%
    dplyr::filter(
      fdr_diff < fdr_threshold,
      z_diff < 0
    ) %>%
    dplyr::mutate(
      from = gsub("[(){}|]", "", ct1),
      to = gsub("[(){}|]", "", ct2),
      weight = abs(z_diff),
      fdr = fdr_diff,
      pval = p_diff
    ) %>%
    dplyr::select(from, to, weight, fdr, pval)
  
  if (nrow(edges) == 0) {
    stop("No significant edges found at FDR < ", fdr_threshold)
  }
  
  # For undirected: create symmetric edges BEFORE making graph
  if (!directed) {
    edges_rev <- edges %>%
      dplyr::rename(from = to, to = from)
    edges <- dplyr::bind_rows(edges, edges_rev)
  }
  
  # Create graph
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  # Simplify: keep max weight, min fdr
  g <- igraph::simplify(
    g, 
    remove.multiple = TRUE, 
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "max", fdr = "min", pval = "min")
  )
  
  # Leiden clustering
  set.seed(123)
  clusters <- igraph::cluster_leiden(
    g, 
    objective_function = "modularity",
    weights = igraph::E(g)$weight,
    resolution_parameter = 1.2
  )
  
  # Calculate modularity
  mod <- igraph::modularity(g, membership = clusters$membership, 
                            weights = igraph::E(g)$weight)
  
  # Add vertex attributes
  igraph::V(g)$cluster <- as.factor(clusters$membership)
  igraph::V(g)$cluster_id <- clusters$membership
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$betweenness <- igraph::betweenness(g, weights = 1/igraph::E(g)$weight)
  igraph::V(g)$strength <- igraph::strength(g, weights = igraph::E(g)$weight)
  
  return(list(
    graph = g,
    clusters = clusters,
    n_clusters = max(clusters$membership),
    modularity = mod
  ))
}

#' Plot PANORAMIC spatial colocalization network
#'
#' Visualize the spatial colocalization network constructed by
#' \code{create_spatial_network()} using \code{ggraph} and \code{tidygraph}.
#' Nodes represent cell types and edges represent significant differential
#' colocalization, with edge width proportional to \eqn{|z|} and node size
#' based on a chosen centrality measure.
#'
#' @param net_result A list returned by \code{create_spatial_network()},
#'   containing at least a \code{graph} element (an \code{igraph} object)
#'   and summary statistics such as \code{n_clusters} and \code{modularity}.
#' @param layout Character string specifying the graph layout passed to
#'   \code{ggraph} (e.g. \code{"fr"}, \code{"kk"}, \code{"stress"}).
#'   Default is \code{"fr"} (Fruchterman–Reingold).
#' @param node_size_by Character name of a vertex attribute used to scale
#'   node sizes (e.g. \code{"degree"}, \code{"strength"}, \code{"betweenness"}).
#'   Default is \code{"degree"}.
#'
#' @return A \code{ggplot} object representing the network plot.
#'
#' @details
#' The function converts the \code{igraph} object into a \code{tidygraph}
#' \code{tbl_graph} and uses \code{ggraph} to draw edges and nodes. Edge
#' width encodes \eqn{|z|}, edge alpha encodes FDR (darker = more significant),
#' node color encodes community assignment (cluster), and node size encodes
#' the chosen centrality metric.
#'
#' The subtitle summarizes the number of clusters, modularity, and graph
#' size (nodes and edges). This visualization is intended for high-level
#' inspection of modules of cell types with coordinated changes in spatial
#' colocalization.
#'
#' @examples
#' \dontrun{
#' net <- create_spatial_network(se_diff)
#' p   <- plot_spatial_network(net, layout = "fr", node_size_by = "strength")
#' print(p)
#' }
#'
#' @export
plot_spatial_network <- function(net_result, layout = "fr", node_size_by = "degree") {
  
  g <- net_result$graph
  
  # Get edge and node data frames
  edge_list <- igraph::as_data_frame(g, what = "edges")
  node_list <- igraph::as_data_frame(g, what = "vertices")
  
  # Create tidygraph from dataframes
  tg <- tidygraph::tbl_graph(nodes = node_list, edges = edge_list, directed = FALSE)
  
  # Plot
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
      subtitle = sprintf("%d clusters | Q = %.3f | %d nodes, %d edges",
                         net_result$n_clusters, net_result$modularity,
                         igraph::vcount(g), igraph::ecount(g))
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}