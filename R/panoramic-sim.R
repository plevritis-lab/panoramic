#' Simulate one SpatialExperiment object for PANORAMIC examples
#'
#' Create a toy \code{SpatialExperiment} with simple spatial patterns that can
#' be used in package examples, vignettes, and tests.
#'
#' @param n_cells Integer number of cells to simulate.
#' @param sample_id Character sample identifier stored in \code{colData}.
#' @param cell_types Character vector of cell-type labels.
#' @param scenario Character string, either \code{"random"} or
#'   \code{"colocalized"}. In \code{"colocalized"}, the first two cell types
#'   are enriched near the center.
#' @param bounds Numeric length-2 vector giving minimum and maximum coordinate.
#' @param center Numeric length-2 vector giving center for colocalized pattern.
#' @param cluster_sd Numeric standard deviation for central clustering.
#' @param n_genes Integer number of toy genes in the counts matrix.
#' @param seed Optional integer for reproducibility.
#'
#' @return A \code{SpatialExperiment} object with simulated coordinates and
#'   cell types.
#'
#' @examples
#' spe <- panoramic_simulate_spe(
#'   n_cells = 120,
#'   sample_id = "sample_1",
#'   scenario = "colocalized",
#'   seed = 1
#' )
#' spe
#'
#' @export
panoramic_simulate_spe <- function(
    n_cells = 200L,
    sample_id = "sample_1",
    cell_types = c("(A)", "(B)", "(C)"),
    scenario = c("random", "colocalized"),
    bounds = c(0, 100),
    center = c(50, 50),
    cluster_sd = 18,
    n_genes = 10L,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  stopifnot(length(bounds) == 2L, length(center) == 2L, length(cell_types) >= 2L)
  
  .simulate <- function() {
    if (identical(scenario, "random")) {
      x <- stats::runif(n_cells, min = bounds[1], max = bounds[2])
      y <- stats::runif(n_cells, min = bounds[1], max = bounds[2])
      ct <- sample(cell_types, size = n_cells, replace = TRUE)
    } else {
      x <- stats::rnorm(n_cells, mean = center[1], sd = cluster_sd)
      y <- stats::rnorm(n_cells, mean = center[2], sd = cluster_sd)
      x <- pmin(pmax(x, bounds[1]), bounds[2])
      y <- pmin(pmax(y, bounds[1]), bounds[2])
      
      d <- sqrt((x - center[1])^2 + (y - center[2])^2)
      p_center <- exp(-(d^2) / (2 * cluster_sd^2))
      p_center <- p_center / max(p_center)
      
      ct <- ifelse(
        stats::runif(n_cells) < p_center,
        sample(cell_types[1:2], size = n_cells, replace = TRUE),
        sample(cell_types, size = n_cells, replace = TRUE)
      )
    }
    
    coords <- cbind(x = x, y = y)
    cell_ids <- paste0(sample_id, "_cell", seq_len(n_cells))
    rownames(coords) <- cell_ids
    
    counts <- matrix(
      stats::rpois(n_genes * n_cells, lambda = 5),
      nrow = n_genes,
      dimnames = list(paste0("gene", seq_len(n_genes)), cell_ids)
    )
    
    SpatialExperiment::SpatialExperiment(
      assays = list(counts = counts),
      colData = S4Vectors::DataFrame(
        cell_type = ct,
        sample_id = sample_id
      ),
      spatialCoords = coords
    )
  }
  
  if (is.null(seed)) {
    .simulate()
  } else {
    withr::with_seed(seed, .simulate())
  }
}


#' Simulate a two-group PANORAMIC example dataset
#'
#' Create a list of \code{SpatialExperiment} objects and matching design table
#' for differential colocalization tutorials.
#'
#' @param n_group1 Integer number of samples in group 1.
#' @param n_group2 Integer number of samples in group 2.
#' @param n_cells_group1 Integer number of cells per group-1 sample.
#' @param n_cells_group2 Integer number of cells per group-2 sample.
#' @param group_labels Character length-2 vector of group names.
#' @param scenario_group1 Scenario passed to \code{panoramic_simulate_spe()} for
#'   group 1.
#' @param scenario_group2 Scenario passed to \code{panoramic_simulate_spe()} for
#'   group 2.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with entries \code{spe_list} and \code{design}.
#'
#' @examples
#' toy <- panoramic_simulate_dataset(seed = 1)
#' names(toy)
#' head(toy$design)
#'
#' @export
panoramic_simulate_dataset <- function(
    n_group1 = 3L,
    n_group2 = 3L,
    n_cells_group1 = 200L,
    n_cells_group2 = 350L,
    group_labels = c("group1", "group2"),
    scenario_group1 = "random",
    scenario_group2 = "colocalized",
    seed = NULL
) {
  stopifnot(length(group_labels) == 2L)
  
  .simulate <- function() {
    n_total <- n_group1 + n_group2
    sample_ids <- paste0("sample_", seq_len(n_total))
    group_vec <- c(rep(group_labels[1], n_group1), rep(group_labels[2], n_group2))
    n_cells_vec <- c(rep(n_cells_group1, n_group1), rep(n_cells_group2, n_group2))
    scenario_vec <- c(
      rep(scenario_group1, n_group1),
      rep(scenario_group2, n_group2)
    )
    
    spe_list <- lapply(seq_len(n_total), function(i) {
      panoramic_simulate_spe(
        n_cells = n_cells_vec[i],
        sample_id = sample_ids[i],
        scenario = scenario_vec[i]
      )
    })
    names(spe_list) <- sample_ids
    
    design <- data.frame(
      sample = sample_ids,
      group = group_vec,
      stringsAsFactors = FALSE
    )
    
    list(spe_list = spe_list, design = design)
  }
  
  if (is.null(seed)) {
    .simulate()
  } else {
    withr::with_seed(seed, .simulate())
  }
}
