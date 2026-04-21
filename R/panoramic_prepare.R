#' Internal: reorder hull vertices counterclockwise
#' Orders a set of 2D coordinates counterclockwise around their centroid.
#' Mainly used for convex spatial window construction to ensure correct polygon orientation for spatstat owin objects. 
#'
#' @param coords Numeric matrix of coordinates (columns: x, y).
#' @return Reordered numeric matrix of coordinates
#' @keywords internal
.reorder_ccw <- function(coords) {
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  cxy <- colMeans(coords)
  ang <- atan2(coords[,2] - cxy[2], coords[,1] - cxy[1])
  coords[order(ang, method = "radix"), , drop = FALSE]
}

#' Prepare PANORAMIC inputs from a list of SpatialExperiment objects
#'
#' Creates SpatialExperiment objects ready for PANORAMIC spatial analyses. 
#' Cell type labels are harmonized, rare cell types (fewer than \code{min_cells}) are dropped per sample, 
#' and a spatial window is computed. Cached spatstat objects are stored within each SpatialExperiment's metadata. 
#' 
#' @param spe_list Named or unnamed list of SpatialExperiment (one per sample)
#' @param design data.frame with at least columns \code{sample}, \code{group} to map samples for meta-analysis. If only one group is used, give all the same group label.
#' @param cell_type Character; name of SpatialExperiment colData column holding cell type labels
#' @param min_cells Integer. Cell types with fewer than this count (per sample) are dropped.
#' @param concavity Numeric passed to concaveman::concaveman(). Controls level of hull detail. 1 is highly detailed, \code{Inf} is a convex hull.
#' @param window one of "concave","convex","rect". Typically use concave.
#' @param BPPARAM BiocParallel param for optional parallel processing.
#'
#' @return List of SpatialExperiment objects with metadata slot \code{panoramic} containing \code{ppp}, cell-type table, spatial window, group/sample info.
#' 
#' @details 
#' This step computes per-sample spatial windows to exclude background, filters
#' rare cell types separately per sample, builds consistent cell-type factor levels, 
#' and caches spatstat objects and type tables for PANORAMIC's spatial statistics. 
#' @examples
#' spe_list <- list(
#'   sample1 = panoramic_simulate_spe(
#'     n_cells = 60,
#'     sample_id = "sample1",
#'     scenario = "random",
#'     seed = 1
#'   )
#' )
#'
#' # Design with a single group ----------------------------------------
#' design <- data.frame(
#'   sample = "sample1",
#'   group  = "group1",
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run panoramic_prepare ---------------------------------------------
#' prepped <- panoramic_prepare(
#'   spe_list,
#'   design      = design,
#'   cell_type   = "cell_type",
#'   min_cells   = 3,
#'   concavity   = 50,
#'   window      = "concave",
#'   BPPARAM     = BiocParallel::SerialParam()
#' )
#'
#' # Inspect cached spatstat objects in metadata
#' names(S4Vectors::metadata(prepped[[1]])$panoramic)
#'
#' @export
panoramic_prepare <- function(
    spe_list, design, cell_type = "cell_type",
    min_cells = 5, concavity = 50, window = c("concave","convex","rect"),
    BPPARAM = BiocParallel::SerialParam()
) {
  window <- match.arg(window)
  if (!is.list(spe_list) || length(spe_list) < 1L) {
    stop("`spe_list` must be a non-empty list of SpatialExperiment objects.")
  }
  if (!is.data.frame(design)) {
    stop("`design` must be a data.frame with columns `sample` and `group`.")
  }
  if (!all(c("sample", "group") %in% colnames(design))) {
    stop("`design` must contain columns `sample` and `group`.")
  }
  design$sample <- as.character(design$sample)
  design$group <- as.character(design$group)
  if (anyNA(design$sample) || any(!nzchar(design$sample))) {
    stop("`design$sample` must be non-missing and non-empty.")
  }
  if (anyNA(design$group) || any(!nzchar(design$group))) {
    stop("`design$group` must be non-missing and non-empty.")
  }
  if (anyDuplicated(design$sample)) {
    dups <- unique(design$sample[duplicated(design$sample)])
    stop(
      "`design$sample` contains duplicates: ",
      paste(dups, collapse = ", "),
      "."
    )
  }
  if (nrow(design) < length(spe_list)) {
    stop("`design` has fewer rows than `spe_list` length.")
  }
  if (is.null(names(spe_list))) {
    # align by order
    names(spe_list) <- as.character(design$sample[seq_along(spe_list)])
  }
  if (anyNA(names(spe_list)) || any(!nzchar(names(spe_list)))) {
    stop("`spe_list` names must be non-missing and non-empty.")
  }
  if (!all(names(spe_list) %in% design$sample))
    stop("All list names must appear in design$sample")
  missing_cell_type <- names(spe_list)[!vapply(spe_list, function(spe) {
    cell_type %in% colnames(SummarizedExperiment::colData(spe))
  }, logical(1))]
  if (length(missing_cell_type) > 0L) {
    stop(
      "Missing `cell_type` column `", cell_type, "` in samples: ",
      paste(missing_cell_type, collapse = ", "),
      "."
    )
  }
  min_cells <- as.integer(min_cells)
  if (!is.finite(min_cells) || length(min_cells) != 1L || min_cells < 1L) {
    stop("`min_cells` must be a positive integer.")
  }
  
  # Harmonize cell-type levels across samples
  all_ct <- unique(unlist(lapply(spe_list, function(spe) {
    as.character(SummarizedExperiment::colData(spe)[[cell_type]])
  })))
  all_ct <- sort(all_ct)
  
  # Per-sample prep function
  .prep_one <- function(sid) {
    spe <- spe_list[[sid]]
    meta <- S4Vectors::metadata(spe)
    meta$panoramic <- meta$panoramic %||% list()
    meta$panoramic$sample_id <- sid
    i_design <- match(sid, design$sample)
    if (is.na(i_design)) {
      stop("Sample `", sid, "` is missing from `design$sample`.", call. = FALSE)
    }
    meta$panoramic$group_id  <- as.character(design$group[i_design])
    
    # coords
    coords <- tryCatch(
      SpatialExperiment::spatialCoords(spe),
      error = function(e) {
        stop("Failed to extract spatial coordinates for sample ", sid, ".", call. = FALSE)
      }
    )
    if (!is.matrix(coords) || ncol(coords) < 2) stop("Bad spatialCoords for ", sid)
    coords <- as.matrix(coords[, seq_len(2L), drop = FALSE])
    if (any(!is.finite(coords))) {
      stop("Non-finite spatial coordinates for sample ", sid, ".", call. = FALSE)
    }
    
    # filter rare cell types (per sample)
    ct <- factor(SummarizedExperiment::colData(spe)[[cell_type]], levels = all_ct)
    keep <- !is.na(ct)
    # drop CT with < min_cells
    tab <- table(ct[keep])
    ok  <- ct %in% names(tab)[tab >= min_cells]
    keep <- keep & ok
    if (!all(keep)) {
      spe <- spe[, keep, drop = FALSE]
      coords <- coords[keep, , drop = FALSE]
      ct <- ct[keep, drop = FALSE]      
    }
    if (nrow(coords) < 1L) {
      stop(
        "No cells remain in sample `", sid, "` after filtering with `min_cells = ",
        min_cells, "`.",
        call. = FALSE
      )
    }
    xr <- range(coords[, 1], na.rm = TRUE)
    yr <- range(coords[, 2], na.rm = TRUE)
    if (!all(is.finite(c(xr, yr))) || diff(xr) <= 0 || diff(yr) <= 0) {
      stop(
        "Sample `", sid, "` has degenerate coordinates (zero-area window) after filtering.",
        call. = FALSE
      )
    }
    
    # window
    make_rect_window <- function() {
      spatstat.geom::owin(xr, yr)
    }
    win <- switch(window,
                  concave = {
                    if (nrow(coords) < 3L) {
                      warning(
                        "Sample `", sid, "` has <3 points after filtering; using rectangular window instead of concave hull.",
                        call. = FALSE
                      )
                      make_rect_window()
                    } else {
                      hull <- try(concaveman::concaveman(coords, concavity = concavity), silent = TRUE)
                      if (inherits(hull, "try-error") || !is.matrix(hull) || nrow(hull) < 3L) {
                        warning(
                          "Concave hull failed for sample `", sid, "`; using rectangular window.",
                          call. = FALSE
                        )
                        make_rect_window()
                      } else {
                        hull <- .reorder_ccw(hull)
                        spatstat.geom::owin(poly = list(x = hull[, 1], y = hull[, 2]))
                      }
                    }
                  },
                  convex = {
                    if (nrow(coords) < 3L) {
                      warning(
                        "Sample `", sid, "` has <3 points after filtering; using rectangular window instead of convex hull.",
                        call. = FALSE
                      )
                      make_rect_window()
                    } else {
                      ch <- grDevices::chull(coords[, 1], coords[, 2])
                      hull <- .reorder_ccw(coords[ch, , drop = FALSE])
                      spatstat.geom::owin(poly = list(x = hull[, 1], y = hull[, 2]))
                    }
                  },
                  rect = {
                    make_rect_window()
                  }
    )
    
    # PPP + cache
    ppp <- spatstat.geom::ppp(coords[,1], coords[,2], marks = ct, window = win)
    meta$panoramic$ppp          <- ppp
    meta$panoramic$ppp_rescaled <- spatstat.geom::rescale(ppp)
    meta$panoramic$marks_tab    <- table(spatstat.geom::marks(ppp))
    meta$panoramic$results      <- list(colocalization_bootstrap = list())
    
    S4Vectors::metadata(spe) <- meta
    spe
  }
  
  ids <- names(spe_list)
  out <- BiocParallel::bplapply(ids, .prep_one, BPPARAM = BPPARAM)
  names(out) <- ids
  out
}

#' Pipe-compatible null operator
#' 
#' Returns b if a is NULL, otherwise a. 
#'
#' @name null_coalesce
#' @aliases %||%
#' 
#' @param a, b Objects to test; if a is NULL b is returned. 
#' 
#' @return The non-NULL object. 
#' 
#' @keywords internal
`%||%` <- function(a,b) if (is.null(a)) b else a
