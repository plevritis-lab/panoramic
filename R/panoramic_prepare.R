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
#' \dontrun{
#' # For multiple samples: 
#' prepped <- panoramic_prepare(spe_list, design, min_cells = 5)
#' }
#' 
#' @export
panoramic_prepare <- function(
    spe_list, design, cell_type = "cell_type",
    min_cells = 5, concavity = 50, window = c("concave","convex","rect"),
    BPPARAM = BiocParallel::SerialParam()
) {
  window <- match.arg(window)
  stopifnot(is.list(spe_list), nrow(design) >= length(spe_list))
  if (is.null(names(spe_list))) {
    # align by order
    names(spe_list) <- as.character(design$sample[seq_along(spe_list)])
  }
  if (!all(names(spe_list) %in% design$sample))
    stop("All list names must appear in design$sample")
  
  # Harmonize cell-type levels across samples
  all_ct <- unique(unlist(lapply(spe_list, function(spe) {
    as.character(SummarizedExperiment::colData(spe)[[cell_type]])
  })))
  all_ct <- sort(all_ct)
  
  # Per-sample prep function
  .prep_one <- function(sid) {
    spe <- spe_list[[sid]]
    meta <- spe@metadata
    meta$panoramic <- meta$panoramic %||% list()
    meta$panoramic$sample_id <- sid
    meta$panoramic$group_id  <- as.character(design$group[match(sid, design$sample)])
    
    # coords
    coords <- tryCatch(SpatialExperiment::spatialCoords(spe),
                       error = function(e) as.matrix(spe@int_colData$spatialCoords))
    if (!is.matrix(coords) || ncol(coords) < 2) stop("Bad spatialCoords for ", sid)
    
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
    
    # window
    win <- switch(window,
                  concave = {
                    hull <- concaveman::concaveman(coords, concavity = concavity)
                    hull <- .reorder_ccw(hull)
                    spatstat.geom::owin(poly = list(x = hull[,1], y = hull[,2]))
                  },
                  convex = {
                    ch <- grDevices::chull(coords[,1], coords[,2])
                    hull <- .reorder_ccw(coords[ch, , drop = FALSE])
                    spatstat.geom::owin(poly = list(x = hull[,1], y = hull[,2]))
                  },
                  rect = {
                    xr <- range(coords[,1]); yr <- range(coords[,2])
                    spatstat.geom::owin(xr, yr)
                  }
    )
    
    # PPP + cache
    ppp <- spatstat.geom::ppp(coords[,1], coords[,2], marks = ct, window = win)
    meta$panoramic$ppp          <- ppp
    meta$panoramic$ppp_rescaled <- spatstat.geom::rescale(ppp)
    meta$panoramic$marks_tab    <- table(spatstat.geom::marks(ppp))
    meta$panoramic$results      <- list(colocalization_bootstrap = list())
    
    spe@metadata <- meta
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
