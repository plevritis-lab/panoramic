#' Internal: reorder hull vertices counterclockwise
#' @keywords internal
.reorder_ccw <- function(coords) {
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  cxy <- colMeans(coords)
  ang <- atan2(coords[,2] - cxy[2], coords[,1] - cxy[1])
  coords[order(ang, method = "radix"), , drop = FALSE]
}

#' Prepare PANORAMIC inputs from a list of SpatialExperiment objects
#'
#' @param spe_list named or unnamed list of SpatialExperiment (one per sample)
#' @param design data.frame/DataFrame with at least \code{sample}, \code{group}
#' @param cell_type character scalar; name of colData column with cell type labels
#' @param min_cells integer; drop cell types with fewer than this per sample
#' @param concavity numeric passed to concaveman::concaveman()
#' @param window one of "concave","convex","rect"
#' @param BPPARAM BiocParallel param for optional parallel prep
#' @return list of prepared SpatialExperiment with cached PPP in metadata
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
      ct <- ct[keep, drop = TRUE]
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

`%||%` <- function(a,b) if (is.null(a)) b else a
