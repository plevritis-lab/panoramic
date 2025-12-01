#' Internal: generate cell-type pairs
#' 
#' Construct the set of cell-type pairs to analyze, optionally including 
#' self-pairs and filtering to cell types present in at least a given 
#' fraction of samples.
#' 
#' @param prep_list List of prepared SpatialExperiment objects, each with \code{metadata(spe)$panoramic$ppp} containing a marked spatstat point pattern. 
#' @param include_self Logical. If \code{TRUE}, include (ct, ct) self-pairs. If \code{FALSE}, only cross-type pairs are returned. 
#' @param min_presence Numeric in \eqn{[0, 1]}. Minimum fraction of samples in which a cell type must appear to be included. \code{0} keeps all observed types, \code{0.5} requires presence in at least half of samples. 
#' 
#' @return A data.frame with columns \code{ct1} and \code{ct2}. 
#' 
#' @keywords internal
.panoramic_pairs <- function(prep_list,
                             include_self = TRUE, min_presence = 0.0) {
  
  stopifnot(min_presence >= 0 && min_presence <= 1)
  
  # Get cell types from each sample
  cts_by_sample <- lapply(prep_list, function(spe) {
    levels(spatstat.geom::marks(spe@metadata$panoramic$ppp))
  })
  
  # This includes cell types present in ANY sample
  all_cts <- Reduce(union, cts_by_sample)
  
  if (length(all_cts) == 0L) stop("No cell types found")
  
  # Optional: Filter by minimum presence across samples
  if (min_presence > 0) {
    ct_counts <- table(unlist(cts_by_sample))
    min_samples <- ceiling(length(prep_list) * min_presence)
    all_cts <- names(ct_counts)[ct_counts >= min_samples]
    
    if (length(all_cts) == 0L) {
      stop("No cell types present in at least ", 
           round(min_presence * 100), "% of samples")
    }
  }
  
  cts <- sort(all_cts)
  
  # Generate pairs
  pairs <- if (include_self) {
    expand.grid(ct1 = cts, ct2 = cts, stringsAsFactors = FALSE)
  } else {
    utils::combn(cts, 2, simplify = FALSE) |> 
      lapply(\(x) data.frame(ct1 = x[1], ct2 = x[2])) |> 
      dplyr::bind_rows()
  }
  
  pairs
}

#' Internal: quiet Loh bootstrap wrapper
#' 
#' Run \code{spatstat.explore::lohboot()} while suppressing console output, 
#' returning the bootstrap object invisibly for further summarization. 
#' 
#' @param ... Arguments passed to \code{spatstat.explore::lohboot()}. 
#' 
#' @return A Loh bootstrap object as returned by \code{spatstat.explore::lohboot()}.
#' 
#' @keywords internal
.lohboot_quiet <- function(...) {
  res <- NULL
  suppressWarnings(utils::capture.output({
    res <- spatstat.explore::lohboot(...)
  }))
  res
}

#' Internal: safe interpolation on a common radius grid
#' 
#' Interpolate values \code{y} at target \code{xout} using \code{stats::aprox()}, 
#' handling missing/non-finite values and low sample sizes. 
#' 
#' @param x Numeric vector of original inputs x.
#' @param y Numeric vector of original y-values.
#' @param xout Numeric vector of target x-values. 
#' 
#' @return Numeric vector of interpolated values, length \code{length(xout)}. 
#' 
#' @keywords internal
.safe_approx <- function(x, y, xout) {
  ok <- is.finite(x) & is.finite(y)
  n  <- sum(ok)
  if (n >= 2L) {
    stats::approx(x[ok], y[ok], xout = xout, rule = 1, ties = "ordered")$y
    # stats::approx(x[ok], y[ok], xout = xout, rule = 2, ties = "ordered")$y
    
  } else if (n == 1L) {
    rep(y[ok][1], length(xout))
  } else {
    rep(NA_real_, length(xout))
  }
}

#' Internal: summarize Loh bootstrap to (r, yi, vi)
#' 
#' Convert a spatstat Loh bootstrap object to aligned vectors of radii, 
#' centered estiamtes, and variance. Variance is taken directly if present, 
#' reconstructed from confidence intervale if needed, or computed from 
#' simulation replicates as a fallback. 
#' 
#' @param loh_obj Loh bootstrap object returned by spatstat's L/K-summary functions with \code{global=TRUE}. 
#' @param center_L Logical. If \code{TRUE}, center the estimate by the theoretical curve (e.g. L(r)-r) when available.
#' @param conf Numeric confidence level used when reconstructing variance from confidence interval half-width. 
#' 
#' @return A data.frame with numeric columns \code{r}, \code{yi}, and \code{vi}, or \code{NULL} if the object cannot be parsed. 
#' 
#' @keywords internal
.summarize_lohboot <- function(loh_obj, center_L = TRUE, conf = 0.95) {
  df <- try(as.data.frame(loh_obj), silent = TRUE)
  if (inherits(df, "try-error") || NROW(df) == 0L) return(NULL)
  
  # r (first fv column or named 'r')
  r <- if ("r" %in% names(df)) as.numeric(df$r) else as.numeric(df[[1]])
  
  # choose estimate column
  est_name <- ({
    cand <- intersect(c("trans","iso","obs","est","border"), names(df))
    if (length(cand)) cand[1] else {
      num <- vapply(df, is.numeric, logical(1))
      num[match(c("r","theo"), names(df), nomatch = 0L)] <- FALSE
      pick <- which(num)
      if (length(pick)) names(df)[pick[length(pick)]] else names(df)[2]
    }
  })
  yi <- as.numeric(df[[est_name]])
  
  # center L by theo (usually r) if present
  if (center_L && "theo" %in% names(df)) yi <- yi - as.numeric(df$theo)
  
  # try variance / sd directly
  vi <- NULL
  for (nm in c("var","variance","sd","SD","Var","VAR")) {
    if (nm %in% names(df)) {
      x  <- df[[nm]]
      vi <- if (grepl("sd", nm, ignore.case = TRUE)) as.numeric(x)^2 else as.numeric(x)
      break
    }
  }
  
  # reconstruct from CI if needed
  if (is.null(vi) || !any(is.finite(vi))) {
    lo_candidates <- c(paste0("lo.", est_name), paste0("lo", est_name), "lo", "lower", "loiso", "lotrans")
    hi_candidates <- c(paste0("hi.", est_name), paste0("hi", est_name), "hi", "upper", "hiiso", "hitrans")
    lo_nm <- intersect(lo_candidates, names(df))
    hi_nm <- intersect(hi_candidates, names(df))
    if (length(lo_nm) >= 1L && length(hi_nm) >= 1L) {
      lo <- as.numeric(df[[lo_nm[1]]])
      hi <- as.numeric(df[[hi_nm[1]]])
      z  <- stats::qnorm((1 + conf)/2)
      halfwidth <- (hi - lo) / 2
      vi <- (halfwidth / z)^2
    }
  }
  
  # fallback: rowVars over sim* columns
  if (is.null(vi) || length(vi) == 0L || !any(is.finite(vi))) {
    simcols <- grep("^sim", names(df), value = TRUE)
    if (length(simcols) >= 2L) {
      vi <- matrixStats::rowVars(as.matrix(df[simcols]))
    }
  }
  if (is.null(vi) || length(vi) == 0L) vi <- rep(NA_real_, length(r))
  
  # enforce equal lengths defensively
  n <- length(r)
  yi <- if (length(yi) == n) yi else rep_len(yi, n)
  vi <- if (length(vi) == n) vi else rep_len(vi, n)
  
  data.frame(r = r, yi = yi, vi = vi)
}

#' Internal: compute spatial stats for one pair in one sample
#' 
#' Compute a spatial summary curve and bootstrap variance for a single
#' cell-type pair within a single prepared sample, evaluated on a common 
#' radius grid. 
#' 
#' @param meta The \code{metadata$panoramic} list for one sample, containing at least a spatstat \code{ppp} object and a \code{marks_tab} table. 
#' @param ct1,ct2 Character. Cell-type labels. 
#' @param stat Character. Summary statistic ("Lcross", "Kcross", "Lest", "Kest"). 
#' @param nsim Integer. Number of Loh bootstrap simulations. 
#' @param correction Character. Edge correction passed to spatstat.
#' @param radii_um Numeric vector of radii (microns) on which to summarize.
#' 
#' @return A data.frame with columns \code{radius_um}, \code{yi}, and \code{vi}, filled with \code{NA} if insufficient cells of either type are present. 
#' 
#' @keywords internal
.one_pair_one_sample <- function(meta, ct1, ct2, stat = "Lcross",
                                 nsim = 100, correction = "translate",
                                 radii_um) {
  tab <- meta$marks_tab
  n1 <- unname(tab[ct1]); n2 <- unname(tab[ct2])
  if (is.na(n1) || n1 < 2L || is.na(n2) || n2 < 2L) {
    return(data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_))
  }
  
  X <- meta$ppp
  
  if (!identical(ct1, ct2) && stat %in% c("Lcross","Kcross")) {
    fun <- if (stat == "Lcross") spatstat.explore::Lcross else spatstat.explore::Kcross
    loh <- .lohboot_quiet(
      spatstat.geom::subset.ppp(X, marks %in% c(ct1, ct2)),
      fun, from = ct1, to = ct2, correction = correction, global = TRUE, nsim = nsim
    )
  } else {
    fun <- if (grepl("^L", stat)) spatstat.explore::Lest else spatstat.explore::Kest
    loh <- .lohboot_quiet(
      spatstat.geom::subset.ppp(X, marks == ct1),
      fun, correction = correction, global = TRUE, nsim = nsim
    )
  }
  
  df <- try(.summarize_lohboot(loh), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df) || sum(is.finite(df$r)) < 1L) {
    return(data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_))
  }
  
  o  <- order(df$r)
  yi <- .safe_approx(df$r[o], df$yi[o], radii_um)
  vi <- .safe_approx(df$r[o], df$vi[o], radii_um)
  
  data.frame(radius_um = radii_um, yi = yi, vi = vi)
}

#' Compute pairwise spatial statistics for PANORAMIC
#' 
#' Compute pairwise spatial summary curves (e.g. Lcross) and bootstrap
#' variances for all requested cell-type pairs and radii across multiple
#' prepared samples, returning a SummarizedExperiment. 
#' 
#' @param prep List of prepared SpatialExperiment objects as returned by \code{panoramic_prepare()}, 
#'  each with \code{metadata(spe)$panoramic}. 
#' @param pairs Either "auto" (generate all cell-type pairs including self-pairs) 
#'  or a data.frame with columns \code{ct1}, \code{ct2}. 
#' @param radii_um Numeric vector of radii (microns) at which to evaluate the colocalization statistic. 
#' @param stat Character. Summary statistic ("Lcross", "Lest", "Kcross", "Kest"). 
#' @param nsim Integer. Number of Loh bootstrap simulations per sample/pair. 
#' @param correction Chatacter. Edge correction method for spatstat ("translate", "border", ...). 
#' @param keep_boot Logical. Reserved for future use (keeping full bootstrap objects). 
#' 
#' @return A SummarizedExperiment with : 
#' \itemize{
#'  \item assay "yi": centered estimates per (ct1, ct2, radius) feature and sample. 
#'  \item assay "vi": variance estimates aligned to "yi".
#' }
#' rowData stores \code{ct1}, \code{ct2}, \code{radius_um}, \code{stat}
#' colData stores \code{sample} and \code{group}. 
#' 
#' @export
panoramic_spatialstats <- function(
    prep, pairs = "auto", radii_um, stat = "Lcross",
    nsim = 100, correction = "translate", keep_boot = FALSE, seed = 123,
    BPPARAM = BiocParallel::SerialParam()
) {
  stopifnot(is.list(prep), length(prep) >= 2, length(radii_um) >= 1)
  if (identical(pairs, "auto")) {
    pairs <- .panoramic_pairs(prep, include_self = TRUE)
  } else {
    stopifnot(all(c("ct1","ct2") %in% colnames(pairs)))
  }
  
  samples <- names(prep)
  
  # Compute per-sample per-pair curves
  withr::with_seed(seed, {
    per_sample <- BiocParallel::bplapply(samples, function(sid) {
      meta <- prep[[sid]]@metadata$panoramic
      
      cur <- lapply(seq_len(nrow(pairs)), function(i) {
        pr <- pairs[i, ]
        df <- .one_pair_one_sample(
          meta,
          pr$ct1, pr$ct2,
          stat       = stat,
          nsim       = nsim,
          correction = correction,
          radii_um   = radii_um 
        )
        names(df)[names(df) == "r"] <- "radius_um"
        df$ct1 <- pr$ct1
        df$ct2 <- pr$ct2
        df$key <- paste(df$ct1, df$ct2, df$radius_um, sep = "|")
        df
      })
      
      do.call(rbind, cur)
    }, BPPARAM = BPPARAM)
  })
  
  # Build feature list from ALL samples
  all_features <- lapply(per_sample, function(df) {
    unique(df[, c("ct1", "ct2", "radius_um")])
  })
  feat <- unique(dplyr::bind_rows(all_features))
  feat$key <- paste(feat$ct1, feat$ct2, feat$radius_um, sep = "|")
  
  # Build matrices
  make_mat <- function(slot) {
    m <- matrix(NA_real_, nrow = nrow(feat), ncol = length(samples),
                dimnames = list(feat$key, samples))
    for (j in seq_along(samples)) {
      df <- per_sample[[j]]
      df$key <- paste(df$ct1, df$ct2, df$radius_um, sep = "|")
      idx <- match(df$key, feat$key)
      m[idx, j] <- df[[slot]]
    }
    m
  }
  yi <- make_mat("yi"); vi <- make_mat("vi")
  
  rowData <- S4Vectors::DataFrame(
    ct1 = feat$ct1, ct2 = feat$ct2, radius_um = feat$radius_um, stat = stat
  )
  colData <- S4Vectors::DataFrame(
    sample = samples,
    group  = vapply(prep, function(s) s@metadata$panoramic$group_id, character(1))
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = yi, vi = vi),
    rowData = rowData, colData = colData,
    metadata = list(panoramic = list(
      radii_um = radii_um,
      pairs = pairs,
      stat = stat
    ))
  )
  se
}


#' PANORAMIC: prepare and compute spatial statistics
#' 
#' One-line convenience wrapper that calles \code{panoramic_prepare()} followed 
#' by \code{panoramic_spatialstats()} to produce the PANORAMIC
#' SummarizedExperiment for donwnstream meta-analysis. 
#' 
#' @inheritParams panoramic_prepare
#' @inheritParams panoramic_spatialstats
#' 
#' @return A SummerizedExperiment as described in \code{panoramic_spatialstats()}. 
#'
#' @examples
#' \dontrun{
#' se <- panoramic(
#'   spe_list, 
#'   design=design, 
#'   radii_um = c(50, 75, 100)
#' )
#' }
#' 
#' @export
panoramic <- function(
    spe_list, design, cell_type = "cell_type",
    pairs = "auto", radii_um, stat = "Lcross", nsim = 100,
    correction = "translate", min_cells = 5L, concavity = 50, window = "concave",
    keep_boot = FALSE, seed = 123,
    BPPARAM = BiocParallel::SerialParam()
) {
  prep <- panoramic_prepare(spe_list, design, cell_type, min_cells, concavity, window, BPPARAM)
  panoramic_spatialstats(prep, pairs, radii_um, stat, nsim, correction, keep_boot, seed, BPPARAM)
}