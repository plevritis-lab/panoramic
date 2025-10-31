#' Internal: generate CT pairs present in all/most samples
#' @keywords internal
.panoramic_pairs <- function(prep_list, cell_type = "cell_type", include_self = TRUE, min_presence = 1.0) {
  stopifnot(min_presence > 0 && min_presence <= 1)
  cts_by_sample <- lapply(prep_list, function(spe) levels(spatstat.geom::marks(spe@metadata$panoramic$ppp)))
  shared <- Reduce(intersect, cts_by_sample)
  if (length(shared) == 0L) stop("No shared cell types found")
  cts <- shared
  pairs <- if (include_self) expand.grid(ct1 = cts, ct2 = cts, stringsAsFactors = FALSE) else {
    utils::combn(cts, 2, simplify = FALSE) |> lapply(\(x) data.frame(ct1 = x[1], ct2 = x[2])) |> 
      dplyr::bind_rows()
  }
  pairs
}

# ---------- helpers ----------------------------------------------------------

.lohboot_quiet <- function(...) {
  res <- NULL
  suppressWarnings(utils::capture.output({
    res <- spatstat.explore::lohboot(...)
  }))
  res
}

.safe_approx <- function(x, y, xout) {
  ok <- is.finite(x) & is.finite(y)
  n  <- sum(ok)
  if (n >= 2L) {
    stats::approx(x[ok], y[ok], xout = xout, rule = 2, ties = "ordered")$y
  } else if (n == 1L) {
    rep(y[ok][1], length(xout))
  } else {
    rep(NA_real_, length(xout))
  }
}

# ---------- corrected lohboot summarizer ------------------------------------
# Derives variance from CI half-width when var/sd are absent; keeps lengths aligned
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

# ---------- per-pair × per-sample worker (uses µm units; never errors) ------
.one_pair_one_sample <- function(meta, ct1, ct2, stat = "Lcross",
                                 nsim = 200, correction = "translate",
                                 radii_um) {
  tab <- meta$marks_tab
  n1 <- unname(tab[ct1]); n2 <- unname(tab[ct2])
  if (is.na(n1) || n1 < 2L || is.na(n2) || n2 < 2L) {
    return(data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_))
  }
  
  X <- meta$ppp  # keep original units (µm); do not rescale
  
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



#' Compute pairwise spatial stats with bootstrap across samples
#'
#' Returns a SummarizedExperiment with assays:
#' \itemize{
#'   \item \code{yi} feature × sample effects
#'   \item \code{vi} feature × sample variances (bootstrap)
#' }
#' RowData stores \code{ct1}, \code{ct2}, \code{radius_um}, \code{stat}.
#'
#' @param prep list from \code{pano_prepare()}
#' @param pairs "auto" or data.frame with ct1, ct2
#' @param radii_um numeric vector of radii at which to evaluate
#' @param stat "Lcross" (default) or "Kcross" / "Lest" / "Kest"
#' @param nsim bootstrap iterations
#' @param correction spatstat edge correction
#' @param keep_boot if TRUE, stores lohboot objects in metadata
#' @param seed RNG seed
#' @param BPPARAM BiocParallel param
#' @return SummarizedExperiment
#' @export
panoramic_spatialstats <- function(
    prep, pairs = "auto", radii_um, stat = "Lcross",
    nsim = 200, correction = "translate", keep_boot = FALSE, seed = 123,
    BPPARAM = BiocParallel::SerialParam()
) {
  stopifnot(is.list(prep), length(prep) >= 2, length(radii_um) >= 1)
  if (identical(pairs, "auto")) {
    pairs <- .panoramic_pairs(prep, include_self = TRUE)
  } else {
    stopifnot(all(c("ct1","ct2") %in% colnames(pairs)))
  }
  
  samples <- names(prep)
  # Compute per-sample per-pair curves, then interpolate at radii_um
  withr::with_seed(seed, {
    # per_sample <- BiocParallel::bplapply(samples, function(sid) {
    #   meta <- prep[[sid]]@metadata$panoramic
    #   cur <- lapply(seq_len(nrow(pairs)), function(i) {
    #     pr <- pairs[i,]
    #     df <- .one_pair_one_sample(meta, pr$ct1, pr$ct2, stat, nsim, correction)
    #     # Interpolate yi/vi at requested radii (df$r in rescaled units; assume microns if inputs were microns)
    #     yi <- stats::approx(df$r, df$yi, xout = radii_um, rule = 2, ties = "ordered")$y
    #     vi <- stats::approx(df$r, df$vi, xout = radii_um, rule = 2, ties = "ordered")$y
    #     data.frame(ct1 = pr$ct1, ct2 = pr$ct2, radius_um = radii_um, yi = yi, vi = vi)
    #   })
    #   dplyr::bind_rows(cur)
    # }, BPPARAM = BPPARAM)
    
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
          radii_um   = radii_um   # <-- pass through
        )
        names(df)[names(df) == "r"] <- "radius_um"   # <-- rename here
        df$ct1 <- pr$ct1
        df$ct2 <- pr$ct2
        df$key <- paste(df$ct1, df$ct2, df$radius_um, sep = "|")
        df
      })
      
      do.call(rbind, cur)
    }, BPPARAM = BPPARAM)
  })
  
  # Assemble SummarizedExperiment: features = (ct1,ct2,r); samples = columns
  feat <- dplyr::bind_rows(per_sample[[1L]])[, c("ct1","ct2","radius_um")]
  feat <- unique(feat)
  feat$key <- paste(feat$ct1, feat$ct2, feat$radius_um, sep = "|")
  
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

#' One-liner convenience: prepare + spatialstats
#' @export
panoramic <- function(
    spe_list, design, cell_type = "cell_type",
    pairs = "auto", radii_um, stat = "Lcross", nsim = 200,
    correction = "translate", min_cells = 5L, concavity = 50, window = "concave",
    keep_boot = FALSE, seed = 123,
    BPPARAM = BiocParallel::SerialParam()
) {
  prep <- panoramic_prepare(spe_list, design, cell_type, min_cells, concavity, window, BPPARAM)
  panoramic_spatialstats(prep, pairs, radii_um, stat, nsim, correction, keep_boot, seed, BPPARAM)
}
