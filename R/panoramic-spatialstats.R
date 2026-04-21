#' Internal: generate cell-type pairs
#' 
#' Construct the set of cell-type pairs to analyze, optionally including 
#' self-pairs and filtering to cell types present in at least a given 
#' fraction of samples.
#' 
#' @param prep_list List of prepared SpatialExperiment objects, each with \code{metadata(spe)$panoramic$ppp} containing a marked spatstat point pattern. 
#' @param include_self Logical. If \code{TRUE}, include (ct, ct) self-pairs. If \code{FALSE}, only cross-type pairs are returned. 
#' @param min_presence Numeric in \eqn{[0, 1]}. Minimum fraction of samples in which a cell type must be observed (at least one cell) to be included. \code{0} keeps all observed types, \code{0.5} requires observed presence in at least half of samples. 
#' 
#' @return A data.frame with columns \code{ct1} and \code{ct2}. 
#' 
#' @keywords internal
.panoramic_pairs <- function(prep_list,
                             include_self = TRUE, min_presence = 0.0) {
  
  stopifnot(min_presence >= 0 && min_presence <= 1)
  
  # Get observed cell types from each sample.
  # Use observed marks instead of factor levels to avoid counting absent types.
  cts_by_sample <- lapply(prep_list, function(spe) {
    marks_chr <- as.character(spatstat.geom::marks(S4Vectors::metadata(spe)$panoramic$ppp))
    sort(unique(marks_chr[!is.na(marks_chr)]), method = "radix")
  })
  
  # This includes cell types present in ANY sample
  all_cts <- Reduce(union, cts_by_sample)
  
  if (length(all_cts) == 0L) stop("No cell types found")
  
  # Optional: Filter by minimum presence across samples
  if (min_presence > 0) {
    ct_counts <- table(unlist(cts_by_sample, use.names = FALSE))
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
    if (length(cts) < 2L) {
      return(data.frame(ct1 = character(0), ct2 = character(0)))
    }
    utils::combn(cts, 2, simplify = FALSE) |> 
      lapply(\(x) data.frame(ct1 = x[1], ct2 = x[2])) |> 
      dplyr::bind_rows()
  }
  
  pairs
}

#' Internal: local function info (spatstat helper)
#'
#' @param key Character scalar naming a supported spatstat summary function.
#'
#' @return A named list describing the global/local function mapping and
#'  options for the requested key, or \code{NULL} if unsupported.
.spatstat_local_function_info <- function(key) {
  TheTable <- list(
    pcf          = list(Global = spatstat.explore::pcf,
                        Local = spatstat.explore::localpcf,
                        L = FALSE, inhom = FALSE, indices = 0),
    Kest         = list(Global = spatstat.explore::Kest,
                        Local = spatstat.explore::localK,
                        L = FALSE, inhom = FALSE, indices = 0),
    Lest         = list(Global = spatstat.explore::Lest,
                        Local = spatstat.explore::localK,
                        L = TRUE,  inhom = FALSE, indices = 0),
    pcfinhom     = list(Global = spatstat.explore::pcfinhom,
                        Local = spatstat.explore::localpcfinhom,
                        L = FALSE, inhom = TRUE,  indices = 0),
    Kinhom       = list(Global = spatstat.explore::Kinhom,
                        Local = spatstat.explore::localKinhom,
                        L = FALSE, inhom = TRUE,  indices = 0),
    Linhom       = list(Global = spatstat.explore::Linhom,
                        Local = spatstat.explore::localKinhom,
                        L = TRUE,  inhom = TRUE,  indices = 0),
    Kcross       = list(Global = spatstat.explore::Kcross,
                        Local = spatstat.explore::localKcross,
                        L = FALSE, inhom = FALSE, indices = 2),
    Lcross       = list(Global = spatstat.explore::Lcross,
                        Local = spatstat.explore::localKcross,
                        L = TRUE,  inhom = FALSE, indices = 2),
    Kdot         = list(Global = spatstat.explore::Kdot,
                        Local = spatstat.explore::localKdot,
                        L = FALSE, inhom = FALSE, indices = 1),
    Ldot         = list(Global = spatstat.explore::Ldot,
                        Local = spatstat.explore::localKdot,
                        L = TRUE,  inhom = FALSE, indices = 1),
    Kcross.inhom = list(Global = spatstat.explore::Kcross.inhom,
                        Local = spatstat.explore::localKcross.inhom,
                        L = FALSE, inhom = TRUE,  indices = 2),
    Lcross.inhom = list(Global = spatstat.explore::Lcross.inhom,
                        Local = spatstat.explore::localKcross.inhom,
                        L = TRUE,  inhom = TRUE,  indices = 2)
  )
  if (length(key) != 1L) {
    stop("Argument must be a single character string or function", call. = FALSE)
  }
  nama <- names(TheTable)
  pos <- if (is.character(key)) {
    match(key, nama)
  } else if (is.function(key)) {
    match(list(key), lapply(TheTable, getElement, name = "Global"))
  } else {
    NULL
  }
  if (is.na(pos)) return(NULL)
  out <- TheTable[[pos]]
  out$GlobalName <- nama[pos]
  out
}

#' Internal: precompute local composition enrichment cache
#'
#' Build per-radius local composition enrichment matrices for a point pattern.
#' For each point and radius, store local neighbor fractions/counts and the
#' overlap area between the radius-disc and tissue window so enrichment can be
#' computed as an edge-corrected proportion difference.
#'
#' @param X A marked \code{ppp} object.
#' @param radii_um Numeric radii at which to precompute local composition.
#'
#' @return A list containing precomputed local composition matrices and indexing
#'  structures used by local composition bootstrap routines.
#'
#' @keywords internal
.precompute_local_comp_cache <- function(X, radii_um) {
  if (!spatstat.geom::is.ppp(X)) {
    stop("`X` must be a spatstat ppp object.", call. = FALSE)
  }
  if (!is.numeric(radii_um) || length(radii_um) < 1L || any(!is.finite(radii_um))) {
    stop("`radii_um` must be a finite numeric vector.", call. = FALSE)
  }

  radii_sorted <- sort(unique(as.numeric(radii_um)))
  n <- spatstat.geom::npoints(X)
  if (n < 1L) {
    stop("Point pattern has no points.", call. = FALSE)
  }

  marks_chr <- as.character(spatstat.geom::marks(X))
  cell_types <- sort(unique(marks_chr))
  type_idx <- match(marks_chr, cell_types)
  n_types <- length(cell_types)
  W <- spatstat.geom::Window(X)
  area_w <- spatstat.geom::area.owin(W)
  if (!is.finite(area_w) || area_w <= 0) {
    stop("Point pattern window has non-positive area.", call. = FALSE)
  }

  global_prop <- tabulate(type_idx, nbins = n_types) / n
  names(global_prop) <- cell_types
  global_density <- tabulate(type_idx, nbins = n_types) / area_w
  names(global_density) <- cell_types

  indices_by_type <- split(seq_len(n), marks_chr)
  for (ct in cell_types) {
    if (is.null(indices_by_type[[ct]])) {
      indices_by_type[[ct]] <- integer(0)
    }
  }

  counts <- matrix(0, nrow = n, ncol = n_types)

  edge_i <- integer(0)
  edge_j <- integer(0)
  edge_d <- numeric(0)
  max_r <- max(radii_sorted)
  max_expected_pairs <- getOption("panoramic.local_comp_max_expected_pairs", 1e8)
  skip_edges <- FALSE
  if (is.finite(max_r) && max_r > 0 && n >= 2L) {
    # Heuristic guardrail: avoid pathological close-pair enumeration that can
    # stall or overflow at very high densities.
    p_in_radius <- min(1, pi * max_r^2 / area_w)
    est_pairs <- (n * (n - 1) / 2) * p_in_radius
    if (is.finite(max_expected_pairs) &&
        max_expected_pairs > 0 &&
        is.finite(est_pairs) &&
        est_pairs > max_expected_pairs) {
      warning(
        sprintf(
          paste0(
            "local_comp_enrichment cache skipped closepairs due to high estimated pair count ",
            "(n=%d, rmax=%.3f, est_pairs=%.3g, threshold=%.3g). ",
            "Returning NA local-composition estimates for this sample; ",
            "increase option 'panoramic.local_comp_max_expected_pairs' to override."
          ),
          n, max_r, est_pairs, max_expected_pairs
        ),
        call. = FALSE
      )
      skip_edges <- TRUE
    }
    if (!skip_edges) {
      cp <- try(
        suppressWarnings(spatstat.geom::closepairs(X, rmax = max_r, what = "ijd")),
        silent = TRUE
      )
      if (!inherits(cp, "try-error") && length(cp$i) > 0L) {
        ord <- order(cp$d)
        edge_i <- cp$i[ord]
        edge_j <- cp$j[ord]
        edge_d <- cp$d[ord]
      } else if (inherits(cp, "try-error")) {
        warning(
          "local_comp_enrichment cache closepairs failed; returning NA local-composition estimates for this sample.",
          call. = FALSE
        )
      }
    }
  }

  local_fraction <- vector("list", length(radii_sorted))
  names(local_fraction) <- as.character(radii_sorted)
  local_total_neighbors <- vector("list", length(radii_sorted))
  names(local_total_neighbors) <- as.character(radii_sorted)
  local_overlap_area <- spatstat.geom::discpartarea(X, r = radii_sorted, W = W)
  if (!is.matrix(local_overlap_area)) {
    local_overlap_area <- matrix(
      local_overlap_area,
      nrow = n,
      ncol = length(radii_sorted)
    )
  }

  edge_ptr <- 0L
  n_edges <- length(edge_i)
  for (k in seq_along(radii_sorted)) {
    r <- radii_sorted[k]
    while (edge_ptr < n_edges && edge_d[edge_ptr + 1L] <= r) {
      ii <- edge_i[edge_ptr + 1L]
      jj <- edge_j[edge_ptr + 1L]
      ti <- type_idx[ii]
      tj <- type_idx[jj]
      if (is.finite(tj)) counts[ii, tj] <- counts[ii, tj] + 1
      if (is.finite(ti)) counts[jj, ti] <- counts[jj, ti] + 1
      edge_ptr <- edge_ptr + 1L
    }

    total_neighbors <- rowSums(counts)
    frac <- matrix(NA_real_, nrow = n, ncol = n_types)
    ok <- total_neighbors > 0
    if (any(ok)) {
      frac[ok, ] <- counts[ok, , drop = FALSE] / total_neighbors[ok]
    }

    local_fraction[[k]] <- frac
    local_total_neighbors[[k]] <- total_neighbors
  }

  list(
    radii_um = radii_sorted,
    cell_types = cell_types,
    global_prop = global_prop,
    global_density = global_density,
    indices_by_type = indices_by_type,
    local_fraction = local_fraction,
    local_total_neighbors = local_total_neighbors,
    local_overlap_area = local_overlap_area,
    area_window = area_w,
    x = X$x,
    y = X$y,
    window = W
  )
}

#' Internal: local composition enrichment score
#'
#' @return Numeric vector of local enrichment scores (percentage-point
#'  difference between observed and expected local composition).
#'
#' @keywords internal
.local_comp_score <- function(
    frac,
    n_neighbors,
    overlap_area,
    global_density
) {
  if (!is.finite(global_density) || global_density < 0) {
    return(rep(NA_real_, length(frac)))
  }
  out <- rep(NA_real_, length(frac))
  ok <- is.finite(frac) & is.finite(n_neighbors) & n_neighbors > 0 &
    is.finite(overlap_area) & overlap_area >= 0
  if (!any(ok)) return(out)
  obs_prop <- pmax(pmin(frac[ok], 1), 0)
  expected_count <- global_density * overlap_area[ok]
  expected_prop <- expected_count / n_neighbors[ok]
  expected_prop <- pmax(pmin(expected_prop, 1), 0)
  out[ok] <- 100 * (obs_prop - expected_prop)
  out
}

#' Internal: Loh-style bootstrap for local composition enrichment
#'
#' @param X A marked \code{ppp} object.
#' @param from,to Character marks defining anchor and target cell types.
#' @param radii_um Numeric radii at which enrichment is evaluated.
#' @param nsim Integer bootstrap replicates.
#' @param confidence Numeric confidence level.
#' @param type Quantile type for bootstrap intervals.
#' @param boot One of \code{"approx"} or \code{"block"}.
#' @param nx,ny Integer block counts (used only for \code{boot = "block"}).
#' @param local_cache Optional cache from \code{.precompute_local_comp_cache()}.
#'
#' @return A data.frame compatible with \code{.summarize_lohboot()}.
#'
#' @keywords internal
.lohboot_local_comp <- function(
    X,
    from,
    to,
    radii_um,
    nsim = 200,
    confidence = 0.95,
    type = 7,
    boot = c("approx", "block"),
    nx = 4,
    ny = nx,
    local_cache = NULL
) {
  # PANORAMIC-specific extension (not provided by spatstat::lohboot):
  # bootstrap for local composition enrichment using cached neighborhood
  # geometry and optional overlap-weighted block resampling.
  boot <- match.arg(boot)
  if (!spatstat.geom::is.ppp(X)) {
    stop("`X` must be a spatstat ppp object.", call. = FALSE)
  }
  if (!is.numeric(radii_um) || length(radii_um) < 1L || any(!is.finite(radii_um))) {
    stop("`radii_um` must be a finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim <= 1L) {
    stop("`nsim` must be an integer > 1.", call. = FALSE)
  }
  if (!(confidence > 0.5 && confidence < 1)) {
    stop("`confidence` must be between 0.5 and 1.", call. = FALSE)
  }

  if (is.null(local_cache)) {
    local_cache <- .precompute_local_comp_cache(X, radii_um)
  }

  anchor_idx <- local_cache$indices_by_type[[from]]
  to_idx <- match(to, local_cache$cell_types)
  radii_match <- match(radii_um, local_cache$radii_um)
  if (anyNA(radii_match)) {
    # Be tolerant to floating-point representation drift when matching radii.
    radii_match <- vapply(radii_um, function(r) {
      tol <- sqrt(.Machine$double.eps) * max(1, abs(r))
      idx <- which(abs(local_cache$radii_um - r) <= tol)
      if (length(idx) == 1L) idx else NA_integer_
    }, integer(1))
  }

  nr <- length(radii_um)
  if (length(anchor_idx) < 1L || is.na(to_idx) || anyNA(radii_match)) {
    return(data.frame(
      r = radii_um,
      theo = rep(0, nr),
      trans = rep(NA_real_, nr),
      lo = rep(NA_real_, nr),
      hi = rep(NA_real_, nr),
      var = rep(NA_real_, nr)
    ))
  }

  n_anchor <- length(anchor_idx)
  global_density <- local_cache$global_density[to_idx]
  if (identical(from, to) &&
      is.finite(local_cache$area_window) &&
      local_cache$area_window > 0) {
    n_to <- length(local_cache$indices_by_type[[to]])
    global_density <- max((n_to - 1) / local_cache$area_window, 0)
  }
  y <- matrix(NA_real_, nrow = nr, ncol = n_anchor)
  for (k in seq_len(nr)) {
    frac <- local_cache$local_fraction[[radii_match[k]]][anchor_idx, to_idx]
    n_neighbors <- local_cache$local_total_neighbors[[radii_match[k]]][anchor_idx]
    overlap_area <- local_cache$local_overlap_area[anchor_idx, radii_match[k]]
    y[k, ] <- .local_comp_score(
      frac = frac,
      n_neighbors = n_neighbors,
      overlap_area = overlap_area,
      global_density = global_density
    )
  }

  if (!any(is.finite(y))) {
    return(data.frame(
      r = radii_um,
      theo = rep(0, nr),
      trans = rep(NA_real_, nr),
      lo = rep(NA_real_, nr),
      hi = rep(NA_real_, nr),
      var = rep(NA_real_, nr)
    ))
  }

  ymean <- rep(NA_real_, nr)
  ystar <- matrix(NA_real_, nrow = nr, ncol = nsim)

  if (boot == "approx") {
    ymean <- rowMeans(y, na.rm = TRUE)
    ymean[is.nan(ymean)] <- NA_real_
    for (b in seq_len(nsim)) {
      ind <- sample.int(n_anchor, size = n_anchor, replace = TRUE)
      ystar[, b] <- rowMeans(y[, ind, drop = FALSE], na.rm = TRUE)
    }
  } else {
    GridTess <- spatstat.geom::quadrats(
      spatstat.geom::boundingbox(local_cache$window),
      nx = nx,
      ny = ny
    )
    BlockIndex <- spatstat.geom::tileindex(
      local_cache$x[anchor_idx],
      local_cache$y[anchor_idx],
      GridTess
    )

    blocks <- spatstat.geom::tiles(GridTess)
    tile_area <- vapply(blocks, spatstat.geom::area.owin, numeric(1))
    overlap_area <- vapply(blocks, function(b) {
      iw <- try(spatstat.geom::intersect.owin(b, local_cache$window), silent = TRUE)
      if (inherits(iw, "try-error")) return(0)
      spatstat.geom::area.owin(iw)
    }, numeric(1))

    keep_blocks <- is.finite(overlap_area) & is.finite(tile_area) &
      overlap_area > 0 & tile_area > 0
    if (sum(keep_blocks) < 2L) {
      warning(
        sprintf(
          "local_comp_enrichment block bootstrap: <2 usable tiles for pair '%s' -> '%s'; returning NA variance/bands.",
          from, to
        ),
        call. = FALSE
      )
      ymean <- rowMeans(y, na.rm = TRUE)
      ymean[is.nan(ymean)] <- NA_real_
    } else {
      indexmap <- cumsum(keep_blocks)
      indexmap[!keep_blocks] <- NA_integer_
      BlockIndex <- indexmap[BlockIndex]
      keep_pts <- !is.na(BlockIndex)
      y <- y[, keep_pts, drop = FALSE]
      BlockIndex <- BlockIndex[keep_pts]
      nmarks <- sum(keep_blocks)
      BlockFactor <- factor(BlockIndex, levels = seq_len(nmarks))

      weights <- overlap_area[keep_blocks] / tile_area[keep_blocks]
      if (any(!is.finite(weights)) || sum(weights) <= 0) {
        warning(
          sprintf(
            "local_comp_enrichment block bootstrap: invalid tile weights for pair '%s' -> '%s'; returning NA variance/bands.",
            from, to
          ),
          call. = FALSE
        )
        ymean <- rowMeans(y, na.rm = TRUE)
        ymean[is.nan(ymean)] <- NA_real_
      } else {
        ymarks <- by(t(y), BlockFactor, colSums, na.rm = TRUE, simplify = FALSE)
        bad <- vapply(ymarks, function(z) is.null(z) || length(z) != nr, logical(1))
        if (any(bad)) {
          ymarks[bad] <- rep(list(numeric(nr)), sum(bad))
        }
        ymarks <- as.matrix(do.call(cbind, ymarks))
        block_counts_kept <- tabulate(BlockIndex, nbins = nmarks)
        yblock_mean <- sweep(ymarks, 2, pmax(block_counts_kept, 1), FUN = "/")

        ymean <- as.numeric((yblock_mean %*% weights) / sum(weights))
        for (b in seq_len(nsim)) {
          ind <- sample(seq_len(nmarks), replace = TRUE, prob = weights)
          ystar[, b] <- rowMeans(yblock_mean[, ind, drop = FALSE], na.rm = TRUE)
        }
      }
    }
  }

  ystar[is.nan(ystar)] <- NA_real_
  ymean[is.nan(ymean)] <- NA_real_

  has_boot <- any(is.finite(ystar))
  # Bootstrap-centered point estimate: mean across bootstrap replicate means.
  # Falls back to deterministic ymean when bootstrap replicates are unavailable.
  ycenter <- rowMeans(ystar, na.rm = TRUE)
  ycenter[is.nan(ycenter)] <- NA_real_
  if (!any(is.finite(ycenter))) ycenter <- ymean

  if (!has_boot) {
    lo <- rep(NA_real_, nr)
    hi <- rep(NA_real_, nr)
    var_y <- rep(NA_real_, nr)
  } else {
    # Use simultaneous (global) bands, aligned with how PANORAMIC consumes
    # Loh-bootstrap curves for downstream interpolation and meta-analysis.
    alpha <- 1 - confidence
    ydev <- apply(abs(sweep(ystar, 1L, ycenter)), 2L, max, na.rm = TRUE)
    ydev[!is.finite(ydev)] <- NA_real_
    crit <- stats::quantile(
      ydev,
      probs = 1 - alpha,
      na.rm = TRUE,
      type = type,
      names = FALSE
    )
    if (!is.finite(crit)) crit <- NA_real_
    lo <- ycenter - crit
    hi <- ycenter + crit

    var_y <- apply(ystar, 1L, stats::var, na.rm = TRUE)
    var_y[is.nan(var_y)] <- NA_real_
  }

  data.frame(
    r = radii_um,
    theo = rep(0, nr),
    trans = ycenter,
    lo = lo,
    hi = hi,
    var = var_y
  )
}

#' Internal: custom Loh bootstrap with overlap-weighted blocks
#'
#' @param X A marked \code{ppp} object.
#' @param fun Character/function passed to spatstat local summary mapping.
#' @param ... Additional arguments forwarded to the selected local function.
#' @param global Logical; if \code{TRUE}, compute global envelopes.
#' @param basicboot Logical; if \code{TRUE}, use basic bootstrap intervals.
#' @param Vcorrection Logical; apply variance correction for K-type functions.
#' @param confidence Numeric confidence level in \code{(0.5, 1)}.
#' @param nx,ny Integer numbers of blocks in x/y for block bootstrap.
#' @param nsim Integer bootstrap replicates.
#' @param type Quantile type.
#'
#' @return A data.frame with columns including \code{r}, theoretical curve,
#'  estimate, and lower/upper interval bounds.
.lohboot_block_weighted <- function(
    X, fun = c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom",
               "Kcross", "Lcross", "Kdot", "Ldot",
               "Kcross.inhom", "Lcross.inhom"),
    ..., global = FALSE, basicboot = FALSE, Vcorrection = FALSE,
    confidence = 0.95, nx = 4, ny = nx, nsim = 200, type = 7
) {
  # Relative to spatstat.explore::lohboot(), this internal variant:
  # 1) weights each block by the fraction of tile area overlapping the
  #    observed tissue window (overlap_area / tile_area),
  # 2) drops blocks with zero overlap or zero points before resampling,
  # 3) reuses this weighted backend consistently across PANORAMIC pipeline
  #    functions to reduce boundary-driven bias under irregular windows.
  stopifnot(spatstat.geom::is.ppp(X))
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim <= 1L) {
    stop("`nsim` must be an integer > 1.", call. = FALSE)
  }
  fun.name <- deparse(substitute(fun))
  if (is.character(fun)) fun <- match.arg(fun)
  info <- .spatstat_local_function_info(fun)
  if (is.null(info)) {
    stop(paste("Loh's bootstrap is not supported for the function", sQuote(fun.name)),
         call. = FALSE)
  }
  fun      <- info$GlobalName
  localfun <- info$Local

  if (!(confidence > 0.5 && confidence < 1)) {
    stop("`confidence` must be between 0.5 and 1.", call. = FALSE)
  }
  alpha <- 1 - confidence
  if (!global) {
    probs <- c(alpha / 2, 1 - alpha / 2)
    rank <- nsim * probs[2L]
  } else {
    probs <- 1 - alpha
    rank <- nsim * probs
  }
  if (abs(rank - round(rank)) > 0.001) {
    warning(
      paste("confidence level", confidence,
            "corresponds to a non-integer rank", rank,
            "so quantiles will be interpolated"),
      call. = FALSE
    )
  }

  f <- localfun(X, ...)
  theo <- f$theo
  correction <- attr(f, "correction")
  switch(correction,
         none = { ckey <- clab <- "un"; cadj <- "uncorrected" },
         border = { ckey <- "border"; clab <- "bord"; cadj <- "border-corrected" },
         translate = { ckey <- clab <- "trans"; cadj <- "translation-corrected" },
         isotropic = { ckey <- clab <- "iso"; cadj <- "Ripley isotropic corrected" })

  types <- levels(spatstat.geom::marks(X))
  from <- spatstat.utils::resolve.1.default(list(from = types[1]), list(...))
  to <- spatstat.utils::resolve.1.default(list(to = types[2]), list(...))
  fromName <- spatstat.utils::make.parseable(paste(from))
  toName <- spatstat.utils::make.parseable(paste(to))
  if (info$indices > 0) {
    X <- attr(f, "Xfrom")
  }

  n <- spatstat.geom::npoints(X)
  y <- as.matrix(as.data.frame(f))[, seq_len(n), drop = FALSE]
  nr <- nrow(y)

  W <- spatstat.geom::Window(X)
  GridTess <- spatstat.geom::quadrats(spatstat.geom::boundingbox(W), nx = nx, ny = ny)
  BlockIndex <- spatstat.geom::tileindex(X$x, X$y, GridTess)

  blocks <- spatstat.geom::tiles(GridTess)
  tile_area <- vapply(blocks, spatstat.geom::area.owin, numeric(1))
  overlap_area <- vapply(blocks, function(b) {
    iw <- try(spatstat.geom::intersect.owin(b, W), silent = TRUE)
    if (inherits(iw, "try-error")) return(0)
    spatstat.geom::area.owin(iw)
  }, numeric(1))

  block_counts <- tabulate(BlockIndex, nbins = length(blocks))
  keep_blocks <- overlap_area > 0 & block_counts > 0
  if (sum(keep_blocks) < 2) {
    stop("Not enough blocks with overlap and at least one point.", call. = FALSE)
  }

  indexmap <- cumsum(keep_blocks)
  indexmap[!keep_blocks] <- NA_integer_
  BlockIndex <- indexmap[BlockIndex]

  keep_pts <- !is.na(BlockIndex)
  y <- y[, keep_pts, drop = FALSE]
  n <- sum(keep_pts)
  BlockIndex <- BlockIndex[keep_pts]
  BlockFactor <- factor(BlockIndex, levels = unique(BlockIndex))

  nmarks <- length(levels(BlockFactor))
  weights <- overlap_area[keep_blocks] / tile_area[keep_blocks]
  weights <- weights[seq_len(nmarks)]
  if (any(!is.finite(weights)) || sum(weights) <= 0) {
    stop("Non-finite or zero overlap weights.", call. = FALSE)
  }

  ymarks <- by(t(y), BlockFactor, colSums, na.rm = TRUE, simplify = FALSE)
  if (any(isempty <- vapply(ymarks, is.null, logical(1)))) {
    ymarks[isempty] <- rep(list(numeric(nr)), sum(isempty))
  }
  ymarks <- as.matrix(do.call(cbind, ymarks))
  block_counts_kept <- as.numeric(table(BlockFactor))
  yblock_mean <- sweep(ymarks, 2, pmax(block_counts_kept, 1), FUN = "/")

  ymean <- as.numeric((yblock_mean %*% weights) / sum(weights))
  ystar <- matrix(, nrow = nr, ncol = nsim)
  for (i in seq_len(nsim)) {
    ind <- sample(seq_len(nmarks), replace = TRUE, prob = weights)
    ystar[, i] <- rowMeans(yblock_mean[, ind, drop = FALSE], na.rm = TRUE)
  }

  if (!global) {
    hilo <- apply(ystar, 1, stats::quantile, probs = probs, na.rm = TRUE, type = type)
    if (Vcorrection && (fun == "Kest" || fun == "Kinhom")) {
      Vcov <- sqrt(1 + 2 * pi * n * (f$r)^2 / spatstat.geom::area.owin(W))
      hilo[1L, ] <- ymean + (ymean - hilo[1L, ]) / Vcov
      hilo[2L, ] <- ymean + (ymean - hilo[2L, ]) / Vcov
      hilo <- hilo[2:1, ]
      basicboot <- FALSE
    }
    if (basicboot) {
      hilo[1L, ] <- 2 * ymean - hilo[1L, ]
      hilo[2L, ] <- 2 * ymean - hilo[2L, ]
      hilo <- hilo[c(2, 1), ]
    }
  } else {
    ydif <- sweep(ystar, 1, ymean)
    ydev <- apply(abs(ydif), 2, max, na.rm = TRUE)
    crit <- stats::quantile(ydev, probs = probs, na.rm = TRUE, type = type)
    hilo <- rbind(ymean - crit, ymean + crit)
  }

  if (info$L) {
    theo <- sqrt(theo / pi)
    ymean <- sqrt(ymean / pi)
    hilo <- sqrt(hilo / pi)
    warning(
      paste(
        "The calculation of confidence intervals for L functions",
        "in lohboot() has changed in spatstat 1.60-0 and later;",
        "they are now computed by transforming the confidence intervals",
        "for the corresponding K functions."
      ),
      call. = FALSE
    )
  }

  df <- data.frame(
    r = f$r,
    theo = theo,
    ymean = ymean,
    lo = hilo[1L, ],
    hi = hilo[2L, ]
  )
  colnames(df)[3L] <- ckey
  df
}

#' Internal: quiet Loh bootstrap wrapper
#' 
#' Run \code{spatstat.explore::lohboot()} while suppressing console output, 
#' returning the bootstrap object invisibly for further summarization. 
#' 
#' @param ... Arguments passed to \code{spatstat.explore::lohboot()}. 
#' @param boot Character. Bootstrap mode: \code{"approx"} (no tiling) or
#'  \code{"block"} (tiled Loh bootstrap with overlap-weighted blocks).
#' @param tile_size Optional numeric scalar giving tile size in the same units
#'  as spatial coordinates. Used only when \code{boot = "block"}; translated
#'  into \code{nx} and \code{ny} for
#'  \code{lohboot()} based on the sample's window bounding box. Cannot be used
#'  with \code{nx}/\code{ny}.
#' @param nx,ny Optional integers giving the number of tiles in x/y directions
#'  for \code{lohboot()} when \code{boot = "block"}. Ignored for
#'  \code{boot = "approx"}.
#' 
#' @return A Loh bootstrap object as returned by \code{spatstat.explore::lohboot()}
#'  or a compatible data.frame for \code{boot = "block"}.
#' 
#' @keywords internal
.lohboot_quiet <- function(X, fun, ..., verbose = FALSE,
                           boot = c("approx", "block"),
                           tile_size = NULL, nx = NULL, ny = NULL) {
  boot <- match.arg(boot)
  is_local_comp <- is.character(fun) &&
    length(fun) == 1L &&
    identical(fun, "local_comp_enrichment")
  if (!is.null(tile_size)) {
    if (!is.numeric(tile_size) || length(tile_size) != 1L ||
        !is.finite(tile_size) || tile_size <= 0) {
      stop("`tile_size` must be a positive numeric scalar.")
    }
  }
  if (boot == "approx") {
    if (!is.null(tile_size) || !is.null(nx) || !is.null(ny)) {
      stop("`tile_size`/`nx`/`ny` are only valid when boot = \"block\".")
    }
    args <- if (is_local_comp) {
      list(X = X, fun = fun, ...)
    } else {
      list(X, fun, ..., block = FALSE)
    }
  } else {
    if (!is.null(tile_size) && (!is.null(nx) || !is.null(ny))) {
      stop("Provide either `tile_size` or `nx`/`ny`, not both.")
    }
    if (!is.null(tile_size)) {
      bb <- spatstat.geom::boundingbox(X)
      width <- diff(bb$xrange)
      height <- diff(bb$yrange)
      if (!is.finite(width) || !is.finite(height)) {
        stop("Non-finite window bounds; cannot derive tile counts.")
      }
      nx <- max(1L, as.integer(ceiling(width / tile_size)))
      ny <- max(1L, as.integer(ceiling(height / tile_size)))
    }
    args <- list(X = X, fun = fun, ..., nx = nx, ny = ny)
  }
  if (isTRUE(verbose)) {
    if (is_local_comp) {
      args$boot <- boot
      args$fun <- NULL
      return(do.call(.lohboot_local_comp, args))
    }
    if (boot == "block") return(do.call(.lohboot_block_weighted, args))
    return(do.call(spatstat.explore::lohboot, args))
  }
  res <- NULL
  suppressWarnings(utils::capture.output({
    if (is_local_comp) {
      args$boot <- boot
      args$fun <- NULL
      res <- do.call(.lohboot_local_comp, args)
    } else if (boot == "block") {
      res <- do.call(.lohboot_block_weighted, args)
    } else {
      res <- do.call(spatstat.explore::lohboot, args)
    }
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
#' @param stat Character. Summary statistic. Default is
#'  \code{"local_comp_enrichment"}. Other options are \code{"Lcross"},
#'  \code{"Kcross"}, \code{"Lest"}, and \code{"Kest"}.
#' @param nsim Integer. Number of Loh bootstrap simulations. 
#' @param correction Character. Edge correction passed to spatstat.
#' @param radii_um Numeric vector of radii (microns) on which to summarize.
#' @param boot Character. Bootstrap mode: \code{"approx"} (no tiling) or
#'  \code{"block"} (tiled Loh bootstrap).
#' @param tile_size Optional numeric scalar giving tile size in the same units
#'  as spatial coordinates. Used only when \code{boot = "block"}.
#' @param nx,ny Optional integers giving the number of tiles in x/y directions
#'  when \code{boot = "block"}.
#' @param local_comp_cache Optional cache object from
#'  \code{.precompute_local_comp_cache()} to avoid repeated neighbor searches
#'  when \code{stat = "local_comp_enrichment"}.
#' 
#' @return A data.frame with columns \code{radius_um}, \code{yi}, and \code{vi}, filled with \code{NA} if insufficient cells of either type are present. 
#' 
#' @keywords internal
.one_pair_one_sample <- function(meta, ct1, ct2, stat = "local_comp_enrichment",
                                 nsim = 100, correction = "translate",
                                 radii_um, verbose = FALSE,
                                 boot = c("approx", "block"),
                                 tile_size = NULL, nx = NULL, ny = NULL,
                                 local_comp_cache = NULL) {
  boot <- match.arg(boot)
  is_L <- grepl("^L", stat)
  stat_use <- if (is_L) sub("^L", "K", stat) else stat

  tab <- meta$marks_tab
  n1 <- unname(tab[ct1]); n2 <- unname(tab[ct2])
  min_to <- if (identical(stat_use, "local_comp_enrichment")) 1L else 2L
  if (is.na(n1) || n1 < 2L || is.na(n2) || n2 < min_to) {
    return(data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_))
  }
  
  X <- meta$ppp
  
  if (identical(stat_use, "local_comp_enrichment")) {
    loh <- .lohboot_quiet(
      X,
      "local_comp_enrichment",
      from = ct1,
      to = ct2,
      radii_um = radii_um,
      nsim = nsim,
      verbose = verbose,
      boot = boot,
      tile_size = tile_size,
      nx = nx,
      ny = ny,
      local_cache = local_comp_cache
    )
  } else if (!identical(ct1, ct2) && identical(stat_use, "Kcross")) {
    fun <- spatstat.explore::Kcross
    loh <- .lohboot_quiet(
      spatstat.geom::subset.ppp(X, marks %in% c(ct1, ct2)),
      fun, from = ct1, to = ct2, correction = correction, global = TRUE, nsim = nsim,
      verbose = verbose, boot = boot, tile_size = tile_size, nx = nx, ny = ny
    )
  } else {
    fun <- if (grepl("^L", stat_use)) spatstat.explore::Lest else spatstat.explore::Kest
    loh <- .lohboot_quiet(
      spatstat.geom::subset.ppp(X, marks == ct1),
      fun, correction = correction, global = TRUE, nsim = nsim,
      verbose = verbose, boot = boot, tile_size = tile_size, nx = nx, ny = ny
    )
  }
  
  center_curve <- !is_L && !identical(stat_use, "local_comp_enrichment")
  df <- try(.summarize_lohboot(loh, center_L = center_curve), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df) || sum(is.finite(df$r)) < 1L) {
    return(data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_))
  }
  
  if (is_L) {
    # Transform K -> L and compute variance via delta method
    k <- df$yi
    v <- df$vi
    l <- sqrt(k / pi)
    v_l <- rep(NA_real_, length(k))
    ok <- is.finite(k) & is.finite(v) & k > 0
    v_l[ok] <- v[ok] / (4 * pi * k[ok])
    # Center L by r (theoretical L)
    l <- l - df$r
    df$yi <- l
    df$vi <- v_l
  }
  
  o  <- order(df$r)
  yi <- .safe_approx(df$r[o], df$yi[o], radii_um)
  vi <- .safe_approx(df$r[o], df$vi[o], radii_um)
  
  data.frame(radius_um = radii_um, yi = yi, vi = vi)
}

#' Compute pairwise spatial statistics for PANORAMIC
#' 
#' Compute pairwise spatial summary curves and bootstrap variances for all
#' requested cell-type pairs and radii across multiple prepared samples,
#' returning a SummarizedExperiment. By default, PANORAMIC computes
#' \code{"local_comp_enrichment"}.
#' 
#' @param prep List of prepared SpatialExperiment objects as returned by \code{panoramic_prepare()}, 
#'  each with \code{metadata(spe)$panoramic}. 
#' @param pairs Either "auto" (generate all cell-type pairs including self-pairs) 
#'  or a data.frame with columns \code{ct1}, \code{ct2}. 
#' @param radii_um Numeric vector of radii (microns) at which to evaluate the colocalization statistic. 
#' @param stat Character. Summary statistic. Default is
#'  \code{"local_comp_enrichment"}. Other supported values are
#'  \code{"Lcross"}, \code{"Lest"}, \code{"Kcross"}, and \code{"Kest"}.
#' @param nsim Integer. Number of Loh bootstrap simulations per sample/pair. 
#' @param boot Character. Bootstrap mode: \code{"approx"} (no tiling) or
#'  \code{"block"} (tiled Loh bootstrap).
#' @param tile_size Optional numeric scalar giving tile size in the same units
#'  as spatial coordinates. Used only when \code{boot = "block"}.
#' @param nx,ny Optional integers giving the number of tiles in x/y directions
#'  when \code{boot = "block"}.
#' @param correction Chatacter. Edge correction method for spatstat ("translate", "border", ...). 
#' @param seed Integer random seed for reproducible bootstrap simulations.
#' @param BPPARAM A BiocParallelParam object controlling parallelisation across samples.
#' @param verbose Logical. If \code{TRUE}, show output from bootstrap calls.
#'  If \code{FALSE} (default), suppress bootstrap console output.
#'
#' @return A SummarizedExperiment with : 
#' \itemize{
#'  \item assay "yi": centered estimates per (ct1, ct2, radius) feature and sample. 
#'  \item assay "vi": variance estimates aligned to "yi".
#' }
#' rowData stores \code{ct1}, \code{ct2}, \code{radius_um}, \code{stat}
#' colData stores \code{sample} and \code{group}. 
#'
#' @details
#' When \code{stat} is an L-function, variance is computed via the delta method
#' from the corresponding K-function estimates and then centered by \code{r}.
#' For \code{stat = "local_comp_enrichment"}, enrichment is computed for each
#' anchor cell as observed minus expected local target proportion within radius
#' \code{r}. The expected proportion is edge-corrected by first computing the
#' expected target count in the overlap area \eqn{|B_r(ct1) \cap W|} and then
#' dividing by the observed local neighbor count:
#' \deqn{100\cdot\left(\hat{p}_{obs}(ct2 \mid ct1, r) - \hat{p}_{exp}^{edge}(ct2 \mid ct1, r)\right)}
#' where \eqn{\hat{p}_{exp}^{edge} = (\lambda_{ct2}\cdot|B_r(ct1)\cap W|)/N_{obs,\cdot}(ct1,r)}.
#' The null expectation is 0 and units are percentage points.
#' Local neighbor composition and overlap areas are precomputed per
#' sample/radius to avoid repeated geometry/neighborhood calculations across
#' cell-type pairs.
#' 
#' @examples
#' library(BiocParallel)
#'
#' spe_list <- list(
#'   sample1 = panoramic_simulate_spe(
#'     n_cells = 40,
#'     sample_id = "sample1",
#'     cell_types = c("A", "B"),
#'     scenario = "random",
#'     seed = 1
#'   ),
#'   sample2 = panoramic_simulate_spe(
#'     n_cells = 40,
#'     sample_id = "sample2",
#'     cell_types = c("A", "B"),
#'     scenario = "random",
#'     seed = 2
#'   )
#' )
#'
#' design <- data.frame(
#'   sample = names(spe_list),
#'   group  = c("group1", "group1"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Prepare inputs for PANORAMIC
#' prep <- panoramic_prepare(
#'   spe_list,
#'   design      = design,
#'   cell_type   = "cell_type",
#'   min_cells   = 2,
#'   window      = "rect",
#'   BPPARAM     = BiocParallel::SerialParam()
#' )
#'
#' # Compute spatial statistics across samples
#' se_stats <- panoramic_spatialstats(
#'   prep      = prep,
#'   pairs     = "auto",
#'   radii_um  = c(10, 20),
#'   nsim      = 5,
#'   correction = "translate",
#'   seed      = 1,
#'   BPPARAM   = BiocParallel::SerialParam()
#' )
#'
#' se_stats
#' @export
panoramic_spatialstats <- function(
    prep, pairs = "auto", radii_um, stat = "local_comp_enrichment",
    nsim = 100, correction = "translate", seed = 123,
    boot = c("approx", "block"), tile_size = NULL, nx = NULL, ny = NULL,
    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
) {
  boot <- match.arg(boot)
  stopifnot(is.list(prep), length(prep) >= 1, length(radii_um) >= 1)
  if (stat %in% c("local_composition_enrichment", "local_comp")) {
    stat <- "local_comp_enrichment"
  }
  valid_stats <- c("local_comp_enrichment", "Lcross", "Kcross", "Lest", "Kest")
  if (!stat %in% valid_stats) {
    stop(
      "`stat` must be one of: ",
      paste(valid_stats, collapse = ", "),
      "."
    )
  }
  if (grepl("^L", stat)) {
    warning(
      "L-function variance is computed via delta method from the corresponding K function; ",
      "estimates are centered by r.",
      call. = FALSE
    )
  }
  if (identical(pairs, "auto")) {
    pairs <- .panoramic_pairs(prep, include_self = TRUE)
  } else {
    stopifnot(all(c("ct1","ct2") %in% colnames(pairs)))
  }
  
  samples <- names(prep)
  
  # Compute per-sample per-pair curves
  compute_one_sample <- function(sid) {
    meta <- S4Vectors::metadata(prep[[sid]])$panoramic
    local_comp_cache <- NULL
    local_comp_cache_error <- NULL
    if (identical(stat, "local_comp_enrichment")) {
      cache_obj <- tryCatch(
        .precompute_local_comp_cache(meta$ppp, radii_um),
        error = function(e) e
      )
      if (inherits(cache_obj, "error")) {
        local_comp_cache_error <- conditionMessage(cache_obj)
      } else {
        local_comp_cache <- cache_obj
      }
      if (!is.null(local_comp_cache_error)) {
        warning(
          sprintf(
            "panoramic_spatialstats: sample '%s' local composition precompute failed: %s. Returning NA for all pairs in this sample.",
            sid, local_comp_cache_error
          ),
          call. = FALSE
        )
      }
    }
    
    cur <- lapply(seq_len(nrow(pairs)), function(i) {
      pr <- pairs[i, ]
      if (!is.null(local_comp_cache_error)) {
        df <- data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_)
        df$ct1 <- pr$ct1
        df$ct2 <- pr$ct2
        df$key <- paste(df$ct1, df$ct2, df$radius_um, sep = "|")
        return(df)
      }
      df <- tryCatch(
        .one_pair_one_sample(
          meta,
          pr$ct1, pr$ct2,
          stat       = stat,
          nsim       = nsim,
          correction = correction,
          radii_um   = radii_um,
          verbose    = verbose,
          boot       = boot,
          tile_size  = tile_size,
          nx         = nx,
          ny         = ny,
          local_comp_cache = local_comp_cache
        ),
        error = function(e) {
          warning(
            sprintf(
              "panoramic_spatialstats: sample '%s', pair '%s'/'%s' failed: %s",
              sid, pr$ct1, pr$ct2, conditionMessage(e)
            ),
            call. = FALSE
          )
          data.frame(radius_um = radii_um, yi = NA_real_, vi = NA_real_)
        }
      )
      names(df)[names(df) == "r"] <- "radius_um"
      df$ct1 <- pr$ct1
      df$ct2 <- pr$ct2
      df$key <- paste(df$ct1, df$ct2, df$radius_um, sep = "|")
      df
    })
    
    do.call(rbind, cur)
  }
  
  # Prefer BiocParallel RNG handling for parallel backends
  if (!is.null(seed) && BiocParallel::bpworkers(BPPARAM) > 1L) {
    BiocParallel::bpRNGseed(BPPARAM) <- seed
    per_sample <- BiocParallel::bplapply(samples, compute_one_sample, BPPARAM = BPPARAM)
  } else if (!is.null(seed)) {
    per_sample <- withr::with_seed(seed, {
      BiocParallel::bplapply(samples, compute_one_sample, BPPARAM = BPPARAM)
    })
  } else {
    per_sample <- BiocParallel::bplapply(samples, compute_one_sample, BPPARAM = BPPARAM)
  }
  
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
    group  = vapply(prep, function(s) S4Vectors::metadata(s)$panoramic$group_id, character(1))
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = yi, vi = vi),
    rowData = rowData, colData = colData,
    metadata = list(panoramic = list(
      radii_um = radii_um,
      pairs = pairs,
      stat = stat,
      boot = boot,
      tile_size = tile_size,
      nx = nx,
      ny = ny
    ))
  )
  se_meta <- S4Vectors::metadata(se)
  se_meta$panoramic$tables <- list(
    spatialstats = .panoramic_extract_spatialstats_table(se, drop_na = FALSE)
  )
  S4Vectors::metadata(se) <- se_meta
  se
}


#' PANORAMIC end-to-end wrapper
#' 
#' One-line convenience wrapper that delegates to \code{panoramic_analyze()}.
#' It runs preparation, spatial statistics, and meta-analysis, and returns
#' both object-level and table-level outputs for downstream work.
#' 
#' @inheritParams panoramic_prepare
#' @inheritParams panoramic_analyze
#' 
#' @return A list with \code{prep}, \code{stats}, \code{pooled}, and
#'  \code{tables} as returned by \code{panoramic_analyze()}. 
#'
#' @examples
#' library(BiocParallel)
#'
#' sample_ids <- paste0("sample", seq_len(4L))
#' group_vec <- c("group1", "group1", "group2", "group2")
#' scenario_vec <- c("random", "random", "colocalized", "colocalized")
#'
#' spe_list <- stats::setNames(
#'   lapply(seq_along(sample_ids), function(i) {
#'     panoramic_simulate_spe(
#'       n_cells = 80,
#'       sample_id = sample_ids[i],
#'       cell_types = c("A", "B"),
#'       scenario = scenario_vec[i],
#'       seed = i
#'     )
#'   }),
#'   sample_ids
#' )
#'
#' design <- data.frame(
#'   sample = sample_ids,
#'   group  = group_vec,
#'   stringsAsFactors = FALSE
#' )
#'
#' out <- panoramic(
#'   spe_list,
#'   design      = design,
#'   cell_type   = "cell_type",
#'   radii_um    = c(10, 20),
#'   nsim        = 5,
#'   min_cells   = 2,
#'   concavity   = 50,
#'   window      = "rect",
#'   seed        = 1,
#'   BPPARAM     = BiocParallel::SerialParam()
#' )
#'
#' names(out)
#' @export
panoramic <- function(
    spe_list, design, cell_type = "cell_type",
    pairs = "auto", radii_um, stat = "local_comp_enrichment", nsim = 100,
    correction = NULL, min_cells = 5L, concavity = 50,
    window = c("concave", "convex", "rect"),
    group_col = "group",
    group_tau2 = c("none", "separate"),
    patient_col = NULL,
    sample_col = NULL,
    tau_structure = c("patient", "patient_sample"),
    method_mv = "REML",
    vi_floor = "group_median",
    seed = 123,
    boot = c("approx", "block"), tile_size = NULL, nx = NULL, ny = NULL,
    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
) {
  panoramic_analyze(
    spe_list = spe_list,
    design = design,
    cell_type = cell_type,
    pairs = pairs,
    radii_um = radii_um,
    stat = stat,
    nsim = nsim,
    correction = correction,
    min_cells = min_cells,
    concavity = concavity,
    window = window,
    group_col = group_col,
    group_tau2 = group_tau2,
    patient_col = patient_col,
    sample_col = sample_col,
    tau_structure = tau_structure,
    method_mv = method_mv,
    vi_floor = vi_floor,
    seed = seed,
    boot = boot,
    tile_size = tile_size,
    nx = nx,
    ny = ny,
    BPPARAM = BPPARAM,
    verbose = verbose
  )
}
