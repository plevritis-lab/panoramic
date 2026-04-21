#' @title Multilevel random-effects meta-analysis of PANORAMIC features
#' @description
#' Perform feature-wise multilevel random-effects meta-analysis on PANORAMIC
#' spatial statistics using \code{metafor::rma.mv()}, accounting for clustering
#' of samples within patients. Optionally estimate group-specific pooled means
#' via fixed-effect moderators.
#'
#' @param se A \code{SummarizedExperiment} produced by
#'  \code{panoramic_spatialstats()} or \code{panoramic()}, containing
#'  assays \code{"yi"} (effect estimates) and \code{"vi"} (within-sample variances).
#' @param patient_col Character scalar giving the name of a column in
#'  \code{colData(se)} that identifies patients.
#' @param group_col Optional character scalar giving the name of a column in
#'  \code{colData(se)} that defines groups (e.g. treatment). If provided,
#'  group-specific pooled means are estimated using fixed-effect moderators
#'  with no intercept. For two-group contrasts, the reported
#'  \code{beta_diff} is computed as \code{group2 - group1}, where group order
#'  is deterministic: factor levels if \code{group_col} is a factor, otherwise
#'  lexicographic order of observed group labels.
#' @param sample_col Optional character scalar giving the sample identifier
#'  column in \code{colData(se)}. If \code{NULL}, \code{colnames(se)} are used.
#' @param tau_structure Character. Random-effects structure to use. One of
#'  \code{"patient"} (random intercept per patient) or \code{"patient_sample"}
#'  (nested random intercepts per patient and sample).
#' @param method Character. Estimator for variance components passed to
#'  \code{metafor::rma.mv()}. Default \code{"REML"}.
#' @param group_tau2 Character or logical. If \code{"separate"} (or \code{TRUE})
#'  and \code{group_col} is provided, fit separate multilevel models within
#'  each group to estimate group-specific heterogeneity (\code{tau2}) values.
#'  Use \code{"none"} (or \code{FALSE}) to skip group-specific \code{tau2}.
#' @param warn_sigma2 Logical. If \code{TRUE}, warn when a large fraction of
#'  features have near-zero variance components (potential boundary fits).
#' @param sigma2_tol Numeric absolute tolerance used to flag near-zero
#'  variance components.
#' @param sigma2_rel Numeric relative tolerance (fraction of the per-feature
#'  median \code{vi}) used to flag near-zero variance components.
#' @param vi_floor Optional character. If \code{"median"}, replace non-positive
#'  \code{vi} values for a feature with the median of its positive \code{vi}
#'  values before fitting. If \code{"group_median"}, compute the median
#'  positive \code{vi} within each group and apply per-group flooring; if any
#'  group has no positive \code{vi} for a feature, that feature is skipped.
#'  If \code{NULL} or \code{NA} (or \code{"none"}), no flooring is applied and
#'  samples with \code{vi <= 0} are dropped.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallelization
#'  across features (rows). Defaults to \code{BiocParallel::SerialParam()}.
#' @param control Optional named list of optimizer control arguments passed to
#'  \code{metafor::rma.mv()} (for example \code{iter.max}, \code{eval.max},
#'  \code{rel.tol}).
#' @param sparse Logical passed to \code{metafor::rma.mv(sparse = ...)}.
#'
#' @return The input \code{SummarizedExperiment}, with additional columns
#'  appended to \code{rowData(se)}:
#'  \itemize{
#'    \item If \code{group_col} is \code{NULL}: \code{mu_hat}, \code{se_mu},
#'      \code{p_mu}, \code{k}, and variance components \code{tau2_patient}
#'      (and \code{tau2_sample} if \code{tau_structure = "patient_sample"}).
#'    \item If \code{group_col} is provided: group-prefixed \code{*_mu_hat},
#'      \code{*_se_mu}, \code{*_p_mu}, \code{*_k} plus variance components
#'      as above (shared across groups). For two groups, contrast columns
#'      \code{beta_diff}, \code{se_diff}, \code{z_diff}, \code{p_diff},
#'      and \code{fdr_diff} are also added.
#'    \item If \code{group_tau2 = TRUE}, group-prefixed heterogeneity
#'      estimates \code{*_tau2_patient} (and \code{*_tau2_sample} when
#'      applicable) are added.
#' }
#' Metadata are stored in \code{metadata(se)$panoramic$meta_mv}.
#'
#' @examples
#' yi <- matrix(c(0.20, 0.15, 0.10, 0.05), nrow = 1)
#' vi <- matrix(c(0.04, 0.05, 0.06, 0.05), nrow = 1)
#' colnames(yi) <- colnames(vi) <- paste0("s", seq_len(4L))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(yi = yi, vi = vi),
#'   rowData = S4Vectors::DataFrame(ct1 = "A", ct2 = "B", radius_um = 10),
#'   colData = S4Vectors::DataFrame(
#'     sample = paste0("s", seq_len(4L)),
#'     patient = paste0("p", seq_len(4L)),
#'     group = c("control", "control", "case", "case")
#'   )
#' )
#'
#' se_mv <- panoramic_meta_mv(
#'   se = se,
#'   patient_col = "patient",
#'   group_col = "group",
#'   sample_col = "sample",
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#' SummarizedExperiment::rowData(se_mv)$beta_diff
#'
#' @export
panoramic_meta_mv <- function(
    se,
    patient_col,
    group_col = NULL,
    sample_col = "sample",
    tau_structure = c("patient", "patient_sample"),
    method = "REML",
    group_tau2 = c("none", "separate"),
    warn_sigma2 = TRUE,
    vi_floor = NULL,
    sigma2_tol = 1e-6,
    sigma2_rel = 1e-4,
    control = NULL,
    sparse = FALSE,
    BPPARAM = BiocParallel::SerialParam()
) {
  tau_structure <- match.arg(tau_structure)
  if (is.null(vi_floor) || (length(vi_floor) == 1L && is.na(vi_floor))) {
    vi_floor <- "none"
  } else {
    vi_floor <- match.arg(vi_floor, c("median", "group_median", "none"))
  }
  stopifnot(is.numeric(sigma2_tol), length(sigma2_tol) == 1L, sigma2_tol >= 0,
            is.numeric(sigma2_rel), length(sigma2_rel) == 1L, sigma2_rel >= 0)
  if (!is.null(control) && !is.list(control)) {
    stop("`control` must be NULL or a named list.")
  }
  sparse <- isTRUE(sparse)
  if (is.logical(group_tau2)) {
    group_tau2 <- if (isTRUE(group_tau2)) "separate" else "none"
  } else {
    group_tau2 <- match.arg(group_tau2)
  }
  group_tau2_flag <- group_tau2 == "separate"
  stopifnot("yi" %in% names(SummarizedExperiment::assays(se)),
            "vi" %in% names(SummarizedExperiment::assays(se)))

  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se))

  if (!patient_col %in% colnames(cd)) {
    stop("`patient_col` not found in colData(se).")
  }
  if (!is.null(group_col) && !group_col %in% colnames(cd)) {
    stop("`group_col` not found in colData(se).")
  }

  if (is.null(sample_col)) {
    sample_id <- colnames(yi)
  } else if (sample_col %in% colnames(cd)) {
    sample_id <- as.character(cd[[sample_col]])
  } else {
    stop("`sample_col` not found in colData(se).")
  }

  patient_id <- as.character(cd[[patient_col]])
  group_id <- if (is.null(group_col)) NULL else as.character(cd[[group_col]])
  all_groups <- if (is.null(group_id)) {
    NULL
  } else {
    group_raw <- cd[[group_col]]
    if (is.factor(group_raw)) {
      lev <- levels(group_raw)
      lev[lev %in% group_id]
    } else {
      sort(unique(group_id[!is.na(group_id)]), method = "radix")
    }
  }
  if (!is.null(group_col) && length(all_groups) == 0L) {
    stop("`group_col` has no non-missing groups in colData(se).")
  }
  group_names <- NULL
  if (!is.null(group_col)) {
    group_names <- make.names(all_groups)
    if (anyDuplicated(group_names)) {
      dup_groups <- unique(group_names[duplicated(group_names)])
      colliding <- vapply(dup_groups, function(dn) {
        paste(all_groups[group_names == dn], collapse = " vs ")
      }, character(1))
      stop(
        "`group_col` labels are ambiguous after sanitization (make.names). ",
        "Please rename groups to unique safe labels. Collisions: ",
        paste(colliding, collapse = "; "),
        "."
      )
    }
  }

  random_formula <- if (tau_structure == "patient") {
    stats::as.formula("~ 1 | patient_id")
  } else {
    stats::as.formula("~ 1 | patient_id/sample_id")
  }

  fit_one <- function(i) {
    y <- yi[i, ]
    v <- vi[i, ]
    if (vi_floor == "median") {
      v_pos <- v[is.finite(v) & v > 0]
      if (length(v_pos) > 0) {
        v_floor <- stats::median(v_pos)
        v[is.finite(v) & v <= 0] <- v_floor
      }
    } else if (vi_floor == "group_median") {
      if (is.null(group_id)) {
        return(c(mu_hat = NA_real_, se_mu = NA_real_, p_mu = NA_real_, k = 0,
                 tau2_patient = NA_real_,
                 tau2_sample = if (tau_structure == "patient_sample") NA_real_ else NULL))
      }
      ok_groups <- TRUE
      for (g in all_groups) {
        idx_g <- group_id == g
        v_pos_g <- v[idx_g & is.finite(v) & v > 0]
        if (length(v_pos_g) == 0L) {
          ok_groups <- FALSE
          break
        }
        v_floor_g <- stats::median(v_pos_g)
        v[idx_g & is.finite(v) & v <= 0] <- v_floor_g
      }
      if (!ok_groups) {
        out <- list()
        for (g in group_names) {
          prefix <- paste0(g, "_")
          out[[paste0(prefix, "mu_hat")]] <- NA_real_
          out[[paste0(prefix, "se_mu")]] <- NA_real_
          out[[paste0(prefix, "p_mu")]] <- NA_real_
          out[[paste0(prefix, "k")]] <- 0
        }
        out$tau2_patient <- NA_real_
        if (tau_structure == "patient_sample") out$tau2_sample <- NA_real_
        return(unlist(out))
      }
    }
    keep <- is.finite(y) & is.finite(v) & v > 0

    if (sum(keep) < 2) {
      if (is.null(group_id)) {
        return(c(mu_hat = NA_real_, se_mu = NA_real_, p_mu = NA_real_, k = sum(keep),
                 tau2_patient = NA_real_,
                 tau2_sample = if (tau_structure == "patient_sample") NA_real_ else NULL))
      }
      out <- list()
      group_id_keep <- make.names(group_id[keep])
      for (g in group_names) {
        prefix <- paste0(g, "_")
        out[[paste0(prefix, "mu_hat")]] <- NA_real_
        out[[paste0(prefix, "se_mu")]] <- NA_real_
        out[[paste0(prefix, "p_mu")]] <- NA_real_
        out[[paste0(prefix, "k")]] <- sum(group_id_keep == g)
      }
      out$tau2_patient <- NA_real_
      if (tau_structure == "patient_sample") out$tau2_sample <- NA_real_
      return(unlist(out))
    }

    df <- data.frame(
      yi = y[keep],
      vi = v[keep],
      patient_id = patient_id[keep],
      sample_id = sample_id[keep],
      stringsAsFactors = FALSE
    )
    if (!is.null(group_id)) df$group_id <- group_id[keep]

    if (is.null(group_id)) {
      fit <- try(metafor::rma.mv(
        yi = y[keep], V = v[keep], random = random_formula,
        data = df, method = method, test = "z",
        control = control, sparse = sparse
      ), silent = TRUE)

      if (inherits(fit, "try-error")) {
        return(c(mu_hat = NA_real_, se_mu = NA_real_, p_mu = NA_real_, k = sum(keep),
                 tau2_patient = NA_real_,
                 tau2_sample = if (tau_structure == "patient_sample") NA_real_ else NULL))
      }

      sigma2 <- unname(fit$sigma2)
      out <- c(
        mu_hat = unname(fit$b[1, 1]),
        se_mu = unname(fit$se[1]),
        p_mu = unname(fit$pval[1]),
        k = sum(keep),
        tau2_patient = if (length(sigma2) >= 1) sigma2[1] else NA_real_
      )
      if (tau_structure == "patient_sample") {
        out <- c(out, tau2_sample = if (length(sigma2) >= 2) sigma2[2] else NA_real_)
      }
      return(out)
    }

    df$group_id <- factor(make.names(df$group_id), levels = group_names)
    fit <- try(metafor::rma.mv(
      yi = y[keep], V = v[keep], mods = ~ 0 + group_id, random = random_formula,
      data = df, method = method, test = "z",
      control = control, sparse = sparse
    ), silent = TRUE)

    if (inherits(fit, "try-error")) {
      out <- list()
      for (g in levels(df$group_id)) {
        prefix <- paste0(g, "_")
        out[[paste0(prefix, "mu_hat")]] <- NA_real_
        out[[paste0(prefix, "se_mu")]] <- NA_real_
        out[[paste0(prefix, "p_mu")]] <- NA_real_
        out[[paste0(prefix, "k")]] <- sum(df$group_id == g)
      }
      out$tau2_patient <- NA_real_
      if (tau_structure == "patient_sample") out$tau2_sample <- NA_real_
      return(unlist(out))
    }

    sigma2 <- unname(fit$sigma2)
    out <- list()
    coefs <- unname(fit$b)
    ses <- unname(fit$se)
    pvals <- unname(fit$pval)

    for (g in group_names) {
      prefix <- paste0(g, "_")
      coef_name <- paste0("group_id", g)
      idx <- match(coef_name, rownames(fit$b))
      if (is.na(idx)) {
        out[[paste0(prefix, "mu_hat")]] <- NA_real_
        out[[paste0(prefix, "se_mu")]] <- NA_real_
        out[[paste0(prefix, "p_mu")]] <- NA_real_
      } else {
        out[[paste0(prefix, "mu_hat")]] <- coefs[idx]
        out[[paste0(prefix, "se_mu")]] <- ses[idx]
        out[[paste0(prefix, "p_mu")]] <- pvals[idx]
      }
      out[[paste0(prefix, "k")]] <- sum(df$group_id == g)
    }
    
    # Optional: contrast for exactly two groups
    if (length(all_groups) == 2L) {
      g1 <- group_names[1]
      g2 <- group_names[2]
      idx1 <- match(paste0("group_id", g1), rownames(fit$b))
      idx2 <- match(paste0("group_id", g2), rownames(fit$b))
      if (!is.na(idx1) && !is.na(idx2)) {
        vb <- fit$vb
        var_diff <- vb[idx1, idx1] + vb[idx2, idx2] - 2 * vb[idx1, idx2]
        se_diff <- if (is.finite(var_diff) && var_diff >= 0) sqrt(var_diff) else NA_real_
        beta_diff <- coefs[idx2] - coefs[idx1]
        z_diff <- beta_diff / se_diff
        p_diff <- 2 * stats::pnorm(-abs(z_diff))
        out[["beta_diff"]] <- beta_diff
        out[["se_diff"]] <- se_diff
        out[["z_diff"]] <- z_diff
        out[["p_diff"]] <- p_diff
      } else {
        out[["beta_diff"]] <- NA_real_
        out[["se_diff"]] <- NA_real_
        out[["z_diff"]] <- NA_real_
        out[["p_diff"]] <- NA_real_
      }
    }
    out$tau2_patient <- if (length(sigma2) >= 1) sigma2[1] else NA_real_
    if (tau_structure == "patient_sample") {
      out$tau2_sample <- if (length(sigma2) >= 2) sigma2[2] else NA_real_
    }
    
    # Optional: group-specific heterogeneity via separate fits
    if (isTRUE(group_tau2_flag)) {
      for (g in group_names) {
        idx_g <- df$group_id == g
        out[[paste0(g, "_tau2_patient")]] <- NA_real_
        if (tau_structure == "patient_sample") {
          out[[paste0(g, "_tau2_sample")]] <- NA_real_
        }
        if (sum(idx_g) >= 2) {
          fit_g <- try(metafor::rma.mv(
            yi = df$yi[idx_g], V = df$vi[idx_g], random = random_formula,
            data = df[idx_g, , drop = FALSE], method = method, test = "z",
            control = control, sparse = sparse
          ), silent = TRUE)
          if (!inherits(fit_g, "try-error")) {
            s2 <- unname(fit_g$sigma2)
            out[[paste0(g, "_tau2_patient")]] <- if (length(s2) >= 1) s2[1] else NA_real_
            if (tau_structure == "patient_sample") {
              out[[paste0(g, "_tau2_sample")]] <- if (length(s2) >= 2) s2[2] else NA_real_
            }
          }
        }
      }
    }
    unlist(out)
  }

  M <- BiocParallel::bplapply(seq_len(nrow(yi)), fit_one, BPPARAM = BPPARAM)
  M <- do.call(rbind, M)

  if ("p_diff" %in% colnames(M)) {
    M <- cbind(M, fdr_diff = p.adjust(M[, "p_diff"], method = "BH"))
  }
  rd_new <- S4Vectors::DataFrame(M)
  rd_old <- S4Vectors::as.data.frame(SummarizedExperiment::rowData(se))
  overlap <- intersect(colnames(rd_old), colnames(rd_new))
  if (length(overlap) > 0L) {
    rd_old[overlap] <- NULL
  }
  SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(cbind(rd_old, as.data.frame(rd_new)))

  if (isTRUE(warn_sigma2)) {
    tau_cols <- intersect(c("tau2_patient", "tau2_sample"), colnames(M))
    if (length(tau_cols) > 0) {
      vi_med <- apply(vi, 1, function(x) {
        x <- x[is.finite(x) & x > 0]
        if (length(x) == 0L) NA_real_ else stats::median(x)
      })
      thresh <- pmax(sigma2_tol, sigma2_rel * vi_med)
      flag <- rep(FALSE, nrow(M))
      for (tc in tau_cols) {
        x <- M[, tc]
        flag <- flag | (is.finite(x) & is.finite(thresh) & x <= thresh)
      }
      frac <- mean(flag, na.rm = TRUE)
      if (is.finite(frac) && frac >= 0.2) {
        warning(
          sprintf(
            "panoramic_meta_mv: %.1f%% of features have near-zero variance components (<= max(%g, %g * median(vi))).",
            100 * frac, sigma2_tol, sigma2_rel
          ),
          call. = FALSE
        )
      }
    }
  }

  m <- S4Vectors::metadata(se)
  if (is.null(m)) m <- list()
  if (is.null(m$panoramic)) m$panoramic <- list()
  m$panoramic$meta_mv <- list(
    method = method,
    patient_col = patient_col,
    sample_col = sample_col,
    group_col = group_col,
    tau_structure = tau_structure,
    group_tau2 = group_tau2,
    vi_floor = vi_floor,
    sparse = sparse,
    control = control
  )
  if (!is.null(group_col) && length(all_groups) == 2L) {
    m$panoramic$contrast <- list(
      control = all_groups[1],
      case = all_groups[2],
      method = "multilevel: single-stage contrast (beta_diff = case - control)"
    )
  }
  m$panoramic$tables <- m$panoramic$tables %||% list()
  m$panoramic$tables$spatialstats <- .panoramic_extract_spatialstats_table(se, drop_na = FALSE)
  m$panoramic$tables$meta <- .panoramic_extract_meta_table(se, drop_na = FALSE)
  m$panoramic$tables$contrast <- tryCatch(
    panoramic_extract_contrast(se),
    error = function(e) NULL
  )
  S4Vectors::metadata(se) <- m

  se
}
