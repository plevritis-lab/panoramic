#' Two-group test (case - control) with standardized outputs
#'
#' Adds to rowData:
#'   beta_diff, se_diff, z_diff, p_diff, fdr_diff, k, tau2,
#'   mu_control, se_control, mu_case, se_case
#' Stores contrast info in S4Vectors::metadata(se)$panoramic$contrast
#'
#' @param se SummarizedExperiment from pano_spatialstats() or pano_meta()
#' @param group_col colData(se) column with group labels (default "group")
#' @param case case-group label (optional; inferred if exactly 2 groups)
#' @param control control-group label (optional; inferred if exactly 2 groups)
#' @param tau2 tau^2 estimator for metafor::rma.uni ("SJ","REML","DL")
#' @param BPPARAM BiocParallel param (default SerialParam())
#' @return SummarizedExperiment with rowData augmented and metadata contrast set
#' @export
pano_test_groups <- function(se, group_col = "group",
                             case = NULL, control = NULL,
                             tau2 = "SJ",
                             BPPARAM = BiocParallel::SerialParam()) {
  # -- validations
  if (!inherits(se, "SummarizedExperiment"))
    stop("`se` must be a SummarizedExperiment.")
  if (!group_col %in% colnames(SummarizedExperiment::colData(se)))
    stop("`group_col` not found in colData(se).")
  if (!all(c("yi","vi") %in% names(SummarizedExperiment::assays(se))))
    stop("Assays 'yi' and 'vi' are required. Run pano_spatialstats() first.")
  
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  g  <- as.character(cd[[group_col]])
  lev <- unique(g)
  
  # infer two groups if needed
  if (is.null(control) || is.null(case)) {
    if (length(lev) != 2L)
      stop("Provide `case` and `control` explicitly, or ensure exactly two groups.")
    if (is.null(control)) control <- lev[1]
    if (is.null(case))    case    <- lev[2]
  }
  grp <- factor(g, levels = c(control, case))
  if (any(is.na(grp)))
    stop("Found samples not in {control, case}: ", paste(setdiff(lev, c(control, case)), collapse = ", "))
  
  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  
  # per-feature worker: returns a fixed-length named numeric vector
  fit_one <- function(i) {
    y <- yi[i, ]; v <- vi[i, ]
    keep <- is.finite(y) & is.finite(v) & v > 0 & !is.na(grp)
    k <- sum(keep)
    
    out <- c(beta_diff = NA_real_, se_diff = NA_real_, z_diff = NA_real_, p_diff = NA_real_,
             k = k, tau2 = NA_real_,
             mu_control = NA_real_, se_control = NA_real_,
             mu_case = NA_real_, se_case = NA_real_)
    
    if (k < 2L) return(out)
    
    gk <- droplevels(grp[keep])
    
    # If only one group contributes for this feature, fit intercept-only and fill that group's mean
    if (nlevels(gk) < 2L) {
      fit0 <- try(metafor::rma.uni(yi = y[keep], vi = v[keep], method = tau2, test = "z"),
                  silent = TRUE)
      if (inherits(fit0, "try-error")) return(out)
      mu  <- unname(fit0$b[1,1]); se0 <- unname(fit0$se[1]); t2 <- unname(fit0$tau2)
      if (levels(gk)[1] == control) {
        out["mu_control"] <- mu; out["se_control"] <- se0
      } else {
        out["mu_case"]    <- mu; out["se_case"]    <- se0
      }
      out["tau2"] <- t2
      return(out)
    }
    
    # Both groups present → meta-regression with intercept = control
    M <- stats::model.matrix(~ gk)  # columns: (Intercept), gkcase
    fit <- try(metafor::rma.uni(yi = y[keep], vi = v[keep], mods = M, method = tau2, test = "z"),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(out)
    
    bnm <- rownames(fit$b); b <- as.numeric(fit$b)
    V   <- try(stats::vcov(fit), silent = TRUE)
    if (inherits(V, "try-error")) V <- matrix(NA_real_, length(b), length(b), dimnames = list(bnm, bnm))
    
    idx_int  <- if (any(bnm %in% c("(Intercept)","intrcpt"))) which(bnm %in% c("(Intercept)","intrcpt"))[1] else 1L
    idx_diff <- setdiff(seq_along(b), idx_int)[1]
    beta   <- if (length(idx_diff)) b[idx_diff] else NA_real_
    se_b   <- if (length(idx_diff)) sqrt(V[idx_diff, idx_diff]) else NA_real_
    z_b    <- beta / se_b
    p_b    <- 2 * stats::pnorm(-abs(z_b))
    
    mu_ctl <- b[idx_int]
    se_ctl <- sqrt(V[idx_int, idx_int])
    cov_ic <- if (length(idx_diff)) V[idx_int, idx_diff] else NA_real_
    mu_cas <- mu_ctl + beta
    se_cas <- sqrt(se_ctl^2 + se_b^2 + 2 * cov_ic)
    
    out[] <- c(beta, se_b, z_b, p_b, k, unname(fit$tau2),
               mu_ctl, se_ctl, mu_cas, se_cas)
    out
  }
  
  res <- BiocParallel::bplapply(seq_len(nrow(yi)), fit_one, BPPARAM = BPPARAM)
  M   <- do.call(rbind, res)
  
  rd_add <- S4Vectors::DataFrame(M)
  rd_add$fdr_diff <- stats::p.adjust(rd_add$p_diff, method = "BH")
  
  # append to rowData
  SummarizedExperiment::rowData(se) <- cbind(SummarizedExperiment::rowData(se), rd_add)
  
  # write standardized contrast info into S4 metadata
  m <- S4Vectors::metadata(se)
  if (is.null(m)) m <- list()
  if (is.null(m$panoramic)) m$panoramic <- list()
  m$panoramic$contrast <- list(group_col = group_col, control = control, case = case)
  S4Vectors::metadata(se) <- m
  
  se
}


