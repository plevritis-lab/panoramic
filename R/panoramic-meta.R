#' Random-effects pooling (feature-wise)
#'
#' Writes pooled estimates into rowData: mu_hat, se_mu, tau2, Q, I2, k, p_mu.
#' @param se SummarizedExperiment from \code{pano_spatialstats()}
#' @param tau2 estimator for between-sample variance ("SJ","REML","DL")
#' @param moderators optional data.frame with per-sample moderators (columns align to colData rows)
#' @param BPPARAM BiocParallel param
#' @export
panoramic_meta <- function(se, tau2 = c("SJ","REML","DL"), moderators = NULL,
                      BPPARAM = BiocParallel::SerialParam()) {
  tau2 <- match.arg(tau2)
  stopifnot("yi" %in% names(SummarizedExperiment::assays(se)),
            "vi" %in% names(SummarizedExperiment::assays(se)))
  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se))
  
  # Build moderator model matrix if provided (e.g., group)
  if (!is.null(moderators)) {
    stopifnot(nrow(moderators) == nrow(cd))
    X <- stats::model.matrix(~ . , data = moderators)
  } else {
    X <- matrix(1, nrow = nrow(cd), ncol = 1) # intercept-only pooling
  }
  
  fit_one <- function(i) {
    y  <- yi[i, ]
    v  <- vi[i, ]
    keep <- is.finite(y) & is.finite(v) & v > 0
    if (sum(keep) < 2) {
      return(c(mu_hat = NA_real_, se_mu = NA_real_, tau2 = NA_real_, Q = NA_real_, I2 = NA_real_, k = sum(keep), p_mu = NA_real_))
    }
    dat <- data.frame(y = y[keep], v = v[keep])
    Xk  <- X[keep, , drop = FALSE]
    # metafor handles moderators; intercept-only if ncol(Xk) == 1
    ctrl <- metafor::rma.uni(yi = dat$y, vi = dat$v, mods = Xk,
                             method = tau2, test = "z")
    c(mu_hat = unname(ctrl$b[1,1]),
      se_mu  = unname(ctrl$se[1]),
      tau2   = unname(ctrl$tau2),
      Q      = unname(ctrl$QE),
      I2     = unname(ctrl$I2),
      k      = sum(keep),
      p_mu   = unname(ctrl$pval[1]))
  }
  
  M <- BiocParallel::bplapply(seq_len(nrow(yi)), fit_one, BPPARAM = BPPARAM)
  M <- do.call(rbind, M)
  SummarizedExperiment::rowData(se)[, c("mu_hat","se_mu","tau2","Q","I2","k","p_mu")] <- S4Vectors::DataFrame(M)
  se
}
