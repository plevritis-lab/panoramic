#' Compare group-specific PANORAMIC meta-analytic estimates
#'
#' Compare pre-computed group-specific random-effects meta-analysis estimates
#' from \code{panoramic_meta()} to test for differential colocalization
#' between two groups. This is a two-stage approach: group-specific pooled
#' effects are treated as independent estimates and compared via a z-test.
#'
#' @param se A \code{SummarizedExperiment} returned by
#'   \code{panoramic_meta()} with \code{group_col} specified, so that
#'   \code{rowData(se)} contains group-specific columns such as
#'   \code{GroupA_mu_hat}, \code{GroupA_se_mu}, \code{GroupB_mu_hat},
#'   \code{GroupB_se_mu}, and corresponding \code{*_tau2}.
#' @param group1 Character scalar giving the label (prefix) for the
#'   control/reference group, exactly as it appears in the column prefixes
#'   (after \code{make.names()}), e.g. \code{"X1"}.
#' @param group2 Character scalar giving the label (prefix) for the
#'   case/comparison group, matching the column prefixes in \code{rowData(se)}.
#'
#' @return The input \code{SummarizedExperiment} with additional
#'   comparison results appended to \code{rowData(se)}:
#'   \itemize{
#'     \item \code{beta_diff}: difference in pooled effects (\code{mu_2 - mu_1})
#'     \item \code{se_diff}: standard error of the difference
#'     \item \code{z_diff}: z-statistic for the difference
#'     \item \code{p_diff}: two-sided p-value
#'     \item \code{fdr_diff}: Benjamini–Hochberg adjusted p-value
#'     \item \code{mu_control}, \code{se_control}, \code{tau2_control}:
#'       pooled effect, standard error, and heterogeneity for \code{group1}
#'     \item \code{mu_case}, \code{se_case}, \code{tau2_case}:
#'       same quantities for \code{group2}
#'   }
#'   Comparison settings are also stored in \code{metadata(se)$panoramic$comparison}.
#'
#' @details
#' This function assumes that the group-specific pooled estimates for
#' \code{group1} and \code{group2} are approximately independent and
#' uses the usual formula for the standard error of a difference of
#' independent estimates. The method is intended as a simple two-stage
#' differential test leveraging the output of \code{panoramic_meta()}.
#'
#' @examples
#' \dontrun{
#' se_meta <- panoramic_meta(se, group_col = "group")
#' se_cmp  <- pano_compare_groups(se_meta, group1 = "Control", group2 = "Treatment")
#' }
#'
#' @export
panoramic_compare_groups <- function(se, group1 = "X1", group2 = "X2") {
  
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  
  # Check for required columns
  required_cols <- c(paste0(group1, "_mu_hat"), paste0(group1, "_se_mu"),
                     paste0(group2, "_mu_hat"), paste0(group2, "_se_mu"),
                     paste0(group1, "_tau2"), paste0(group2, "_tau2"))
  
  if (!all(required_cols %in% colnames(rd))) {
    stop("Required columns not found. Run panoramic_meta() with group_col first.\n",
         "Expected: ", paste(required_cols, collapse = ", "))
  }
  
  # Extract group-specific estimates
  mu_1 <- rd[[paste0(group1, "_mu_hat")]]
  se_1 <- rd[[paste0(group1, "_se_mu")]]
  tau2_1 <- rd[[paste0(group1, "_tau2")]]
  
  mu_2 <- rd[[paste0(group2, "_mu_hat")]]
  se_2 <- rd[[paste0(group2, "_se_mu")]]
  tau2_2 <- rd[[paste0(group2, "_tau2")]]
  
  # Compute difference
  beta_diff <- mu_2 - mu_1
  
  # Standard error of difference (independent estimates)
  se_diff <- sqrt(se_1^2 + se_2^2)
  
  # Z-test
  z_diff <- beta_diff / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))
  
  # FDR correction
  fdr_diff <- p.adjust(p_diff, method = "BH")
  
  # Create results data frame
  results <- S4Vectors::DataFrame(
    beta_diff = beta_diff,
    se_diff = se_diff,
    z_diff = z_diff,
    p_diff = p_diff,
    fdr_diff = fdr_diff,
    mu_control = mu_1,
    se_control = se_1,
    tau2_control = tau2_1,
    mu_case = mu_2,
    se_case = se_2,
    tau2_case = tau2_2
  )
  
  # Add to rowData
  SummarizedExperiment::rowData(se) <- cbind(SummarizedExperiment::rowData(se), results)
  
  # Store metadata
  m <- S4Vectors::metadata(se)
  if (is.null(m)) m <- list()
  if (is.null(m$panoramic)) m$panoramic <- list()
  m$panoramic$comparison <- list(
    control = group1,
    case = group2,
    method = "two-stage: independent group meta-analyses compared via z-test"
  )
  S4Vectors::metadata(se) <- m
  
  se
}


# #' Two-group test (case - control) with standardized outputs
# #'
# #' Adds to rowData:
# #'   beta_diff, se_diff, z_diff, p_diff, fdr_diff, k, tau2,
# #'   mu_control, se_control, mu_case, se_case
# #' Stores contrast info in S4Vectors::metadata(se)$panoramic$contrast
# #'
# #' @param se SummarizedExperiment from pano_spatialstats() or pano_meta()
# #' @param group_col colData(se) column with group labels (default "group")
# #' @param case case-group label (optional; inferred if exactly 2 groups)
# #' @param control control-group label (optional; inferred if exactly 2 groups)
# #' @param tau2 tau^2 estimator for metafor::rma.uni ("SJ","REML","DL")
# #' @param BPPARAM BiocParallel param (default SerialParam())
# #' @return SummarizedExperiment with rowData augmented and metadata contrast set
# pano_test_groups <- function(se, group_col = "group",
#                              case = NULL, control = NULL,
#                              tau2 = "SJ",
#                              BPPARAM = BiocParallel::SerialParam()) {
#   # -- validations
#   if (!inherits(se, "SummarizedExperiment"))
#     stop("`se` must be a SummarizedExperiment.")
#   if (!group_col %in% colnames(SummarizedExperiment::colData(se)))
#     stop("`group_col` not found in colData(se).")
#   if (!all(c("yi","vi") %in% names(SummarizedExperiment::assays(se))))
#     stop("Assays 'yi' and 'vi' are required. Run pano_spatialstats() first.")
  
#   cd <- as.data.frame(SummarizedExperiment::colData(se))
#   g  <- as.character(cd[[group_col]])
#   lev <- unique(g)
  
#   # infer two groups if needed
#   if (is.null(control) || is.null(case)) {
#     if (length(lev) != 2L)
#       stop("Provide `case` and `control` explicitly, or ensure exactly two groups.")
#     if (is.null(control)) control <- lev[1]
#     if (is.null(case))    case    <- lev[2]
#   }
#   grp <- factor(g, levels = c(control, case))
#   if (any(is.na(grp)))
#     stop("Found samples not in {control, case}: ", paste(setdiff(lev, c(control, case)), collapse = ", "))
  
#   yi <- SummarizedExperiment::assay(se, "yi")
#   vi <- SummarizedExperiment::assay(se, "vi")
  
#   # per-feature worker: returns a fixed-length named numeric vector
#   fit_one <- function(i) {
#     y <- yi[i, ]; v <- vi[i, ]
#     keep <- is.finite(y) & is.finite(v) & v > 0 & !is.na(grp)
#     k <- sum(keep)
    
#     out <- c(beta_diff = NA_real_, se_diff = NA_real_, z_diff = NA_real_, p_diff = NA_real_,
#              k = k, tau2 = NA_real_,
#              mu_control = NA_real_, se_control = NA_real_,
#              mu_case = NA_real_, se_case = NA_real_)
    
#     if (k < 2L) return(out)
    
#     gk <- droplevels(grp[keep])
    
#     # If only one group contributes for this feature, fit intercept-only and fill that group's mean
#     if (nlevels(gk) < 2L) {
#       fit0 <- try(metafor::rma.uni(yi = y[keep], vi = v[keep], method = tau2, test = "z"),
#                   silent = TRUE)
#       if (inherits(fit0, "try-error")) return(out)
#       mu  <- unname(fit0$b[1,1]); se0 <- unname(fit0$se[1]); t2 <- unname(fit0$tau2)
#       if (levels(gk)[1] == control) {
#         out["mu_control"] <- mu; out["se_control"] <- se0
#       } else {
#         out["mu_case"]    <- mu; out["se_case"]    <- se0
#       }
#       out["tau2"] <- t2
#       return(out)
#     }
    
#     # Both groups present → meta-regression with intercept = control
#     M <- stats::model.matrix(~ gk)  # columns: (Intercept), gkcase
#     fit <- try(metafor::rma.uni(yi = y[keep], vi = v[keep], mods = M, method = tau2, test = "z"),
#                silent = TRUE)
#     if (inherits(fit, "try-error")) return(out)
    
#     bnm <- rownames(fit$b); b <- as.numeric(fit$b)
#     V   <- try(stats::vcov(fit), silent = TRUE)
#     if (inherits(V, "try-error")) V <- matrix(NA_real_, length(b), length(b), dimnames = list(bnm, bnm))
    
#     idx_int  <- if (any(bnm %in% c("(Intercept)","intrcpt"))) which(bnm %in% c("(Intercept)","intrcpt"))[1] else 1L
#     idx_diff <- setdiff(seq_along(b), idx_int)[1]
#     beta   <- if (length(idx_diff)) b[idx_diff] else NA_real_
#     se_b   <- if (length(idx_diff)) sqrt(V[idx_diff, idx_diff]) else NA_real_
#     z_b    <- beta / se_b
#     p_b    <- 2 * stats::pnorm(-abs(z_b))
    
#     mu_ctl <- b[idx_int]
#     se_ctl <- sqrt(V[idx_int, idx_int])
#     cov_ic <- if (length(idx_diff)) V[idx_int, idx_diff] else NA_real_
#     mu_cas <- mu_ctl + beta
#     se_cas <- sqrt(se_ctl^2 + se_b^2 + 2 * cov_ic)
    
#     out[] <- c(beta, se_b, z_b, p_b, k, unname(fit$tau2),
#                mu_ctl, se_ctl, mu_cas, se_cas)
#     out
#   }
  
#   res <- BiocParallel::bplapply(seq_len(nrow(yi)), fit_one, BPPARAM = BPPARAM)
#   M   <- do.call(rbind, res)
  
#   rd_add <- S4Vectors::DataFrame(M)
#   rd_add$fdr_diff <- stats::p.adjust(rd_add$p_diff, method = "BH")
  
#   # append to rowData
#   SummarizedExperiment::rowData(se) <- cbind(SummarizedExperiment::rowData(se), rd_add)
  
#   # write standardized contrast info into S4 metadata
#   m <- S4Vectors::metadata(se)
#   if (is.null(m)) m <- list()
#   if (is.null(m$panoramic)) m$panoramic <- list()
#   m$panoramic$contrast <- list(group_col = group_col, control = control, case = case)
#   S4Vectors::metadata(se) <- m
  
#   se
# }


