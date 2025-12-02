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
#' library(SpatialExperiment)
#' library(S4Vectors)
#' library(BiocParallel)
#'
#' set.seed(1)
#'
#' make_spe <- function(n_cells, sample_id, group_label) {
#'   coords <- cbind(
#'     x = runif(n_cells, 0, 100),
#'     y = runif(n_cells, 0, 100)
#'   )
#'   ct <- sample(c("A", "B"), size = n_cells, replace = TRUE)
#'
#'   counts <- matrix(
#'     rpois(5 * n_cells, lambda = 5),
#'     nrow = 5,
#'     dimnames = list(
#'       paste0("gene", seq_len(5)),
#'       paste0(sample_id, "_cell", seq_len(n_cells))
#'     )
#'   )
#'
#'   SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'       cell_type = ct,
#'       sample_id = sample_id,
#'       group     = group_label
#'     ),
#'     spatialCoords = coords
#'   )
#' }
#'
#' # Two control and two case samples -------------------------------
#' spe1 <- make_spe(20, "sample1", "control")
#' spe2 <- make_spe(22, "sample2", "control")
#' spe3 <- make_spe(24, "sample3", "case")
#' spe4 <- make_spe(26, "sample4", "case")
#'
#' spe_list <- list(
#'   sample1 = spe1,
#'   sample2 = spe2,
#'   sample3 = spe3,
#'   sample4 = spe4
#' )
#'
#' design <- data.frame(
#'   sample = paste0("sample", 1:4),
#'   group  = c("control", "control", "case", "case"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # PANORAMIC spatial stats ---------------------------------------
#' se_stats <- panoramic(
#'   spe_list,
#'   design      = design,
#'   cell_type   = "cell_type",
#'   radii_um    = c(10, 20),
#'   stat        = "Lcross",
#'   nsim        = 5,
#'   correction  = "translate",
#'   min_cells   = 2,
#'   concavity   = 50,
#'   window      = "rect",
#'   seed        = 1,
#'   BPPARAM     = BiocParallel::SerialParam()
#' )
#'
#' # Group-specific meta-analysis and comparison -------------------
#' se_meta <- panoramic_meta(se_stats, tau2 = "SJ", group_col = "group")
#'
#' # Column prefixes are make.names("control") and make.names("case"),
#' # which are "control" and "case" here.
#' se_diff <- panoramic_compare_groups(
#'   se_meta,
#'   group1 = "control",
#'   group2 = "case"
#' )
#'
#' se_diff
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