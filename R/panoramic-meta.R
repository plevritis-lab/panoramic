#' Random-effects meta-analysis of PANORAMIC features 
#' Perform feature-wise random-effects meta-analysis on PANORAMIC spatial 
#' statistics, optionally estimating separate between-sample variances
#' (\eqn{\tau^2}) per group, which is recommended. Pooled estimates and diagnostics are written 
#' into \code{rowData(se)}. 
#' 
#' @param se A \code{SummarizedExperiment} produced by 
#'  \code{panoramic_spatialstats()} or \code{panoramic()}, containing
#'  assays \code{"yi"} (effect estimates) and \code{"vi"} (within-sample variances) for each feature (ct1, ct2, radius) and sample. 
#' @param tau2 Character. Estimator for the between-sample variance \eqn{\tau^2}. One of \code{"SJ"} (Sidik-Jonkman), \code{"REML"}, 
#'  or \code{"DL"} (DerSimonian-Laird), passed to \code{metafor::rma.uni()}. 
#' @param group_col Optional character scalar giving the name of a column 
#'  in \code{colData(se)} that defines groups (e.g. treatment) for 
#'  group-specific pooling. If \code{NULL}, a single meta-analysis is 
#'  fitted across all samples for each feature. 
#' @param BPPARAM A \code{BiocParallelParam} object controlling 
#'  parallelization across features (rows). Defaults to
#'  \code{BiocParallel::SerialParam()}. 
#' 
#' @return The input \code{SummarizedExperiment}, with additional columns appended to \code{rowData(se)}: 
#'  \itemize{
#'    \item If \code{group_col} is \code{NULL}, columns: 
#'      \code{mu_hat}, \code{se_mu}, \code{tau2}, \code{Q}, \code{I2}, \code{k}, \code{p_mu}. 
#'    \item If \code{group_col} is provided, the same set of columns 
#'      is created per group, prefixed by a group label, e.g. \code{GroupA_mu_hat}, etc. 
#' }
#' PANORAMIC meta-analysis settings are also stored in 
#' \code{metadata(se)$panoramix$meta}. 
#' 
#' @details
#' For each feature (cell-type pair and radius), the function fits a
#' univariate random-effects model via \code{metafor::rma.uni()}, using
#' the supplied \code{yi} and \code{vi} across samples. Features with
#' fewer than 2 finite observations (and positive variances) receive
#' \code{NA} for all meta-analytic quantities.
#'
#' When \code{group_col} is supplied, separate random-effects models are
#' fitted within each group level, producing group-specific pooled
#' effects and heterogeneity statistics. This supports differential
#' colocalization analysis across biological conditions.
#'
#' @examples
#' library(SpatialExperiment)
#' library(S4Vectors)
#' library(BiocParallel)
#'
#' set.seed(1)
#'
#' make_spe <- function(n_cells, sample_id) {
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
#'       sample_id = sample_id
#'     ),
#'     spatialCoords = coords
#'   )
#' }
#'
#' spe1 <- make_spe(20, "sample1")
#' spe2 <- make_spe(22, "sample2")
#'
#' spe_list <- list(sample1 = spe1, sample2 = spe2)
#'
#' design <- data.frame(
#'   sample = c("sample1", "sample2"),
#'   group  = c("group1", "group2"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Compute spatial statistics across samples
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
#' # Global random-effects meta-analysis
#' se_meta_global <- panoramic_meta(se_stats, tau2 = "SJ")
#'
#' # Group-specific meta-analysis by 'group'
#' se_meta_group <- panoramic_meta(se_stats, tau2 = "SJ", group_col = "group")
#'
#' @export 

panoramic_meta <- function(se, tau2 = c("SJ","REML","DL"), 
                           group_col = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
  tau2 <- match.arg(tau2)
  stopifnot("yi" %in% names(SummarizedExperiment::assays(se)),
            "vi" %in% names(SummarizedExperiment::assays(se)))
  
  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se))
  
  # If group_col provided, fit separate models per group
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(cd)) {
      stop("`group_col` not found in colData(se).")
    }
    
    grp <- as.character(cd[[group_col]])
    groups <- unique(grp)
    
    fit_one <- function(i) {
      y <- yi[i, ]
      v <- vi[i, ]
      
      # Initialize output for all groups
      out <- list()
      
      for (g in groups) {
        idx <- grp == g
        yg <- y[idx]
        vg <- v[idx]
        keep <- is.finite(yg) & is.finite(vg) & vg > 0
        
        prefix <- paste0(make.names(g), "_")
        
        if (sum(keep) < 2) {
          out[[paste0(prefix, "mu_hat")]] <- NA_real_
          out[[paste0(prefix, "se_mu")]] <- NA_real_
          out[[paste0(prefix, "tau2")]] <- NA_real_
          out[[paste0(prefix, "Q")]] <- NA_real_
          out[[paste0(prefix, "I2")]] <- NA_real_
          out[[paste0(prefix, "k")]] <- sum(keep)
          out[[paste0(prefix, "p_mu")]] <- NA_real_
        } else {
          fit_g <- try(metafor::rma.uni(yi = yg[keep], vi = vg[keep], 
                                        method = tau2, test = "z"),
                       silent = TRUE)
          
          if (inherits(fit_g, "try-error")) {
            out[[paste0(prefix, "mu_hat")]] <- NA_real_
            out[[paste0(prefix, "se_mu")]] <- NA_real_
            out[[paste0(prefix, "tau2")]] <- NA_real_
            out[[paste0(prefix, "Q")]] <- NA_real_
            out[[paste0(prefix, "I2")]] <- NA_real_
            out[[paste0(prefix, "k")]] <- sum(keep)
            out[[paste0(prefix, "p_mu")]] <- NA_real_
          } else {
            out[[paste0(prefix, "mu_hat")]] <- unname(fit_g$b[1,1])
            out[[paste0(prefix, "se_mu")]] <- unname(fit_g$se[1])
            out[[paste0(prefix, "tau2")]] <- unname(fit_g$tau2)
            out[[paste0(prefix, "Q")]] <- unname(fit_g$QE)
            out[[paste0(prefix, "I2")]] <- unname(fit_g$I2)
            out[[paste0(prefix, "k")]] <- sum(keep)
            out[[paste0(prefix, "p_mu")]] <- unname(fit_g$pval[1])
          }
        }
      }
      
      unlist(out)
    }
    
  } else {
    # Original behavior: single pooling across all samples
    fit_one <- function(i) {
      y <- yi[i, ]
      v <- vi[i, ]
      keep <- is.finite(y) & is.finite(v) & v > 0
      
      if (sum(keep) < 2) {
        return(c(mu_hat = NA_real_, se_mu = NA_real_, tau2 = NA_real_, 
                 Q = NA_real_, I2 = NA_real_, k = sum(keep), p_mu = NA_real_))
      }
      
      fit <- try(metafor::rma.uni(yi = y[keep], vi = v[keep], 
                                  method = tau2, test = "z"),
                 silent = TRUE)
      
      if (inherits(fit, "try-error")) {
        return(c(mu_hat = NA_real_, se_mu = NA_real_, tau2 = NA_real_, 
                 Q = NA_real_, I2 = NA_real_, k = sum(keep), p_mu = NA_real_))
      }
      
      c(mu_hat = unname(fit$b[1,1]),
        se_mu = unname(fit$se[1]),
        tau2 = unname(fit$tau2),
        Q = unname(fit$QE),
        I2 = unname(fit$I2),
        k = sum(keep),
        p_mu = unname(fit$pval[1]))
    }
  }
  
  M <- BiocParallel::bplapply(seq_len(nrow(yi)), fit_one, BPPARAM = BPPARAM)
  M <- do.call(rbind, M)
  
  # Add results to rowData
  rd_new <- S4Vectors::DataFrame(M)
  SummarizedExperiment::rowData(se) <- cbind(SummarizedExperiment::rowData(se), rd_new)
  
  # Store metadata
  m <- S4Vectors::metadata(se)
  if (is.null(m)) m <- list()
  if (is.null(m$panoramic)) m$panoramic <- list()
  m$panoramic$meta <- list(tau2_method = tau2, 
                           group_specific = !is.null(group_col),
                           group_col = group_col)
  S4Vectors::metadata(se) <- m
  
  se
}