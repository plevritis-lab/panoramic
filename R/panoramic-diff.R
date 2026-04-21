#' Extract contrast results from PANORAMIC meta-analysis
#'
#' Convenience helper to retrieve group contrast results (e.g., \code{beta_diff},
#' \code{p_diff}) from a \code{SummarizedExperiment} produced by
#' \code{panoramic_meta_mv(...)} when contrast columns are present.
#'
#' @param se A \code{SummarizedExperiment} with contrast columns in
#'  \code{rowData(se)}.
#' @param feature_cols Optional character vector of feature columns to include
#'  (e.g., \code{c("ct1","ct2","radius_um")}). Only columns that exist will
#'  be returned.
#' @param what One of \code{"contrast"}, \code{"meta"}, or
#'  \code{"spatialstats"}.
#' @param drop_na Logical; if \code{TRUE}, drop rows with missing extracted
#'  statistics.
#'
#' @return A \code{data.frame} with requested feature columns and contrast
#'  statistics (\code{beta_diff}, \code{se_diff}, \code{z_diff},
#'  \code{p_diff}, \code{fdr_diff}) when \code{what = "contrast"}.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(dummy = matrix(0, nrow = 1, ncol = 1)),
#'   rowData = S4Vectors::DataFrame(
#'     ct1 = "A",
#'     ct2 = "B",
#'     radius_um = 10,
#'     beta_diff = 0.25,
#'     se_diff = 0.10,
#'     z_diff = 2.5,
#'     p_diff = 0.012,
#'     fdr_diff = 0.02
#'   )
#' )
#' panoramic_extract_contrast(se)
#'
#' @export
panoramic_extract_contrast <- function(
    se,
    feature_cols = c("ct1", "ct2", "radius_um"),
    what = c("contrast", "meta", "spatialstats"),
    drop_na = FALSE
) {
  what <- match.arg(what)

  if (identical(what, "spatialstats")) {
    return(.panoramic_extract_spatialstats_table(se, drop_na = drop_na))
  }

  if (identical(what, "meta")) {
    return(.panoramic_extract_meta_table(
      se = se,
      feature_cols = unique(c(feature_cols, "stat")),
      drop_na = drop_na
    ))
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  contrast_cols <- c("beta_diff", "se_diff", "z_diff", "p_diff", "fdr_diff")
  missing <- setdiff(contrast_cols, colnames(rd))
  if (length(missing) > 0) {
    stop("Contrast columns not found in rowData(se): ",
         paste(missing, collapse = ", "),
         ". Run panoramic_meta_mv() first.")
  }
  feature_cols <- intersect(feature_cols, colnames(rd))
  out <- rd[, c(feature_cols, contrast_cols), drop = FALSE]
  if (isTRUE(drop_na)) {
    out <- out[stats::complete.cases(out[, contrast_cols, drop = FALSE]), , drop = FALSE]
  }
  out <- .panoramic_add_coloc_direction(out)
  rownames(out) <- NULL
  out
}
