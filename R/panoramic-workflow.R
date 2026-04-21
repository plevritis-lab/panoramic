#' Internal: extract one metadata value per sample from a SpatialExperiment list
#'
#' @param spe_list Named list of SpatialExperiment objects.
#' @param field Character scalar naming a colData column in each sample.
#'
#' @return Character vector aligned to names(spe_list), one value per sample.
#' @keywords internal
.panoramic_one_value_per_sample <- function(spe_list, field) {
  vals <- vapply(names(spe_list), function(sid) {
    spe <- spe_list[[sid]]
    cd <- SummarizedExperiment::colData(spe)
    if (!field %in% colnames(cd)) {
      stop(
        "`", field, "` not found in colData for sample `", sid, "`.",
        call. = FALSE
      )
    }
    x <- unique(as.character(cd[[field]]))
    x <- x[!is.na(x)]
    if (length(x) != 1L) {
      stop(
        "Expected exactly one non-missing `", field, "` value in sample `", sid,
        "`. Found: ", paste(x, collapse = ", "),
        call. = FALSE
      )
    }
    x
  }, character(1))
  vals
}

#' Internal: attach sample-level metadata to PANORAMIC colData
#'
#' @param se_stats SummarizedExperiment from panoramic_spatialstats().
#' @param spe_list Named list of SpatialExperiment objects used to produce se_stats.
#' @param fields Character vector of colData fields to attach.
#'
#' @return \code{se_stats} with additional columns in \code{colData(se_stats)}.
#' @keywords internal
.panoramic_attach_sample_metadata <- function(se_stats, spe_list, fields) {
  stopifnot("sample" %in% colnames(SummarizedExperiment::colData(se_stats)))
  if (length(fields) == 0L) return(se_stats)

  if (is.null(names(spe_list))) {
    stop("`spe_list` must be named to attach sample metadata.", call. = FALSE)
  }

  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se_stats))
  sid <- as.character(cd$sample)
  if (!all(sid %in% names(spe_list))) {
    miss <- setdiff(unique(sid), names(spe_list))
    stop(
      "Some se_stats samples are not present in names(spe_list): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  for (field in unique(fields)) {
    if (field %in% colnames(cd)) next
    map <- .panoramic_one_value_per_sample(spe_list, field)
    cd[[field]] <- map[match(sid, names(map))]
  }

  SummarizedExperiment::colData(se_stats) <- S4Vectors::DataFrame(cd)
  se_stats
}

#' Run PANORAMIC end-to-end, including pooling/meta-analysis
#'
#' \code{panoramic_analyze()} provides a streamlined API for the most common
#' workflow: prepare data, compute spatial statistics, and run
#' pooling/meta-analysis in one call.
#'
#' @param spe_list Named list of SpatialExperiment objects (one per sample).
#' @param design data.frame with at least columns \code{sample} and \code{group}.
#' @param cell_type Character; colData column containing cell-type labels.
#' @param pairs Either \code{"auto"} or a data.frame with columns \code{ct1}, \code{ct2}.
#' @param radii_um Numeric vector of radii in microns.
#' @param stat Character spatial statistic passed to \code{panoramic_spatialstats()}.
#' @param nsim Integer bootstrap replicates.
#' @param correction Optional edge-correction method for spatstat-based
#'  statistics. Ignored when \code{stat = "local_comp_enrichment"}.
#' @param min_cells Minimum cells per type per sample.
#' @param concavity Concavity for concave hull windows.
#' @param window One of \code{"concave"}, \code{"convex"}, or \code{"rect"}.
#' @param group_col Character group column for meta-analysis (defaults to \code{"group"}).
#' @param group_tau2 Controls whether PANORAMIC additionally computes
#'  group-specific heterogeneity (\code{tau2}) by fitting per-group multilevel
#'  models. Use \code{"none"} for faster fitting or \code{"separate"} for more
#'  detailed heterogeneity summaries.
#' @param patient_col Optional patient-id column for multilevel pooling. If
#'  \code{NULL}, PANORAMIC uses \code{"patient"} from \code{design} when
#'  available, otherwise falls back to \code{sample_col}.
#' @param sample_col Optional sample-id column used in
#'  \code{panoramic_meta_mv()}. If \code{NULL}, PANORAMIC uses
#'  \code{"sample"} (aligned to \code{design$sample}).
#' @param tau_structure Random-effects structure passed to
#'  \code{panoramic_meta_mv()}.
#' @param method_mv Method passed to \code{metafor::rma.mv()} in \code{panoramic_meta_mv()}.
#' @param vi_floor Variance-flooring mode in \code{panoramic_meta_mv()}.
#'  Default is \code{"group_median"}. Use \code{"median"} to replace
#'  non-positive \code{vi} with per-feature median positive \code{vi}, or
#'  \code{"none"} to disable flooring.
#' @param seed Optional random seed.
#' @param boot Bootstrap mode \code{"approx"} or \code{"block"}.
#' @param tile_size Optional tile size for block bootstrap.
#' @param nx,ny Optional tile counts for block bootstrap.
#' @param BPPARAM BiocParallel backend.
#' @param verbose Logical verbosity passed to \code{panoramic_spatialstats()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{prep}: output of \code{panoramic_prepare()}
#'   \item \code{stats}: output of \code{panoramic_spatialstats()}
#'   \item \code{pooled}: output of \code{panoramic_meta_mv()}
#'   \item \code{tables}: pre-flattened data.frames for convenient result extraction
#' }
#'
#' @examples
#' toy <- panoramic_simulate_dataset(seed = 1)
#' out <- panoramic_analyze(
#'   spe_list = toy$spe_list,
#'   design = toy$design,
#'   cell_type = "cell_type",
#'   radii_um = c(10, 20),
#'   stat = "local_comp_enrichment",
#'   nsim = 5,
#'   min_cells = 2,
#'   window = "rect",
#'   group_col = "group",
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#' names(out)
#'
#' @export
panoramic_analyze <- function(
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
  window <- match.arg(window)
  group_tau2 <- match.arg(group_tau2)
  tau_structure <- match.arg(tau_structure)
  if (is.null(vi_floor) || (length(vi_floor) == 1L && is.na(vi_floor))) {
    vi_floor <- "group_median"
  } else {
    vi_floor <- match.arg(vi_floor, c("group_median", "median", "none"))
  }

  if (!is.data.frame(design) || !"sample" %in% colnames(design)) {
    stop("`design` must be a data.frame containing a `sample` column.", call. = FALSE)
  }
  design$sample <- as.character(design$sample)
  if (!is.null(sample_col) && length(sample_col) != 1L) {
    stop("`sample_col` must be NULL or a single column name.", call. = FALSE)
  }
  if (is.null(sample_col) || is.na(sample_col) || !nzchar(sample_col)) {
    sample_col <- "sample"
  }
  if (!is.null(patient_col) && length(patient_col) != 1L) {
    stop("`patient_col` must be NULL or a single column name.", call. = FALSE)
  }
  if (is.null(patient_col) || is.na(patient_col) || !nzchar(patient_col)) {
    patient_col <- if ("patient" %in% colnames(design)) "patient" else sample_col
  }

  if (identical(stat, "local_comp_enrichment")) {
    if (!is.null(correction) &&
        !(length(correction) == 1L &&
          is.character(correction) &&
          identical(correction, "translate"))) {
      warning(
        "`correction` is ignored when `stat = \"local_comp_enrichment\"`.",
        call. = FALSE
      )
    }
    correction_use <- "translate"
  } else if (is.null(correction) || (length(correction) == 1L && !nzchar(correction))) {
    correction_use <- "translate"
  } else {
    correction_use <- correction
  }

  prep <- panoramic_prepare(
    spe_list = spe_list,
    design = design,
    cell_type = cell_type,
    min_cells = min_cells,
    concavity = concavity,
    window = window,
    BPPARAM = BPPARAM
  )

  se_stats <- panoramic_spatialstats(
    prep = prep,
    pairs = pairs,
    radii_um = radii_um,
    stat = stat,
    nsim = nsim,
    correction = correction_use,
    seed = seed,
    boot = boot,
    tile_size = tile_size,
    nx = nx,
    ny = ny,
    BPPARAM = BPPARAM,
    verbose = verbose
  )

  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se_stats))
  idx <- match(as.character(cd$sample), design$sample)
  if (anyNA(idx)) {
    miss <- unique(as.character(cd$sample)[is.na(idx)])
    stop(
      "Some se_stats samples are missing from `design$sample`: ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  design_match <- design[idx, , drop = FALSE]
  fields <- unique(c(patient_col, sample_col, group_col))
  for (field in fields) {
    if (is.null(field) || !nzchar(field) || field %in% colnames(cd)) next
    if (field %in% colnames(design_match)) {
      cd[[field]] <- as.character(design_match[[field]])
    }
  }
  SummarizedExperiment::colData(se_stats) <- S4Vectors::DataFrame(cd)

  se_stats_mv <- .panoramic_attach_sample_metadata(se_stats, prep, fields = fields)
  se_pooled <- panoramic_meta_mv(
    se = se_stats_mv,
    patient_col = patient_col,
    group_col = group_col,
    sample_col = sample_col,
    tau_structure = tau_structure,
    method = method_mv,
    group_tau2 = group_tau2,
    vi_floor = vi_floor,
    BPPARAM = BPPARAM
  )

  contrast_tbl <- tryCatch(
    panoramic_extract_contrast(se_pooled),
    error = function(e) NULL
  )

  list(
    prep = prep,
    stats = se_stats,
    pooled = se_pooled,
    tables = list(
      spatialstats = .panoramic_extract_spatialstats_table(se_stats, drop_na = FALSE),
      meta = .panoramic_extract_meta_table(se_pooled, drop_na = FALSE),
      contrast = contrast_tbl
    )
  )
}
