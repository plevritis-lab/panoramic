test_that("simulation helper returns expected structure", {
  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 80,
    n_cells_group2 = 100,
    seed = 1
  )

  expect_true(is.list(toy))
  expect_true(all(c("spe_list", "design") %in% names(toy)))
  expect_length(toy$spe_list, 4)
  expect_equal(nrow(toy$design), 4)
  expect_true(all(c("sample", "group") %in% colnames(toy$design)))
})


test_that("stepwise workflow runs and exposes spatialstats table", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocParallel")

  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 90,
    n_cells_group2 = 110,
    seed = 1
  )

  prep <- panoramic_prepare(
    spe_list = toy$spe_list,
    design = toy$design,
    cell_type = "cell_type",
    min_cells = 2,
    window = "rect",
    BPPARAM = BiocParallel::SerialParam()
  )

  se <- panoramic_spatialstats(
    prep = prep,
    radii_um = c(10, 20),
    stat = "Lcross",
    nsim = 5,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_true(all(c("yi", "vi") %in% names(SummarizedExperiment::assays(se))))
  expect_true(nrow(se) > 0)
  expect_equal(ncol(se), 4)

  stats_tbl <- panoramic_extract_contrast(se, what = "spatialstats")
  expect_true(is.data.frame(stats_tbl))
  expect_true(all(c("ct1", "ct2", "radius_um", "sample", "group", "yi", "vi") %in% colnames(stats_tbl)))
})

test_that("panoramic wrapper runs end-to-end and returns tables", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocParallel")

  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 80,
    n_cells_group2 = 100,
    seed = 1
  )

  out <- panoramic(
    spe_list = toy$spe_list,
    design = toy$design,
    cell_type = "cell_type",
    radii_um = c(10, 20),
    stat = "Lcross",
    nsim = 5,
    min_cells = 2,
    window = "rect",
    group_col = "group",
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(is.list(out))
  expect_true(all(c("prep", "stats", "pooled", "tables") %in% names(out)))
  expect_s4_class(out$stats, "SummarizedExperiment")
  expect_s4_class(out$pooled, "SummarizedExperiment")

  expect_true(is.list(out$tables))
  expect_true(all(c("spatialstats", "meta", "contrast") %in% names(out$tables)))
  expect_true(is.data.frame(out$tables$spatialstats))
  expect_true(is.data.frame(out$tables$meta))

  meta_tbl <- panoramic_extract_contrast(out$pooled, what = "meta")
  expect_true(is.data.frame(meta_tbl))
  expect_true(all(c("ct1", "ct2", "radius_um") %in% colnames(meta_tbl)))

  contrast_tbl <- panoramic_extract_contrast(out$pooled, what = "contrast")
  expect_true(is.data.frame(contrast_tbl))
  expect_true(all(c("beta_diff", "se_diff", "z_diff", "p_diff", "fdr_diff") %in% colnames(contrast_tbl)))
})

test_that("retired helpers are not available", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocParallel")

  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 50,
    n_cells_group2 = 50,
    seed = 1
  )

  prep <- panoramic_prepare(
    spe_list = toy$spe_list,
    design = toy$design,
    cell_type = "cell_type",
    min_cells = 2,
    window = "rect",
    BPPARAM = BiocParallel::SerialParam()
  )
  se <- panoramic_spatialstats(
    prep = prep,
    radii_um = c(10, 20),
    stat = "Lcross",
    nsim = 5,
    BPPARAM = BiocParallel::SerialParam()
  )
  se_meta <- panoramic_meta_mv(
    se = se,
    patient_col = "sample",
    group_col = "group",
    sample_col = "sample",
    tau_structure = "patient",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_s4_class(se_meta, "SummarizedExperiment")
  expect_false(exists("panoramic_compare_groups", where = asNamespace("panoramic"), inherits = FALSE))
})

test_that("volcano and representative plotting helpers run", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocParallel")

  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 80,
    n_cells_group2 = 100,
    seed = 1
  )

  out <- panoramic(
    spe_list = toy$spe_list,
    design = toy$design,
    cell_type = "cell_type",
    radii_um = 10,
    stat = "Lcross",
    nsim = 5,
    min_cells = 2,
    window = "rect",
    group_col = "group",
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  p <- plot_volcano(
    se_diff = out$pooled,
    sig_col = "p_diff",
    alpha = 0.05,
    label_top = 3
  )
  expect_s3_class(p, "ggplot")

  reps <- plot_representative_samples(
    se_stats = out$stats,
    se_meta = out$pooled,
    spe_list = out$prep,
    sig_col = "p_diff",
    alpha = 1,
    top_n = 2,
    group_col = "group",
    cell_type_col = "cell_type",
    sample_col = "sample"
  )
  expect_true(is.list(reps))
  expect_true(all(c("plots", "index") %in% names(reps)))
  expect_true(is.list(reps$plots))
  expect_true(is.data.frame(reps$index))
})

test_that("auto pair generation uses observed marks, not unused factor levels", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  make_prepped <- function(mark_values) {
    ppp <- spatstat.geom::ppp(
      x = seq_along(mark_values),
      y = seq_along(mark_values),
      window = spatstat.geom::owin(c(0, 5), c(0, 5)),
      marks = mark_values
    )
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(dummy = matrix(0, nrow = 1, ncol = length(mark_values)))
    )
    meta <- S4Vectors::metadata(se)
    meta$panoramic <- list(ppp = ppp)
    S4Vectors::metadata(se) <- meta
    se
  }

  prep <- list(
    s1 = make_prepped(factor(c("A", "B"), levels = c("A", "B", "C"))),
    s2 = make_prepped(factor(c("A", "A"), levels = c("A", "C")))
  )

  pairs <- panoramic:::.panoramic_pairs(prep, include_self = TRUE, min_presence = 1)
  expect_equal(nrow(pairs), 1L)
  expect_identical(pairs$ct1, "A")
  expect_identical(pairs$ct2, "A")
})

test_that("panoramic_prepare validates design schema and duplicate sample ids", {
  toy_spe <- panoramic_simulate_spe(n_cells = 25, sample_id = "sample1", seed = 1)
  spe_list <- list(sample1 = toy_spe)

  bad_design_missing <- data.frame(sample = "sample1", stringsAsFactors = FALSE)
  expect_error(
    panoramic_prepare(
      spe_list = spe_list,
      design = bad_design_missing,
      cell_type = "cell_type",
      window = "rect"
    ),
    "must contain columns `sample` and `group`"
  )

  bad_design_dup <- data.frame(
    sample = c("sample1", "sample1"),
    group = c("g1", "g1"),
    stringsAsFactors = FALSE
  )
  expect_error(
    panoramic_prepare(
      spe_list = spe_list,
      design = bad_design_dup,
      cell_type = "cell_type",
      window = "rect"
    ),
    "contains duplicates"
  )
})

test_that("panoramic_meta_mv rejects ambiguous group labels after sanitization", {
  yi <- matrix(c(0.20, 0.15, 0.10, 0.05), nrow = 1)
  vi <- matrix(c(0.04, 0.05, 0.06, 0.05), nrow = 1)
  colnames(yi) <- colnames(vi) <- paste0("s", seq_len(4L))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = yi, vi = vi),
    rowData = S4Vectors::DataFrame(ct1 = "A", ct2 = "B", radius_um = 10),
    colData = S4Vectors::DataFrame(
      sample = paste0("s", seq_len(4L)),
      patient = paste0("p", seq_len(4L)),
      group = c("A-B", "A B", "C", "C")
    )
  )

  expect_error(
    panoramic_meta_mv(
      se = se,
      patient_col = "patient",
      group_col = "group",
      sample_col = "sample",
      BPPARAM = BiocParallel::SerialParam()
    ),
    "ambiguous after sanitization"
  )
})

test_that("create_spatial_network validates required columns and flags inverted thresholds", {
  se_missing <- SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = matrix(0, nrow = 1, ncol = 1), vi = matrix(1, nrow = 1, ncol = 1)),
    rowData = S4Vectors::DataFrame(ct1 = "A", ct2 = "B", p_diff = 0.01, fdr_diff = 0.01)
  )
  expect_error(
    create_spatial_network(se_missing),
    "missing required columns"
  )

  se_warn <- SummarizedExperiment::SummarizedExperiment(
    assays = list(yi = matrix(0, nrow = 2, ncol = 1), vi = matrix(1, nrow = 2, ncol = 1)),
    rowData = S4Vectors::DataFrame(
      ct1 = c("A", "B"),
      ct2 = c("B", "C"),
      z_diff = c(2.0, -1.2),
      p_diff = c(0.20, 0.30),
      fdr_diff = c(0.20, 0.30)
    )
  )
  expect_warning(
    create_spatial_network(se_warn, sig_operator = "gt", fdr_threshold = 0.05),
    "inverts the usual FDR rule"
  )
})
