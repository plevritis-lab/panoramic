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


test_that("panoramic workflow runs on simulated data", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocParallel")

  toy <- panoramic_simulate_dataset(
    n_group1 = 2,
    n_group2 = 2,
    n_cells_group1 = 90,
    n_cells_group2 = 110,
    seed = 1
  )

  se <- panoramic(
    spe_list = toy$spe_list,
    design = toy$design,
    cell_type = "cell_type",
    radii_um = c(10, 20),
    stat = "Lcross",
    nsim = 5,
    min_cells = 2,
    window = "rect",
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_true(all(c("yi", "vi") %in% names(SummarizedExperiment::assays(se))))
  expect_true(nrow(se) > 0)
  expect_equal(ncol(se), 4)
})
