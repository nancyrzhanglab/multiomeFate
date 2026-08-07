context("Test the package's public API")

## X2: cheap insurance against an accidental `@export`, which became a live risk
## once every function in `R/` acquired a roxygen block -- flipping `@noRd` to
## `@export` is a one-character edit, and nothing else would notice.
test_that("the exported names are exactly the agreed set (X2)", {
  expected_vec <- sort(c(
    # Estimation core
    "cyfer", "cyfer_finalize", "lineage_imputation",
    "lineage_imputation_sequence",
    # Plotting
    "plot_anova", "plot_cellGrowthUmap", "plot_lineageScatterplot",
    "plot_simplex", "plot_trainTest",
    # Entropy
    "compute_entropy",
    # Simulators
    "generate_simulation", "generate_simulation_attachFuture",
    "generate_simulation_plastic",
    # Barcoding pipeline
    "barcode_clustering", "barcode_combine", "barcoding_assignment",
    "barcoding_posterior"))

  expect_true(identical(sort(getNamespaceExports("multiomeFate")),
                        expected_vec))
  expect_true(length(expected_vec) == 17)
})

## These two were considered for export in the 2026-08-06 review and
## deliberately kept private: `construct_folds()` is an implementation detail of
## `cyfer()`, and `evaluate_nll()` is a scoring helper whose sign convention
## makes it easy to misuse. `data_loader()` was *un*exported in 1.0.2.001
## because it reads a hardcoded lab path.
test_that("the deliberately-private functions stay private (X2)", {
  for(function_name in c("construct_folds", "evaluate_nll", "data_loader")){
    label <- paste0("function ", function_name)
    expect_true(!(function_name %in% getNamespaceExports("multiomeFate")),
                info = label)
    expect_true(is.function(get(function_name,
                                envir = asNamespace("multiomeFate"))),
                info = label)
  }
})

## Every exported function must have a help page, or the export is not really an
## API. `man/` is roxygen-generated, so this catches an `@export` added without
## the block being filled in.
test_that("every exported function has a help page", {
  rd_vec <- sub("\\.Rd$", "", list.files(
    system.file("man", package = "multiomeFate")))
  if(length(rd_vec) == 0){
    # Under devtools::load_all() there is no installed man/ directory to read.
    skip("man/ is not available from an uninstalled package")
  }

  for(function_name in getNamespaceExports("multiomeFate")){
    expect_true(function_name %in% rd_vec, info = function_name)
  }
})
