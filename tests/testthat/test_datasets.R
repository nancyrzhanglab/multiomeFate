context("Test the bundled datasets")

## CRAN-prep 2026-09-06. The datasets shipped with 15 lineages until 1.0.3.000,
## which cannot identify a 30-feature fit: `cyfer()` requires `p + 1` lineages
## in every training fold, and the vignette had to fall back to 5 features.
## They now carry all 50 lineages of the original simulations. Pin the contract
## the vignette and the examples depend on, so a future "make the example
## faster" subset cannot silently reintroduce the problem.
test_that("the bundled datasets clear the identifiability floor at 5 folds", {
  for(dataset_name in c("priming_simulation", "plastic_simulation")){
    dataset <- get(utils::data(list = dataset_name, package = "multiomeFate",
                               envir = environment()))
    num_lineages <- length(dataset$lineage_future_count)
    num_features <- ncol(dataset$cell_features)
    # the largest of 5 round-robin folds holds ceiling(L / 5) lineages
    num_train <- num_lineages - ceiling(num_lineages / 5)
    expect_true(num_train >= num_features + 1, info = dataset_name)

    expect_true(is.matrix(dataset$cell_features), info = dataset_name)
    expect_true(!is.null(rownames(dataset$cell_features)), info = dataset_name)
    expect_true(!is.null(colnames(dataset$cell_features)), info = dataset_name)
    expect_true(!"Intercept" %in% colnames(dataset$cell_features), info = dataset_name)
    expect_true(nrow(dataset$cell_features) == length(dataset$cell_lineage),
                info = dataset_name)
    expect_true(setequal(names(dataset$lineage_future_count),
                         unique(dataset$cell_lineage)), info = dataset_name)
    expect_true(all(dataset$lineage_future_count > 0), info = dataset_name)
    expect_true(identical(rownames(dataset$tab_mat),
                          names(dataset$lineage_future_count)), info = dataset_name)
    expect_true(all(dataset$tab_mat[, "future"] == dataset$lineage_future_count),
                info = dataset_name)
    expect_true(all(dataset$tab_mat[, "now"] ==
                      as.numeric(table(dataset$cell_lineage)[rownames(dataset$tab_mat)])),
                info = dataset_name)
  }
})
