context("Test cyfer")

# .construct_lineage_data() calls .lineage_cleanup() internally, which prepends
# an Intercept column. Strip it so tests match what end users would pass.
.raw_features <- function(res) {
  res$cell_features[, setdiff(colnames(res$cell_features), "Intercept"), drop = FALSE]
}

test_that("cyfer errors if cell_features has no row names", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  rownames(cell_features) <- NULL

  expect_error(
    cyfer(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      lineage_future_count = res$lineage_future_count,
      lambda_initial = NA,
      lambda_sequence_length = 5,
      num_folds = 3
    ),
    "row names"
  )
})

test_that("cyfer errors if cell_features has no column names", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  colnames(cell_features) <- NULL

  expect_error(
    cyfer(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      lineage_future_count = res$lineage_future_count,
      lambda_initial = NA,
      lambda_sequence_length = 5,
      num_folds = 3
    ),
    "column names"
  )
})

test_that("cyfer runs end-to-end and returns correct structure", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count,
    lambda_initial = NA,
    lambda_sequence_length = 10,
    num_folds = 3,
    seed_number = 1,
    verbose = 0
  )

  expect_s3_class(cv, "cyfer")
  expect_equal(length(cv), 3)
  expect_true(all(c("test_loglik", "train_loglik", "train_fit") %in% names(cv[[1]])))
  expect_equal(length(cv[[1]]$test_loglik), 10)
})

test_that("cyfer_finalize runs after cyfer", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count,
    lambda_initial = NA,
    lambda_sequence_length = 10,
    num_folds = 3,
    seed_number = 1,
    verbose = 0
  )

  fit <- cyfer_finalize(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    fit_res = cv,
    lineage_future_count = res$lineage_future_count
  )

  expect_true(is.list(fit))
  expect_true(all(c("cell_imputed_score", "coefficient_vec", "lambda",
                     "lineage_imputed_count") %in% names(fit)))
  expect_equal(length(fit$cell_imputed_score), nrow(cell_features))
  expect_true(is.numeric(fit$lambda) && length(fit$lambda) == 1)
})
