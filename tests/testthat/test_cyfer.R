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

######################################
## Factor-valued cell_lineage end-to-end
##
## With the gradient bug present, optim() never moved inside a fold, so
## test_loglik was bit-identical at every lambda and which.min() on the tie
## returned index 1 -- always the largest lambda in the (decreasing) sequence.

test_that("cyfer gives the same result for factor and character cell_lineage", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  run <- function(cell_lineage) {
    cyfer(cell_features = cell_features,
          cell_lineage = cell_lineage,
          lineage_future_count = res$lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 8,
          num_folds = 3,
          seed_number = 1,
          verbose = 0)
  }

  cv_chr <- run(res$cell_lineage)
  cv_fac <- run(factor(res$cell_lineage))

  expect_equal(sapply(cv_fac, function(x) x$test_loglik),
               sapply(cv_chr, function(x) x$test_loglik))
})

test_that("cyfer held-out loglik varies across the lambda path", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = factor(res$cell_lineage),
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 8,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  for (fold in seq_along(cv)) {
    expect_gt(diff(range(cv[[fold]]$test_loglik)), 0)
  }
})

test_that("cyfer does not pin the selected lambda to the top of the path", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  selected_lambda <- function(lambda_initial) {
    cv <- cyfer(cell_features = cell_features,
                cell_lineage = factor(res$cell_lineage),
                lineage_future_count = res$lineage_future_count,
                lambda_initial = lambda_initial,
                lambda_sequence_length = 10,
                num_folds = 3,
                seed_number = 1,
                verbose = 0)
    lambda_sequence <- cv[[1]]$train_fit$lambda_sequence
    test_mat <- sapply(cv, function(x) x$test_loglik)
    lambda_sequence[which.min(apply(test_mat, 1, stats::median))]
  }

  # Raising the ceiling must not drag the selection along with it.
  expect_lt(selected_lambda(50), 50)
  expect_lt(selected_lambda(500), 50)
})

test_that("cyfer folds are reproducible for a given seed_number", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  run <- function() {
    cyfer(cell_features = cell_features,
          cell_lineage = res$cell_lineage,
          lineage_future_count = res$lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 4,
          num_folds = 3,
          seed_number = 5,
          verbose = 0)
  }

  # Perturb the ambient RNG state between the two calls: seed_number is
  # documented as controlling reproducibility, so it must dominate.
  stats::runif(1)
  a <- run()
  stats::runif(7)
  b <- run()

  expect_equal(sapply(b, function(x) x$test_loglik),
               sapply(a, function(x) x$test_loglik))
})

test_that("cyfer rejects more folds than lineages", {
  set.seed(10)
  res <- .construct_lineage_data()   # 10 lineages
  cell_features <- .raw_features(res)

  expect_error(
    cyfer(cell_features = cell_features,
          cell_lineage = res$cell_lineage,
          lineage_future_count = res$lineage_future_count,
          lambda_initial = 1,
          lambda_sequence_length = 3,
          num_folds = 13,
          verbose = 0),
    "num_folds"
  )
})
