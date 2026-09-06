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

##############################################################################
## Test cyfer and cyfer_finalize (test-plan 2026-08-06)
##############################################################################


.raw_features <- function(res){
  res$cell_features[, setdiff(colnames(res$cell_features), "Intercept"),
                    drop = FALSE]
}

# See plan item P2 -- lineage names whose sort order reverses their appearance
# order, so positional indexing cannot masquerade as name indexing.
.scramble_lineage_names <- function(res){
  uniq_lineages <- unique(res$cell_lineage)
  new_names <- paste0("z", sprintf("%03d", rev(seq_along(uniq_lineages))))
  name_map <- stats::setNames(new_names, uniq_lineages)

  res$cell_lineage <- unname(name_map[res$cell_lineage])
  names(res$lineage_future_count) <-
    unname(name_map[names(res$lineage_future_count)])
  res
}

.test_loglik_mat <- function(cv){
  sapply(cv, function(lis){lis$test_loglik})
}

######################################
## Section 1 -- fold construction inside cyfer (plan items D1-D6)

## Plan item D1. `.lineage_cleanup()` exists to reconcile lineages named in
## `lineage_future_count` that have no cells at the first time point -- but
## `cyfer()` builds its folds *before* any cleanup runs, so such a lineage is
## still dealt into a fold. It contributes no cells, which silently unbalances
## the folds. Kevin: run the cleanup first, dropping cell-less lineages before
## making folds.
test_that("cyfer ignores lineages that have no cells at the first time point", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  lineage_future_count_phantom <- c(res$lineage_future_count,
                                    "lin:901" = 5,
                                    "lin:902" = 3,
                                    "lin:903" = 8)

  run <- function(lineage_future_count){
    cyfer(cell_features = cell_features,
          cell_lineage = res$cell_lineage,
          lineage_future_count = lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 5,
          num_folds = 3,
          seed_number = 1,
          verbose = 0)
  }

  expect_equal(.test_loglik_mat(run(lineage_future_count_phantom)),
               .test_loglik_mat(run(res$lineage_future_count)))
})

## Plan item D1, the destructive case. If *every* lineage in a fold has no
## cells, `cv_cell_list[[fold]]` is NULL and `cell_features[-NULL,,drop=F]`
## returns zero rows -- R's `x[-integer(0)]` trap -- so `cyfer()` trains on an
## empty matrix. With the cleanup running first there is only one real lineage
## left here, so the correct outcome is the `num_folds` error instead.
test_that("cyfer counts only cell-bearing lineages against num_folds", {
  set.seed(10)
  res <- .construct_lineage_data()

  keep_lineage <- "lin:1"
  idx_vec <- which(res$cell_lineage == keep_lineage)
  cell_features <- .raw_features(res)[idx_vec, , drop = FALSE]
  cell_lineage <- res$cell_lineage[idx_vec]

  lineage_future_count <- c(res$lineage_future_count[keep_lineage],
                            "lin:901" = 5,
                            "lin:902" = 3,
                            "lin:903" = 8)

  expect_error(
    cyfer(cell_features = cell_features,
          cell_lineage = cell_lineage,
          lineage_future_count = lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 3,
          num_folds = 2,
          seed_number = 1,
          verbose = 0),
    "num_folds"
  )
})

## Plan item D2. `seed_number = NULL` is documented as "do not touch the
## stream", and every `set.seed()` in `cyfer()` sits behind that guard.
test_that("cyfer runs with seed_number = NULL", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  set.seed(4)
  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 4,
              num_folds = 3,
              seed_number = NULL,
              verbose = 0)

  expect_s3_class(cv, "cyfer")
  expect_equal(length(cv), 3)
  expect_true(all(is.finite(.test_loglik_mat(cv))))
})

## Plan item D3. The paper's real runs all pass `savefile_tmp`, since a 20-fold
## fit takes long enough that losing it to a crash matters.
test_that("cyfer writes a resumable checkpoint to savefile_tmp", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  savefile_tmp <- tempfile(fileext = ".RData")

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 4,
              num_folds = 3,
              savefile_tmp = savefile_tmp,
              seed_number = 1,
              verbose = 0)

  expect_true(file.exists(savefile_tmp))

  load_env <- new.env()
  load(savefile_tmp, envir = load_env)
  expect_true(all(c("cv_fit_list", "date_of_run") %in% ls(load_env)))
  expect_equal(length(load_env$cv_fit_list), length(cv))
  expect_equal(.test_loglik_mat(load_env$cv_fit_list), .test_loglik_mat(cv))
})

## Plan item D4. Both boundaries are tested in `construct_folds()` but neither
## has ever been run end-to-end, where an empty training or test split would
## actually bite.
test_that("cyfer runs at both ends of the admissible num_folds range", {
  set.seed(10)
  res <- .construct_lineage_data()   # 10 lineages
  cell_features <- .raw_features(res)

  for(num_folds in c(2, 10)){
    label <- paste0("num_folds=", num_folds)

    cv <- cyfer(cell_features = cell_features,
                cell_lineage = res$cell_lineage,
                lineage_future_count = res$lineage_future_count,
                lambda_initial = 10,
                lambda_sequence_length = 4,
                num_folds = num_folds,
                seed_number = 1,
                verbose = 0)

    expect_equal(length(cv), num_folds, info = label)
    expect_true(all(is.finite(.test_loglik_mat(cv))), info = label)
    expect_true(all(is.finite(sapply(cv, function(lis){lis$train_loglik}))),
                info = label)
  }
})

## Plan item D5. `cyfer()` returns only the three loglik/fit elements, so which
## lineages a fold actually trained on is unobservable from the outside. Mock
## the inner fit and read the training split directly, then check it against the
## folds the same seed produces.
test_that("cyfer trains each fold on exactly the complement of its held-out lineages", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  num_folds <- 3
  seed_number <- 1
  lambda_sequence_length <- 4

  train_list <- list()
  fake_sequence <- function(cell_features,
                            cell_lineage,
                            lineage_future_count,
                            lambda_initial = NA,
                            lambda_max = 101,
                            lambda_min = 0.01,
                            lambda_sequence_length = 50,
                            multipler = 1e4,
                            verbose = 1){
    train_list[[length(train_list) + 1]] <<- list(
      cell_lineage = cell_lineage,
      cell_names = rownames(cell_features),
      lineage_future_count = lineage_future_count
    )

    coefficient_vec <- rep(0, ncol(cell_features) + 1)
    names(coefficient_vec) <- c("Intercept", colnames(cell_features))
    list(fit_list = replicate(lambda_sequence_length,
                              list(coefficient_vec = coefficient_vec),
                              simplify = FALSE),
         lambda_sequence = seq(lambda_initial, 0,
                               length.out = lambda_sequence_length))
  }

  with_mocked_bindings(
    cv <- cyfer(cell_features = cell_features,
                cell_lineage = res$cell_lineage,
                lineage_future_count = res$lineage_future_count,
                lambda_initial = 10,
                lambda_sequence_length = lambda_sequence_length,
                num_folds = num_folds,
                seed_number = seed_number,
                verbose = 0),
    lineage_imputation_sequence = fake_sequence,
    .package = "multiomeFate"
  )

  # `cyfer()` seeds immediately before `construct_folds()`, so this reproduces
  # the very folds it used.
  set.seed(seed_number)
  folds <- construct_folds(cell_lineage = res$cell_lineage,
                           lineage_future_count = res$lineage_future_count,
                           num_folds = num_folds)

  expect_equal(length(train_list), num_folds)

  for(i in seq_len(num_folds)){
    label <- paste0("fold ", i)
    held_out <- folds$fold_lineage_list[[i]]
    expected_train <- setdiff(names(res$lineage_future_count), held_out)
    test_cell_names <- rownames(cell_features)[folds$cv_cell_list[[i]]]

    expect_setequal(names(train_list[[i]]$lineage_future_count), expected_train)
    expect_setequal(unique(train_list[[i]]$cell_lineage), expected_train)

    # train and test must partition the cells, and neither may be empty
    expect_gt(length(train_list[[i]]$cell_names), 0)
    expect_gt(length(test_cell_names), 0)
    expect_equal(length(intersect(train_list[[i]]$cell_names, test_cell_names)),
                 0, info = label)
    expect_setequal(c(train_list[[i]]$cell_names, test_cell_names),
                    rownames(cell_features))
  }
})

## Plan item D6. The paper's own call sites take `lineage_future_count` straight
## from `table()` without the documented `> 0` filter (see
## multiomeFate_analysis/kevin/Writeup15_percentage/), so zero counts are part
## of the supported input surface whether or not the docs say so.
test_that("cyfer accepts lineages with a zero future count", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lineage_future_count <- res$lineage_future_count
  lineage_future_count[c(1, 3, 5)] <- 0

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = lineage_future_count,
              lambda_initial = 5,
              lambda_sequence_length = 5,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  expect_true(all(is.finite(.test_loglik_mat(cv))))
  expect_true(all(is.finite(sapply(cv, function(lis){lis$train_loglik}))))
})

######################################
## Section 2 -- cyfer_finalize (plan items E1-E6)

## Plan item E1. `sapply()` over a length-1 `test_loglik` returns a vector, not
## a 1-row matrix, so `apply(test_vec, 1, ...)` dies with "dim(X) must have a
## positive length" -- a message that names nothing the caller controls.
## Kevin: support a length-1 path.
test_that("cyfer_finalize supports a length-1 lambda path", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lambda_initial <- 5

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = lambda_initial,
              lambda_sequence_length = 1,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  fit <- cyfer_finalize(cell_features = cell_features,
                        cell_lineage = res$cell_lineage,
                        fit_res = cv,
                        lineage_future_count = res$lineage_future_count)

  expect_equal(fit$lambda, lambda_initial)
  expect_equal(length(fit$cell_imputed_score), nrow(cell_features))
})

## Plan item E2. The selection rule -- minimize the median across folds of the
## held-out objective -- is stated in the roxygen and in CLAUDE_kevin.md but has
## never been checked against what the function returns.
test_that("cyfer_finalize picks the lambda minimizing the median held-out objective", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 8,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  lambda_sequence <- cv[[1]]$train_fit$lambda_sequence
  median_vec <- apply(.test_loglik_mat(cv), 1, stats::median)
  lambda_expected <- lambda_sequence[which.min(median_vec)]

  fit <- cyfer_finalize(cell_features = cell_features,
                        cell_lineage = res$cell_lineage,
                        fit_res = cv,
                        lineage_future_count = res$lineage_future_count)

  expect_equal(fit$lambda, lambda_expected)
})

## Plan item E2, tie behavior. A flat CV curve means the fit never moved.
## `which.min()` on a tie returns the first index and `lambda_sequence` is
## decreasing, so the reported lambda is the ceiling the caller supplied.
## CLAUDE_kevin.md calls this out as a diagnostic symptom, so pin it rather than
## leaving it to be rediscovered.
test_that("cyfer_finalize resolves a tied CV curve to the largest lambda", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 6,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  cv_flat <- cv
  for(i in seq_along(cv_flat)){
    cv_flat[[i]]$test_loglik[] <- 1
  }

  fit <- cyfer_finalize(cell_features = cell_features,
                        cell_lineage = res$cell_lineage,
                        fit_res = cv_flat,
                        lineage_future_count = res$lineage_future_count)

  expect_equal(fit$lambda, cv[[1]]$train_fit$lambda_sequence[1])
})

## Plan item E3. `cell_imputed_score` is returned on the log10 scale while
## `lineage_imputed_count` is accumulated on the natural scale. Kevin: confirm
## there is no math error; if the two agree, the behavior is intended and gets
## documented rather than changed.
test_that("cyfer_finalize returns cell scores in log10 and counts on the natural scale", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 6,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  fit <- cyfer_finalize(cell_features = cell_features,
                        cell_lineage = res$cell_lineage,
                        fit_res = cv,
                        lineage_future_count = res$lineage_future_count)

  cell_features_intercept <- cbind(Intercept = 1, cell_features)
  score_expected <- log10(exp(as.numeric(cell_features_intercept %*%
                                           fit$coefficient_vec)))
  expect_equal(unname(fit$cell_imputed_score), score_expected)

  for(lineage in names(fit$lineage_imputed_count)){
    idx_vec <- which(res$cell_lineage == lineage)
    expect_equal(unname(fit$lineage_imputed_count[lineage]),
                 sum(10^fit$cell_imputed_score[idx_vec]),
                 label = paste0("lineage ", lineage))
  }
})

## Plan item E4. `cyfer()` takes a `seed_number` and its folds are reproducible;
## `cyfer_finalize()` then refits with ten unseeded random restarts, so the
## coefficients the paper reports are not pinned to a stream. Kevin: add a seed.
test_that("cyfer_finalize is reproducible for a given seed_number", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 6,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  run <- function(){
    cyfer_finalize(cell_features = cell_features,
                   cell_lineage = res$cell_lineage,
                   fit_res = cv,
                   lineage_future_count = res$lineage_future_count,
                   seed_number = 5)
  }

  # perturb the ambient stream between the calls: seed_number must dominate it
  stats::runif(1)
  fit_a <- run()
  stats::runif(7)
  fit_b <- run()

  expect_equal(fit_b$coefficient_vec, fit_a$coefficient_vec)
  expect_equal(fit_b$cell_imputed_score, fit_a$cell_imputed_score)
})

## Plan item E5. Cells whose lineage is absent from `lineage_future_count` are
## dropped by `.lineage_cleanup()` during the refit but still scored, because
## the full `cell_features` is multiplied through. Kevin: intended -- and every
## scored cell must be named so the caller can tell which cells came back.
test_that("cyfer_finalize names every scored cell, including unfitted ones", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  # hold one lineage out of the counts entirely; its cells are never fitted
  unfitted_lineage <- "lin:1"
  lineage_future_count <-
    res$lineage_future_count[names(res$lineage_future_count) !=
                               unfitted_lineage]
  idx_unfitted <- which(res$cell_lineage == unfitted_lineage)
  expect_gt(length(idx_unfitted), 0)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 5,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  fit <- cyfer_finalize(cell_features = cell_features,
                        cell_lineage = res$cell_lineage,
                        fit_res = cv,
                        lineage_future_count = lineage_future_count)

  expect_equal(names(fit$cell_imputed_score), rownames(cell_features))
  expect_false(anyNA(fit$cell_imputed_score))
  expect_true(all(is.finite(fit$cell_imputed_score)))
  expect_true(all(rownames(cell_features)[idx_unfitted] %in%
                    names(fit$cell_imputed_score)))
})

## Plan item E6. `cyfer_finalize()` prepends its own intercept with
## `cbind(1, cell_features)` unconditionally, so a caller who already supplied a
## constant column silently gets a collinear design -- and one named "Intercept"
## produces two columns of that name, caught only by a bare `stopifnot`.
## Kevin: detect a constant column and error.
test_that("cyfer_finalize rejects a constant column in cell_features", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- cbind(.raw_features(res), constant_col = 1)

  cv <- cyfer(cell_features = cell_features,
              cell_lineage = res$cell_lineage,
              lineage_future_count = res$lineage_future_count,
              lambda_initial = 10,
              lambda_sequence_length = 4,
              num_folds = 3,
              seed_number = 1,
              verbose = 0)

  expect_error(
    cyfer_finalize(cell_features = cell_features,
                   cell_lineage = res$cell_lineage,
                   fit_res = cv,
                   lineage_future_count = res$lineage_future_count),
    "constant"
  )
})

######################################
## Section 3 -- invariances through the whole CV pipeline (plan item F3)
##
## `lineage_imputation()` alone cannot exercise the fold machinery, which is
## where the ordering assumptions actually live: `construct_folds()` sorts by
## count and deals by position, then `cyfer()` subsets cells by index and counts
## by name. Both historical bugs lived in exactly that seam.

test_that("cyfer is invariant to the order of the cells", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  run <- function(cell_features, cell_lineage){
    cyfer(cell_features = cell_features,
          cell_lineage = cell_lineage,
          lineage_future_count = res$lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 5,
          num_folds = 3,
          seed_number = 1,
          verbose = 0)
  }

  set.seed(2)
  idx_vec <- sample(nrow(cell_features))

  expect_equal(.test_loglik_mat(run(cell_features[idx_vec, , drop = FALSE],
                                    res$cell_lineage[idx_vec])),
               .test_loglik_mat(run(cell_features, res$cell_lineage)))
})

test_that("cyfer is invariant to relabelling the lineages", {
  set.seed(10)
  res <- .construct_lineage_data()
  res_scrambled <- .scramble_lineage_names(res)
  cell_features <- .raw_features(res)

  run <- function(cell_lineage, lineage_future_count){
    cyfer(cell_features = cell_features,
          cell_lineage = cell_lineage,
          lineage_future_count = lineage_future_count,
          lambda_initial = 10,
          lambda_sequence_length = 5,
          num_folds = 3,
          seed_number = 1,
          verbose = 0)
  }

  expect_equal(.test_loglik_mat(run(res_scrambled$cell_lineage,
                                    res_scrambled$lineage_future_count)),
               .test_loglik_mat(run(res$cell_lineage,
                                    res$lineage_future_count)))
})


test_that("cyfer resolves lambda_initial = NA once, so every fold shares one lambda_sequence", {
  # `lineage_imputation_sequence()` used to re-derive `lambda_initial` inside
  # each fold, from that fold's TRAINING lineages. Every input to
  # `.compute_initial_parameters()` (`future_total`, `current_total`, `term2`,
  # `num_lineages`) is a sum or count over the lineages it is handed, so each fold
  # built a different grid -- while `cyfer_finalize()` medians the folds'
  # `test_loglik` vectors by POSITION and reads the winning lambda off fold 1's
  # grid alone. Row k was then an average of held-out error measured at k
  # different penalties. Nothing errored, because the paths share a LENGTH.
  #
  # The fix resolves `lambda_initial` once, on the full data, before the folds
  # are built (R/lineage_cv.R).
  #
  # This test needs the BOTTLENECK regime -- one lineage carrying nearly all the
  # future mass. On benign data the per-fold heuristic lands in the same place
  # whichever lineages are held out, so the buggy and fixed trees agree and a
  # test built on the default fixture would pass against both. Here, whether the
  # dominant lineage is in the training set moves `future_total` by two orders of
  # magnitude, which is what makes the grids diverge.
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  # One clone carries nearly all the mass, but the rest stay non-zero on purpose.
  # Zeroing them makes some fold's TRAINING response almost entirely zero, which
  # is a degenerate fit -- BFGS wanders and trips the `maxit` warning, adding
  # noise unrelated to what is being tested here.
  lineage_future_count <- res$lineage_future_count
  lineage_future_count[] <- c(80, 3, 2, 2, 1, 2, 1, 3, 2, 1)

  # The fixture has teeth only if the heuristic actually MOVES when the dominant
  # lineage is withheld -- otherwise "all folds agree" is vacuous. Assert that
  # first, so this test cannot silently degrade into a tautology.
  lambda_full <- .compute_initial_parameters(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = lineage_future_count
  )$lambda_initial

  keep <- res$cell_lineage != names(lineage_future_count)[1]
  lambda_without_dominant <- .compute_initial_parameters(
    cell_features = cell_features[keep, , drop = FALSE],
    cell_lineage = res$cell_lineage[keep],
    lineage_future_count = lineage_future_count[-1]
  )$lambda_initial

  expect_false(isTRUE(all.equal(lambda_full, lambda_without_dominant)))

  cv <- cyfer(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = lineage_future_count,
    lambda_initial = NA,
    lambda_sequence_length = 6,
    num_folds = 5,
    seed_number = 1,
    verbose = 0
  )

  path_list <- lapply(cv, function(x) x$train_fit$lambda_sequence)

  # every fold fit on the SAME grid, elementwise -- not merely the same length
  for (kk in seq_along(path_list)[-1]) {
    expect_equal(path_list[[kk]], path_list[[1]])
  }

  # and that shared grid is the one derived from the FULL data, which is what
  # distinguishes "resolved once up front" from "every fold happened to agree"
  expect_equal(path_list[[1]][1], lambda_full)
})

test_that("cyfer refuses a fit whose unpenalized endpoint is not identified", {
  # H5. CYFER's effective sample size is the LINEAGE count, not the cell count:
  # the Fisher information has rank at most min(L, p+1) however many cells were
  # sequenced. The lambda path always ends at exactly 0, so the last fit on every
  # path is unpenalized -- and underdetermined whenever the training folds hold
  # fewer than p+1 lineages. Nothing errored before the guard; `optim()` returned
  # whatever point it reached on a flat ridge.
  #
  # The guard is in cyfer(), after the folds are built, and compares p+1 against
  # the SMALLEST training set, i.e. L minus the largest fold.
  make_data <- function(L, p, n_per = 8, seed = 1) {
    set.seed(seed)
    cell_features <- matrix(stats::rnorm(L * n_per * p), ncol = p)
    rownames(cell_features) <- paste0("cell", seq_len(nrow(cell_features)))
    colnames(cell_features) <- paste0("f", seq_len(p))
    cell_lineage <- rep(paste0("clone", seq_len(L)), each = n_per)
    names(cell_lineage) <- rownames(cell_features)
    list(cell_features = cell_features,
         cell_lineage = cell_lineage,
         lineage_future_count = stats::setNames(stats::rpois(L, 5) + 1,
                                                paste0("clone", seq_len(L))))
  }

  run <- function(L, p, num_folds, n_per = 8) {
    d <- make_data(L, p, n_per)
    cyfer(cell_features = d$cell_features,
          cell_lineage = d$cell_lineage,
          lineage_future_count = d$lineage_future_count,
          lambda_initial = 3,
          lambda_sequence_length = 2,
          num_folds = num_folds,
          seed_number = 10,
          verbose = 0)
  }

  # L = 6 over 3 folds leaves 4 training lineages, against 9 coefficients
  expect_error(run(L = 6, p = 8, num_folds = 3),
               "training folds have 4 lineages but 9 coefficients")
  expect_error(run(L = 6, p = 8, num_folds = 3), "not identified")

  # The boundary is exact. Same L and same folds, so the training set is 4
  # lineages either way and only the coefficient count moves: p+1 = 4 is
  # admissible, p+1 = 5 is not. n_per = 200 here only to keep the output clean --
  # at the floor the unpenalized fit is marginally determined and warns about
  # convergence with few cells. That is conditioning, not identifiability, and it
  # does not move the guard either way.
  expect_s3_class(run(L = 6, p = 3, num_folds = 3, n_per = 200), "cyfer")
  expect_error(run(L = 6, p = 4, num_folds = 3),
               "training folds have 4 lineages but 5 coefficients")

  # It is the LINEAGE count that binds. Sequencing 25x the cells
  # adds no rank to the Fisher information, so the verdict must not move.
  expect_error(run(L = 6, p = 4, num_folds = 3, n_per = 8),
               "training folds have 4 lineages but 5 coefficients")
  expect_error(run(L = 6, p = 4, num_folds = 3, n_per = 200),
               "training folds have 4 lineages but 5 coefficients")

  # More lineages lift the same feature count over the floor
  expect_s3_class(run(L = 8, p = 5, num_folds = 4, n_per = 200), "cyfer")
})

test_that("cyfer catches a cell_lineage that is not row-aligned with cell_features", {
  # `cell_lineage` is matched to `cell_features` by POSITION --
  # `.lineage_cleanup()` builds `cell_lineage_idx_list` with
  # `which(cell_lineage == lineage)`, and those integers index ROWS of
  # `cell_features`. A permuted `cell_lineage` therefore fits the wrong cells to
  # every lineage, and does it silently: nothing is missing, nothing is NA, and
  # the fit runs to completion on a scrambled design.
  #
  # The guard is at the top of `.lineage_cleanup()`: a length check always, plus
  # a names-vs-rownames check WHEN THE CALLER KEPT THE NAMES. It reaches the
  # cyfer() path only because `cyfer()` coerces with
  # `setNames(as.character(x), names(x))`; a plain `as.character()` drops the
  # names and leaves nothing to check against.
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  cell_lineage <- stats::setNames(res$cell_lineage, rownames(cell_features))

  run <- function(cell_features, cell_lineage, lambda_initial = NA) {
    cyfer(cell_features = cell_features,
          cell_lineage = cell_lineage,
          lineage_future_count = res$lineage_future_count,
          lambda_initial = lambda_initial,
          lambda_sequence_length = 3,
          num_folds = 3,
          seed_number = 10,
          verbose = 0)
  }

  # correctly aligned and named: must still run, or the guard is useless
  expect_s3_class(run(cell_features, cell_lineage), "cyfer")

  # Same cells, permuted. Names travel with the values, so the
  # vector is self-consistent and only the comparison against rownames sees it.
  set.seed(2)
  permuted <- cell_lineage[sample(length(cell_lineage))]
  expect_error(run(cell_features, permuted), "DIFFERENT ORDER")

  # The same misalignment is caught on the numeric-lambda path too. 
  # With `lambda_initial = NA` the guard first sees the FULL data, 
  # where the two name the same cells and the mismatch is a permutation. 
  # With a numeric lambda that full-data call is skipped, so the first 
  # check happens inside a fold -- and subsetting the same POSITIONS 
  # out of a permuted vector and its matrix yields two genuinely
  # different cell sets.
  expect_error(run(cell_features, permuted, lambda_initial = 3),
               "name different cells")

  # A permutation is detectable only while the names survive. This is the
  # documented escape hatch: an UNNAMED cell_lineage asserts "already aligned",
  # and the same scrambled vector then fits silently. Pinned deliberately, so
  # that dropping the name-preserving coercion shows up as a failure here.
  expect_s3_class(run(cell_features, unname(permuted)), "cyfer")

  # names that refer to different cells entirely
  relabelled <- stats::setNames(cell_lineage, paste0("other", seq_along(cell_lineage)))
  expect_error(run(cell_features, relabelled), "name different cells")

  # named lineage against a feature matrix with no row names: nothing to check
  # against, so refuse rather than assume
  unnamed_rows <- cell_features
  rownames(unnamed_rows) <- NULL
  expect_error(run(unnamed_rows, cell_lineage), "row names")

  # length mismatch fires whether or not the names survive
  expect_error(run(cell_features, cell_lineage[-1]),
               "but `cell_features` has")
  expect_error(run(cell_features, unname(cell_lineage)[-1]),
               "but `cell_features` has")
})

## CRAN-prep 2026-09-06 (D1). The paper's real-data scripts pass
## `assigned_lineage` straight from the Seurat object, where every cell the
## barcode assignment left unassigned is `NA`, and they rely on those cells
## being excluded from the fit but still scored by `cyfer_finalize()` -- the
## same contract as a cell whose lineage NAME is absent from
## `lineage_future_count` (plan item E5). A guard that errors on `NA` while
## passing an absent name is a trap, so the two cases are pinned together.
test_that("cyfer excludes NA-lineage cells from the fit, scores them, and says so", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  cell_lineage_na <- res$cell_lineage
  na_idx <- c(1, 5, 9)
  cell_lineage_na[na_idx] <- NA
  cell_lineage_absent <- res$cell_lineage
  cell_lineage_absent[na_idx] <- "unassigned"

  run <- function(cell_lineage){
    cv <- cyfer(cell_features = cell_features,
                cell_lineage = cell_lineage,
                lineage_future_count = res$lineage_future_count,
                lambda_initial = 3,
                lambda_sequence_length = 3,
                num_folds = 3,
                seed_number = 10,
                verbose = 0)
    cyfer_finalize(cell_features = cell_features,
                   cell_lineage = cell_lineage,
                   fit_res = cv,
                   lineage_future_count = res$lineage_future_count)
  }

  expect_message(fit_na <- run(cell_lineage_na), "3 of .* cells have an `NA` lineage")
  fit_absent <- run(cell_lineage_absent)

  # every cell is scored and named, the NA ones included
  expect_equal(names(fit_na$cell_imputed_score), rownames(cell_features))
  expect_true(all(is.finite(fit_na$cell_imputed_score[na_idx])))

  # NA and an absent name are the same case: identical fits
  expect_equal(fit_na$coefficient_vec, fit_absent$coefficient_vec)
  expect_equal(fit_na$cell_imputed_score, fit_absent$cell_imputed_score)
  expect_equal(fit_na$lambda, fit_absent$lambda)
})

test_that("cyfer rejects a non-positive lambda_initial", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  for(lambda_initial in c(0, -1)){
    expect_error(cyfer(cell_features = cell_features,
                       cell_lineage = res$cell_lineage,
                       lineage_future_count = res$lineage_future_count,
                       lambda_initial = lambda_initial,
                       lambda_sequence_length = 3,
                       num_folds = 3,
                       verbose = 0),
                 "single positive number",
                 info = paste0("lambda_initial = ", lambda_initial))
  }
})
