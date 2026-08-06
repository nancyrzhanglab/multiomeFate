context("Test cyfer and cyfer_finalize (test-plan 2026-08-06)")

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

