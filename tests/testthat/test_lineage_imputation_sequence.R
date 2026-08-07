context("Test lineage imputation sequence")

## .compute_initial_parameters is correct

test_that(".compute_initial_parameters works", {
  set.seed(10)
  n_each <- 30
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  
  res <- .compute_initial_parameters(cell_features = cell_features,
                                     cell_lineage = cell_lineage,
                                     lineage_future_count = lineage_future_count,
                                     multipler = 10)
  
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("coefficient_initial", "lambda_initial")))
  expect_true(res$lambda_initial > 0)
  expect_true(all(names(res$coefficient_initial) == colnames(cell_features)))
})

## lineage_imputation_sequence is correct

test_that("lineage_imputation_sequence works", {
  trials <- 10
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    
    res <- .construct_lineage_data()
    cell_features <- res$cell_features
    cell_lineage <- res$cell_lineage
    cell_lineage_idx_list <- res$cell_lineage_idx_list
    lineage_future_count <- res$lineage_future_count
    
    res <- lineage_imputation_sequence(cell_features = cell_features,
                                       cell_lineage = cell_lineage,
                                       lineage_future_count = lineage_future_count,
                                       lambda_sequence_length = 25,
                                       verbose = 0)
    
    coef_mat <- sapply(res$fit_list, function(x){x$coefficient_vec})
    
    bool1 <- any(!is.na(coef_mat))
    bool2 <- all(dim(coef_mat) == c(3, 25))
    bool3 <- is.list(res)
    
    all(c(bool1, bool2, bool3))
  })
  
  expect_true(all(bool_vec))
})

######################################
## Auto-computed lambda_initial
##
## lambda_max and lambda_min both defaulted to 101, so
## `min(max(lambda_initial, lambda_max), lambda_min)` collapsed to the constant
## 101 and discarded the data-driven heuristic on the line above.

# Characterization test only -- this one also passes against the pre-fix code,
# since the old min(max(x, 8), 2) == 2 happens to land inside [2, 8]. The
# "does not collapse to a constant" test below is the actual regression guard.
test_that(".compute_initial_parameters respects lambda_min and lambda_max", {
  set.seed(10)
  res <- .construct_lineage_data()

  out <- .compute_initial_parameters(cell_features = res$cell_features,
                                     cell_lineage = res$cell_lineage,
                                     lineage_future_count = res$lineage_future_count,
                                     lambda_min = 2,
                                     lambda_max = 8,
                                     multipler = 10)

  expect_gte(out$lambda_initial, 2)
  expect_lte(out$lambda_initial, 8)
})

test_that(".compute_initial_parameters does not collapse to a constant", {
  set.seed(10)
  res <- .construct_lineage_data()

  # A wide bracket must let the data-driven value through rather than pinning
  # to an endpoint.
  low <- .compute_initial_parameters(cell_features = res$cell_features,
                                     cell_lineage = res$cell_lineage,
                                     lineage_future_count = res$lineage_future_count,
                                     lambda_min = 0, lambda_max = 1e6,
                                     multipler = 1)
  high <- .compute_initial_parameters(cell_features = res$cell_features,
                                      cell_lineage = res$cell_lineage,
                                      lineage_future_count = res$lineage_future_count,
                                      lambda_min = 0, lambda_max = 1e6,
                                      multipler = 100)

  expect_false(isTRUE(all.equal(low$lambda_initial, high$lambda_initial)))
})

##############################################################################
## Test lineage_imputation_sequence (test-plan 2026-08-06)
##############################################################################


.raw_features <- function(res){
  res$cell_features[, setdiff(colnames(res$cell_features), "Intercept"),
                    drop = FALSE]
}

## Plan item C1. `cyfer_finalize()` resolves a tied CV curve with `which.min()`,
## so the *direction* of `lambda_sequence` decides which lambda a degenerate fit
## reports -- and "selected lambda == the ceiling I set" is documented in
## CLAUDE_kevin.md as the signature of a broken fit. That contract has never
## been pinned, so pin it here rather than leaving it implicit in the caller.
test_that("lineage_imputation_sequence runs lambda_initial down to zero", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  grid <- expand.grid(lambda_initial = c(1, 10, 100),
                      lambda_sequence_length = c(2, 5, 12))

  for(i in seq_len(nrow(grid))){
    lambda_initial <- grid$lambda_initial[i]
    lambda_sequence_length <- grid$lambda_sequence_length[i]
    label <- paste0("lambda_initial=", lambda_initial,
                    ", length=", lambda_sequence_length)

    fit <- lineage_imputation_sequence(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      lineage_future_count = res$lineage_future_count,
      lambda_initial = lambda_initial,
      lambda_sequence_length = lambda_sequence_length,
      verbose = 0
    )

    lambda_sequence <- fit$lambda_sequence
    expect_equal(length(lambda_sequence), lambda_sequence_length, info = label)
    expect_equal(length(fit$fit_list), lambda_sequence_length, info = label)
    expect_equal(lambda_sequence[1], lambda_initial, info = label)
    expect_equal(lambda_sequence[lambda_sequence_length], 0, info = label)
    expect_true(all(diff(lambda_sequence) < 0), info = label)
  }
})

## Plan item C2. Each fit is meant to warm-start from the previous one. Nothing
## observable in the return value proves it: `fit_list[[i]]` is the best of
## eleven restarts, so the warm start need not be the one that won. Mock the
## inner call and read the initialization it was handed.
test_that("lineage_imputation_sequence warm-starts from the previous fit", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lambda_sequence_length <- 5

  init_list <- list()
  fake_imputation <- function(cell_features,
                              cell_lineage,
                              coefficient_initial_list,
                              lineage_future_count,
                              lambda = 0,
                              random_initializations = 10,
                              upper_randomness = 5,
                              verbose = 1){
    init_list[[length(init_list) + 1]] <<- coefficient_initial_list

    coefficient_vec <- rep(length(init_list), ncol(cell_features) + 1)
    names(coefficient_vec) <- c("Intercept", colnames(cell_features))
    list(fit = list(coefficient_initial = coefficient_initial_list,
                    coefficient_vec = coefficient_vec,
                    convergence = 0,
                    lambda = lambda,
                    objective_val = 0))
  }

  with_mocked_bindings(
    fit <- lineage_imputation_sequence(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      lineage_future_count = res$lineage_future_count,
      lambda_initial = 10,
      lambda_sequence_length = lambda_sequence_length,
      verbose = 0
    ),
    lineage_imputation = fake_imputation,
    .package = "multiomeFate"
  )

  expect_equal(length(init_list), lambda_sequence_length)

  for(i in seq.int(2, lambda_sequence_length)){
    expect_equal(init_list[[i]],
                 fit$fit_list[[i-1]]$coefficient_vec,
                 label = paste0("warm start into fit ", i))
  }
})

## Plan item C3. With every future count zero, future_total/current_total is 0,
## so the intercept initialization is log(0) = -Inf and lambda_initial is
## 0 * Inf = NaN. Both are handed straight to `optim` without comment.
##
## NOT the same fix as B6, despite the shared symptom: B5/B6 change locals
## inside `lineage_imputation()`, whereas these two quantities are computed here
## in `.compute_initial_parameters()`. Verified by applying the B5/B6 change in
## isolation -- B6 goes green and C3 stays red.
##
## Kevin's fix is +1 smoothing throughout:
##   coefficient_initial["Intercept"] <- log((future_total+1)/(current_total+1))
##   term1 <- future_total*(1 - log((future_total+1)/(current_total+1)))
##   term2 <- sum(lineage_future_count * log1p(lineage_current_count))
test_that(".compute_initial_parameters is finite on all-zero future counts", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lineage_future_count <- res$lineage_future_count
  lineage_future_count[] <- 0
  lambda_min <- 0.01
  lambda_max <- 101

  out <- .compute_initial_parameters(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = lineage_future_count,
    lambda_min = lambda_min,
    lambda_max = lambda_max
  )

  expect_true(is.finite(out$lambda_initial))
  expect_true(all(is.finite(out$coefficient_initial)))

  # NaN slips through `min(max(x, lambda_min), lambda_max)` untouched, so the
  # clamp cannot be trusted to have bounded anything unless it is checked.
  expect_gte(out$lambda_initial, lambda_min)
  expect_lte(out$lambda_initial, lambda_max)

  # a total wipeout must start the intercept below zero -- predicting decay --
  # rather than at -Inf
  expect_lt(unname(out$coefficient_initial["Intercept"]), 0)
})

## Plan item C3, the arithmetic. The smoothing moves `lambda_initial` and the
## intercept on *every* input, not only the degenerate one, so pin the
## expressions rather than only the finiteness property -- otherwise a fix that
## smoothed just one of the three terms would still look correct here.
## Re-stating the formula in the test is deliberate: the formula is the
## decision, and this is the only place it gets checked.
test_that(".compute_initial_parameters smooths the growth ratio by one count", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  # A wide bracket is essential here: with the shipped `multipler = 1e4` the
  # heuristic saturates at `lambda_max = 101` on ordinary data, so the smoothed
  # and unsmoothed values both clamp to 101 and the comparison below would pass
  # against either implementation.
  multipler <- 1
  lambda_min <- 0
  lambda_max <- 1e12

  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = res$lineage_future_count)
  lineage_current_count <- sapply(tmp$cell_lineage_idx_list, length)
  future_total <- sum(tmp$lineage_future_count)
  current_total <- sum(lineage_current_count)
  num_lineages <- length(tmp$lineage_future_count)

  out <- .compute_initial_parameters(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    multipler = multipler
  )

  intercept_expected <- log((future_total + 1)/(current_total + 1))
  expect_equal(unname(out$coefficient_initial["Intercept"]), intercept_expected)

  # every non-intercept coefficient starts at zero
  idx_notintercept <- which(names(out$coefficient_initial) != "Intercept")
  expect_true(all(out$coefficient_initial[idx_notintercept] == 0))

  term1 <- future_total*(1 - intercept_expected)
  term2 <- sum(tmp$lineage_future_count * log1p(lineage_current_count))
  lambda_expected <- -multipler*(term1 - term2)/num_lineages
  lambda_expected <- min(max(lambda_expected, lambda_min), lambda_max)

  expect_equal(out$lambda_initial, lambda_expected)
})

## Plan item C3, the property the smoothing must not destroy. The intercept
## initialization is meant to be the log growth ratio, so it has to stay
## monotone in the future total and keep its sign at the break-even point.
test_that(".compute_initial_parameters keeps the intercept tracking growth", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lineage_names <- names(res$lineage_future_count)
  num_cells_vec <- table(res$cell_lineage)[lineage_names]

  intercept_of <- function(count_ratio){
    lineage_future_count <- stats::setNames(
      as.numeric(num_cells_vec)*count_ratio,
      lineage_names
    )
    out <- .compute_initial_parameters(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      lineage_future_count = lineage_future_count
    )
    unname(out$coefficient_initial["Intercept"])
  }

  ratio_vec <- c(0, 0.25, 0.5, 1, 2, 10)
  intercept_vec <- sapply(ratio_vec, intercept_of)

  expect_true(all(is.finite(intercept_vec)))
  expect_true(all(diff(intercept_vec) > 0))
  # shrinking clones start negative, growing clones start positive
  expect_true(all(intercept_vec[ratio_vec < 1] < 0))
  expect_true(all(intercept_vec[ratio_vec > 1] > 0))
})

## Plan item C4. The default here is 10 while the only caller passes 1e4, so the
## default has never described the shipped behavior. Kevin: set the default to
## 1e4.
test_that(".compute_initial_parameters defaults multipler to 1e4", {
  expect_equal(eval(formals(.compute_initial_parameters)$multipler), 1e4)

  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  out_default <- .compute_initial_parameters(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count
  )
  out_explicit <- .compute_initial_parameters(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count,
    multipler = 1e4
  )

  expect_equal(out_default$lambda_initial, out_explicit$lambda_initial)
})

## Plan item C5. A length-1 path is the degenerate end of the lambda sweep and
## is what breaks `cyfer_finalize()` (plan item E1); pin the producer side here.
test_that("lineage_imputation_sequence supports a length-1 lambda path", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  fit <- lineage_imputation_sequence(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    lineage_future_count = res$lineage_future_count,
    lambda_initial = 5,
    lambda_sequence_length = 1,
    verbose = 0
  )

  expect_equal(length(fit$lambda_sequence), 1)
  expect_equal(fit$lambda_sequence, 5)
  expect_equal(length(fit$fit_list), 1)
  expect_true(all(is.finite(fit$fit_list[[1]]$coefficient_vec)))
})

