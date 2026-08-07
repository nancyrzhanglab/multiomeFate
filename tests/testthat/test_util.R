context("Testing utils")

## This file previously covered `.log_sum_exp()`, `.exp_ratio()`,
## `.mult_vec_mat()` and `.mult_mat_vec()`. All four were dead code -- called
## from nowhere in `R/`, only from this file -- and were deleted in 1.0.2.001
## (defect D5). Their tests went with them. `.nonzero_col()` is the only
## survivor in `R/util.R`.

## UT5: the two modes of `.nonzero_col()` return element-aligned vectors, so
## zipping them reconstructs the stored column. The function reads the
## compressed-column slots directly rather than subsetting, so a mismatch
## between the `@i` and `@x` reads would be invisible in either mode alone.

test_that(".nonzero_col returns aligned indices and values (UT5)", {
  set.seed(10)
  dense_mat <- matrix(stats::rbinom(60, size = 1, prob = 0.4) * stats::runif(60),
                      nrow = 10,
                      ncol = 6)
  # A column with no structural non-zeros at all, to hit the early return.
  dense_mat[, 3] <- 0
  sparse_mat <- Matrix::Matrix(dense_mat, sparse = TRUE)

  for(col_idx in seq_len(ncol(sparse_mat))){
    label <- paste0("column ", col_idx)
    idx_vec <- .nonzero_col(sparse_mat, col_idx = col_idx, bool_value = FALSE)
    value_vec <- .nonzero_col(sparse_mat, col_idx = col_idx, bool_value = TRUE)

    expect_true(length(idx_vec) == length(value_vec), info = label)

    reconstructed_vec <- rep(0, nrow(sparse_mat))
    reconstructed_vec[idx_vec] <- value_vec
    expect_true(sum(abs(reconstructed_vec - dense_mat[, col_idx])) <= 1e-10,
                info = label)
  }
})

test_that(".nonzero_col returns numeric(0) on an all-zero column (UT5)", {
  sparse_mat <- Matrix::Matrix(matrix(0, nrow = 4, ncol = 3), sparse = TRUE)

  expect_true(length(.nonzero_col(sparse_mat, col_idx = 2,
                                  bool_value = FALSE)) == 0)
  expect_true(length(.nonzero_col(sparse_mat, col_idx = 2,
                                  bool_value = TRUE)) == 0)
})

test_that(".nonzero_col rejects a dense matrix and an out-of-range column", {
  dense_mat <- matrix(1:12, nrow = 4, ncol = 3)
  sparse_mat <- Matrix::Matrix(dense_mat, sparse = TRUE)

  expect_error(.nonzero_col(dense_mat, col_idx = 1, bool_value = TRUE))
  expect_error(.nonzero_col(sparse_mat, col_idx = 4, bool_value = TRUE))
  expect_error(.nonzero_col(sparse_mat, col_idx = 0, bool_value = TRUE))
})

##############################################################################
## Estimation internals (test-plan-full 2026-08-06, section 7)
##############################################################################

## UT6: `.lineage_cleanup()` is the single entry point through which every
## estimation function normalizes its inputs, and when `lineage_future_count`
## and `cell_lineage` disagree it **drops cells silently** at the default
## `verbose = 0`. This is the most surprising behaviour in the estimation core:
## a caller who mistypes one lineage name loses every cell of it and gets a
## perfectly ordinary-looking fit back. Assert both the drop and the silence,
## so the silence is a documented choice rather than an oversight.
test_that(".lineage_cleanup drops mismatched cells silently (UT6)", {
  set.seed(10)
  cell_features <- matrix(stats::rnorm(20 * 3), nrow = 20, ncol = 3,
                          dimnames = list(paste0("cell", seq_len(20)),
                                          paste0("feat", seq_len(3))))
  cell_lineage <- rep(paste0("L", seq_len(4)), each = 5)
  names(cell_lineage) <- rownames(cell_features)
  # L4 is absent from the counts; L9 has counts but no cells.
  lineage_future_count <- stats::setNames(c(10, 20, 30, 40),
                                          c("L1", "L2", "L3", "L9"))

  expect_silent(
    res <- .lineage_cleanup(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            lineage_future_count = lineage_future_count))

  # The five L4 cells are gone, without a word.
  expect_true(nrow(res$cell_features) == 15)
  expect_true(length(res$cell_lineage) == 15)
  expect_true(!("L4" %in% res$cell_lineage))
  expect_true(all(res$uniq_lineages == c("L1", "L2", "L3")))
  expect_true(!("L9" %in% names(res$lineage_future_count)))

  # At verbose > 0 it does warn, so the information exists -- it is just off by
  # default, and every caller in the package uses the default.
  expect_warning(.lineage_cleanup(cell_features = cell_features,
                                  cell_lineage = cell_lineage,
                                  lineage_future_count = lineage_future_count,
                                  verbose = 1))
})

## UT7: every returned object must agree on lineage ordering, which is
## `sort()` order. `.lineage_objective()` and `.lineage_gradient()` index
## different ones of these (see UT8), so a divergence here would misalign the
## objective from its own gradient.
test_that(".lineage_cleanup returns everything in sort order (UT7)", {
  set.seed(10)
  cell_features <- matrix(stats::rnorm(30 * 3), nrow = 30, ncol = 3,
                          dimnames = list(paste0("cell", seq_len(30)),
                                          paste0("feat", seq_len(3))))
  # Deliberately unsorted lineage names, including ones whose sort order differs
  # from their first appearance.
  cell_lineage <- rep(c("zeta", "alpha", "Mu", "beta", "gamma", "Delta"),
                      each = 5)
  names(cell_lineage) <- rownames(cell_features)
  lineage_future_count <- stats::setNames(c(1, 2, 3, 4, 5, 6),
                                          c("gamma", "zeta", "Delta", "Mu",
                                            "beta", "alpha"))

  res <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)

  expect_true(identical(res$uniq_lineages, sort(unique(cell_lineage))))
  expect_true(identical(names(res$lineage_future_count), res$uniq_lineages))
  expect_true(identical(names(res$cell_lineage_idx_list), res$uniq_lineages))

  # And the index list really does point at the cells of that lineage.
  for(lineage in res$uniq_lineages){
    idx_vec <- res$cell_lineage_idx_list[[lineage]]
    expect_true(all(res$cell_lineage[idx_vec] == lineage), info = lineage)
    expect_true(length(idx_vec) == sum(res$cell_lineage == lineage),
                info = lineage)
  }
})

## UT8: the objective walks lineages by `names(lineage_future_count)` while the
## gradient walks them by `names(cell_lineage_idx_list)`. `.lineage_cleanup()`
## makes the two identical, so this is fine today -- but a caller who builds the
## arguments by hand can hand the pair a mismatched ordering and get an
## objective and a gradient describing different models, with no error. Assert
## the identity so that arrangement fails loudly instead.
test_that(".lineage_objective and .lineage_gradient share an ordering (UT8)", {
  set.seed(10)
  cell_features <- matrix(stats::rnorm(30 * 2), nrow = 30, ncol = 2,
                          dimnames = list(paste0("cell", seq_len(30)),
                                          paste0("feat", seq_len(2))))
  cell_lineage <- rep(c("c", "a", "b"), each = 10)
  names(cell_lineage) <- rownames(cell_features)
  lineage_future_count <- stats::setNames(c(12, 7, 20), c("b", "c", "a"))

  res <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)

  expect_true(identical(names(res$lineage_future_count),
                        names(res$cell_lineage_idx_list)))
})

## UT9: the function returns the **negative** penalized log-likelihood -- lower
## is better -- which is what the 1.0.2.001 rename from
## `evaluate_loglikelihood()` to `evaluate_nll()` records. Assert the sign
## against a hand computation, because a value that is merely "a number that
## goes down as the fit improves" would pass any invariant test under either
## sign convention.
test_that("evaluate_nll returns the negative log-likelihood (UT9)", {
  set.seed(10)
  cell_features <- matrix(stats::rnorm(20 * 2, sd = 0.3), nrow = 20, ncol = 2,
                          dimnames = list(paste0("cell", seq_len(20)),
                                          paste0("feat", seq_len(2))))
  cell_lineage <- rep(c("a", "b"), each = 10)
  names(cell_lineage) <- rownames(cell_features)
  lineage_future_count <- stats::setNames(c(15, 25), c("a", "b"))
  coefficient_vec <- stats::setNames(c(0.5, 0.2, -0.1),
                                     c("Intercept", "feat1", "feat2"))

  res_val <- evaluate_nll(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          coefficient_vec = coefficient_vec,
                          lineage_future_count = lineage_future_count,
                          lambda = 0)

  # Hand computation of the same objective: per lineage, the summed rate minus
  # the observed count times its log, averaged over lineages.
  design_mat <- cbind(Intercept = 1, cell_features)
  rate_vec <- exp(as.numeric(design_mat %*% coefficient_vec))
  expected_val <- mean(sapply(c("a", "b"), function(lineage){
    idx <- which(cell_lineage == lineage)
    total_val <- sum(rate_vec[idx])
    total_val - lineage_future_count[lineage] * log(total_val)
  }))

  expect_true(abs(res_val - expected_val) <= 1e-8)

  ## Lower really is better. Generate counts from a known coefficient vector and
  ## check the truth scores below a deliberately wrong vector -- if the sign
  ## were flipped this comparison would reverse, whereas any "the number goes
  ## down as the fit improves" invariant would pass under either convention.
  truth_vec <- stats::setNames(c(1, 0.8, -0.4),
                               c("Intercept", "feat1", "feat2"))
  truth_rate_vec <- exp(as.numeric(design_mat %*% truth_vec))
  exact_count_vec <- stats::setNames(
    sapply(c("a", "b"), function(lineage){
      sum(truth_rate_vec[which(cell_lineage == lineage)])
    }), c("a", "b"))

  wrong_vec <- stats::setNames(c(-2, -1.5, 2),
                               c("Intercept", "feat1", "feat2"))
  score_truth <- evaluate_nll(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              coefficient_vec = truth_vec,
                              lineage_future_count = exact_count_vec,
                              lambda = 0)
  score_wrong <- evaluate_nll(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              coefficient_vec = wrong_vec,
                              lineage_future_count = exact_count_vec,
                              lambda = 0)
  expect_true(score_truth < score_wrong)

  # The old name is gone, deliberately, with no deprecated alias.
  expect_true(!exists("evaluate_loglikelihood",
                      where = asNamespace("multiomeFate"),
                      inherits = FALSE))
})
