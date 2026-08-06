context("Test compute_entropy (test-plan 2026-08-06)")

## Plan item H1. `R/compute_entropy.R` has no tests at all. `compute_entropy()`
## itself needs a Seurat object, but `.shannon_entropy()` is pure and carries
## the only arithmetic in the file -- including the `y[x == 0] <- 0` guard that
## exists solely to keep log2(0) = -Inf out of the sum.

test_that(".shannon_entropy is log2(k) on a uniform distribution", {
  for(k in c(2, 3, 4, 8, 16)){
    expect_equal(.shannon_entropy(rep(1, k)), log2(k),
                 label = paste0("k=", k))
    # entropy depends on the proportions, not on the scale of the counts
    expect_equal(.shannon_entropy(rep(37, k)), log2(k),
                 label = paste0("k=", k, ", rescaled"))
  }
})

test_that(".shannon_entropy is zero for a degenerate distribution", {
  # a single category, which is the early-return branch
  expect_equal(.shannon_entropy(5), 0)
  # all the mass on one of several categories
  expect_equal(.shannon_entropy(c(9, 0, 0, 0)), 0)
})

test_that(".shannon_entropy ignores zero-count categories", {
  # zero-count categories contribute nothing, so padding must not change the
  # answer and must not introduce NaN through log2(0)
  expect_equal(.shannon_entropy(c(1, 1, 0, 0)), 1)
  expect_equal(.shannon_entropy(c(3, 1)), .shannon_entropy(c(3, 1, 0, 0, 0)))
  expect_false(is.na(.shannon_entropy(c(1, 1, 0, 0))))
})

test_that(".shannon_entropy accepts the table it is called with", {
  # `compute_entropy()` only ever passes a `table`, never a bare vector
  tab_vec <- table(c("a", "a", "b", "b"))
  expect_equal(.shannon_entropy(tab_vec), 1)

  tab_vec <- table(c("a", "a", "a", "b"))
  expect_equal(.shannon_entropy(tab_vec),
               -(0.75*log2(0.75) + 0.25*log2(0.25)))
})

test_that(".shannon_entropy is maximized by the uniform distribution", {
  set.seed(10)
  trials <- 100
  k <- 5

  bool_vec <- sapply(seq_len(trials), function(trial){
    set.seed(trial)
    count_vec <- stats::rpois(k, lambda = 10) + 1
    .shannon_entropy(count_vec) <= log2(k) + 1e-9
  })

  expect_true(all(bool_vec))
})
