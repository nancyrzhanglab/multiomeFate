context("Test plastic simulations")

## .reorder_by_contribution is correct

test_that(".reorder_by_contribution works", {
  set.seed(10)
  n <- 101
  vec <- abs(stats::rnorm(n))
  res <- .reorder_by_contribution(vec)
  
  expect_true(length(res) == n)
  expect_true(any(!duplicated(res)))
  expect_true(all(res > 0))
  expect_true(all(res <= n))
  expect_true(all(res %% 1 == 0))
  
  expect_true(all(diff(vec[res[seq(1, n/2, by = 2)]]) >= 0))
  expect_true(all(diff(vec[res[seq(2, n/2, by = 2)]]) <= 0))
})

############

## .compute_plastic_probabilities is correct

test_that(".compute_plastic_probabilities works", {
  set.seed(10)
  n <- 101
  cell_contribution_truth <- stats::rnorm(n)
  names(cell_contribution_truth) <- paste0("cell:", 1:n)
  num_lineages <- 10
  
  res <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    rho = NA
  )
  
  # image(t(res$prob_mat))
  
  expect_true(all(abs(rowSums(res$prob_mat) - 1) <= 1e-6))
  expect_true(all(dim(res$prob_mat) == c(n, num_lineages)))
})

############

## .assign_plastic_lineages is correct

test_that(".assign_plastic_lineages works", {
  set.seed(10)
  n <- 101
  cell_contribution_truth <- stats::rnorm(n)
  names(cell_contribution_truth) <- paste0("cell:", 1:n)
  num_lineages <- 10
  
  tmp <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    rho = NA
  )
  prob_mat <- tmp$prob_mat
  
  lineage_assignment <- .assign_plastic_lineages(prob_mat)
  
  expect_true(all(!is.na(lineage_assignment)))
  expect_true(length(lineage_assignment) == n)
  expect_true(is.factor(lineage_assignment))
  expect_true(all(table(lineage_assignment) <= ceiling(n/num_lineages)))
})

##########################

## generate_simulation_plastic is correct

test_that("generate_simulation_plastic works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)

  res <- generate_simulation_plastic(embedding_mat)
  
  expect_true(is.list(res))
  # plot(res$cell_fate_potential, res$lineage_future_size[res$lineage_assignment])
  # mean_val = sapply(levels(res$lineage_assignment), function(lineage){mean(res$cell_fate_potential[res$lineage_assignment == lineage])}); plot(mean_val, res$lineage_future_size)
})
