context("Test simulation attach-future")

## .adjust_coefficient_intercept is correct

test_that(".adjust_coefficient_intercept works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  coefficient_intercept <- 0
  embedding_coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + coefficient_intercept))
  
  simulation_res <- generate_simulation(
    embedding_mat = previous_cell_embedding_mat,
    coefficient_intercept = coefficient_intercept,
    embedding_coefficient_vec = embedding_coefficient_vec,
    lineage_spread = 1,
    lineage_prior = lineage_prior,
    num_lineages = K
  )
  
  future_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  num_future_cells <- nrow(future_cell_embedding_mat)
  
  res <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  
  stopifnot(length(res) == 1)
  
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + res)
  expect_true(sum(cell_contribution) >= num_future_cells)
})

test_that(".adjust_coefficient_intercept ensures enough future cells", {
  trials <- 1000
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    n <- 1000
    cell_contribution <- stats::rnorm(n, mean = 1)
    coefficient_intercept <- stats::rnorm(1)
    num_future_cells <- round(stats::runif(1, min = 100, max = 2000))
    
    res <- .adjust_coefficient_intercept(
      cell_contribution = cell_contribution,
      coefficient_intercept = coefficient_intercept,
      num_future_cells = num_future_cells
    )
    
    new_sum <- sum(round(exp(cell_contribution)*exp(-coefficient_intercept)*exp(res)))
    new_sum >= num_future_cells
  })
  
  expect_true(all(bool_vec))
})

############

## .pushforward_func_constructor is correct

test_that(".pushforward_func_constructor works", {
  a <- 1.5
  b <- 1:5
  res <- .pushforward_func_constructor(a = a,
                                       b = b)
  vec <- rep(2, 5)
  new_vec <- res(vec)
  expect_true(sum(abs((a*vec + b) - new_vec)) <= 1e-6)
})

## .compute_pushforward_fit is correct

test_that(".compute_pushforward_fit works", {
  set.seed(10)
  n <- 100; d <- 5
  future_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  previous_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  
  res <- .compute_pushforward_fit(
    future_mat = future_mat,
    previous_mat = previous_mat
  )
  
  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == sort(c("a", "b", "pushforward_func"))))
  
  prediction1 <- t(sapply(1:n, function(i){
    res$pushforward_func(previous_mat[i,])
  }))
  prediction2 <- t(sapply(1:n, function(i){
    res$a * previous_mat[i,] + res$b
  }))
  
  expect_true(sum(abs(prediction1 - prediction2)) <= 1e-6)
})

test_that(".compute_pushforward_fit predicts zero", {
  set.seed(10)
  n <- 100; d <- 5
  previous_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_mat <- matrix(0, nrow = n, ncol = d)
  
  res <- .compute_pushforward_fit(
    future_mat = future_mat,
    previous_mat = previous_mat
  )
  
  expect_true(abs(res$a) <= 1e-6)
  expect_true(sum(abs(res$b)) <= 1e-6)
})

test_that(".compute_pushforward_fit predicts itself", {
  set.seed(10)
  n <- 100; d <- 5
  previous_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  
  res <- .compute_pushforward_fit(
    future_mat = previous_mat,
    previous_mat = previous_mat
  )
  
  expect_true(abs(res$a - 1) <= 1e-6)
  expect_true(sum(abs(res$b)) <= 1e-6)
})

test_that(".compute_pushforward_fit is reasonable", {
  set.seed(10)
  n <- 100; d <- 5
  previous_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_mat1 <- sweep(previous_mat,
                       MARGIN = 2,
                       STATS = 1:5,
                       FUN = "+")
  future_mat2 <- sweep(10*future_mat1, 
                       MARGIN = 2,
                       STATS = 11:15,
                       FUN = "+")
  
  res1 <- .compute_pushforward_fit(
    future_mat = future_mat1,
    previous_mat = previous_mat
  )
  res2 <- .compute_pushforward_fit(
    future_mat = future_mat2,
    previous_mat = previous_mat
  )
  
  prediction1 <- t(sapply(1:n, function(i){
    res1$pushforward_func(previous_mat[i,])
  }))
  prediction2 <- t(sapply(1:n, function(i){
    res2$pushforward_func(previous_mat[i,])
  }))
  
  one_to_one <- sum(sapply(1:n, function(i){
    .l2norm(prediction1[i,] - future_mat1[i,])^2
  }))
  one_to_two <- sum(sapply(1:n, function(i){
    .l2norm(prediction1[i,] - future_mat2[i,])^2
  }))
  two_to_one <- sum(sapply(1:n, function(i){
    .l2norm(prediction2[i,] - future_mat1[i,])^2
  }))
  two_to_two <- sum(sapply(1:n, function(i){
    .l2norm(prediction2[i,] - future_mat2[i,])^2
  }))
  
  expect_true(one_to_one < one_to_two)
  expect_true(one_to_one < two_to_one)
  expect_true(two_to_two < one_to_two)
  expect_true(two_to_two < two_to_one)
})

#####

## .compute_pushforward is correct

test_that(".compute_pushforward works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  coefficient_intercept <- 0
  embedding_coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + new_coefficient_intercept)
  cell_contribution_rounded <- round(cell_contribution)
  
  res <- .compute_pushforward(
    cell_contribution = cell_contribution_rounded,
    future_cell_embedding_mat = future_cell_embedding_mat,
    num_pushforward_training_iter = 10,
    num_subsamples = 50,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    verbose = 0
  )
  
  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == sort(c("a", "b", "pushforward_func"))))
})

##################

## .compute_previous_to_future_mapping is correct

test_that(".compute_previous_to_future_mapping works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_cell_embedding_mat <- matrix(stats::rnorm(2*n*d), nrow = 2*n, ncol = d)
  coefficient_intercept <- 0
  embedding_coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + new_coefficient_intercept)
  cell_contribution_rounded <- round(cell_contribution)
  
  pushforward_res <- .compute_pushforward(
    cell_contribution = cell_contribution_rounded,
    future_cell_embedding_mat = future_cell_embedding_mat,
    num_pushforward_training_iter = 10,
    num_subsamples = 50,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    verbose = 0
  )
  
  sd_vec <- apply(future_cell_embedding_mat, 2, stats::sd)
  res <- .compute_previous_to_future_mapping(
    future_cell_embedding_mat = future_cell_embedding_mat,
    lineage_spread = 1,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    pushforward_func = pushforward_res$pushforward_func,
    sd_vec = sd_vec
  )
  
  expect_true(is.matrix(res))
  expect_true(nrow(res) == nrow(previous_cell_embedding_mat))
  expect_true(ncol(res) == nrow(future_cell_embedding_mat))
})

###################

## .dmvnorm_log_many_samples is correct

test_that(".dmvnorm_log_many_samples works", {
  set.seed(10)
  p <- 4; n <- 20
  mean_vec <- rep(0, p)
  sigma <- diag(p)
  x_mat <- MASS::mvrnorm(n, mu = mean_vec, Sigma = sigma)
  res1 <- .dmvnorm_log_many_samples(mean = mean_vec,
                                sigma = sigma,
                                x_mat = x_mat)
  
  res2 <- sapply(1:n, function(i){
    .dmvnorm(x = x_mat[i,],
             mean = mean_vec,
             sigma = sigma,
             log = TRUE)
  })
  
  expect_true(sum(abs(res1 - res2)) <= 1e-5)
})

test_that(".dmvnorm_log_many_samples is correct", {
  trials <- 200
  p <- 20
  n <- 200
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    mean_vec <- stats::runif(p)
    tmp <- matrix(stats::runif(p^2), nrow = p)
    tmp <- tmp + t(tmp)
    sigma <- crossprod(tmp)
    x_mat <- MASS::mvrnorm(n, mu = mean_vec, Sigma = sigma)
    res1 <- .dmvnorm_log_many_samples(mean = mean_vec,
                                      sigma = sigma,
                                      x_mat = x_mat)
    
    res2 <- sapply(1:n, function(i){
      .dmvnorm(x = x_mat[i,],
               mean = mean_vec,
               sigma = sigma,
               log = TRUE)
    })
    
    sum(abs(res1 - res2)) <= 1e-4
  })
 
  expect_true(all(bool_vec))
})

###################

## .assign_future_to_previous is correct

test_that(".assign_future_to_previous works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_cell_embedding_mat <- matrix(stats::rnorm(2*n*d), nrow = 2*n, ncol = d)
  rownames(previous_cell_embedding_mat) <- paste0("prev:", 1:nrow(previous_cell_embedding_mat))
  rownames(future_cell_embedding_mat) <- paste0("fut:", 1:nrow(future_cell_embedding_mat))
  
  coefficient_intercept <- 0
  embedding_coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- ceiling(exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) + new_coefficient_intercept)
  names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  cell_contribution_rounded <- ceiling(cell_contribution)
  
  pushforward_res <- .compute_pushforward(
    cell_contribution = cell_contribution_rounded,
    future_cell_embedding_mat = future_cell_embedding_mat,
    num_pushforward_training_iter = 10,
    num_subsamples = 50,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    verbose = 0
  )
  
  sd_vec <- apply(future_cell_embedding_mat, 2, stats::sd)
  mapping_mat <- .compute_previous_to_future_mapping(
    future_cell_embedding_mat = future_cell_embedding_mat,
    lineage_spread = 1,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    pushforward_func = pushforward_res$pushforward_func,
    sd_vec = sd_vec
  )
  
  res <- .assign_future_to_previous(
    mapping_mat = mapping_mat,
    previous_cell_contribution = cell_contribution_rounded
  )
  
  expect_true(length(res$prev_cell_num_progenitor) == length(cell_contribution_rounded))
  expect_true(all(cell_contribution_rounded >= res$prev_cell_num_progenitor))
  expect_true(length(res$future_cell_assignment) == num_future_cells)
  expect_true(all(res$future_cell_assignment %in% rownames(previous_cell_embedding_mat)))
  expect_true(all(names(res$future_cell_assignment) %in% rownames(future_cell_embedding_mat)))
  
  tab_vec <- table(res$future_cell_assignment)
  vec <- rep(0, length(res$prev_cell_num_progenitor))
  names(vec) <- names(res$prev_cell_num_progenitor)
  vec[names(tab_vec)] <- tab_vec
  expect_true(sum(abs(vec - res$prev_cell_num_progenitor)) <= 1e-6)
})

##################

## generate_simulation_attachFuture is correct

test_that("generate_simulation_attachFuture works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  future_cell_embedding_mat <- matrix(stats::rnorm(2*n*d), nrow = 2*n, ncol = d)
  rownames(previous_cell_embedding_mat) <- paste0("prev:", 1:nrow(previous_cell_embedding_mat))
  rownames(future_cell_embedding_mat) <- paste0("fut:", 1:nrow(future_cell_embedding_mat))
  
  coefficient_intercept <- 0
  embedding_coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  lineage_assignment <- factor(sample(paste0("lineage:", 1:K), size = n, replace = TRUE))
  names(lineage_assignment) <- rownames(previous_cell_embedding_mat)
  
  res <- generate_simulation_attachFuture(coefficient_intercept = coefficient_intercept,
                                          embedding_coefficient_vec = embedding_coefficient_vec,
                                          future_cell_embedding_mat = future_cell_embedding_mat,
                                          lineage_assignment = lineage_assignment,
                                          previous_cell_embedding_mat = previous_cell_embedding_mat,
                                          verbose = 0)
  
  expect_true(is.list(res))
  
  cell_contribution <- round(exp(res$coefficient_intercept + as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec)))
  names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  tmp <- table(res$future_cell_assignment)
  res_tabulate <- rep(0, length = length(cell_contribution))
  names(res_tabulate) <- names(cell_contribution)
  res_tabulate[names(tmp)] <- tmp
  
  expect_true(all(res_tabulate <= cell_contribution))
  
  # non-trivial differences in lineages
  tmp <- apply(res$mapping_mat, 1, function(x){diff(range(x))})
  expect_true(any(tmp > 0))
})

##############################################################################
## generate_simulation_attachFuture() -- top level (test-plan-full 2026-08-06)
##############################################################################

.attachFuture_fixture <- function(d = 5,
                                  num_lineages = 5,
                                  num_previous = 100,
                                  seed_number = 10){
  set.seed(seed_number)
  previous_cell_embedding_mat <- matrix(stats::rnorm(num_previous * d),
                                        nrow = num_previous,
                                        ncol = d)
  future_cell_embedding_mat <- matrix(stats::rnorm(2 * num_previous * d),
                                      nrow = 2 * num_previous,
                                      ncol = d)
  rownames(previous_cell_embedding_mat) <- paste0("prev:",
                                                  seq_len(num_previous))
  rownames(future_cell_embedding_mat) <- paste0("fut:",
                                                seq_len(2 * num_previous))
  lineage_assignment <- factor(sample(paste0("lineage:",
                                             seq_len(num_lineages)),
                                      size = num_previous, replace = TRUE))
  names(lineage_assignment) <- rownames(previous_cell_embedding_mat)

  list(embedding_coefficient_vec = rep(1, d),
       future_cell_embedding_mat = future_cell_embedding_mat,
       lineage_assignment = lineage_assignment,
       previous_cell_embedding_mat = previous_cell_embedding_mat)
}

.run_attachFuture <- function(fixture, ...){
  generate_simulation_attachFuture(
    coefficient_intercept = 0,
    embedding_coefficient_vec = fixture$embedding_coefficient_vec,
    future_cell_embedding_mat = fixture$future_cell_embedding_mat,
    lineage_assignment = fixture$lineage_assignment,
    previous_cell_embedding_mat = fixture$previous_cell_embedding_mat,
    verbose = 0,
    ...)
}

## AF1
test_that("num_pushforward_training_iter changes the fit (AF1)", {
  fixture <- .attachFuture_fixture()

  set.seed(10)
  res_one <- .run_attachFuture(fixture, num_pushforward_training_iter = 1)
  set.seed(10)
  res_many <- .run_attachFuture(fixture, num_pushforward_training_iter = 20)

  expect_true(!identical(res_one$mapping_mat, res_many$mapping_mat))
})
## AF2: conservation. Every future cell is handed to exactly one parent, and
## the per-parent tally must add back up to the number of future cells.
test_that("every future cell is assigned to exactly one parent (AF2)", {
  fixture <- .attachFuture_fixture()
  set.seed(10)
  res <- .run_attachFuture(fixture)

  num_future <- nrow(fixture$future_cell_embedding_mat)
  expect_true(length(res$future_cell_assignment) == num_future)
  expect_true(!any(is.na(res$future_cell_assignment)))
  # Named by future cell, but *not* in input order: the mapping matrix's
  # columns are reordered by column sum inside
  # `.compute_previous_to_future_mapping()`, and that order carries through.
  expect_true(setequal(names(res$future_cell_assignment),
                       rownames(fixture$future_cell_embedding_mat)))
  expect_true(anyDuplicated(names(res$future_cell_assignment)) == 0)
  expect_true(!identical(names(res$future_cell_assignment),
                         rownames(fixture$future_cell_embedding_mat)))
  expect_true(all(res$future_cell_assignment %in%
                    rownames(fixture$previous_cell_embedding_mat)))
  expect_true(sum(res$prev_cell_num_progenitor) == num_future)

  # And the tally really is the tabulation of the assignment.
  tally_vec <- table(res$future_cell_assignment)
  expect_true(all(res$prev_cell_num_progenitor[names(tally_vec)] ==
                    as.numeric(tally_vec)))
})

## AF3: the quota is what keeps the simulated expansion consistent with the fate
## potentials that generated it -- a parent may not produce more progeny than
## its rounded contribution allows.
test_that("no parent exceeds its rounded contribution (AF3)", {
  fixture <- .attachFuture_fixture()
  set.seed(10)
  res <- .run_attachFuture(fixture)

  contribution_vec <- round(exp(
    res$coefficient_intercept +
      as.numeric(fixture$previous_cell_embedding_mat %*%
                   fixture$embedding_coefficient_vec)))
  names(contribution_vec) <- rownames(fixture$previous_cell_embedding_mat)

  shared_names <- names(res$prev_cell_num_progenitor)
  expect_true(all(res$prev_cell_num_progenitor <=
                    contribution_vec[shared_names]))
})

## AF4
test_that("future_lineage_size matches the recomputed assignment (AF4)", {
  fixture <- .attachFuture_fixture()
  set.seed(10)
  res <- .run_attachFuture(fixture)

  lineage_vec <- fixture$lineage_assignment[
    names(res$prev_cell_num_progenitor)]
  expected_vec <- sapply(levels(droplevels(lineage_vec)), function(lev){
    sum(res$prev_cell_num_progenitor[which(lineage_vec == lev)])
  })

  expect_true(all(names(res$future_lineage_size) %in% names(expected_vec)) ||
                all(names(expected_vec) %in% names(res$future_lineage_size)))
  shared_names <- intersect(names(expected_vec), names(res$future_lineage_size))
  expect_true(length(shared_names) > 0)
  expect_true(all(res$future_lineage_size[shared_names] ==
                    expected_vec[shared_names]))
  expect_true(sum(res$future_lineage_size) ==
                nrow(fixture$future_cell_embedding_mat))
})

## AF5: the loop exists because rounding can push the total below the target
## even after the exact shift, so the guarantee is on the *rounded* sum.
test_that(".adjust_coefficient_intercept guarantees the rounded total (AF5)", {
  set.seed(10)
  for(num_future_cells in c(50, 200, 1000)){
    label <- paste0("num_future_cells = ", num_future_cells)
    cell_contribution <- exp(stats::rnorm(100))

    new_intercept <- .adjust_coefficient_intercept(
      cell_contribution = cell_contribution,
      coefficient_intercept = 0,
      num_future_cells = num_future_cells)

    adjusted_vec <- cell_contribution * exp(new_intercept)
    expect_true(sum(round(adjusted_vec)) >= num_future_cells, info = label)
  }

  # An input that cannot converge inside max_iter must error rather than
  # silently return an intercept that does not meet the guarantee.
  expect_error(.adjust_coefficient_intercept(
    cell_contribution = exp(stats::rnorm(100)),
    coefficient_intercept = 0,
    num_future_cells = 1e6,
    interval_add = 1e-8,
    max_iter = 5))
})

## AF6: `mapping_mat` comes back scaled by 1e3 and rounded, so a caller who
## forgets that gets numbers that look like counts but are per-mille weights.
test_that("mapping_mat is returned scaled by 1e3 and rounded (AF6)", {
  fixture <- .attachFuture_fixture()
  set.seed(10)
  res <- .run_attachFuture(fixture)

  expect_true(all(res$mapping_mat == round(res$mapping_mat)))
  expect_true(max(res$mapping_mat) > 1)
  # Columns were probability vectors before scaling, so they now sum to ~1e3.
  column_sum_vec <- colSums(res$mapping_mat)
  expect_true(max(abs(column_sum_vec - 1e3)) <= 5)
})

## AF7: two independent implementations of the same density live in one package
## -- `.dmvnorm()` in R/simulation.R and `.dmvnorm_log_many_samples()` here.
## That is exactly the situation where a cross-check earns its place.
test_that(".dmvnorm_log_many_samples agrees with .dmvnorm (AF7)", {
  set.seed(10)
  for(d in c(2, 3, 5)){
    label <- paste0("d = ", d)
    x_mat <- matrix(stats::rnorm(40 * d), nrow = 40, ncol = d)
    mean_vec <- stats::rnorm(d)
    sigma_mat <- diag(stats::runif(d, min = 0.5, max = 2))

    res_many <- .dmvnorm_log_many_samples(x_mat = x_mat,
                                          mean = mean_vec,
                                          sigma = sigma_mat)
    res_single <- .dmvnorm(x = x_mat,
                           mean = mean_vec,
                           sigma = sigma_mat,
                           log = TRUE)

    expect_true(length(res_many) == nrow(x_mat), info = label)
    expect_true(max(abs(res_many - res_single)) <= 1e-8, info = label)
  }
})

## AF8: the kernel covariance must be built from *variances*
## (`lineage_spread * diag(sd_vec^2)`), matching the convention
## `.form_gaussian_distribution()` uses in R/simulation.R. It previously used
## unsquared standard deviations, which made every parent-child kernel too tight
## and gave `lineage_spread` a different meaning in the two files. Assert the
## algebra directly against `.dmvnorm()` rather than the shape of the result,
## since an unsquared diagonal also produces a perfectly plausible matrix.
test_that(".compute_previous_to_future_mapping uses variances (AF8)", {
  set.seed(10)
  d <- 3
  previous_mat <- matrix(stats::rnorm(6 * d), nrow = 6, ncol = d,
                         dimnames = list(paste0("prev:", seq_len(6)), NULL))
  future_mat <- matrix(stats::rnorm(8 * d), nrow = 8, ncol = d,
                       dimnames = list(paste0("fut:", seq_len(8)), NULL))
  # sd of 2 is the case that separates the two forms: 2 versus 4 on the
  # diagonal.
  sd_vec <- rep(2, d)
  lineage_spread <- 1.5
  pushforward_func <- .pushforward_func_constructor(a = 1, b = rep(0, d))

  res <- .compute_previous_to_future_mapping(
    future_cell_embedding_mat = future_mat,
    lineage_spread = lineage_spread,
    previous_cell_embedding_mat = previous_mat,
    pushforward_func = pushforward_func,
    sd_vec = sd_vec)

  # Rebuild the same matrix by hand with the variance form.
  expected_mat <- sapply(seq_len(nrow(previous_mat)), function(i){
    .dmvnorm(x = future_mat,
             mean = pushforward_func(previous_mat[i, ]),
             sigma = lineage_spread * diag(sd_vec^2),
             log = TRUE)
  })
  expected_mat <- t(expected_mat)
  expected_mat <- exp(expected_mat - max(expected_mat))
  expected_mat <- sweep(expected_mat, MARGIN = 2,
                        STATS = colSums(expected_mat), FUN = "/")
  expected_mat <- expected_mat[, order(colSums(expected_mat),
                                       decreasing = TRUE), drop = FALSE]

  expect_true(all(dim(res) == dim(expected_mat)))
  expect_true(max(abs(sort(as.numeric(res)) -
                        sort(as.numeric(expected_mat)))) <= 1e-8)

  # And the unsquared form gives a materially different answer, so this is not
  # a test that would pass either way.
  wrong_mat <- sapply(seq_len(nrow(previous_mat)), function(i){
    .dmvnorm(x = future_mat,
             mean = pushforward_func(previous_mat[i, ]),
             sigma = lineage_spread * diag(sd_vec),
             log = TRUE)
  })
  wrong_mat <- exp(t(wrong_mat) - max(wrong_mat))
  wrong_mat <- sweep(wrong_mat, MARGIN = 2, STATS = colSums(wrong_mat),
                     FUN = "/")
  expect_true(max(abs(sort(as.numeric(res)) -
                        sort(as.numeric(wrong_mat)))) > 1e-4)
})

## AF9
test_that(".compute_pushforward has no unused previous_cell_potential (AF9)", {
  expect_true(!("previous_cell_potential" %in%
                  names(formals(.compute_pushforward))))
})