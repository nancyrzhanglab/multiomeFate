context("Test simulation attach-future")

## .adjust_coefficient_intercept is correct

test_that(".adjust_coefficient_intercept works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  previous_cell_embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  coefficient_intercept <- 0
  coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + coefficient_intercept))
  
  simulation_res <- generate_simulation(
    embedding_mat = previous_cell_embedding_mat,
    coefficient_intercept = coefficient_intercept,
    coefficient_vec = coefficient_vec,
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
  
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + res)
  expect_true(abs(sum(cell_contribution) - num_future_cells) <= 2)
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
  future_mat2 <- sweep(10*future_mat, 
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
  coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + new_coefficient_intercept)
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
  coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- round(exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + new_coefficient_intercept)
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
  coefficient_vec <- rep(1, ncol(previous_cell_embedding_mat))
  cell_contribution <- ceiling(exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + coefficient_intercept))
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(previous_cell_embedding_mat %*% coefficient_vec) + new_coefficient_intercept)
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
  
  expect_true(length(res$prev_lineage_size) == length(cell_contribution_rounded))
  expect_true(all(cell_contribution_rounded >= res$prev_lineage_size))
  expect_true(all(res$prev_lineage_assignment %in% rownames(previous_cell_embedding_mat)))
  expect_true(all(names(res$prev_lineage_assignment) %in% rownames(future_cell_embedding_mat)))
  
  tab_vec <- table(res$prev_lineage_assignment)
  vec <- rep(0, length(res$prev_lineage_size))
  names(vec) <- names(res$prev_lineage_size)
  vec[names(tab_vec)] <- tab_vec
  expect_true(sum(abs(vec - res$prev_lineage_size)) <= 1e-6)
})

