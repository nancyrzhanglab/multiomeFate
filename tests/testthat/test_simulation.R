context("Test simulation")

## .kmeans_seed is correct

test_that(".kmeans_seed works", {
  set.seed(10)
  n <- 100; d <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  K <- 5
  res <- .kmeans_seed(
    embedding_mat = embedding_mat,
    K = K
  )
  
  expect_true(length(res) == K)
})

## .form_gaussian_distribution is correct

test_that(".form_gaussian_distribution works", {
  set.seed(10)
  n <- 100; d <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  K <- 5
  cluster_idx <- .kmeans_seed(
    embedding_mat = embedding_mat,
    K = K
  )
  
  res <- .form_gaussian_distribution(
    cluster_idx = cluster_idx[1],
    embedding_mat = embedding_mat,
    rho = 1
  )
  
  expect_true(class(res) == "gaussian")
})

## .form_gaussian_distributions is correct

test_that(".form_gaussian_distribution works", {
  set.seed(10)
  n <- 100; d <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  K <- 5
  cluster_idx <- .kmeans_seed(
    embedding_mat = embedding_mat,
    K = K
  )
  
  res <- .form_gaussian_distributions(
    cluster_idx = cluster_idx,
    embedding_mat = embedding_mat,
    rho = 1
  )
  
  expect_true(is.list(res))
})

## .dmvnorm is correct

test_that(".dmvnorm works", {
  set.seed(10)
  n <- 100; d <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  K <- 5
  cluster_idx <- .kmeans_seed(
    embedding_mat = embedding_mat,
    K = K
  )
  gaussian_list <- .form_gaussian_distributions(
    cluster_idx = cluster_idx,
    embedding_mat = embedding_mat,
    rho = 1
  )
  
  res <- .dmvnorm(x = embedding_mat[1,],
                  mean = gaussian_list[[1]]$mean, 
                  sigma = gaussian_list[[1]]$cov, 
                  log = TRUE, 
                  checkSymmetry = FALSE)
  
  expect_true(is.numeric(res))
})

############

## .log_sum_exp_normalization is correct

test_that(".log_sum_exp_normalization", {
  set.seed(10)
  p <- abs(stats::rnorm(5))
  p <- p/sum(p)
  
  x <- log(10*p)
  res <- .log_sum_exp_normalization(x)
  
  expect_true(sum(abs(res - p)) <= 1e-6)
  expect_true(sum(abs(res - exp(x)/sum(exp(x)))) <= 1e-6)
})

## .compute_posteriors is correct

test_that(".compute_posteriors works", {
  set.seed(10)
  n <- 100; d <- 5
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  K <- 5
  cluster_idx <- .kmeans_seed(
    embedding_mat = embedding_mat,
    K = K
  )
  gaussian_list <- .form_gaussian_distributions(
    cluster_idx = cluster_idx,
    embedding_mat = embedding_mat,
    rho = 1
  )
  
  res <- .compute_posteriors(
    embedding_mat = embedding_mat,
    gaussian_list = gaussian_list,
    lineage_prior = rep(1/K, length = K)
  )
  
  expect_true(all(abs(rowSums(res)-1) <= 1e-6))
})

## generate_simulation is correct

test_that("generate_simulation works", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)

  res <- generate_simulation(
    embedding_mat = embedding_mat,
    coefficient_intercept = 0,
    coefficient_vec = rep(1, ncol(embedding_mat)),
    lineage_concentration = 1,
    lineage_prior = lineage_prior,
    num_lineages = K
  )
  
  expect_true(is.list(res))
})