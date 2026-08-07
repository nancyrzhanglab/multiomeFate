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
    embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
    lineage_spread = 1,
    lineage_prior = lineage_prior,
    num_lineages = K
  )
  
  expect_true(is.list(res))
})

test_that("generate_simulation can take in fate embeddings as well", {
  set.seed(10)
  n <- 100; d <- 5
  K <- 5
  lineage_prior <- rep(1/K, length = K)
  embedding_mat <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  fatefeatures_mat <- matrix(stats::rnorm(n*2), nrow = n, ncol = 2)
  fatefeatures_coefficient_vec <- c(2,-2)
  
  res <- generate_simulation(
    embedding_mat = embedding_mat,
    coefficient_intercept = 0,
    embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
    fatefeatures_coefficient_vec = fatefeatures_coefficient_vec,
    fatefeatures_mat = fatefeatures_mat, 
    lineage_spread = 1,
    lineage_prior = lineage_prior,
    num_lineages = K
  )
  
  expect_true(is.list(res))
})

##############################################################################
## generate_simulation() -- top-level behaviour (test-plan-full 2026-08-06)
##############################################################################

## A clustered embedding, so that the k-means seeding has real structure to
## find and the priming simulator can actually separate lineages by mean fate
## potential. A spherical Gaussian blob would make SM1 vacuous.
.clustered_embedding_mat <- function(d = 4,
                                     num_clusters = 5,
                                     num_per_cluster = 40,
                                     seed_number = 10){
  set.seed(seed_number)
  center_mat <- matrix(stats::rnorm(num_clusters * d, sd = 3),
                       nrow = num_clusters,
                       ncol = d)
  embedding_mat <- do.call(rbind, lapply(seq_len(num_clusters), function(k){
    matrix(stats::rnorm(num_per_cluster * d, sd = 0.5),
           nrow = num_per_cluster,
           ncol = d) + rep(center_mat[k, ], each = num_per_cluster)
  }))
  rownames(embedding_mat) <- paste0("cell:", seq_len(nrow(embedding_mat)))
  colnames(embedding_mat) <- paste0("dim", seq_len(d))

  embedding_mat
}

## SM2: pins the generative model against the roxygen and against
## `cyfer_finalize()`'s scale convention. If this drifts, every simulation-based
## characterization number in the package is measuring something else.
test_that("cell_fate_potential_truth is log10(exp(intercept + X b)) (SM2)", {
  embedding_mat <- .clustered_embedding_mat()
  d <- ncol(embedding_mat)

  param_grid <- expand.grid(coefficient_intercept = c(-1, 0, 2),
                            seed_number = 1:2)

  for(i in seq_len(nrow(param_grid))){
    label <- paste0("intercept = ", param_grid[i, "coefficient_intercept"],
                    ", seed = ", param_grid[i, "seed_number"])
    set.seed(param_grid[i, "seed_number"])
    embedding_coefficient_vec <- stats::runif(d, min = -0.5, max = 0.5)

    set.seed(param_grid[i, "seed_number"])
    res <- generate_simulation(
      embedding_mat = embedding_mat,
      coefficient_intercept = param_grid[i, "coefficient_intercept"],
      embedding_coefficient_vec = embedding_coefficient_vec,
      num_lineages = 5)

    expected_vec <- log10(exp(
      as.numeric(embedding_mat %*% embedding_coefficient_vec) +
        param_grid[i, "coefficient_intercept"]))
    expect_true(max(abs(res$cell_fate_potential_truth - expected_vec)) <= 1e-10,
                info = label)
  }
})

test_that("fatefeatures_mat enters cell_fate_potential_truth additively (SM2)", {
  embedding_mat <- .clustered_embedding_mat()
  d <- ncol(embedding_mat)
  set.seed(10)
  fatefeatures_mat <- matrix(stats::rnorm(nrow(embedding_mat) * 2),
                             nrow = nrow(embedding_mat),
                             ncol = 2)
  embedding_coefficient_vec <- rep(0.3, d)
  fatefeatures_coefficient_vec <- c(0.5, -0.2)

  set.seed(10)
  res <- generate_simulation(
    embedding_mat = embedding_mat,
    embedding_coefficient_vec = embedding_coefficient_vec,
    fatefeatures_coefficient_vec = fatefeatures_coefficient_vec,
    fatefeatures_mat = fatefeatures_mat,
    num_lineages = 5)

  expected_vec <- log10(exp(
    as.numeric(embedding_mat %*% embedding_coefficient_vec) +
      as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)))
  expect_true(max(abs(res$cell_fate_potential_truth - expected_vec)) <= 1e-10)
})

## SM3: the two fate-potential fields are deliberately on different footings --
## `cell_fate_potential` is log10(realized + 1) and `cell_fate_potential_truth`
## is log10(expected) with no `+1`. They are not two estimates of one quantity,
## and "fixing" one to match the other would silently change what every
## simulation reports. Same shape of trap as the exp()/10^ note in
## ?cyfer_finalize.
test_that("the +1 asymmetry between the two fate fields is real (SM3)", {
  embedding_mat <- .clustered_embedding_mat()

  set.seed(10)
  res <- generate_simulation(embedding_mat = embedding_mat,
                             bool_add_randomness = FALSE,
                             embedding_coefficient_vec =
                               rep(0.2, ncol(embedding_mat)),
                             num_lineages = 5)

  # Even with no Poisson draw at all, the two differ by exactly the +1.
  expect_true(max(abs(res$cell_fate_potential -
                        log10(10^res$cell_fate_potential_truth + 1))) <= 1e-10)
  expect_true(max(abs(res$cell_fate_potential -
                        res$cell_fate_potential_truth)) > 1e-3)
})

## SM4
test_that("bool_add_randomness = FALSE is exact, TRUE is Poisson (SM4)", {
  embedding_mat <- .clustered_embedding_mat()
  embedding_coefficient_vec <- rep(0.2, ncol(embedding_mat))

  set.seed(10)
  res_a <- generate_simulation(embedding_mat = embedding_mat,
                               bool_add_randomness = FALSE,
                               embedding_coefficient_vec =
                                 embedding_coefficient_vec,
                               num_lineages = 5)
  set.seed(10)
  res_b <- generate_simulation(embedding_mat = embedding_mat,
                               bool_add_randomness = FALSE,
                               embedding_coefficient_vec =
                                 embedding_coefficient_vec,
                               num_lineages = 5)
  expect_true(identical(res_a$cell_fate_potential, res_b$cell_fate_potential))

  # With randomness on, the realized counts must average to the expected ones.
  expected_vec <- 10^res_a$cell_fate_potential_truth
  num_trials <- 200
  realized_mat <- sapply(seq_len(num_trials), function(trial){
    set.seed(trial)
    res <- generate_simulation(embedding_mat = embedding_mat,
                               bool_add_randomness = TRUE,
                               embedding_coefficient_vec =
                                 embedding_coefficient_vec,
                               num_lineages = 5)
    10^res$cell_fate_potential - 1
  })
  observed_vec <- rowMeans(realized_mat)
  expect_true(stats::cor(observed_vec, expected_vec) >= 0.99)
  expect_true(abs(mean(observed_vec) - mean(expected_vec)) /
                mean(expected_vec) <= 0.05)
})

## SM5: ties the two returned objects together. Note the rounding is applied to
## the lineage *sum*, not to each cell, so summing rounded per-cell counts is
## not the same thing.
test_that("lineage_future_size is the rounded within-lineage sum (SM5)", {
  embedding_mat <- .clustered_embedding_mat()

  set.seed(10)
  res <- generate_simulation(embedding_mat = embedding_mat,
                             embedding_coefficient_vec =
                               rep(0.2, ncol(embedding_mat)),
                             num_lineages = 5)

  realized_vec <- 10^res$cell_fate_potential - 1
  expected_vec <- sapply(levels(res$lineage_assignment), function(lev){
    round(sum(realized_vec[which(res$lineage_assignment == lev)]))
  })

  expect_true(all(names(res$lineage_future_size) ==
                    levels(res$lineage_assignment)))
  expect_true(max(abs(res$lineage_future_size - expected_vec)) <= 1e-8)
})

## SM6: the rows of `summary_mat` are assembled positionally by `rbind()` inside
## `.compute_summary_lineages()` and named afterwards, so a reordering there
## would relabel every row without any other symptom.
test_that("summary_mat rows are what their names claim (SM6)", {
  embedding_mat <- .clustered_embedding_mat()

  set.seed(10)
  res <- generate_simulation(embedding_mat = embedding_mat,
                             embedding_coefficient_vec =
                               rep(0.2, ncol(embedding_mat)),
                             num_lineages = 5)

  expect_true(all(rownames(res$summary_mat) ==
                    c("mean", "median", "sd", "range", "future_size")))
  expect_true(all(colnames(res$summary_mat) ==
                    levels(res$lineage_assignment)))

  for(lineage in levels(res$lineage_assignment)){
    idx <- which(res$lineage_assignment == lineage)
    truth_vec <- res$cell_fate_potential_truth[idx]

    expect_true(abs(res$summary_mat["mean", lineage] -
                      mean(truth_vec)) <= 1e-10, info = lineage)
    expect_true(abs(res$summary_mat["median", lineage] -
                      stats::median(truth_vec)) <= 1e-10, info = lineage)
    expect_true(abs(res$summary_mat["sd", lineage] -
                      stats::sd(truth_vec)) <= 1e-10, info = lineage)
    expect_true(abs(res$summary_mat["range", lineage] -
                      diff(range(truth_vec))) <= 1e-10, info = lineage)
    expect_true(res$summary_mat["future_size", lineage] ==
                  res$lineage_future_size[lineage], info = lineage)
  }
})

## SM9: the simulators take no `seed_number` -- seeding is the caller's job --
## so this is the guard that would catch someone adding an un-seeded draw that
## bypasses the ambient stream.
test_that("generate_simulation is reproducible under set.seed (SM9)", {
  embedding_mat <- .clustered_embedding_mat()

  for(seed_val in c(1, 5, 10)){
    label <- paste0("seed = ", seed_val)
    set.seed(seed_val)
    res_a <- generate_simulation(embedding_mat = embedding_mat,
                                 num_lineages = 5)
    set.seed(seed_val)
    res_b <- generate_simulation(embedding_mat = embedding_mat,
                                 num_lineages = 5)

    expect_true(identical(res_a$lineage_assignment, res_b$lineage_assignment),
                info = label)
    expect_true(identical(res_a$cell_fate_potential, res_b$cell_fate_potential),
                info = label)
    expect_true(identical(res_a$lineage_future_size, res_b$lineage_future_size),
                info = label)
  }
})

## SM11: `.kmeans_seed()` over-clusters into 2K and keeps clusters 1..K. Which
## half that is depends on k-means' arbitrary numbering, so it is effectively a
## random half -- deliberate, per Kevin, and pinned here so nobody "corrects" it
## to K centres or to sampling K of the 2K.
test_that(".kmeans_seed over-clusters into 2K and keeps K (SM11)", {
  embedding_mat <- .clustered_embedding_mat()

  for(K in c(3, 5, 8)){
    label <- paste0("K = ", K)
    set.seed(10)
    res <- .kmeans_seed(embedding_mat = embedding_mat, K = K)

    expect_true(length(res) == K, info = label)
    expect_true(all(res %in% seq_len(nrow(embedding_mat))), info = label)
    expect_true(anyDuplicated(res) == 0, info = label)

    # The seeds must come from the first K of 2K k-means clusters.
    set.seed(10)
    kmeans_res <- suppressWarnings(stats::kmeans(embedding_mat, centers = 2 * K))
    expect_true(all(sort(unique(kmeans_res$cluster[res])) == seq_len(K)),
                info = label)
  }
})
