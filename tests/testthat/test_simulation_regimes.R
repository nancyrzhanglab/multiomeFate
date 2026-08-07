context("Test that the simulators simulate the regimes they are named for")

## The claims in this file are the ones the paper's simulation figures rest on:
## that `generate_simulation()` and `generate_simulation_plastic()` produce two
## *different* regimes, and that CYFER recovers truth from the first of them.
## Nothing else in the suite checks either.

## A clustered embedding, so the k-means seeding has real structure to find.
## A single spherical blob would make the regime test vacuous, since every
## lineage would then have the same mean fate potential by construction.
.regime_embedding_mat <- function(d = 4,
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

## The preprocessing every caller of `cyfer()` repeats: drop lineages with no
## future cells, drop lineages with a single cell (CV needs at least two), and
## scale the features. Lives in the analysis scripts rather than the package,
## which is why it is spelled out here rather than called.
.fit_cyfer_on_simulation <- function(sim_res,
                                     embedding_mat,
                                     num_folds = 3,
                                     seed_number = 1){
  cell_lineage <- as.character(sim_res$lineage_assignment)
  names(cell_lineage) <- names(sim_res$lineage_assignment)

  keep_vec <- names(sim_res$lineage_future_size)[
    sim_res$lineage_future_size > 0]
  size_table <- table(cell_lineage)
  keep_vec <- intersect(keep_vec, names(size_table)[size_table >= 2])
  cell_idx <- which(cell_lineage %in% keep_vec)

  cell_features <- scale(embedding_mat[cell_idx, , drop = FALSE])
  cell_lineage <- cell_lineage[cell_idx]
  lineage_future_count <- sim_res$lineage_future_size[keep_vec]

  cv_fit_list <- cyfer(cell_features = cell_features,
                       cell_lineage = cell_lineage,
                       lineage_future_count = lineage_future_count,
                       lambda_initial = NA,
                       lambda_sequence_length = 10,
                       num_folds = num_folds,
                       seed_number = seed_number,
                       verbose = 0)
  fit_res <- cyfer_finalize(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            fit_res = cv_fit_list,
                            lineage_future_count = lineage_future_count)

  list(cell_idx = cell_idx,
       fit_res = fit_res)
}

## SM1: the single most important test in this file. "Priming" means the
## variation lives *between* lineages -- pre-existing high-potential clones --
## and "plastic" means it lives *within* them. If the two simulators ever stop
## separating on these two statistics, the simulation figures no longer
## illustrate what their captions say.
test_that("priming spreads means, plastic spreads within-lineage sd (SM1)", {
  embedding_mat <- .regime_embedding_mat()
  embedding_coefficient_vec <- rep(0.2, ncol(embedding_mat))

  for(num_lineages in c(5, 10)){
    label <- paste0("num_lineages = ", num_lineages)

    set.seed(10)
    priming_res <- generate_simulation(
      embedding_mat = embedding_mat,
      embedding_coefficient_vec = embedding_coefficient_vec,
      num_lineages = num_lineages)
    set.seed(10)
    plastic_res <- generate_simulation_plastic(
      embedding_mat = embedding_mat,
      embedding_coefficient_vec = embedding_coefficient_vec,
      num_lineages = num_lineages)

    priming_mean_spread <- stats::var(priming_res$summary_mat["mean", ])
    plastic_mean_spread <- stats::var(plastic_res$summary_mat["mean", ])
    priming_sd_val <- mean(priming_res$summary_mat["sd", ])
    plastic_sd_val <- mean(plastic_res$summary_mat["sd", ])

    expect_true(priming_mean_spread > plastic_mean_spread, info = label)
    expect_true(plastic_sd_val > priming_sd_val, info = label)
  }
})

## SM7: the only end-to-end test in the package, and exactly what the vignette
## claims works -- simulate from the CYFER model, fit CYFER, recover the truth.
test_that("cyfer recovers the simulated fate potential end to end (SM7)", {
  embedding_mat <- .regime_embedding_mat()

  set.seed(10)
  sim_res <- generate_simulation(
    embedding_mat = embedding_mat,
    embedding_coefficient_vec = rep(0.3, ncol(embedding_mat)),
    num_lineages = 10)

  fit_list <- .fit_cyfer_on_simulation(sim_res = sim_res,
                                       embedding_mat = embedding_mat)
  truth_vec <- sim_res$cell_fate_potential_truth[fit_list$cell_idx]
  correlation_val <- stats::cor(fit_list$fit_res$cell_imputed_score, truth_vec)

  expect_true(correlation_val >= 0.9)
  expect_true(length(fit_list$fit_res$cell_imputed_score) ==
                length(fit_list$cell_idx))
})

## SM8: `fatefeatures_mat` is signal that drives the true fate potential but is
## never handed to the estimator, so it is the misspecification knob the
## revision simulations turn. If withholding it did not degrade recovery, those
## simulations would not be measuring what they claim.
test_that("withheld fatefeatures signal degrades recovery (SM8)", {
  embedding_mat <- .regime_embedding_mat()
  set.seed(99)
  fatefeatures_mat <- matrix(stats::rnorm(nrow(embedding_mat) * 3),
                             nrow = nrow(embedding_mat),
                             ncol = 3)
  colnames(fatefeatures_mat) <- paste0("ff", seq_len(3))

  coefficient_vec <- c(0, 0.5, 1, 2)
  correlation_vec <- sapply(coefficient_vec, function(coefficient_val){
    set.seed(10)
    sim_res <- generate_simulation(
      embedding_mat = embedding_mat,
      embedding_coefficient_vec = rep(0.3, ncol(embedding_mat)),
      fatefeatures_coefficient_vec = rep(coefficient_val, 3),
      fatefeatures_mat = fatefeatures_mat,
      num_lineages = 10)
    fit_list <- .fit_cyfer_on_simulation(sim_res = sim_res,
                                         embedding_mat = embedding_mat)
    stats::cor(fit_list$fit_res$cell_imputed_score,
               sim_res$cell_fate_potential_truth[fit_list$cell_idx])
  })

  # Recovery must be near-perfect with no withheld signal, and must degrade
  # monotonically as the withheld coefficient grows.
  expect_true(correlation_vec[1] >= 0.9)
  for(i in seq_len(length(correlation_vec) - 1)){
    label <- paste0("fatefeatures coefficient ", coefficient_vec[i], " -> ",
                    coefficient_vec[i + 1])
    expect_true(correlation_vec[i] > correlation_vec[i + 1], info = label)
  }
  expect_true(correlation_vec[length(correlation_vec)] < 0.6)
})
