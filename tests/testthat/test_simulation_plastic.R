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
  cell_contribution_truth <- abs(stats::rnorm(n))
  names(cell_contribution_truth) <- paste0("cell:", 1:n)
  num_lineages <- 10
  
  res <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    gamma = 1,
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
  cell_contribution_truth <- abs(stats::rnorm(n))
  names(cell_contribution_truth) <- paste0("cell:", 1:n)
  num_lineages <- 10
  
  tmp <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    gamma = 1,
    rho = NA
  )
  prob_mat <- tmp$prob_mat
  
  lineage_assignment <- .assign_plastic_lineages(enforce_equal_size = TRUE,
                                                 prob_mat = prob_mat)
  
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
  expect_true(all(names(res$cell_fate_potential) == names(res$lineage_assignment)))
})

test_that("generate_simulation_plastic generates lineages with different variances", {
  set.seed(10)
  n <- 1000; d <- 5
  K <- 5
  embedding_mat <- matrix(stats::rt(n*d, df = 3), nrow = n, ncol = d)
  
  res1 <- generate_simulation_plastic(embedding_mat, 
                                      num_lineages = 50)
  sd1 <- sd(res1$summary_mat["sd",])
  
  res2 <- generate_simulation_plastic(embedding_mat, 
                                      num_lineages = 50,
                                      lineage_sd_spread = 1)
  sd2 <- sd(res2$summary_mat["sd",])
  
  expect_true(sd1 > sd2)
})


test_that("generate_simulation_plastic generates lineages with different means", {
  set.seed(10)
  n <- 1000; d <- 5
  K <- 5
  embedding_mat <- matrix(stats::rt(n*d, df = 3), nrow = n, ncol = d)
  
  res1 <- generate_simulation_plastic(embedding_mat, 
                                      num_lineages = 20,
                                      lineage_mean_spread = NA, 
                                      lineage_sd_spread = 1)
  sd1 <- sd(res1$summary_mat["mean",])
  
  res2 <- generate_simulation_plastic(embedding_mat, 
                                      num_lineages = 20,
                                      lineage_mean_spread = 1, 
                                      lineage_sd_spread = 1)
  sd2 <- sd(res2$summary_mat["mean",])
  
  expect_true(sd1 > sd2)
})


##############################################################################
## generate_simulation_plastic() specifics (test-plan-full 2026-08-06)
##############################################################################

.plastic_embedding_mat <- function(d = 4,
                                   n = 200,
                                   seed_number = 10){
  set.seed(seed_number)
  embedding_mat <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  rownames(embedding_mat) <- paste0("cell:", seq_len(n))
  colnames(embedding_mat) <- paste0("dim", seq_len(d))

  embedding_mat
}

## SM12: `lineage_mean_spread` accepts only 1 and NA, and anything else is
## silently coerced to 1 after a warning. That warning is the only thing between
## a caller and a simulation that is not the one they asked for, so pin both the
## warning and the coercion.
test_that("lineage_mean_spread honours only 1 and NA (SM12)", {
  embedding_mat <- .plastic_embedding_mat()

  for(spread_val in c(2, 5, 0.5)){
    label <- paste0("lineage_mean_spread = ", spread_val)
    set.seed(10)
    expect_warning(generate_simulation_plastic(
      embedding_mat = embedding_mat,
      lineage_mean_spread = spread_val,
      num_lineages = 5),
      "can only handle NA or 1", info = label)

    set.seed(10)
    res_other <- suppressWarnings(generate_simulation_plastic(
      embedding_mat = embedding_mat,
      lineage_mean_spread = spread_val,
      num_lineages = 5))
    set.seed(10)
    res_one <- generate_simulation_plastic(embedding_mat = embedding_mat,
                                           lineage_mean_spread = 1,
                                           num_lineages = 5)

    expect_true(identical(res_other$lineage_assignment,
                          res_one$lineage_assignment), info = label)
  }

  # 1 and NA are the two supported values and neither warns on its own.
  set.seed(10)
  expect_silent(generate_simulation_plastic(embedding_mat = embedding_mat,
                                            lineage_mean_spread = 1,
                                            num_lineages = 5))
  set.seed(10)
  expect_silent(generate_simulation_plastic(embedding_mat = embedding_mat,
                                            lineage_mean_spread = NA,
                                            lineage_sd_spread = 2,
                                            num_lineages = 5))
})

## SM13: the size cap exists so that lineage size does not itself carry the
## signal -- in the plastic regime lineages differ in composition, not in how
## many cells they hold.
test_that("lineage_mean_spread = 1 caps lineage sizes (SM13)", {
  embedding_mat <- .plastic_embedding_mat()
  n <- nrow(embedding_mat)

  for(num_lineages in c(4, 5, 8)){
    label <- paste0("num_lineages = ", num_lineages)
    set.seed(10)
    res <- generate_simulation_plastic(embedding_mat = embedding_mat,
                                       lineage_mean_spread = 1,
                                       num_lineages = num_lineages)
    size_vec <- table(res$lineage_assignment)
    expect_true(all(size_vec <= ceiling(n / num_lineages)), info = label)
    expect_true(sum(size_vec) == n, info = label)
  }
})

## SM14: the spread ladder is the mechanism of the plastic regime -- lineage 1
## is the widest and lineage K the tightest. `lineage_sd_vec` is constructed
## with `seq()` so it is monotone by construction; the realized `sd` row only
## has to be monotone on average, since the assignment is a draw.
test_that("the lineage spread ladder decreases (SM14)", {
  cell_contribution_truth <- exp(stats::rnorm(400))
  names(cell_contribution_truth) <- paste0("cell:", seq_len(400))

  res <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = 8,
    gamma = 1,
    rho = 3)

  expect_true(all(diff(res$lineage_sd_vec) < 0))
  expect_true(res$lineage_sd_vec[1] > res$lineage_sd_vec[8])

  # The realized spread must trend the same way. Compare the first and last
  # thirds rather than asserting strict monotonicity of a random assignment.
  embedding_mat <- .plastic_embedding_mat(n = 400)
  set.seed(10)
  sim_res <- generate_simulation_plastic(embedding_mat = embedding_mat,
                                         lineage_sd_spread = 3,
                                         num_lineages = 8)
  sd_vec <- sim_res$summary_mat["sd", ]
  expect_true(mean(sd_vec[1:3]) > mean(sd_vec[6:8]))
})

## SM15: the de-duplication at the meeting point is the fiddly part -- an
## off-by-one there would drop or repeat the middle element, which would then
## silently drop a cell from the assignment loop.
test_that(".reorder_by_contribution returns a permutation (SM15)", {
  for(n in c(1, 2, 3, 4, 5, 10, 11, 50, 51)){
    label <- paste0("n = ", n)
    set.seed(n)
    vec <- stats::runif(n)
    res <- .reorder_by_contribution(vec)

    expect_true(length(res) == n, info = label)
    expect_true(anyDuplicated(res) == 0, info = label)
    expect_true(all(sort(res) == seq_len(n)), info = label)
  }

  # It alternates tails: smallest, largest, second smallest, second largest.
  vec <- c(5, 1, 4, 2, 3)
  res <- .reorder_by_contribution(vec)
  expect_true(all(res == c(2, 1, 4, 3, 5)))
})

## SM16: every cell must be placed. The uniform fallback for a row whose
## probabilities are all at or below 1e-6 is what makes that true for a cell in
## the far tail of every lineage's Gaussian, and it is reachable in practice.
test_that(".assign_plastic_lineages places every cell (SM16)", {
  set.seed(10)
  prob_mat <- matrix(stats::runif(200 * 5), nrow = 200, ncol = 5)
  prob_mat <- prob_mat / rowSums(prob_mat)
  rownames(prob_mat) <- paste0("cell:", seq_len(200))
  colnames(prob_mat) <- paste0("lineage:", seq_len(5))

  for(enforce_equal_size in c(TRUE, FALSE)){
    label <- paste0("enforce_equal_size = ", enforce_equal_size)
    set.seed(10)
    res <- .assign_plastic_lineages(enforce_equal_size = enforce_equal_size,
                                    prob_mat = prob_mat)

    expect_true(length(res) == nrow(prob_mat), info = label)
    expect_true(!any(is.na(res)), info = label)
    expect_true(all(levels(res) == colnames(prob_mat)), info = label)
    expect_true(all(names(res) == rownames(prob_mat)), info = label)
  }

  # All levels are retained even when a lineage receives no cells at all.
  degenerate_mat <- prob_mat
  degenerate_mat[, 2:5] <- 0
  degenerate_mat[, 1] <- 1
  set.seed(10)
  res <- .assign_plastic_lineages(enforce_equal_size = FALSE,
                                  prob_mat = degenerate_mat)
  expect_true(all(levels(res) == colnames(prob_mat)))
  expect_true(length(unique(as.character(res))) == 1)

  # The all-below-1e-6 uniform fallback places the cell rather than erroring.
  tiny_mat <- matrix(1e-9, nrow = 20, ncol = 4)
  rownames(tiny_mat) <- paste0("cell:", seq_len(20))
  colnames(tiny_mat) <- paste0("lineage:", seq_len(4))
  set.seed(10)
  res <- .assign_plastic_lineages(enforce_equal_size = FALSE,
                                  prob_mat = tiny_mat)
  expect_true(!any(is.na(res)))
  expect_true(length(res) == 20)
})

## SM17
test_that("both spreads NA is an error (SM17)", {
  embedding_mat <- .plastic_embedding_mat()
  expect_error(generate_simulation_plastic(embedding_mat = embedding_mat,
                                           lineage_mean_spread = NA,
                                           lineage_sd_spread = NA,
                                           num_lineages = 5))
})