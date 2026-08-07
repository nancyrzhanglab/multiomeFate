context("Testing barcoding assignments")

## Fixtures -------------------------------------------------------------------

## Counts drawn from the model `barcoding_posterior()` assumes: a cell's true
## barcode is read at `signal_rate`, every other barcode at the much lower
## ambient `background_rate`. Returns the truth vector so recovery can be
## scored, which is the only way to catch an inversion in `gamma`.
.barcode_count_mat <- function(background_rate = 2,
                               num_cells = 300,
                               num_lineages = 10,
                               signal_rate = 200,
                               seed_number = 10){
  set.seed(seed_number)
  truth_idx <- rep(seq_len(num_lineages), length.out = num_cells)
  lin_mat <- matrix(stats::rpois(num_lineages * num_cells,
                                 lambda = background_rate),
                    nrow = num_lineages,
                    ncol = num_cells)
  for(i in seq_len(num_cells)){
    lin_mat[truth_idx[i], i] <- stats::rpois(1, lambda = signal_rate)
  }
  rownames(lin_mat) <- paste0("bc", seq_len(num_lineages))
  colnames(lin_mat) <- paste0("cell", seq_len(num_cells))

  list(lin_mat = lin_mat,
       truth_vec = rownames(lin_mat)[truth_idx])
}

## Six barcodes over two disjoint cell blocks: barcodes 1-3 share one per-cell
## rate profile, 4-6 share another. Within-block correlation runs about 0.94 and
## across-block about -0.66, so `cor_threshold = 0.55` separates them cleanly.
.planted_cluster_mat <- function(num_cells = 200,
                                 seed_number = 10){
  set.seed(seed_number)
  half_val <- num_cells / 2
  rate_a_vec <- stats::runif(half_val, min = 5, max = 50)
  rate_b_vec <- stats::runif(half_val, min = 5, max = 50)
  zero_vec <- rep(0, half_val)

  lin_mat <- rbind(c(stats::rpois(half_val, rate_a_vec), zero_vec),
                   c(stats::rpois(half_val, rate_a_vec), zero_vec),
                   c(stats::rpois(half_val, rate_a_vec), zero_vec),
                   c(zero_vec, stats::rpois(half_val, rate_b_vec)),
                   c(zero_vec, stats::rpois(half_val, rate_b_vec)),
                   c(zero_vec, stats::rpois(half_val, rate_b_vec)))
  rownames(lin_mat) <- paste0("bc", seq_len(6))
  colnames(lin_mat) <- paste0("cell", seq_len(num_cells))

  Matrix::Matrix(lin_mat, sparse = TRUE)
}

## Five barcodes built from three independent latent rate profiles, so that
## bc1-bc2, bc3-bc4, bc2-bc5 and bc4-bc5 clear 0.55 while every other pair does
## not. Crucially the bridging pair (bc4, bc5) is *last* in the column-major
## order `which(..., arr.ind = TRUE)` produces, so by the time it is processed
## bc1/bc2/bc5 sit in one cluster and bc3/bc4 in another -- which is what makes
## the two clusters merge and `warn_merging` fire.
.transitive_cluster_mat <- function(num_cells = 400,
                                    seed_number = 10){
  set.seed(seed_number)
  rate_a_vec <- stats::runif(num_cells, min = 20, max = 200)
  rate_b_vec <- stats::runif(num_cells, min = 20, max = 200)
  rate_c_vec <- stats::runif(num_cells, min = 20, max = 200)

  lin_mat <- rbind(stats::rpois(num_cells, rate_a_vec),
                   stats::rpois(num_cells, rate_a_vec + rate_c_vec),
                   stats::rpois(num_cells, rate_b_vec),
                   stats::rpois(num_cells, rate_b_vec + rate_c_vec),
                   stats::rpois(num_cells, rate_c_vec))
  rownames(lin_mat) <- paste0("bc", seq_len(5))
  colnames(lin_mat) <- paste0("cell", seq_len(num_cells))

  Matrix::Matrix(lin_mat, sparse = TRUE)
}

##############################################################################
## barcoding_posterior()
##############################################################################

## BC1: the one test that would catch a sign or inversion error in `gamma`.
## Everything else about this function is an invariant that a flipped ratio
## would still satisfy.
test_that("barcoding_posterior recovers the true barcode (BC1)", {
  param_grid <- expand.grid(num_lineages = c(5, 10, 20),
                            signal_rate = c(50, 200),
                            seed_number = 1:3)

  for(i in seq_len(nrow(param_grid))){
    label <- paste0("num_lineages = ", param_grid[i, "num_lineages"],
                    ", signal_rate = ", param_grid[i, "signal_rate"],
                    ", seed = ", param_grid[i, "seed_number"])
    fixture <- .barcode_count_mat(num_lineages = param_grid[i, "num_lineages"],
                                  signal_rate = param_grid[i, "signal_rate"],
                                  seed_number = param_grid[i, "seed_number"])
    res <- barcoding_posterior(lin_mat = fixture$lin_mat)

    called_vec <- rownames(res$posterior_mat)[
      apply(res$posterior_mat, 2, which.max)]
    expect_true(mean(called_vec == fixture$truth_vec) >= 0.95, info = label)
  }
})

## BC2
test_that("barcoding_posterior columns are probability vectors (BC2)", {
  fixture <- .barcode_count_mat()
  res <- barcoding_posterior(lin_mat = fixture$lin_mat)

  expect_true(all(res$posterior_mat >= 0))
  expect_true(all(res$posterior_mat <= 1))
  expect_true(max(abs(colSums(res$posterior_mat) - 1)) <= 1e-10)
  expect_true(all(dim(res$posterior_mat) == dim(fixture$lin_mat)))
  expect_true(identical(dimnames(res$posterior_mat),
                        dimnames(fixture$lin_mat)))
})

## BC3: the overflow guard is an algebraic identity, not an approximation, so
## the two branches must agree exactly wherever the plain branch is usable.
test_that(".multinomial_posterior_vector branches agree below the shift (BC3)", {
  lgamma_vec <- log(c(0.5, 1, 2, 3, 4))

  for(count_scale in c(0, 1, 2)){
    label <- paste0("count_scale = ", count_scale)
    lin_count <- c(0, 1, 2, 1, 0) * count_scale
    # Precondition: the plain branch is the one that would be taken.
    expect_true(max(lin_count * lgamma_vec) <= 10, info = label)

    res_plain <- .multinomial_posterior_vector(bool_force_rebase = FALSE,
                                               lgamma = lgamma_vec,
                                               lin_count = lin_count)
    res_shift <- .multinomial_posterior_vector(bool_force_rebase = TRUE,
                                               lgamma = lgamma_vec,
                                               lin_count = lin_count)
    res_direct <- exp(lgamma_vec)^lin_count
    res_direct <- res_direct / sum(res_direct)

    expect_true(max(abs(res_plain - res_shift)) <= 1e-12, info = label)
    expect_true(max(abs(res_plain - res_direct)) <= 1e-12, info = label)
  }
})

## BC4: counts large enough that `exp(lgammaX)` overflows to `Inf`, which the
## plain branch would turn into `Inf/Inf = NaN`. The default
## `bool_force_rebase = FALSE` must still take the safe branch, because the
## `max(lgammaX) > 10` test fires long before the overflow does.
test_that(".multinomial_posterior_vector survives overflow (BC4)", {
  lgamma_vec <- log(c(2, 3, 4))
  lin_count <- c(1200, 10, 5)

  # Precondition: the naive form really does overflow here.
  expect_true(any(is.infinite(exp(lin_count * lgamma_vec))))

  res <- .multinomial_posterior_vector(bool_force_rebase = FALSE,
                                       lgamma = lgamma_vec,
                                       lin_count = lin_count)
  expect_true(all(is.finite(res)))
  expect_true(abs(sum(res) - 1) <= 1e-12)

  # Just above and just below the max(lgammaX) > 10 switch, the two branches
  # must still describe the same distribution.
  for(target_val in c(9.9, 10.1)){
    label <- paste0("max(lgammaX) = ", target_val)
    lin_count <- c(target_val / lgamma_vec[1], 0, 0)
    res_default <- .multinomial_posterior_vector(bool_force_rebase = FALSE,
                                                 lgamma = lgamma_vec,
                                                 lin_count = lin_count)
    res_forced <- .multinomial_posterior_vector(bool_force_rebase = TRUE,
                                                lgamma = lgamma_vec,
                                                lin_count = lin_count)
    expect_true(max(abs(res_default - res_forced)) <= 1e-12, info = label)
  }
})

## BC5: `beta1_mean` uses `na.rm = TRUE`, so a barcode nobody maximizes at
## cannot poison the posterior. Removing that `na.rm` would blank the whole
## matrix silently, so the invariant is worth pinning even though it holds by
## construction today.
test_that("a barcode no cell maximizes at gives NA beta1, not NA posterior", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  lin_mat["bc1", ] <- 0

  res <- barcoding_posterior(lin_mat = lin_mat)

  expect_true(res$lineage_num_winner[1] == 0)
  expect_true(is.na(res$beta1[1]))
  expect_true(all(!is.na(res$beta1[-1])))
  expect_true(!is.na(res$beta1_mean))
  expect_true(all(!is.na(res$posterior_mat)))
  expect_true(max(abs(colSums(res$posterior_mat) - 1)) <= 1e-10)
})

## BC6: the winsorization exists to stop `gamma = beta1_mean / beta0` exploding
## on a barcode whose background is essentially zero. Assert the bound rather
## than a value, since the quantiles move with the fixture.
test_that("beta0 winsorization keeps gamma finite and bounded (BC6)", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  # bc3 has signal in the cells it wins and no background anywhere else, so its
  # raw beta0 is exactly 0 and its unwinsorized gamma would be Inf.
  won_idx <- which(fixture$truth_vec == "bc3")
  lin_mat["bc3", -won_idx] <- 0

  res <- barcoding_posterior(lin_mat = lin_mat)

  expect_true(all(is.finite(res$gamma)))
  expect_true(res$beta0[3] == 0)

  lower_val <- stats::quantile(res$beta0[res$beta0 > 1e-8], 0.02)
  expect_true(abs(res$gamma[3] - res$beta1_mean / lower_val) <= 1e-10)
  expect_true(all(res$gamma <= res$beta1_mean / lower_val + 1e-10))
})

## BC7: `names(lin_mat)` is NULL for a matrix, so assigning it strips the names
## rather than setting them. Use identical() rather than all(x == y) here --
## `all(NULL == character(n))` is `all(logical(0))`, which is TRUE, so the
## obvious form of this test passes against unnamed output.
test_that("beta0 and beta1 are named by barcode (BC7)", {
  fixture <- .barcode_count_mat()
  res <- barcoding_posterior(lin_mat = fixture$lin_mat)

  expect_true(identical(names(res$beta0), rownames(fixture$lin_mat)))
  expect_true(identical(names(res$beta1), rownames(fixture$lin_mat)))
  expect_true(identical(names(res$gamma), rownames(fixture$lin_mat)))
  expect_true(identical(names(res$lineage_num_winner),
                        rownames(fixture$lin_mat)))
  # Indexable by barcode name, which is the whole point.
  expect_true(res$beta0["bc3"] == res$beta0[3])
})

## BC8: `gamma[b]` divides the *global* mean signal rate by the barcode's own
## background, never `beta1[b]`. This is deliberate -- `beta1[b]` is estimated
## from however many cells happened to maximize at `b`, which for a rare barcode
## is one cell or none -- so pin it, because switching to the per-barcode rate
## would leave BC1 passing on easy fixtures and fail on rare barcodes.
test_that("gamma uses the global mean signal rate, not beta1[b] (BC8)", {
  fixture <- .barcode_count_mat()
  res <- barcoding_posterior(lin_mat = fixture$lin_mat)

  implied_beta0_vec <- res$beta1_mean / res$gamma
  # The numerator is constant across barcodes, so the implied denominators must
  # reproduce the winsorized beta0 exactly.
  lower_val <- stats::quantile(res$beta0[res$beta0 > 1e-8], 0.02)
  upper_val <- stats::quantile(res$beta0[res$beta0 > 1e-8], 0.98)
  expected_vec <- pmax(pmin(res$beta0, upper_val), lower_val)
  expect_true(max(abs(implied_beta0_vec - expected_vec)) <= 1e-10)

  # And it is genuinely not the per-barcode rate: beta1 varies across barcodes.
  expect_true(stats::sd(res$beta1) > 0)
})
## The whole pipeline works in `dgCMatrix`, so `barcoding_posterior()` must
## accept one and give bit-identical answers to the dense equivalent. Callers
## used to have to bridge the two with a bare `as.matrix()`.
test_that("barcoding_posterior accepts sparse and dense alike", {
  fixture <- .barcode_count_mat()
  sparse_mat <- Matrix::Matrix(fixture$lin_mat, sparse = TRUE)
  expect_true(inherits(sparse_mat, "dgCMatrix"))

  res_dense <- barcoding_posterior(lin_mat = fixture$lin_mat)
  res_sparse <- barcoding_posterior(lin_mat = sparse_mat)

  expect_true(max(abs(as.matrix(res_sparse$posterior_mat) -
                        res_dense$posterior_mat)) <= 1e-12)
  expect_true(max(abs(res_sparse$gamma - res_dense$gamma)) <= 1e-12)
  expect_true(max(abs(res_sparse$beta0 - res_dense$beta0)) <= 1e-12)
  expect_true(identical(names(res_sparse$beta0), names(res_dense$beta0)))

  # Row names are what become the lineage names, so they are required.
  unnamed_mat <- fixture$lin_mat
  rownames(unnamed_mat) <- NULL
  expect_error(barcoding_posterior(lin_mat = unnamed_mat))
})

##############################################################################
## barcode_clustering()
##############################################################################

## BC9
test_that("barcode_clustering recovers planted clusters (BC9)", {
  lin_mat <- .planted_cluster_mat()
  res <- barcode_clustering(lin_mat = lin_mat,
                            cell_lower_limit = 50,
                            cor_threshold = 0.55,
                            warn_merging = FALSE)

  expect_true(length(res$lineage_clusters) == 2)
  cluster_list <- lapply(res$lineage_clusters, sort)
  expect_true(any(sapply(cluster_list, function(vec){
    identical(vec, c("bc1", "bc2", "bc3"))
  })))
  expect_true(any(sapply(cluster_list, function(vec){
    identical(vec, c("bc4", "bc5", "bc6"))
  })))
  expect_true(all(res$minimum_correlation >= 0.55))
})

## BC10
test_that("barcode_clustering returns all-NULL when nothing correlates (BC10)", {
  set.seed(10)
  lin_mat <- matrix(stats::rpois(6 * 200, lambda = 20), nrow = 6, ncol = 200)
  rownames(lin_mat) <- paste0("bc", seq_len(6))
  colnames(lin_mat) <- paste0("cell", seq_len(200))
  lin_mat <- Matrix::Matrix(lin_mat, sparse = TRUE)

  res <- barcode_clustering(lin_mat = lin_mat,
                            cell_lower_limit = 50,
                            cor_threshold = 0.55,
                            warn_merging = FALSE)

  expect_true(is.null(res$arr_idx))
  expect_true(is.null(res$lineage_clusters))
  expect_true(is.null(res$uniq_lineage))
  expect_true(is.null(res$minimum_correlation))
})

## BC11: a correlation over a handful of cells is not evidence, so
## `cell_lower_limit` must veto the merge however perfect the correlation.
test_that("barcodes below cell_lower_limit are never merged (BC11)", {
  lin_mat <- as.matrix(.planted_cluster_mat())
  # Two perfectly correlated barcodes, supported on 5 cells each.
  small_vec <- rep(0, ncol(lin_mat))
  small_vec[seq_len(5)] <- c(10, 20, 30, 40, 50)
  lin_mat <- rbind(lin_mat, bc7 = small_vec, bc8 = small_vec)
  lin_mat <- Matrix::Matrix(lin_mat, sparse = TRUE)

  res <- barcode_clustering(lin_mat = lin_mat,
                            cell_lower_limit = 100,
                            cor_threshold = 0.55,
                            warn_merging = FALSE)

  expect_true(!any(c("bc7", "bc8") %in% unlist(res$lineage_clusters)))
})

## BC12: single-linkage is transitive, so a chain of merges can pull together
## two barcodes that were never directly correlated. That is intended behaviour
## rather than a bug, but `minimum_correlation` is the diagnostic that reveals
## how far a cluster was stretched -- assert it drops below `cor_threshold`, or
## a reader will assume every pair in a cluster cleared the bar.
test_that("transitive merging chains clusters and warns (BC12)", {
  lin_mat <- .transitive_cluster_mat()

  expect_warning(barcode_clustering(lin_mat = lin_mat,
                                    cell_lower_limit = 50,
                                    cor_threshold = 0.55,
                                    warn_merging = TRUE),
                 "Merging happening")

  res <- barcode_clustering(lin_mat = lin_mat,
                            cell_lower_limit = 50,
                            cor_threshold = 0.55,
                            warn_merging = FALSE)

  expect_true(length(res$lineage_clusters) == 1)
  expect_true(identical(sort(res$lineage_clusters[[1]]),
                        paste0("bc", seq_len(5))))
  expect_true(res$minimum_correlation[1] < 0.55)
})

## BC13: the merge path asserts `length(val) == 2`, i.e. that a single pair can
## bridge at most two existing clusters. A pair has two members and each member
## belongs to at most one cluster, so it holds by construction. Pin it over the
## fixture that exercises the merge, since the assertion's message is unhelpful
## if it ever does fire.
test_that("a bridging pair joins exactly two clusters (BC13)", {
  lin_mat <- .transitive_cluster_mat()
  res <- barcode_clustering(lin_mat = lin_mat,
                            cell_lower_limit = 50,
                            cor_threshold = 0.55,
                            warn_merging = FALSE)

  expect_true(all(!is.na(res$uniq_lineage[, "Cluster"])))
  expect_true(length(unique(res$uniq_lineage[, "Cluster"])) == 1)
})

##############################################################################
## barcode_combine()
##############################################################################

## BC14: the strongest invariant available for this function -- merging rows
## must move counts between rows, never create or destroy them.
test_that("barcode_combine conserves counts exactly (BC14)", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  cluster_list <- list(c1 = c("bc1", "bc2"), c2 = c("bc4", "bc5", "bc6"))

  res <- barcode_combine(lin_mat = lin_mat, lineage_clusters = cluster_list)

  expect_true(all(colSums(res) == colSums(lin_mat)))
  expect_true(sum(res) == sum(lin_mat))
  expect_true(ncol(res) == ncol(lin_mat))
  expect_true(all(colnames(res) == colnames(lin_mat)))
})

## BC15
test_that("barcode_combine row count and naming are as documented (BC15)", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  cluster_list <- list(c1 = c("bc2", "bc1"), c2 = c("bc6", "bc4", "bc5"))

  res <- barcode_combine(lin_mat = lin_mat, lineage_clusters = cluster_list)

  expect_true(nrow(res) == nrow(lin_mat) -
                sum(lengths(cluster_list)) + length(cluster_list))
  # The surviving row takes the first barcode of the cluster in *sorted* order,
  # not the order the caller listed them in.
  expect_true("bc1" %in% rownames(res))
  expect_true("bc4" %in% rownames(res))
  expect_true(!any(c("bc2", "bc5", "bc6") %in% rownames(res)))
  expect_true(all(res["bc1", ] == lin_mat["bc1", ] + lin_mat["bc2", ]))
})

## BC16: row order is deliberately not preserved -- untouched barcodes first,
## then one row per cluster. Downstream code must index by name.
test_that("barcode_combine does not preserve row order (BC16)", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  cluster_list <- list(c1 = c("bc1", "bc2"))

  res <- barcode_combine(lin_mat = lin_mat, lineage_clusters = cluster_list)

  expect_true(rownames(res)[nrow(res)] == "bc1")
  expect_true(all(rownames(res)[seq_len(nrow(res) - 1)] ==
                    paste0("bc", 3:10)))
})

## BC17
test_that("barcode_combine skips NA tombstone entries (BC17)", {
  fixture <- .barcode_count_mat()
  lin_mat <- fixture$lin_mat
  cluster_list <- list(c1 = c("bc1", "bc2"), c2 = NA)

  res <- barcode_combine(lin_mat = lin_mat, lineage_clusters = cluster_list)

  expect_true(nrow(res) == nrow(lin_mat) - 1)
  expect_true(all(colSums(res) == colSums(lin_mat)))
  expect_true(all(res["bc1", ] == lin_mat["bc1", ] + lin_mat["bc2", ]))
})

## BC18
test_that("barcode_combine gives the same answer dense and sparse (BC18)", {
  fixture <- .barcode_count_mat()
  cluster_list <- list(c1 = c("bc1", "bc2"), c2 = c("bc4", "bc5"))

  res_dense <- barcode_combine(lin_mat = fixture$lin_mat,
                               lineage_clusters = cluster_list)
  res_sparse <- barcode_combine(
    lin_mat = Matrix::Matrix(fixture$lin_mat, sparse = TRUE),
    lineage_clusters = cluster_list)

  expect_true(all(rownames(res_dense) == rownames(res_sparse)))
  expect_true(max(abs(as.matrix(res_sparse) - as.matrix(res_dense))) == 0)
})

##############################################################################
## barcoding_assignment()
##############################################################################

## BC19: the criterion is a *margin* between the top two posteriors, not a level.
## This is the function's whole contract and the thing most likely to be
## misremembered as "assign when the top posterior is confident enough".
test_that("barcoding_assignment thresholds the margin, not the level (BC19)", {
  posterior_mat <- cbind(near_tie = c(0.95, 0.90, 0.01),
                         clear_winner = c(0.40, 0.10, 0.05))
  rownames(posterior_mat) <- c("bcA", "bcB", "bcC")

  res <- barcoding_assignment(posterior_mat = posterior_mat,
                              difference_val = 0.2)

  # High posterior, tiny margin -> unassigned.
  expect_true(is.na(res["near_tie"]))
  # Much lower posterior, comfortable margin -> assigned.
  expect_true(res["clear_winner"] == "bcA")
})

## BC20
test_that("barcoding_assignment is monotone in difference_val (BC20)", {
  fixture <- .barcode_count_mat(signal_rate = 8)
  res <- barcoding_posterior(lin_mat = fixture$lin_mat)
  posterior_mat <- res$posterior_mat

  assigned_all <- barcoding_assignment(posterior_mat = posterior_mat,
                                       difference_val = 0)
  assigned_none <- barcoding_assignment(posterior_mat = posterior_mat,
                                        difference_val = 1.1)
  expect_true(all(!is.na(assigned_all)))
  expect_true(all(is.na(assigned_none)))

  # The assigned set must shrink monotonically as the margin tightens.
  threshold_vec <- c(0, 0.2, 0.5, 0.8, 1.1)
  keep_list <- lapply(threshold_vec, function(threshold_val){
    which(!is.na(barcoding_assignment(posterior_mat = posterior_mat,
                                      difference_val = threshold_val)))
  })
  for(i in seq_len(length(keep_list) - 1)){
    label <- paste0("difference_val ", threshold_vec[i], " -> ",
                    threshold_vec[i + 1])
    expect_true(all(keep_list[[i + 1]] %in% keep_list[[i]]), info = label)
  }
})

## BC21
test_that("barcoding_assignment output is cell-named and cell-length (BC21)", {
  fixture <- .barcode_count_mat()
  res <- barcoding_posterior(lin_mat = fixture$lin_mat)
  assigned_vec <- barcoding_assignment(posterior_mat = res$posterior_mat)

  expect_true(length(assigned_vec) == ncol(res$posterior_mat))
  expect_true(all(names(assigned_vec) == colnames(res$posterior_mat)))
  expect_true(all(stats::na.omit(assigned_vec) %in%
                    rownames(res$posterior_mat)))
})

## BC22
test_that("exact ties are left unassigned at any positive margin (BC22)", {
  posterior_mat <- cbind(tied = c(0.5, 0.5, 0))
  rownames(posterior_mat) <- c("bcA", "bcB", "bcC")

  for(threshold_val in c(1e-8, 0.01, 0.2)){
    label <- paste0("difference_val = ", threshold_val)
    res <- barcoding_assignment(posterior_mat = posterior_mat,
                                difference_val = threshold_val)
    expect_true(is.na(res["tied"]), info = label)
  }

  # At exactly zero the comparison is `>=`, so a tie is assigned.
  res <- barcoding_assignment(posterior_mat = posterior_mat,
                              difference_val = 0)
  expect_true(res["tied"] == "bcA")
})

##############################################################################
## Defects found while writing these tests -- see the plan file, section 0
##############################################################################

## D6: `barcode_clustering()` returns `lineage_clusters = NULL` whenever no pair
## clears `cor_threshold`, and `barcode_combine()` documents NULL as "nothing to
## do" -- but `stopifnot(is.list(lineage_clusters))` rejects it first, because
## `is.list(NULL)` is FALSE. The `if(all(is.null(lineage_clusters)))` branch
## immediately below is therefore unreachable, and piping the two functions
## together on a dataset with no correlated barcodes crashes.
test_that("barcode_combine passes NULL through unchanged (D6)", {
  fixture <- .barcode_count_mat()
  res <- barcode_combine(lin_mat = fixture$lin_mat, lineage_clusters = NULL)
  expect_true(identical(res, fixture$lin_mat))
})

## D7: an all-tombstone list survives the NULL check, is emptied by the NA
## filter, and then hits `for(i in 1:len)` with `len == 0` -- the `1:0` hazard,
## which iterates `c(1, 0)` and subscripts an empty list.
test_that("an all-tombstone lineage_clusters behaves like NULL (D7)", {
  fixture <- .barcode_count_mat()
  res <- barcode_combine(lin_mat = fixture$lin_mat,
                         lineage_clusters = list(c1 = NA, c2 = NA))
  expect_true(identical(res, fixture$lin_mat))
})

## D8: `lin_mat[lineage_excluded_idx,]` has no `drop = FALSE`, so when exactly
## one barcode is left unclustered the row collapses to a vector and `rbind()`
## names the surviving row after the *variable* -- "lin_untouched" -- rather
## than the barcode. A real lineage name is silently replaced.
test_that("a single untouched barcode keeps its own name (D8)", {
  fixture <- .barcode_count_mat(num_lineages = 3)
  res <- barcode_combine(lin_mat = fixture$lin_mat,
                         lineage_clusters = list(c1 = c("bc1", "bc2")))

  expect_true(all(sort(rownames(res)) == c("bc1", "bc3")))
})

## D9: same missing `drop = FALSE`, on the cluster side. A cluster naming a
## single barcode collapses to a vector and `Matrix::colSums()` fails. Not
## reachable from `barcode_clustering()`, whose clusters are built from pairs
## and so always hold at least two barcodes, but reachable from a hand-built
## list.
test_that("a cluster of one barcode is a no-op for that barcode (D9)", {
  fixture <- .barcode_count_mat()
  res <- barcode_combine(lin_mat = fixture$lin_mat,
                         lineage_clusters = list(c1 = "bc1"))

  expect_true(nrow(res) == nrow(fixture$lin_mat))
  expect_true(all(colSums(res) == colSums(fixture$lin_mat)))
})

##############################################################################
## Pre-existing coverage
##############################################################################

test_that(".multinomial_posterior works", {
  set.seed(10)
  nlineages <- 20
  n <- 100
  gamma <- seq(0.1,10,length.out=nlineages)
  lin_mat <- matrix(sample(0:5, n*nlineages, replace = T, prob = 6:1),
                    nrow = nlineages, ncol = n)
  res <- .multinomial_posterior(
    bool_force_rebase = T,
    gamma = gamma,
    lin_mat = lin_mat
  )

  expect_true(is.matrix(res))
  expect_true(sum(abs(colSums(res) - 1)) <= 1e-3)
  expect_true(all(dim(res) == dim(lin_mat)))
})

##############################

## .multinomial_posterior_vector is correct

test_that(".multinomial_posterior_vector is correct", {
  trials <- 2000

  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    nlineages <- 20
    gamma <- seq(0.1,10,length.out=nlineages)
    lgamma <- log(gamma)
    lin_count <- sample(0:20, nlineages, replace = T, prob = 21:1)
    res <- .multinomial_posterior_vector(
      bool_force_rebase = F,
      lgamma = lgamma,
      lin_count = lin_count
    )
    res2 <- .multinomial_posterior_vector(
      bool_force_rebase = T,
      lgamma = lgamma,
      lin_count = lin_count
    )

    res3 <- gamma^(lin_count)
    res3 <- res3/sum(res3)

    bstarc <-  which.max(lin_count)
    ldeltabc <- lgamma - lgamma[bstarc]
    diff_vec <- lin_count - lin_count[bstarc]
    tmp <- exp(lin_count*ldeltabc + diff_vec*lgamma[bstarc])
    res4 <- tmp/sum(tmp)

    sum(abs(res - res2)) <= 1e-5 & sum(abs(res - res3)) <= 1e-5 & sum(abs(res - res4)) <= 1e-5
  })

  expect_true(all(bool_vec))
})
