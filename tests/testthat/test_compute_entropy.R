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

##############################################################################
## compute_entropy() -- the whole function (test-plan-full 2026-08-06)
##############################################################################

## `test_compute_entropy.R` covered `.shannon_entropy()` well and
## `compute_entropy()` not at all, because the latter needs a Seurat object.
## It has more branching than anything else in this section, so build one.

## `cell_imputation_mat` is one `cyfer_finalize()` fit per candidate fate,
## column-bound -- that is what makes it feed `plot_simplex()`. Scores are on
## the log10 scale, which is why `bool_10_power` defaults to TRUE.
.entropy_fixture <- function(num_cells = 60,
                             num_fates = 3,
                             seed_number = 10){
  set.seed(seed_number)
  cell_names <- paste0("cell", seq_len(num_cells))
  fate_names <- c("Monocyte", "Neutrophil", "Undifferentiated")[
    seq_len(num_fates)]

  cell_imputation_mat <- matrix(stats::runif(num_cells * num_fates,
                                             min = -0.5, max = 1.5),
                                nrow = num_cells,
                                ncol = num_fates,
                                dimnames = list(cell_names, fate_names))

  count_mat <- matrix(stats::rpois(30 * num_cells, lambda = 3),
                      nrow = 30,
                      dimnames = list(paste0("gene", seq_len(30)), cell_names))
  seurat_object <- suppressWarnings(
    Seurat::CreateSeuratObject(counts = count_mat))
  seurat_object$celltype <- rep(fate_names, length.out = num_cells)
  seurat_object$lineage <- rep(paste0("L", seq_len(6)), length.out = num_cells)
  seurat_object$timepoint <- rep(c("day0", "day7"), each = num_cells / 2)

  list(cell_imputation_mat = cell_imputation_mat,
       seurat_object = seurat_object)
}

.run_compute_entropy <- function(fixture, ...){
  compute_entropy(cell_imputation_mat = fixture$cell_imputation_mat,
                  later_timepoint = "day7",
                  seurat_object = fixture$seurat_object,
                  variable_celltype = "celltype",
                  variable_lineage = "lineage",
                  variable_timepoint = "timepoint",
                  ...)
}

## EN1
test_that("composition columns sum to 1 after jitter (EN1)", {
  fixture <- .entropy_fixture()
  fate_names <- colnames(fixture$cell_imputation_mat)

  for(bool_jitter in c(TRUE, FALSE)){
    label <- paste0("bool_jitter = ", bool_jitter)
    set.seed(10)
    res <- .run_compute_entropy(fixture, bool_jitter = bool_jitter)

    composition_mat <- as.matrix(res[, fate_names])
    expect_true(max(abs(rowSums(composition_mat) - 1)) <= 1e-10, info = label)
    expect_true(all(composition_mat >= 0), info = label)
  }
})

## EN2: `cellsize` is captured *before* normalization, so it is the only place
## the magnitude survives -- everything else in the frame is compositional.
test_that("cellsize is the pre-normalization row sum (EN2)", {
  fixture <- .entropy_fixture()

  set.seed(10)
  res <- .run_compute_entropy(fixture, bool_10_power = TRUE)
  expected_vec <- rowSums(10^fixture$cell_imputation_mat)
  expect_true(max(abs(res$cellsize - expected_vec[rownames(res)])) <= 1e-10)

  set.seed(10)
  res_raw <- .run_compute_entropy(fixture, bool_10_power = FALSE)
  expected_raw_vec <- rowSums(fixture$cell_imputation_mat)
  expect_true(max(abs(res_raw$cellsize -
                        expected_raw_vec[rownames(res_raw)])) <= 1e-10)
})

## EN3: the same exp()/10^ trap the "Scales" section of ?cyfer_finalize
## documents, now in a second place. `bool_10_power = TRUE` is correct for a
## `cell_imputed_score` (log10); `FALSE` is correct for an already-exponentiated
## matrix. Pre-exponentiating and passing FALSE must give the same answer as
## passing the log10 scores with TRUE.
test_that("bool_10_power agrees with pre-exponentiating (EN3)", {
  fixture <- .entropy_fixture()

  set.seed(10)
  res_power <- .run_compute_entropy(fixture, bool_10_power = TRUE,
                                    bool_jitter = FALSE)

  fixture_pre <- fixture
  fixture_pre$cell_imputation_mat <- 10^fixture$cell_imputation_mat
  set.seed(10)
  res_pre <- .run_compute_entropy(fixture_pre, bool_10_power = FALSE,
                                  bool_jitter = FALSE)

  fate_names <- colnames(fixture$cell_imputation_mat)
  expect_true(identical(rownames(res_power), rownames(res_pre)))
  expect_true(max(abs(as.matrix(res_power[, fate_names]) -
                        as.matrix(res_pre[, fate_names]))) <= 1e-10)
  expect_true(max(abs(res_power$cellsize - res_pre$cellsize)) <= 1e-10)
})

## EN4: the threshold is applied *after* `bool_10_power`, so it is on the count
## scale, not the log10 scale -- easy to misread as "drop cells whose score is
## below 0.01".
test_that("min_imputation drops rows on the count scale (EN4)", {
  fixture <- .entropy_fixture()
  # Force one cell far below any plausible count threshold.
  fixture$cell_imputation_mat["cell1", ] <- -10

  set.seed(10)
  res <- .run_compute_entropy(fixture, min_imputation = 0.01)
  expect_true(!("cell1" %in% rownames(res)))

  # Raising the threshold onto the count scale removes more cells; the same
  # numeric value read as a log10 score would remove none.
  count_vec <- rowSums(10^fixture$cell_imputation_mat)
  threshold_val <- stats::median(count_vec)
  set.seed(10)
  res_strict <- .run_compute_entropy(fixture, min_imputation = threshold_val)

  expect_true(nrow(res_strict) < nrow(res))
  expect_true(all(count_vec[rownames(res_strict)] > threshold_val))
})

## EN5: `entropy` and `dominant_fate` are **lineage**-level despite sitting in a
## per-cell frame -- every cell of a lineage carries the same value, computed
## from the observed cell types at `later_timepoint`. The column names suggest
## per-cell quantities and they are not.
test_that("entropy and dominant_fate are lineage-level (EN5)", {
  fixture <- .entropy_fixture()
  set.seed(10)
  res <- .run_compute_entropy(fixture)

  for(lineage in unique(res$lineage)){
    idx <- which(res$lineage == lineage)
    label <- paste0("lineage ", lineage)
    expect_true(length(unique(res$entropy[idx])) == 1, info = label)
    expect_true(length(unique(as.character(res$dominant_fate[idx]))) == 1,
                info = label)
  }

  # And the value is the entropy of the later-timepoint cell types, bumped.
  metadata <- fixture$seurat_object@meta.data
  for(lineage in unique(res$lineage)){
    seurat_idx <- intersect(which(metadata[, "timepoint"] == "day7"),
                            which(metadata[, "lineage"] == lineage))
    tab_vec <- table(metadata[seurat_idx, "celltype"])
    if(length(tab_vec) == 0) next()
    idx <- which(res$lineage == lineage)[1]
    expect_true(abs(res$entropy[idx] -
                      (.shannon_entropy(tab_vec) + 0.01)) <= 1e-10,
                info = lineage)
    expect_true(as.character(res$dominant_fate[idx]) ==
                  names(tab_vec)[which.max(tab_vec)], info = lineage)
  }
})

## EN6: the two columns disagree about how they represent "no cells at the later
## time point" -- `dominant_fate` gets the literal string "NA" coerced to a
## factor level, `entropy` gets a real NA (plus the bump, so it stays NA). That
## inconsistency is load-bearing for the plot; pin it before someone unifies it.
test_that("a lineage absent at later_timepoint splits NA two ways (EN6)", {
  fixture <- .entropy_fixture()
  # L1 exists only at day0.
  lineage_vec <- fixture$seurat_object$lineage
  timepoint_vec <- fixture$seurat_object$timepoint
  lineage_vec[timepoint_vec == "day7" & lineage_vec == "L1"] <- "L2"
  fixture$seurat_object$lineage <- lineage_vec

  set.seed(10)
  res <- .run_compute_entropy(fixture)
  idx <- which(res$lineage == "L1")

  expect_true(length(idx) > 0)
  expect_true(all(as.character(res$dominant_fate[idx]) == "NA"))
  expect_true(!any(is.na(as.character(res$dominant_fate[idx]))))
  expect_true("NA" %in% levels(res$dominant_fate))
  expect_true(all(is.na(res$entropy[idx])))
})

## EN7: `entropy_bump` is added to every entropy, so a single-fate lineage
## reports 0.01 rather than 0. Cosmetic -- entropy drives point size -- but it
## means the column is not a Shannon entropy, and a downstream analysis that
## treated it as one would be slightly wrong everywhere.
test_that("entropy_bump means the column is not a Shannon entropy (EN7)", {
  fixture <- .entropy_fixture()
  # Make L1 single-fate at day7, so its true Shannon entropy is exactly 0.
  celltype_vec <- fixture$seurat_object$celltype
  lineage_vec <- fixture$seurat_object$lineage
  timepoint_vec <- fixture$seurat_object$timepoint
  celltype_vec[lineage_vec == "L1" & timepoint_vec == "day7"] <- "Monocyte"
  fixture$seurat_object$celltype <- celltype_vec

  set.seed(10)
  res <- .run_compute_entropy(fixture, entropy_bump = 0.01)
  idx <- which(res$lineage == "L1")
  expect_true(abs(res$entropy[idx[1]] - 0.01) <= 1e-10)

  set.seed(10)
  res_nobump <- .run_compute_entropy(fixture, entropy_bump = 0)
  idx <- which(res_nobump$lineage == "L1")
  expect_true(abs(res_nobump$entropy[idx[1]]) <= 1e-10)
})

## EN8
test_that("bool_jitter = FALSE is deterministic, TRUE is not (EN8)", {
  fixture <- .entropy_fixture()
  fate_names <- colnames(fixture$cell_imputation_mat)

  # No seed at all: the un-jittered path must still be reproducible.
  res_a <- .run_compute_entropy(fixture, bool_jitter = FALSE)
  res_b <- .run_compute_entropy(fixture, bool_jitter = FALSE)
  expect_true(identical(res_a, res_b))

  set.seed(1)
  res_c <- .run_compute_entropy(fixture, bool_jitter = TRUE)
  set.seed(2)
  res_d <- .run_compute_entropy(fixture, bool_jitter = TRUE)
  expect_true(max(abs(as.matrix(res_c[, fate_names]) -
                        as.matrix(res_d[, fate_names]))) > 1e-6)

  # There is no seed_number argument -- seeding is the caller's job.
  set.seed(1)
  res_e <- .run_compute_entropy(fixture, bool_jitter = TRUE)
  expect_true(max(abs(as.matrix(res_c[, fate_names]) -
                        as.matrix(res_e[, fate_names]))) <= 1e-10)
  expect_true(!("seed_number" %in% names(formals(compute_entropy))))
})

## EN10
test_that("a single surviving cell stays a matrix (EN10)", {
  fixture <- .entropy_fixture()
  # Only cell1 clears the threshold.
  fixture$cell_imputation_mat[, ] <- -10
  fixture$cell_imputation_mat["cell1", ] <- 1

  set.seed(10)
  res <- .run_compute_entropy(fixture, min_imputation = 0.5)

  expect_true(nrow(res) == 1)
  expect_true(rownames(res) == "cell1")
  expect_true(!is.na(res$celltype[1]))
})
