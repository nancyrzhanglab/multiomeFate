context("Test construct_folds")

# Helper: build a lineage_future_count / cell_lineage pair with `L` lineages.
.fold_fixture <- function(L, cells_each = 2, counts = NULL) {
  lineage_names <- sprintf("lin%04d", 1:L)
  if (is.null(counts)) counts <- seq(L, 1)  # distinct, descending-sortable
  lineage_future_count <- stats::setNames(counts, lineage_names)
  cell_lineage <- rep(lineage_names, each = cells_each)
  list(lineage_future_count = lineage_future_count,
       cell_lineage = cell_lineage,
       lineage_names = lineage_names)
}

## The central invariant: the folds must be a *partition* of the lineages.
## Bug B1 (loop bound `1:num_per_fold`) drops a lineage and leaves NA entries;
## bug B2 (`unique(pmin(idx, num_lineages))` clamping instead of dropping) puts the
## last lineage into every fold. Both violate this invariant.

test_that("construct_folds assigns every lineage to exactly one fold", {
  fx <- .fold_fixture(L = 100)

  set.seed(10)
  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 10)

  assigned <- unlist(res$fold_lineage_list, use.names = FALSE)

  expect_false(any(is.na(assigned)))
  expect_equal(sort(assigned), sort(fx$lineage_names))
  expect_equal(anyDuplicated(assigned), 0L)
})

test_that("construct_folds partitions lineages across a grid of sizes and fold counts", {
  grid <- expand.grid(L = c(7, 10, 12, 20, 21, 37, 40, 53, 100, 200, 416),
                      num_folds = c(2, 3, 5, 10, 20))
  # num_folds == num_lineages is the boundary case (one lineage per fold) and is
  # deliberately included; beyond it, construct_folds() errors by design.
  grid <- grid[grid$num_folds <= grid$L, ]

  for (i in seq_len(nrow(grid))) {
    L <- grid$L[i]
    k <- grid$num_folds[i]
    fx <- .fold_fixture(L = L)

    set.seed(42)
    res <- construct_folds(cell_lineage = fx$cell_lineage,
                           lineage_future_count = fx$lineage_future_count,
                           num_folds = k)
    assigned <- unlist(res$fold_lineage_list, use.names = FALSE)

    label <- paste0("L=", L, ", num_folds=", k)
    expect_false(any(is.na(assigned)), info = label)
    expect_equal(sort(assigned), sort(fx$lineage_names), info = label)
    expect_equal(anyDuplicated(assigned), 0L, info = label)
    expect_equal(length(res$fold_lineage_list), k, info = label)
  }
})

## Bug B3: `sample(x)` on a length-1 vector returns a permutation of `1:x`, not `x`.
## This fires when the final block holds exactly one lineage, and corrupts the
## ordering by recycling one lineage name across many slots.
test_that("construct_folds is correct when a shuffle block holds a single lineage", {
  # L = 7, num_folds = 3 -> num_per_fold = 3, blocks are 1-3, 4-6, 7-7
  fx <- .fold_fixture(L = 7)

  set.seed(1)
  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 3)
  assigned <- unlist(res$fold_lineage_list, use.names = FALSE)

  expect_equal(sort(assigned), sort(fx$lineage_names))
  expect_equal(anyDuplicated(assigned), 0L)
})

test_that("construct_folds keeps the smallest-count lineage in a fold", {
  # Lineages are sorted by descending future count, so the tail of the ordering
  # holds the smallest counts -- exactly where bugs B1/B2 do their damage.
  fx <- .fold_fixture(L = 40)
  smallest <- names(which.min(fx$lineage_future_count))

  set.seed(3)
  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 5)
  assigned <- unlist(res$fold_lineage_list, use.names = FALSE)

  expect_true(smallest %in% assigned)
  expect_equal(sum(assigned == smallest), 1L)
})

test_that("construct_folds partitions the cells, not just the lineages", {
  fx <- .fold_fixture(L = 40, cells_each = 3)

  set.seed(7)
  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 5)

  cell_idx <- unlist(res$cv_cell_list, use.names = FALSE)
  expect_equal(sort(cell_idx), seq_along(fx$cell_lineage))
  expect_equal(anyDuplicated(cell_idx), 0L)
})

test_that("construct_folds balances the number of lineages per fold", {
  fx <- .fold_fixture(L = 100)

  set.seed(11)
  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 10)

  fold_sizes <- sapply(res$fold_lineage_list, length)
  expect_lte(max(fold_sizes) - min(fold_sizes), 1)
})

test_that("construct_folds accepts a factor cell_lineage", {
  fx <- .fold_fixture(L = 40, cells_each = 3)

  set.seed(5)
  res_chr <- construct_folds(cell_lineage = fx$cell_lineage,
                             lineage_future_count = fx$lineage_future_count,
                             num_folds = 5)
  set.seed(5)
  res_fac <- construct_folds(cell_lineage = factor(fx$cell_lineage),
                             lineage_future_count = fx$lineage_future_count,
                             num_folds = 5)

  expect_equal(res_fac$fold_lineage_list, res_chr$fold_lineage_list)
  expect_equal(res_fac$cv_cell_list, res_chr$cv_cell_list)
})

######################################
## Input validation
##
## Removing the `pmin` clamp was correct, but that clamp was also what
## accidentally kept every fold non-empty. `num_folds > num_lineages` now
## produces empty folds, which `cyfer()` cannot consume -- so reject it here,
## where the message can name the offending argument.

test_that("construct_folds rejects more folds than lineages", {
  fx <- .fold_fixture(L = 10)

  expect_error(
    construct_folds(cell_lineage = fx$cell_lineage,
                    lineage_future_count = fx$lineage_future_count,
                    num_folds = 13),
    "num_folds"
  )
})

test_that("construct_folds allows exactly one lineage per fold", {
  fx <- .fold_fixture(L = 10)

  res <- construct_folds(cell_lineage = fx$cell_lineage,
                         lineage_future_count = fx$lineage_future_count,
                         num_folds = 10)

  expect_equal(sort(unlist(res$fold_lineage_list, use.names = FALSE)),
               sort(fx$lineage_names))
  expect_true(all(sapply(res$fold_lineage_list, length) == 1))
})

test_that("construct_folds rejects degenerate fold counts", {
  fx <- .fold_fixture(L = 10)

  for (bad in list(1, 0, -1, NA, c(2, 3))) {
    expect_error(
      construct_folds(cell_lineage = fx$cell_lineage,
                      lineage_future_count = fx$lineage_future_count,
                      num_folds = bad),
      "num_folds"
    )
  }
})

test_that("construct_folds rejects an empty or unnamed lineage_future_count", {
  expect_error(
    construct_folds(cell_lineage = character(0),
                    lineage_future_count = stats::setNames(numeric(0), character(0)),
                    num_folds = 3),
    "lineage_future_count"
  )
  expect_error(
    construct_folds(cell_lineage = rep("a", 4),
                    lineage_future_count = c(4, 2),
                    num_folds = 2),
    "lineage_future_count"
  )
})

test_that("construct_folds rejects duplicated lineage names", {
  fx <- .fold_fixture(L = 10)
  dup <- c(fx$lineage_future_count, fx$lineage_future_count[1])

  expect_error(
    construct_folds(cell_lineage = fx$cell_lineage,
                    lineage_future_count = dup,
                    num_folds = 3),
    "duplicat"
  )
})
