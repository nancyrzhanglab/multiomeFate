context("Test lineage imputation sequence")

## .compute_initial_parameters is correct

test_that(".compute_initial_parameters works", {
  set.seed(10)
  n_each <- 30
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  
  res <- .compute_initial_parameters(cell_features = cell_features,
                                     cell_lineage = cell_lineage,
                                     lineage_future_count = lineage_future_count,
                                     multipler = 10)
  
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("coefficient_initial", "lambda_initial")))
  expect_true(res$lambda_initial > 0)
  expect_true(all(names(res$coefficient_initial) == colnames(cell_features)))
})

## lineage_imputation_sequence is correct

test_that("lineage_imputation_sequence works", {
  trials <- 10
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    
    res <- .construct_lineage_data()
    cell_features <- res$cell_features
    cell_lineage <- res$cell_lineage
    cell_lineage_idx_list <- res$cell_lineage_idx_list
    lineage_future_count <- res$lineage_future_count
    
    res <- lineage_imputation_sequence(cell_features = cell_features,
                                       cell_lineage = cell_lineage,
                                       lineage_future_count = lineage_future_count,
                                       lambda_sequence_length = 25,
                                       verbose = 0)
    
    coef_mat <- sapply(res$fit_list, function(x){x$coefficient_vec})
    
    bool1 <- any(!is.na(coef_mat))
    bool2 <- all(dim(coef_mat) == c(3, 25))
    bool3 <- is.list(res)
    
    all(c(bool1, bool2, bool3))
  })
  
  expect_true(all(bool_vec))
})

######################################
## Auto-computed lambda_initial
##
## lambda_max and lambda_min both defaulted to 101, so
## `min(max(lambda_initial, lambda_max), lambda_min)` collapsed to the constant
## 101 and discarded the data-driven heuristic on the line above.

# Characterization test only -- this one also passes against the pre-fix code,
# since the old min(max(x, 8), 2) == 2 happens to land inside [2, 8]. The
# "does not collapse to a constant" test below is the actual regression guard.
test_that(".compute_initial_parameters respects lambda_min and lambda_max", {
  set.seed(10)
  res <- .construct_lineage_data()

  out <- .compute_initial_parameters(cell_features = res$cell_features,
                                     cell_lineage = res$cell_lineage,
                                     lineage_future_count = res$lineage_future_count,
                                     lambda_min = 2,
                                     lambda_max = 8,
                                     multipler = 10)

  expect_gte(out$lambda_initial, 2)
  expect_lte(out$lambda_initial, 8)
})

test_that(".compute_initial_parameters does not collapse to a constant", {
  set.seed(10)
  res <- .construct_lineage_data()

  # A wide bracket must let the data-driven value through rather than pinning
  # to an endpoint.
  low <- .compute_initial_parameters(cell_features = res$cell_features,
                                     cell_lineage = res$cell_lineage,
                                     lineage_future_count = res$lineage_future_count,
                                     lambda_min = 0, lambda_max = 1e6,
                                     multipler = 1)
  high <- .compute_initial_parameters(cell_features = res$cell_features,
                                      cell_lineage = res$cell_lineage,
                                      lineage_future_count = res$lineage_future_count,
                                      lambda_min = 0, lambda_max = 1e6,
                                      multipler = 100)

  expect_false(isTRUE(all.equal(low$lambda_initial, high$lambda_initial)))
})
