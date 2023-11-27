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
