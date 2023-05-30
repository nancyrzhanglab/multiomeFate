context("Test lineage imputation")

## .lineage_objective is correct

test_that(".lineage_objective works", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  coefficient_vec <- res$coefficient_vec
  lineage_future_count <- res$lineage_future_count
  
  res <- .lineage_objective(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            cell_lineage_idx_list = cell_lineage_idx_list,
                            coefficient_vec = coefficient_vec,
                            lineage_future_count = lineage_future_count)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  
  # ensure we are indeed minimizing
  res2 <- .lineage_objective(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            cell_lineage_idx_list = cell_lineage_idx_list,
                            coefficient_vec = c(3,0),
                            lineage_future_count = lineage_future_count)
  expect_true(res2 >= res) 
})

test_that(".lineage_objective is equivalent to the long-form calculation", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    coef_vec <- runif(2)
    
    obj1 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec,
                               lineage_future_count = lineage_future_count)
    
    obj2 <- 0
    uniq_lineage <- sort(unique(cell_lineage))
    for(i in 1:length(uniq_lineage)){
      lineage <- uniq_lineage[i]
      cell_idx <- which(cell_lineage == lineage)
      exp_vec <- sapply(cell_idx, function(j){
        as.numeric(exp(cell_features[j,,drop = F] %*% coef_vec))
      })
      obj2 <- obj2 + sum(exp_vec) - lineage_future_count[lineage] * log(sum(exp_vec))
    }
    
    abs(obj1 - obj2) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

######################

## .lineage_gradient is correct

test_that(".lineage_gradient works", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  coefficient_vec <- res$coefficient_vec
  lineage_future_count <- res$lineage_future_count
  
  res <- .lineage_gradient(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            cell_lineage_idx_list = cell_lineage_idx_list,
                            coefficient_vec = coefficient_vec,
                            lineage_future_count = lineage_future_count)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == colnames(cell_features)))
})

## tests inspired by https://github.com/linnykos/permanent_notes/blob/master/convex_optimization/nonconvex-scribed.pdf

test_that(".lineage_gradient works", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    coef_vec1 <- runif(2)
    coef_vec2 <- runif(2)
    
    obj1 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec1,
                               lineage_future_count = lineage_future_count)
    obj2 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec2,
                               lineage_future_count = lineage_future_count)
    grad_vec1 <- .lineage_gradient(cell_features = cell_features,
                                   cell_lineage = cell_lineage,
                                   cell_lineage_idx_list = cell_lineage_idx_list,
                                   coefficient_vec = coef_vec1,
                                   lineage_future_count = lineage_future_count)
    
    obj2 >= obj1 + as.numeric(grad_vec1 %*% (coef_vec2 - coef_vec1))
  })
  
  expect_true(all(bool_vec))
})
