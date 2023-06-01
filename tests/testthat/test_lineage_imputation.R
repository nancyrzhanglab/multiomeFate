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

test_that(".lineage_gradient is equivalent to the long-form calculation", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    coef_vec <- runif(2)
    
    res1 <- .lineage_gradient(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              cell_lineage_idx_list = cell_lineage_idx_list,
                              coefficient_vec = coef_vec,
                              lineage_future_count = lineage_future_count)
    
    res2 <- c(0,0)
    uniq_lineage <- sort(unique(cell_lineage))
    exp_vec <- sapply(1:nrow(cell_features), function(i){
      as.numeric(exp(cell_features[i,,drop = F] %*% coef_vec))
    })
    
    for(lineage in uniq_lineage){
      cell_idx_vec <- which(cell_lineage == lineage)
      
      for(cell_idx in cell_idx_vec){
        res2 <- res2 + exp_vec[cell_idx] * cell_features[cell_idx,,drop = F]
      }
      
      n_future <- lineage_future_count[lineage]
      for(cell_idx in cell_idx_vec){
        res2 <- res2 - (n_future * exp_vec[cell_idx] / sum(exp_vec[cell_idx_vec])) * cell_features[cell_idx,,drop = F]
      }
    }
    
    abs(sum(res1 - res2) <= 1e-6)
    
  })
  
  expect_true(all(bool_vec))
})

## tests inspired by https://github.com/linnykos/permanent_notes/blob/master/convex_optimization/nonconvex-scribed.pdf
test_that(".lineage_gradient has the correct mathematical property (for 1-cell-per-lineage, where it's convex)", {
  res <- .construct_lineage_data(n_each = 1)
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  trials <- 1000
  
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
    
    obj2 + 1e-6 >= obj1 + as.numeric(grad_vec1 %*% (coef_vec2 - coef_vec1))
  })
  
  expect_true(all(bool_vec))
})

test_that(".lineage_gradient matches the automatic differentiator", {
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  lineage_future_count <- res$lineage_future_count
  trials <- 1000
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    coef_vec <- runif(2)
    
    grad_vec1 <- .lineage_gradient(cell_features = cell_features,
                                   cell_lineage = cell_lineage,
                                   cell_lineage_idx_list = cell_lineage_idx_list,
                                   coefficient_vec = coef_vec,
                                   lineage_future_count = lineage_future_count)
    grad_vec2 <- numDeriv::grad(.lineage_objective, 
                                coef_vec, 
                                side = NULL,
                                cell_features = cell_features,
                                cell_lineage = cell_lineage,
                                cell_lineage_idx_list = cell_lineage_idx_list,
                                lineage_future_count = lineage_future_count)
    
    sum(abs(grad_vec1 - grad_vec2)) <= 1e-3
  })
  
  expect_true(all(bool_vec))
})


test_that(".lineage_gradient seems sensible in 1-dimension", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    res <- .construct_lineage_data(p = 1, seed = trial)
    cell_features <- res$cell_features
    cell_lineage <- res$cell_lineage
    cell_lineage_idx_list <- res$cell_lineage_idx_list
    lineage_future_count <- res$lineage_future_count
    
    set.seed(trial)
    coef_target <- runif(1)
    
    obj_val <- .lineage_objective(cell_features = cell_features,
                                  cell_lineage = cell_lineage,
                                  cell_lineage_idx_list = cell_lineage_idx_list,
                                  coefficient_vec = coef_target,
                                  lineage_future_count = lineage_future_count)
    grad_val <- .lineage_gradient(cell_features = cell_features,
                                  cell_lineage = cell_lineage,
                                  cell_lineage_idx_list = cell_lineage_idx_list,
                                  coefficient_vec = coef_target,
                                  lineage_future_count = lineage_future_count)
    
    coef_jitter <- coef_target + seq(-0.01,0.01,by=0.001)
    obj_vec <- sapply(coef_jitter, function(coef_val){
      .lineage_objective(cell_features = cell_features,
                         cell_lineage = cell_lineage,
                         cell_lineage_idx_list = cell_lineage_idx_list,
                         coefficient_vec = coef_val,
                         lineage_future_count = lineage_future_count)
    })
    
    lower_bound_vec <- sapply(coef_jitter, function(coef_val){
      obj_val - grad_val*(coef_target-coef_val)
    })
    
    ## you want to uncomment this line and make coef_jitter wider to see the non-convexity
    # plot(coef_jitter, obj_vec); points(coef_target, obj_val, pch = 16, col = "red"); lines(coef_jitter, lower_bound_vec, col = "red")
    all(lower_bound_vec <= obj_vec + 1) # just for jitter
  })
  
  expect_true(all(bool_vec))
})

############################

## lineage_imputation is correct

test_that("lineage_imputation works", {
  set.seed(10)
  n_each <- 30
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  true_coefficient <- res$coefficient_vec
  coefficient_initial <- true_coefficient/2
  lineage_future_count <- res$lineage_future_count
  
  res <- lineage_imputation(cell_features,
                            cell_lineage,
                            coefficient_initial,
                            lineage_future_count,
                            verbose = 0)
})