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
                            lambda = 0,
                            lineage_future_count = lineage_future_count)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  
  # ensure we are indeed minimizing
  res2 <- .lineage_objective(cell_features = cell_features,
                             cell_lineage = cell_lineage,
                             cell_lineage_idx_list = cell_lineage_idx_list,
                             coefficient_vec = c(0,3,0),
                             lambda = 0,
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
    coef_vec <- runif(3)
    names(coef_vec) <- colnames(cell_features)
    lambda <- runif(1, min = 0, max = 100)
    
    obj1 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec,
                               lambda = lambda,
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
    obj2 <- obj2/length(uniq_lineage) + lambda*.l2norm(coef_vec[-1])^2
    
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
                           lambda = 0,
                           lineage_future_count = lineage_future_count)
  
  expect_true(is.numeric(res))
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
    coef_vec <- runif(3)
    names(coef_vec) <- colnames(cell_features)
    lambda <- runif(1, min = 0, max = 100)
    
    res1 <- .lineage_gradient(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              cell_lineage_idx_list = cell_lineage_idx_list,
                              coefficient_vec = coef_vec,
                              lambda = lambda,
                              lineage_future_count = lineage_future_count)
    
    res2 <- c(0,0,0)
    names(res2) <- colnames(cell_features)
    uniq_lineage <- sort(unique(cell_lineage))
    exp_vec <- sapply(1:nrow(cell_features), function(i){
      as.numeric(exp(cell_features[i,,drop = F] %*% coef_vec))
    })
    
    colname_vec <- colnames(cell_features)
    colname_vec <- colname_vec[colname_vec != "Intercept"]
    intercept_idx <- which(colnames(cell_features) == "Intercept")
    
    for(lineage in uniq_lineage){
      cell_idx_vec <- which(cell_lineage == lineage)
      
      for(cell_idx in cell_idx_vec){
        res2[intercept_idx] <- res2[intercept_idx] + exp_vec[cell_idx]
        res2[colname_vec] <- res2[colname_vec] + 
          exp_vec[cell_idx] * cell_features[cell_idx,colname_vec,drop = F]
      }
      
      n_future <- lineage_future_count[lineage]
      for(cell_idx in cell_idx_vec){
        res2[intercept_idx] <- res2[intercept_idx] - (n_future * exp_vec[cell_idx] / sum(exp_vec[cell_idx_vec]))
        res2[colname_vec] <- res2[colname_vec] - 
          (n_future * exp_vec[cell_idx] / sum(exp_vec[cell_idx_vec])) * cell_features[cell_idx,colname_vec,drop = F]
      }
    }
    
    res2 <- res2/length(uniq_lineage)
    res2[colname_vec] <- res2[colname_vec] + 2*lambda*coef_vec[colname_vec]
    
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
    coef_vec1 <- runif(3)
    coef_vec2 <- runif(3)
    names(coef_vec1) <- colnames(cell_features)
    names(coef_vec2) <- colnames(cell_features)
    
    obj1 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec1,
                               lambda = 0,
                               lineage_future_count = lineage_future_count)
    obj2 <- .lineage_objective(cell_features = cell_features,
                               cell_lineage = cell_lineage,
                               cell_lineage_idx_list = cell_lineage_idx_list,
                               coefficient_vec = coef_vec2,
                               lambda = 0,
                               lineage_future_count = lineage_future_count)
    grad_vec1 <- .lineage_gradient(cell_features = cell_features,
                                   cell_lineage = cell_lineage,
                                   cell_lineage_idx_list = cell_lineage_idx_list,
                                   coefficient_vec = coef_vec1,
                                   lambda = 0,
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
    coef_vec <- runif(3)
    names(coef_vec) <- colnames(cell_features)
    lambda <- runif(1, min = 0, max = 100)
    
    grad_vec1 <- .lineage_gradient(cell_features = cell_features,
                                   cell_lineage = cell_lineage,
                                   cell_lineage_idx_list = cell_lineage_idx_list,
                                   coefficient_vec = coef_vec,
                                   lambda = lambda,
                                   lineage_future_count = lineage_future_count)
    grad_vec2 <- numDeriv::grad(.lineage_objective, 
                                coef_vec, 
                                side = NULL,
                                cell_features = cell_features,
                                cell_lineage = cell_lineage,
                                cell_lineage_idx_list = cell_lineage_idx_list,
                                lambda = lambda,
                                lineage_future_count = lineage_future_count)
    
    sum(abs(grad_vec1 - grad_vec2)) <= 1e-3
  })
  
  expect_true(all(bool_vec))
})


test_that(".lineage_gradient seems sensible in 1-dimension", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    res <- .construct_lineage_data(coefficient_vec = 1, p = 1, seed = trial)
    cell_features <- res$cell_features
    cell_lineage <- res$cell_lineage
    cell_lineage_idx_list <- res$cell_lineage_idx_list
    lineage_future_count <- res$lineage_future_count
    
    set.seed(trial)
    coef_target <- runif(2)
    names(coef_target) <- colnames(cell_features)
    
    obj_val <- .lineage_objective(cell_features = cell_features,
                                  cell_lineage = cell_lineage,
                                  cell_lineage_idx_list = cell_lineage_idx_list,
                                  coefficient_vec = coef_target,
                                  lambda = 0,
                                  lineage_future_count = lineage_future_count)
    grad_val <- .lineage_gradient(cell_features = cell_features,
                                  cell_lineage = cell_lineage,
                                  cell_lineage_idx_list = cell_lineage_idx_list,
                                  coefficient_vec = coef_target,
                                  lambda = 0,
                                  lineage_future_count = lineage_future_count)
    
    coef_jitter <- cbind(coef_target[1], coef_target[2] + seq(-0.01,0.01,by=0.001))
    obj_vec <- sapply(1:nrow(coef_jitter), function(kk){
      .lineage_objective(cell_features = cell_features,
                         cell_lineage = cell_lineage,
                         cell_lineage_idx_list = cell_lineage_idx_list,
                         coefficient_vec = coef_jitter[kk,],
                         lambda = 0,
                         lineage_future_count = lineage_future_count)
    })
    
    lower_bound_vec <- sapply(1:nrow(coef_jitter), function(kk){
      obj_val - grad_val %*% (coef_target-coef_jitter[kk,])
    })
    
    ## you want to uncomment this line and make coef_jitter wider to see the non-convexity
    # plot(coef_jitter[,2], obj_vec); points(coef_target, obj_val, pch = 16, col = "red"); lines(coef_jitter[,2], lower_bound_vec, col = "red")
    all(lower_bound_vec <= obj_vec + 1) # just for jitter
  })
  
  expect_true(all(bool_vec))
})

############################

## lineage_imputation is correct

test_that("lineage_imputation works", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  true_coefficient <- res$coefficient_vec
  coefficient_initial <- true_coefficient/2
  lineage_future_count <- res$lineage_future_count
  
  res <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            coefficient_initial_list = coefficient_initial,
                            lineage_future_count = lineage_future_count,
                            lambda = 0,
                            verbose = 0)
  
  expect_true(is.list(res))
})

test_that("lineage_imputation works can set all non-intercept terms to 0 for large enough lambda", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  cell_lineage_idx_list <- res$cell_lineage_idx_list
  true_coefficient <- res$coefficient_vec
  coefficient_initial <- true_coefficient/2
  lineage_future_count <- res$lineage_future_count
  
  res <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = cell_lineage,
                            coefficient_initial_list = coefficient_initial,
                            lineage_future_count = lineage_future_count,
                            lambda = 1e8,
                            verbose = 0)
  
  expect_true(all(abs(res$fit$coefficient_vec[-1]) <= 1e-4))
})
