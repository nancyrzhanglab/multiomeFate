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

######################################
## Factor-valued cell_lineage
##
## `cell_lineage` is documented as "character or factor". Indexing a *named*
## vector with a factor (`lineage_future_count[cell_lineage]`) silently uses the
## factor's integer codes rather than its labels. On a full-data fit the level
## order happens to coincide with the sorted count vector so the result is
## correct by accident; inside a CV fold the factor retains its unused levels
## while lineage_future_count has been trimmed, and the lookup returns NA.

test_that(".lineage_gradient gives the same answer for factor and character cell_lineage", {
  set.seed(10)
  res <- .construct_lineage_data()

  grad_of <- function(cell_lineage) {
    tmp <- .lineage_cleanup(cell_features = res$cell_features,
                            cell_lineage = cell_lineage,
                            lineage_future_count = res$lineage_future_count)
    .lineage_gradient(cell_features = tmp$cell_features,
                      cell_lineage = tmp$cell_lineage,
                      cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                      coefficient_vec = res$coefficient_vec,
                      lambda = 0,
                      lineage_future_count = tmp$lineage_future_count)
  }

  expect_equal(grad_of(factor(res$cell_lineage)), grad_of(res$cell_lineage))
})

test_that(".lineage_gradient is finite when cell_lineage is a factor with unused levels", {
  # This is the cross-validation situation: subsetting a factor never drops
  # levels, so the training factor still carries the held-out lineages.
  set.seed(10)
  res <- .construct_lineage_data()

  held_out <- c("lin:1", "lin:2", "lin:3")
  keep_idx <- which(!res$cell_lineage %in% held_out)
  cell_lineage_train <- factor(res$cell_lineage)[keep_idx]
  lineage_future_count_train <-
    res$lineage_future_count[!names(res$lineage_future_count) %in% held_out]

  expect_true(all(held_out %in% levels(cell_lineage_train)))  # levels really are retained

  tmp <- .lineage_cleanup(cell_features = res$cell_features[keep_idx, , drop = FALSE],
                          cell_lineage = cell_lineage_train,
                          lineage_future_count = lineage_future_count_train)
  grad <- .lineage_gradient(cell_features = tmp$cell_features,
                            cell_lineage = tmp$cell_lineage,
                            cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                            coefficient_vec = res$coefficient_vec,
                            lambda = 0,
                            lineage_future_count = tmp$lineage_future_count)

  expect_false(anyNA(grad))
  expect_true(all(is.finite(grad)))

  # and it must equal the character-input answer
  tmp_chr <- .lineage_cleanup(cell_features = res$cell_features[keep_idx, , drop = FALSE],
                              cell_lineage = res$cell_lineage[keep_idx],
                              lineage_future_count = lineage_future_count_train)
  grad_chr <- .lineage_gradient(cell_features = tmp_chr$cell_features,
                                cell_lineage = tmp_chr$cell_lineage,
                                cell_lineage_idx_list = tmp_chr$cell_lineage_idx_list,
                                coefficient_vec = res$coefficient_vec,
                                lambda = 0,
                                lineage_future_count = tmp_chr$lineage_future_count)
  expect_equal(grad, grad_chr)
})

test_that("lineage_imputation actually moves off its starting point with a factor cell_lineage", {
  set.seed(10)
  res <- .construct_lineage_data()

  held_out <- c("lin:1", "lin:2", "lin:3")
  keep_idx <- which(!res$cell_lineage %in% held_out)
  lineage_future_count_train <-
    res$lineage_future_count[!names(res$lineage_future_count) %in% held_out]
  coefficient_initial <- res$coefficient_vec/2

  fit <- lineage_imputation(cell_features = res$cell_features[keep_idx, , drop = FALSE],
                            cell_lineage = factor(res$cell_lineage)[keep_idx],
                            coefficient_initial_list = coefficient_initial,
                            lineage_future_count = lineage_future_count_train,
                            lambda = 0,
                            random_initializations = 0,
                            verbose = 0)

  expect_false(isTRUE(all.equal(unname(fit$fit$coefficient_vec),
                                unname(fit$fit$coefficient_initial))))
})

test_that("lineage_imputation gives the same fit for factor and character cell_lineage", {
  set.seed(10)
  res <- .construct_lineage_data()
  coefficient_initial <- res$coefficient_vec/2

  fit_of <- function(cell_lineage) {
    set.seed(1)
    lineage_imputation(cell_features = res$cell_features,
                       cell_lineage = cell_lineage,
                       coefficient_initial_list = coefficient_initial,
                       lineage_future_count = res$lineage_future_count,
                       lambda = 0.1,
                       random_initializations = 2,
                       verbose = 0)$fit$coefficient_vec
  }

  expect_equal(fit_of(factor(res$cell_lineage)), fit_of(res$cell_lineage))
})

######################################
## .lineage_cleanup reconciles partial mismatches
##
## The guard was `all(a != b)`, which is TRUE only when *every* element differs.
## A partial mismatch -- the case the branch exists to handle -- slipped through.

test_that(".lineage_cleanup drops lineages missing from lineage_future_count", {
  set.seed(10)
  res <- .construct_lineage_data()

  # remove one lineage from the counts; its cells should be dropped
  lineage_future_count <- res$lineage_future_count[names(res$lineage_future_count) != "lin:1"]
  n_dropped <- sum(res$cell_lineage == "lin:1")
  expect_gt(n_dropped, 0)

  tmp <- .lineage_cleanup(cell_features = res$cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = lineage_future_count)

  expect_false("lin:1" %in% tmp$cell_lineage)
  expect_equal(length(tmp$cell_lineage), length(res$cell_lineage) - n_dropped)
  expect_equal(nrow(tmp$cell_features), length(res$cell_lineage) - n_dropped)
  expect_setequal(names(tmp$lineage_future_count), unique(tmp$cell_lineage))
})

test_that(".lineage_cleanup drops counts for lineages absent from cell_lineage", {
  set.seed(10)
  res <- .construct_lineage_data()

  lineage_future_count <- c(res$lineage_future_count, "lin:999" = 5)

  tmp <- .lineage_cleanup(cell_features = res$cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = lineage_future_count)

  expect_false("lin:999" %in% names(tmp$lineage_future_count))
  expect_setequal(names(tmp$lineage_future_count), unique(tmp$cell_lineage))
})

test_that(".lineage_gradient errors rather than returning an NA gradient", {
  # A misalignment between cell_lineage and lineage_future_count used to yield a
  # NaN gradient, which BFGS cannot use: optim() returned its starting value
  # while still reporting convergence = 0.
  set.seed(10)
  res <- .construct_lineage_data()

  tmp <- .lineage_cleanup(cell_features = res$cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = res$lineage_future_count)
  # drop one lineage's count, leaving its cells behind
  broken_count <- tmp$lineage_future_count[names(tmp$lineage_future_count) != "lin:1"]

  expect_error(
    .lineage_gradient(cell_features = tmp$cell_features,
                      cell_lineage = tmp$cell_lineage,
                      cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                      coefficient_vec = res$coefficient_vec,
                      lambda = 0,
                      lineage_future_count = broken_count),
    "misaligned"
  )
})

test_that(".lineage_gradient distinguishes numerical degeneracy from misalignment", {
  # Underflow makes every scalar1 zero, so denom_vec is zero and 0/0 is NaN --
  # with inputs that are perfectly aligned. The error must not blame the inputs.
  set.seed(10)
  res <- .construct_lineage_data()

  tmp <- .lineage_cleanup(cell_features = res$cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = res$lineage_future_count)
  coefficient_vec <- res$coefficient_vec
  coefficient_vec[] <- -1000

  expect_true(all(tmp$cell_lineage %in% names(tmp$lineage_future_count)))
  expect_error(
    .lineage_gradient(cell_features = tmp$cell_features,
                      cell_lineage = tmp$cell_lineage,
                      cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                      coefficient_vec = coefficient_vec,
                      lambda = 0,
                      lineage_future_count = tmp$lineage_future_count),
    "overflow|underflow|non-finite"
  )
})
