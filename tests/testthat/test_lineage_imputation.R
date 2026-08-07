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

##############################################################################
## Test lineage_imputation (test-plan 2026-08-06)
##############################################################################


# `.construct_lineage_data()` prepends an Intercept column; strip it so the
# fixtures match what a user actually passes.
.raw_features <- function(res){
  res$cell_features[, setdiff(colnames(res$cell_features), "Intercept"),
                    drop = FALSE]
}

.zero_init <- function(cell_features){
  stats::setNames(rep(0, ncol(cell_features)), colnames(cell_features))
}

.random_inits <- function(fit, num_supplied){
  idx_vec <- seq.int(num_supplied + 1, length(fit$res_list))
  do.call(rbind, lapply(fit$res_list[idx_vec],
                        function(lis){lis$coefficient_initial}))
}

# See plan item P2 -- lineage names whose sort order is the reverse of their
# appearance order, so positional indexing cannot masquerade as name indexing.
.scramble_lineage_names <- function(res){
  uniq_lineages <- unique(res$cell_lineage)
  new_names <- paste0("z", sprintf("%03d", rev(seq_along(uniq_lineages))))
  name_map <- stats::setNames(new_names, uniq_lineages)

  res$cell_lineage <- unname(name_map[res$cell_lineage])
  names(res$lineage_future_count) <-
    unname(name_map[names(res$lineage_future_count)])
  res
}

######################################
## Section 1 -- the optimizer does what it claims

## Plan item A3. `lineage_imputation()` has only ever been asserted to return a
## list. The 1.0.2.000 symptom was an optimizer that never moved while still
## reporting convergence = 0, so the load-bearing property is that each restart
## ends no worse than where it started, measured on the *penalized* objective
## that `optim` actually minimizes.
test_that("lineage_imputation improves on every supplied initialization", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  colname_vec <- colnames(cell_features)
  lambda <- 0.5

  init_list <- list(.zero_init(cell_features),
                    res$coefficient_vec[colname_vec]/2,
                    -res$coefficient_vec[colname_vec])

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = init_list,
                            lineage_future_count = res$lineage_future_count,
                            lambda = lambda,
                            random_initializations = 0,
                            verbose = 0)

  for(i in seq_along(init_list)){
    objective_initial <- evaluate_nll(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      coefficient_vec = fit$res_list[[i]]$coefficient_initial,
      lineage_future_count = res$lineage_future_count,
      lambda = lambda
    )

    expect_lte(fit$res_list[[i]]$objective_val,
               objective_initial + 1e-8,
               label = paste0("initialization ", i))
  }
})

## Plan item A4.
test_that("lineage_imputation reports the best restart and keeps them all", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  num_random <- 4

  init_list <- list(.zero_init(cell_features),
                    res$coefficient_vec[colnames(cell_features)]/2)

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = init_list,
                            lineage_future_count = res$lineage_future_count,
                            lambda = 0.1,
                            random_initializations = num_random,
                            verbose = 0)

  expect_equal(length(fit$res_list), length(init_list) + num_random)

  objective_vec <- sapply(fit$res_list, function(lis){lis$objective_val})
  expect_equal(fit$fit$objective_val, min(objective_vec))
  expect_equal(fit$fit, fit$res_list[[which.min(objective_vec)]])
})

## Plan item A5. The ridge penalty is on the non-intercept coefficients only, so
## their norm must be weakly decreasing in lambda. A path-direction inversion --
## the class of bug that made `cyfer()` always select the largest lambda --
## shows up here as a norm that grows.
test_that("lineage_imputation shrinks non-intercept coefficients as lambda grows", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lambda_vec <- c(0, 0.01, 0.1, 1, 10, 100)

  norm_vec <- sapply(lambda_vec, function(lambda){
    set.seed(1)
    fit <- lineage_imputation(cell_features = cell_features,
                              cell_lineage = res$cell_lineage,
                              coefficient_initial_list =
                                .zero_init(cell_features),
                              lineage_future_count = res$lineage_future_count,
                              lambda = lambda,
                              random_initializations = 5,
                              verbose = 0)
    coefficient_vec <- fit$fit$coefficient_vec
    idx_notintercept <- which(names(coefficient_vec) != "Intercept")
    .l2norm(coefficient_vec[idx_notintercept])
  })

  expect_true(all(diff(norm_vec) <= 1e-6),
              info = paste0("norms: ", paste0(round(norm_vec, 6),
                                              collapse = ", ")))
})

######################################
## Section 2 -- API surface

## Plan item B1. The list-of-vectors form of `coefficient_initial_list` is
## documented but never exercised; only bare vectors appear in the suite.
test_that("lineage_imputation accepts a list of initializations", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  colname_vec <- colnames(cell_features)

  init_list <- list(.zero_init(cell_features),
                    res$coefficient_vec[colname_vec]/2,
                    -res$coefficient_vec[colname_vec])

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = init_list,
                            lineage_future_count = res$lineage_future_count,
                            lambda = 0,
                            random_initializations = 0,
                            verbose = 0)

  expect_equal(length(fit$res_list), length(init_list))

  for(i in seq_along(init_list)){
    expect_equal(fit$res_list[[i]]$coefficient_initial[colname_vec],
                 init_list[[i]],
                 label = paste0("initialization ", i))
    expect_equal(unname(fit$res_list[[i]]$coefficient_initial["Intercept"]), 0)
  }

  objective_vec <- sapply(fit$res_list, function(lis){lis$objective_val})
  expect_equal(fit$fit$objective_val, min(objective_vec))
})

## Plan item B2.
test_that(".append_intercept_term prepends an Intercept only when absent", {
  vec_without <- c("p:1" = 1, "p:2" = 2)
  vec_with <- c("Intercept" = 5, "p:1" = 1, "p:2" = 2)

  res_list <- .append_intercept_term(list(vec_without, vec_with))

  expect_equal(length(res_list), 2)
  expect_equal(names(res_list[[1]]), c("Intercept", "p:1", "p:2"))
  expect_equal(unname(res_list[[1]]["Intercept"]), 0)
  expect_equal(res_list[[1]][c("p:1", "p:2")], vec_without)

  # an existing Intercept is left exactly as supplied
  expect_equal(res_list[[2]], vec_with)

  expect_error(.append_intercept_term(vec_without))
})

## Plan item B3.
test_that("lineage_imputation returns an object of class lineage_imputation", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = .zero_init(cell_features),
                            lineage_future_count = res$lineage_future_count,
                            lambda = 0,
                            random_initializations = 0,
                            verbose = 0)

  expect_s3_class(fit, "lineage_imputation")
  expect_true(all(c("fit", "res_list") %in% names(fit)))
})

## Plan item B4. The overflow guard is only ever tested by calling
## `.lineage_gradient()` directly. Users reach it through `lineage_imputation()`
## with unscaled features, which is the documented failure mode ("Scale
## `cell_features`, or use a larger lambda") -- so the diagnostic has to survive
## the trip through `optim()`.
test_that("lineage_imputation surfaces the overflow diagnostic to the caller", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res) * 1e3
  init_vec <- stats::setNames(rep(50, ncol(cell_features)),
                              colnames(cell_features))

  expect_error(
    lineage_imputation(cell_features = cell_features,
                       cell_lineage = res$cell_lineage,
                       coefficient_initial_list = init_vec,
                       lineage_future_count = res$lineage_future_count,
                       lambda = 0,
                       random_initializations = 0,
                       verbose = 0),
    "overflow|underflow|non-finite"
  )
})

######################################
## Section 3 -- the random-restart range (lineage_imputation.R:110-115)

## Plan item B5. When every lineage has as many future cells as current cells,
## max_count_ratio == 1, so max_limit == 0 and min_value == 0. Every "random"
## restart is then runif(p, 0, 0) -- the zero vector -- and the ten restarts
## that exist to escape local optima silently collapse into one.
test_that("lineage_imputation draws distinct restarts when the count ratio is 1", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  num_cells_vec <- table(res$cell_lineage)
  lineage_names <- names(res$lineage_future_count)
  lineage_future_count <- stats::setNames(
    as.numeric(num_cells_vec[lineage_names]),
    lineage_names
  )
  expect_equal(max(lineage_future_count/as.numeric(num_cells_vec[lineage_names])),
               1)

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = .zero_init(cell_features),
                            lineage_future_count = lineage_future_count,
                            lambda = 0,
                            random_initializations = 5,
                            verbose = 0)

  init_mat <- .random_inits(fit, num_supplied = 1)
  expect_gt(sum(apply(init_mat, 2, stats::sd)), 0)
})

## Plan item B5, the threshold. The test above pins `max_count_ratio == 1`
## exactly, where the restart interval has width 0 and the "are the draws
## distinct" check bites. It does NOT cover the band the `<= 1e-6` threshold was
## introduced for: at `max_count_ratio = 1 + 1e-9` the interval is width 2.9e-10
## -- numerically degenerate, but not zero, so every draw is technically
## distinct and a distinctness check passes vacuously. Assert the interval has
## real width instead, which is the property the rule actually buys.
##
## Stated as a spread, not as `min_value`/`max_limit`: both are locals inside
## `lineage_imputation()`, and pinning the arithmetic would re-implement the
## formula in the test rather than test it. The intended width is ~2, so 0.5 is
## a wide margin over 20 restarts.
.restart_spread <- function(count_ratio, num_random = 20){
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  num_cells_vec <- table(res$cell_lineage)
  lineage_names <- names(res$lineage_future_count)
  lineage_future_count <- stats::setNames(
    as.numeric(num_cells_vec[lineage_names])*count_ratio,
    lineage_names
  )

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = .zero_init(cell_features),
                            lineage_future_count = lineage_future_count,
                            lambda = 0,
                            random_initializations = num_random,
                            verbose = 0)

  .random_inits(fit, num_supplied = 1)
}

test_that("lineage_imputation explores a real range when growth is flat", {
  # exactly flat, and just inside the threshold band above it
  for(count_ratio in c(1, 1 + 1e-9)){
    label <- paste0("count ratio = ", format(count_ratio, digits = 12))
    init_mat <- .restart_spread(count_ratio)

    expect_gt(diff(range(init_mat)), 0.5)
    expect_true(all(init_mat <= 1e-6), info = label)
  }
})

test_that("lineage_imputation leaves the restart floor at zero when clones grow", {
  # The growth branch must keep min_value = 0: widening it here would start
  # seeding the optimizer with shrinkage coefficients on data that only grows.
  for(count_ratio in c(1 + 1e-4, 2, 100)){
    label <- paste0("count ratio = ", count_ratio)
    init_mat <- .restart_spread(count_ratio)

    expect_true(all(init_mat >= 0), info = label)
    expect_gt(diff(range(init_mat)), 0)
  }
})

## Plan item B6. With every future count zero, max_count_ratio == 0 and
## log(0) = -Inf propagates into runif(), which warns "NAs produced" and hands
## `optim` an NA start. The user then sees "non-finite value supplied by optim",
## which names neither the argument nor the cause.
test_that("lineage_imputation fails informatively on all-zero future counts", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  lineage_future_count <- res$lineage_future_count
  lineage_future_count[] <- 0

  warning_vec <- character(0)
  error_msg <- NA_character_

  withCallingHandlers(
    tryCatch(
      lineage_imputation(cell_features = cell_features,
                         cell_lineage = res$cell_lineage,
                         coefficient_initial_list = .zero_init(cell_features),
                         lineage_future_count = lineage_future_count,
                         lambda = 0,
                         random_initializations = 5,
                         verbose = 0),
      error = function(e){error_msg <<- conditionMessage(e)}
    ),
    warning = function(w){
      warning_vec <<- c(warning_vec, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(any(grepl("NAs produced", warning_vec)),
               label = "the random-restart draw emits an NA warning")
  expect_false(isTRUE(grepl("non-finite value supplied by optim", error_msg)),
               label = "the failure is reported by optim rather than by us")
})

## Plan item B7. Shrinking clones (max_count_ratio < 1) is the drug-bottleneck
## regime and covers most of the paper's real data. Kevin: this is the intended
## behavior -- the restart range should be entirely negative there.
test_that("lineage_imputation draws negative restarts when clones shrink", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  num_cells_vec <- table(res$cell_lineage)
  lineage_names <- names(res$lineage_future_count)
  # every clone loses half its cells, so the largest ratio is 0.5 < 1
  lineage_future_count <- stats::setNames(
    as.numeric(num_cells_vec[lineage_names])/2,
    lineage_names
  )

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = .zero_init(cell_features),
                            lineage_future_count = lineage_future_count,
                            lambda = 0,
                            random_initializations = 5,
                            verbose = 0)

  init_mat <- .random_inits(fit, num_supplied = 1)
  expect_true(all(init_mat < 0))
})

## Plan item B8. `upper_randomness` caps the restart draw, but no test has ever
## made the cap bind, so an inverted or dropped `pmin()` would go unnoticed.
test_that("lineage_imputation honours the upper_randomness cap", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  lineage_names <- names(res$lineage_future_count)
  # explosive growth, so max_limit is comfortably above the cap below
  lineage_future_count <- stats::setNames(rep(1e6, length(lineage_names)),
                                          lineage_names)
  upper_randomness <- 1e-3

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = .zero_init(cell_features),
                            lineage_future_count = lineage_future_count,
                            lambda = 0,
                            random_initializations = 5,
                            upper_randomness = upper_randomness,
                            verbose = 0)

  init_mat <- .random_inits(fit, num_supplied = 1)
  expect_true(all(init_mat <= upper_randomness))
  expect_true(any(init_mat == upper_randomness))
})

######################################
## Section 4 -- invariances (plan items F1, F2)
##
## Both historical bugs -- the factor indexing in `.lineage_gradient()` and the
## mis-dealt lineage in `construct_folds()` -- were index-versus-name failures.
## These two invariances are the cheapest way to keep that whole class out.

## Plan item F1.
test_that("lineage_imputation is invariant to the order of the cells", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  fit_of <- function(cell_features, cell_lineage){
    set.seed(3)
    lineage_imputation(cell_features = cell_features,
                       cell_lineage = cell_lineage,
                       coefficient_initial_list = .zero_init(cell_features),
                       lineage_future_count = res$lineage_future_count,
                       lambda = 0.1,
                       random_initializations = 2,
                       verbose = 0)$fit$coefficient_vec
  }

  set.seed(2)
  idx_vec <- sample(nrow(cell_features))

  expect_equal(fit_of(cell_features[idx_vec, , drop = FALSE],
                      res$cell_lineage[idx_vec]),
               fit_of(cell_features, res$cell_lineage))
})

## Plan item F2. Relabelling the lineages so that their sort order reverses must
## not move the fit: the objective is a sum over lineages, and nothing in it may
## depend on where a lineage lands once `.lineage_cleanup()` sorts the names.
test_that("lineage_imputation is invariant to relabelling the lineages", {
  set.seed(10)
  res <- .construct_lineage_data()
  res_scrambled <- .scramble_lineage_names(res)
  cell_features <- .raw_features(res)

  fit_of <- function(cell_lineage, lineage_future_count){
    set.seed(3)
    lineage_imputation(cell_features = cell_features,
                       cell_lineage = cell_lineage,
                       coefficient_initial_list = .zero_init(cell_features),
                       lineage_future_count = lineage_future_count,
                       lambda = 0.1,
                       random_initializations = 2,
                       verbose = 0)$fit$coefficient_vec
  }

  expect_equal(fit_of(res_scrambled$cell_lineage,
                      res_scrambled$lineage_future_count),
               fit_of(res$cell_lineage, res$lineage_future_count))
})

######################################
## Section 5 -- evaluate_nll (plan items G1-G4)
##
## Every `train_loglik` and `test_loglik` value that `cyfer()` reports, and
## therefore every lambda `cyfer_finalize()` selects, comes out of this
## unexported function. It has never been tested.

## Plan item G1.
test_that("evaluate_nll agrees with .lineage_objective", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  trials <- 20

  bool_vec <- sapply(seq_len(trials), function(trial){
    set.seed(trial)
    coefficient_vec <- stats::runif(ncol(res$cell_features))
    names(coefficient_vec) <- colnames(res$cell_features)
    lambda <- stats::runif(1, min = 0, max = 10)

    val1 <- evaluate_nll(
      cell_features = cell_features,
      cell_lineage = res$cell_lineage,
      coefficient_vec = coefficient_vec,
      lineage_future_count = res$lineage_future_count,
      lambda = lambda
    )

    tmp <- .lineage_cleanup(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            lineage_future_count = res$lineage_future_count)
    val2 <- .lineage_objective(cell_features = tmp$cell_features,
                               cell_lineage = tmp$cell_lineage,
                               cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                               coefficient_vec = coefficient_vec,
                               lambda = lambda,
                               lineage_future_count = tmp$lineage_future_count)

    abs(val1 - val2) <= 1e-9
  })

  expect_true(all(bool_vec))
})

## Plan item G2. `cyfer()` evaluates both the training and the held-out curves
## at lambda = 0 on purpose -- the reported value must be the unpenalized
## objective, so that curves fit at different lambdas stay comparable.
test_that("evaluate_nll at lambda = 0 is the unpenalized objective", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)

  coefficient_vec <- res$coefficient_vec
  lambda <- 3
  idx_notintercept <- which(names(coefficient_vec) != "Intercept")

  val_unpenalized <- evaluate_nll(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    coefficient_vec = coefficient_vec,
    lineage_future_count = res$lineage_future_count,
    lambda = 0
  )
  val_penalized <- evaluate_nll(
    cell_features = cell_features,
    cell_lineage = res$cell_lineage,
    coefficient_vec = coefficient_vec,
    lineage_future_count = res$lineage_future_count,
    lambda = lambda
  )

  expect_equal(val_penalized - val_unpenalized,
               lambda*.l2norm(coefficient_vec[idx_notintercept])^2)
})

## Plan item G3.
test_that("evaluate_nll is lower at the fit than at its start", {
  set.seed(10)
  res <- .construct_lineage_data()
  cell_features <- .raw_features(res)
  coefficient_initial <- .zero_init(cell_features)

  fit <- lineage_imputation(cell_features = cell_features,
                            cell_lineage = res$cell_lineage,
                            coefficient_initial_list = coefficient_initial,
                            lineage_future_count = res$lineage_future_count,
                            lambda = 0,
                            random_initializations = 0,
                            verbose = 0)

  loglik_of <- function(coefficient_vec){
    evaluate_nll(cell_features = cell_features,
                           cell_lineage = res$cell_lineage,
                           coefficient_vec = coefficient_vec,
                           lineage_future_count = res$lineage_future_count,
                           lambda = 0)
  }

  expect_lte(loglik_of(fit$fit$coefficient_vec),
             loglik_of(fit$fit$coefficient_initial) + 1e-8)
})

## Plan item G4. `.lineage_gradient()` errors on a non-finite value but
## `.lineage_objective()` does not: it returns Inf (underflow) or NaN
## (overflow). On a held-out fold that value flows straight into `test_loglik`,
## where `which.min()` skips it silently -- so a degenerate fold is
## indistinguishable from a merely poor one. Kevin: the objective should guard
## too.
test_that(".lineage_objective errors rather than returning a non-finite value", {
  set.seed(10)
  res <- .construct_lineage_data()

  tmp <- .lineage_cleanup(cell_features = res$cell_features,
                          cell_lineage = res$cell_lineage,
                          lineage_future_count = res$lineage_future_count)

  objective_at <- function(value){
    coefficient_vec <- res$coefficient_vec
    coefficient_vec[] <- value
    .lineage_objective(cell_features = tmp$cell_features,
                       cell_lineage = tmp$cell_lineage,
                       cell_lineage_idx_list = tmp$cell_lineage_idx_list,
                       coefficient_vec = coefficient_vec,
                       lambda = 0,
                       lineage_future_count = tmp$lineage_future_count)
  }

  # underflow: every exp() is 0, so log(0) = -Inf and the objective is +Inf
  expect_error(objective_at(-1000), "overflow|underflow|non-finite")
  # overflow: every exp() is Inf, so the objective is Inf - Inf = NaN
  expect_error(objective_at(1000), "overflow|underflow|non-finite")
})
