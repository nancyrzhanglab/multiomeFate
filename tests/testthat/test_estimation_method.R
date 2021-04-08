context("Test estimation method")

## .glmnet_fancy is correct

test_that(".glmnet_fancy works", {
  set.seed(10)
  n <- 100; p <- 10
  x <- matrix(abs(rnorm(n*p)), n, p)
  beta <- runif(p)
  y <- rpois(n, lambda = x %*% beta)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, family = "poisson", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(all(sort(names(res)) == c("val_int", "vec_coef")))
  expect_true(all(sapply(res, is.numeric)))
  expect_true(length(res$val_int) == 1)
  expect_true(length(res$vec_coef) == p)
})

test_that(".glmnet_fancy respects intercept", {
  set.seed(10)
  n <- 100; p <- 10
  x <- matrix(abs(rnorm(n*p)), n, p)
  beta <- runif(p)
  y <- rpois(n, lambda = x %*% beta + 1)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, family = "poisson", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(abs(res$val_int) <= 1e-6)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, family = "poisson", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = T, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(abs(res$val_int) >= 1e-6)
})

test_that(".glmnet_fancy respects switch", {
  set.seed(10)
  n <- 1000; p <- 10
  x <- matrix(abs(rnorm(n*p)), n, p)
  beta <- runif(p); beta[1:round(p/2)] <- 0
  y <- rpois(n, lambda = x %*% beta)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, family = "poisson", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(any(abs(res$vec_coef) <= 1e-6))
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, family = "poisson", switch = T, 
                       switch_cutoff = 2, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(all(abs(res$vec_coef) >= 1e-6))
})

########################

## .transform_est_matrix is correct

test_that(".transform_est_matrix works", {
  set.seed(10)
  p1 <- 20; p2 <- 5; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly", rec_method = "nn_yonly",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  
  list_res <- lapply(1:p2, function(i){
    len <- length(est_options$ht_map[[as.character(i)]])
    list(val_int = runif(1), vec_coef = runif(len))
  })
  
  res <- .transform_est_matrix(list_res, est_options, p1)
  
  expect_true(all(sort(names(res)) == sort(c("mat_g", "vec_g"))))
  expect_true(all(dim(res$mat_g) == c(p1, p2)))
  expect_true(length(res$vec_g) == p2)
  for(i in 1:p2){
    idx <- est_options$ht_map[[as.character(i)]]
    expect_true(sum(abs(res$mat_g[-idx,i])) <= 1e-6)
  }
})

#333333################3

## .estimate_g_glmnet is correct

test_that(".estimate_g_glmnet works", {
  set.seed(11)
  p1 <- 25; p2 <- 6; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  n <- 100
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly", rec_method = "nn_yonly",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  
  mat_x1 <- matrix(sample(c(0,1), n*p1, replace = T), n, p1)
  mat_y2 <- matrix(rpois(n*p2, lambda = 3), n, p2)
  
  res <- .estimate_g_glmnet(mat_x1, mat_y2, df$df_y, est_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_g", "vec_g"))))
  expect_true(all(dim(res$mat_g) == c(p1, p2)))
  expect_true(length(res$vec_g) == p2)
  for(i in 1:p2){
    idx <- est_options$ht_map[[as.character(i)]]
    expect_true(sum(abs(res$mat_g[-idx,i])) <= 1e-6)
  }
})
