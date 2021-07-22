context("Test estimation method")

.generate_df_test <- function(seed = 10){
  set.seed(seed)
  g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                   directed = F)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
  idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
  combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                            num_per_wave = num_per_wave, 
                                            distinct_waves = distinct_waves)
  
  df <- simulate_data_input(combn_wave_mat)
}

## .glmnet_fancy is correct

test_that(".glmnet_fancy works", {
  set.seed(10)
  n <- 100; p <- 10
  x <- matrix(abs(rnorm(n*p)), n, p)
  beta <- runif(p)
  y <- rnorm(n, mean = x %*% beta)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = F, 
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
  y <- rnorm(n, mean = x %*% beta + 1)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(abs(res$val_int) <= 1e-6)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = F, 
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
  y <- rnorm(n, mean = x %*% beta)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(any(abs(res$vec_coef) <= 1e-6))
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = T, 
                       switch_cutoff = 2, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(all(abs(res$vec_coef) >= 1e-6))
})

test_that(".glmnet_fancy can handle when there's only one variable", {
  set.seed(10)
  n <- 1000; p <- 1
  x <- matrix(abs(rnorm(n*p)), n, p)
  beta <- runif(p); beta[1:round(p/2)] <- 0
  y <- rnorm(n, mean = x %*% beta)
  
  set.seed(10)
  res <- .glmnet_fancy(x, y, weights = rep(1,n), family = "gaussian", switch = F, 
                       switch_cutoff = 10, alpha = 1, standardize = F,
                       intercept = F, cv = T, nfolds = 3, cv_choice = "lambda.1se", 
                       bool_round = T)
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("val_int", "vec_coef"))))
  expect_true(length(res$vec_coef) == 1)
})

########################

## .transform_est_matrix is correct

test_that(".transform_est_matrix works", {
  df <- .generate_df_test()
  p1 <- nrow(df$df_x); p2 <- nrow(df$df_y)
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
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
  df <- .generate_df_test()
  p1 <- nrow(df$df_x); p2 <- nrow(df$df_y)
  n <- 100
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  
  mat_x1 <- matrix(sample(c(0,1), n*p1, replace = T), n, p1)
  mat_y2 <- matrix(rpois(n*p2, lambda = 3), n, p2)
  
  res <- .estimate_g_glmnet(mat_x1, mat_y2, weights = rep(1,n), est_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_g", "vec_g"))))
  expect_true(all(dim(res$mat_g) == c(p1, p2)))
  expect_true(length(res$vec_g) == p2)
  for(i in 1:p2){
    idx <- est_options$ht_map[[as.character(i)]]
    expect_true(sum(abs(res$mat_g[-idx,i])) <= 1e-6)
  }
})
