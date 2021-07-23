context("Test candidate method")

## .candidate_set_nn_any is correct

test_that(".candidate_set_nn_any works", {
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- matrix(runif(n*p1), n, p1)
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = sample(c(1, NA), size = n, prob = c(1,2), replace = T))
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet",
                            cand_method = "nn_any", rec_method = "distant_cor", 
                            options = list(nn_nn = 3))
  nn_obj <- nearest_neighbor(mat_x, options$nn_options)
  nn_mat <- .query_nn(nn_obj, options$nn_options)
  
  res <- .candidate_set_nn_any(df_res, nn_mat, options$cand_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("vec_cand", "diagnostic"))))
  expect_true(is.numeric(res$vec_cand))
  expect_true(all(res$vec_cand %% 1 == 0))
  expect_true(all(res$vec_cand > 0))
  expect_true(all(res$vec_cand <= n))
  expect_true(all(is.na(df_res$order_rec[res$vec_cand])))
})

################################

## .candidate_set_nn_freq is correct

test_that(".candidate_set_nn_freq works", {
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- matrix(runif(n*p1), n, p1)
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = sample(c(1, NA), size = n, prob = c(1,2), replace = T))
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_freq", rec_method = "distant_cor",
                            options = list(nn_nn = 3, cand_num_cand = 20))
  nn_obj <- nearest_neighbor(mat_x, options$nn_options)
  nn_mat <- .query_nn(nn_obj, options$nn_options)
  
  res <- .candidate_set_nn_freq(df_res, nn_mat, options$cand_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("vec_cand", "diagnostic"))))
  expect_true(is.numeric(res$vec_cand))
  expect_true(all(res$vec_cand %% 1 == 0))
  expect_true(all(res$vec_cand > 0))
  expect_true(all(res$vec_cand <= n))
  expect_true(all(is.na(df_res$order_rec[res$vec_cand])))
  expect_true(length(res$vec_cand) == options$cand_options$num_cand)
})

test_that(".candidate_set_nn_freq works even if cand_options$num_cand is too large", {
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- matrix(runif(n*p1), n, p1)
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = c(rep(1, n-10), rep(NA, 10)))
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_freq", rec_method = "distant_cor",
                            options = list(nn_nn = 3, cand_num_cand = 50))
  nn_obj <- nearest_neighbor(mat_x, options$nn_options)
  nn_mat <- .query_nn(nn_obj, options$nn_options)
  
  res <- .candidate_set_nn_freq(df_res, nn_mat, options$cand_options)
  
  expect_true(length(res$vec_cand) == 10)
  expect_true(length(unique(res$vec_cand)) == 10)
})
