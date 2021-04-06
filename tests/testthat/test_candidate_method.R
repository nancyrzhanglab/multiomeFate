context("Test candidate method")

## .candidate_set_nn is correct

test_that(".candidate_set_nn works", {
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- matrix(runif(n*p1), n, p1)
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = sample(c(1, NA), size = n, prob = c(1,2), replace = T))
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn", rec_method = "nn",
                            options = list())
  
  res <- .candidate_set_nn(mat_x, df_res, options$cand_options)
  
  expect_true(is.numeric(res))
  expect_true(all(res %% 1 == 0))
  expect_true(all(res > 0))
  expect_true(all(res <= n))
  expect_true(all(is.na(df_res$order_rec[res])))
  
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- t(sapply(1:n, function(x){runif(p1)+x}))
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = c(rep(1, n/2), rep(NA, n/2)))
  res2 <- .candidate_set_nn(mat_x, df_res, options$cand_options)
  
  expect_true(is.numeric(res))
  expect_true(all(res %% 1 == 0))
  expect_true(all(res > 0))
  expect_true(all(res <= n))
  expect_true(all(is.na(df_res$order_rec[res2])))
  expect_true(length(res2) <= length(res))
})
