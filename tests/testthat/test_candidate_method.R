context("Test candidate method")

##[[note to self: Need more tests here. for example, add tests of correct output format
## when candidate function early terminates]]

## .candidate_set_nn_xonly_avg is correct

test_that(".candidate_set_nn_xonly_avg works", {
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- matrix(runif(n*p1), n, p1)
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = sample(c(1, NA), size = n, prob = c(1,2), replace = T))
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list())
  
  res <- .candidate_set_nn_xonly_avg(mat_x, df_res, options$cand_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("vec_cand", "diagnostic"))))
  expect_true(is.numeric(res$vec_cand))
  expect_true(all(res$vec_cand %% 1 == 0))
  expect_true(all(res$vec_cand > 0))
  expect_true(all(res$vec_cand <= n))
  expect_true(all(is.na(df_res$order_rec[res$vec_cand])))
  
  set.seed(10)
  n <- 100; p1 <- 20
  mat_x <- t(sapply(1:n, function(x){runif(p1)+x}))
  df_res <- data.frame(init_state = rep(NA, n), 
                       order_rec = c(rep(1, n/2), rep(NA, n/2)))
  res2 <- .candidate_set_nn_xonly_avg(mat_x, df_res, options$cand_options)
  
  expect_true(is.numeric(res2$vec_cand))
  expect_true(all(res2$vec_cand %% 1 == 0))
  expect_true(all(res2$vec_cand > 0))
  expect_true(all(res2$vec_cand <= n))
  expect_true(all(is.na(df_res$order_rec[res2$vec_cand])))
  expect_true(length(res2$vec_cand) <= length(res$vec_cand))
})
