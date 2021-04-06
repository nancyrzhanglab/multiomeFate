context("Test forming method")

## .init_est_matrices is correct

test_that(".init_est_matrices works", {
  set.seed(10)
  n <- 100; p1 <- 10; p2 <- 5
  mat_x <- matrix(runif(n*p1), n, p1)
  mat_y <- matrix(runif(n*p2), n, p2)
  vec_start <- 1:10
  list_end <- list(80:90, 91:100)
  
  res <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y1", "mat_y2", "idx1"))))
  len <- length(c(vec_start, unlist(list_end)))
  expect_true(all(dim(res$mat_x1) == c(len, p1)))
  expect_true(all(dim(res$mat_y2) == c(len, p2)))
  expect_true(all(dim(res$mat_y1) == c(length(unlist(list_end)), p2)))
  
  expect_true(sum(abs(mat_y[res$idx1,] - res$mat_y1)) <= 1e-6)
})

#######

## .update_estimation_literal is correct

test_that(".update_estimation_literal works", {
  set.seed(10)
  options <- .chrom_options(form_method = "literal", est_method = "glmnet_yonly", 
                            cand_method = "nn", rec_method = "nn",
                            options = list())
  n <- 100; p1 <- 10; p2 <- 5
  mat_x <- matrix(runif(n*p1), n, p1)
  mat_y <- matrix(runif(n*p2), n, p2)
  mat_x1 <- mat_x[c(11:20, 1:10),]
  mat_y2 <- mat_y[c(1:10, 1:10),]
  mat_y1 <- mat_y[11:20,]
  idx1 <- 11:20
  rec <- list(vec_from = 21:30, list_to = lapply(11:20,function(x){x}))
  
  res <- .update_estimation_literal(mat_x, mat_y,
                                    mat_x1, mat_y1, mat_y2, idx1,
                                    rec, options$form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y1", "mat_y2", "idx1"))))
  expect_true(length(res$idx1) == nrow(res$mat_y1))
  expect_true(all(dim(res$mat_x1) == c(nrow(mat_x1) + length(rec$vec_from), ncol(mat_x1))))
  expect_true(all(dim(res$mat_y2) == c(nrow(mat_y2) + length(rec$vec_from), ncol(mat_y2))))
  expect_true(all(dim(res$mat_y1) == c(nrow(mat_y1) + length(rec$vec_from), ncol(mat_y1))))
})

## .update_estimation_average is correct

test_that(".update_estimation_average works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  idx1 <- 81:100; mat_y1 <- mat_y[idx1,]; mat_x1 <- mat_x[idx1,]
  mat_y2 <- mat_y[c(91:100, 91:100),]
  vec_cand <- 60:80
  options <- .chrom_options(form_method = "average", est_method = "glmnet_yonly", 
                            cand_method = "nn", rec_method = "nn",
                            options = list())
  rec <- .recruit_next_nn(mat_x, vec_cand, mat_y1, idx1, res_g, 
                          options$rec_options)
  
  res <- .update_estimation_average(mat_x, mat_y,
                                    mat_x1, mat_y1, mat_y2, idx1,
                                    rec, options$form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y1", "mat_y2", "idx1"))))
  expect_true(length(res$idx1) == nrow(res$mat_y1))
  expect_true(all(dim(res$mat_x1) == c(nrow(mat_x1) + length(rec$vec_from), ncol(mat_x1))))
  expect_true(all(dim(res$mat_y2) == c(nrow(mat_y2) + length(rec$vec_from), ncol(mat_y2))))
  expect_true(all(dim(res$mat_y1) == c(nrow(mat_y1) + length(rec$vec_from), ncol(mat_y1))))
})
