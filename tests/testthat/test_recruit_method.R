context("Test recruiting method")

## .predict_yfromx is correct

test_that(".predict_yfromx works", {
  set.seed(12)
  p1 <- 25; p2 <- 5; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  n <- 100
  options <- .chrom_options(form_method = "literal", est_method = "glmnet_yonly", 
                            cand_method = "nn", rec_method = "nn",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  mat_x1 <- matrix(sample(c(0,1), n*p1, replace = T), n, p1)
  mat_y2 <- matrix(rpois(n*p2, lambda = 3), n, p2)
  res_g <- .estimate_g_glmnet_yonly(mat_x1, mat_y2, df$df_y, est_options)
  
  res <- .predict_yfromx(mat_x1, res_g)
  
  expect_true(is.matrix(res))
  expect_true(all(!is.na(res)))
  expect_true(all(dim(res) == c(n,p2)))
})

#############

## .recruit_next_nn is correct

test_that(".recruit_next_nn works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  idx1 <- 90:100; mat_y1 <- mat_y[idx1,]
  vec_cand <- 80:89
  options <- .chrom_options(form_method = "literal", est_method = "glmnet_yonly", 
                            cand_method = "nn", rec_method = "nn",
                            options = list())
  
  res <- .recruit_next_nn(mat_x, vec_cand, mat_y1, idx1, res_g, 
                                 options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$vec_from) == length(res$list_to))
  expect_true(is.list(res$list_to))
})
