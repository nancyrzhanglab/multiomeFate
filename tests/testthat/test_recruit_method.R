context("Test recruiting method")

## .predict_yfromx is correct

test_that(".predict_yfromx works", {
  set.seed(12)
  p1 <- 25; p2 <- 5; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  n <- 100
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  mat_x1 <- matrix(sample(c(0,1), n*p1, replace = T), n, p1)
  mat_y2 <- matrix(rpois(n*p2, lambda = 3), n, p2)
  res_g <- .estimate_g_glmnet(mat_x1, mat_y2, est_options)
  
  res <- .predict_yfromx(mat_x1, res_g)
  
  expect_true(is.matrix(res))
  expect_true(all(!is.na(res)))
  expect_true(all(dim(res) == c(n,p2)))
})

#############

## .recruit_next_nn_yonly is correct

test_that(".recruit_next_nn_yonly works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100; n <- timepoints
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  vec_start <- 1:10
  list_end <- list(90:100)
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
  df_res <- .init_chrom_df(n, vec_start, list_end, paste0("n", 1:n))
  vec_cand <- 80:89
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list())
  
  res <- .recruit_next_nn_yonly(mat_x, mat_y, vec_cand, res_g, df_res,
                                 options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec$list_to))
})

#################

## .recruit_next is correct

test_that(".recruit_next works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100; n <- timepoints
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  vec_start <- 1:10
  list_end <- list(90:100)
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
  df_res <- .init_chrom_df(n, vec_start, list_end, paste0("n", 1:n))
  vec_cand <- 80:89
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list(rec_run_diagnostic = T))
  
  res <- .recruit_next(mat_x, mat_y, vec_cand, res_g, df_res,
                                options$rec_options)
  
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true("postprocess" %in% names(res$diagnostic))
  expect_true(is.data.frame(res$diagnostic$postprocess$df_diag))
  expect_true(nrow(res$diagnostic$postprocess$df_diag) == length(vec_cand))
  expect_true(all(sort(colnames(res$diagnostic$postprocess$df_diag)) == sort(c("idx", "forward_num", "current_num", "backward_num", "selected"))))
  expect_true(length(res$diagnostic$postprocess$lis_nn) == length(vec_cand))
})

###########################

## .recruit_diagnostic_global is correct

test_that(".recruit_diagnostic_global works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100; n <- timepoints
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  vec_start <- 1:10
  list_end <- list(90:100)
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
  df_res <- .init_chrom_df(n, vec_start, list_end, paste0("n", 1:n))
  vec_cand <- 80:89
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list(rec_run_diagnostic = T))
  rec_options <- options$rec_options
  vec_matched <- which(!is.na(df_res$order_rec))
  res_rec <- .recruit_next_nn_yonly(mat_x, mat_y, vec_cand, res_g, df_res, rec_options)
  
  res <- .recruit_diagnostic_global(mat_x, mat_y, vec_cand, res_g, 
                                               df_res, res_rec, rec_options)
  
  expect_true(is.data.frame(res$df_diag))
  expect_true(length(res$lis_nn) == length(vec_cand))
  expect_true(nrow(res$df_diag) == length(vec_cand))
  expect_true(all(sort(colnames(res$df_diag)) == sort(c("idx", "forward_num", "current_num", "backward_num", "selected"))))
})

test_that(".recruit_diagnostic_global gives a sensible output", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  timepoints <- 100; n <- timepoints
  mat_x <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  mat_y <- .predict_yfromx(mat_x, res_g)
  vec_start <- 1:10
  list_end <- list(90:100)
  df_res <- .init_chrom_df(n, vec_start, list_end, paste0("n", 1:n))
  df_res$order_rec <- c(rep(NA,50), rep(1, 50))
  vec_cand <- c(40:50)
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                            options = list(rec_run_diagnostic = T, rec_num_rec = 2,
                                           est_switch = F))
  rec_options <- options$rec_options
  vec_matched <- which(!is.na(df_res$order_rec))
  res_rec <- .recruit_next_nn_yonly(mat_x, mat_y, vec_cand, res_g, df_res, rec_options)
  
  res <- .recruit_diagnostic_global(mat_x, mat_y, vec_cand, res_g,
                                    df_res, res_rec, rec_options)
  
  expect_true(all(res$df_diag$backward_num[which(res$df_diag$selected)] == 0))
  
  ## now try with an aggressive g_function
  est_options <- options$est_options
  est_options <- .gene_peak_map(df$df_x, df$df_y, est_options)
  res_g2 <- .estimate_g_glmnet(mat_x[1:90,], mat_y[11:100,], est_options)
  res_rec2 <- .recruit_next_nn_yonly(mat_x, mat_y, vec_cand, res_g2, df_res, rec_options)
  res2 <- .recruit_diagnostic_global(mat_x, mat_y, vec_cand, res_g2,
                                    df_res, res_rec2, rec_options)
  expect_true(all(res2$df_diag$backward_num == 0))
})

