context("Test recruiting method")

## .predict_yfromx is correct

test_that(".predict_yfromx works", {
  set.seed(12)
  p1 <- 25; p2 <- 5; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  n <- 100
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet",
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list())
  est_options <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  mat_x1 <- matrix(sample(c(0,1), n*p1, replace = T), n, p1)
  mat_y2 <- matrix(rpois(n*p2, lambda = 3), n, p2)
  res_g <- .estimate_g_glmnet(mat_x1, mat_y2, est_options)
  
  res <- .predict_yfromx(mat_x1, res_g, family = options$est_options$family)
  
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
  timepoints <- 20
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)
  
  vec_start <- which(dat$df_info$time <= 0.1)
  list_end <- list(which(dat$df_info$time >= 0.9))
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$obs_x, dat$obs_y, df$df_x, df$df_y,
                                     vec_start, list_end, rec_method = "nn",
                                     options = list(rec_num_rec = 5))
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res <- .recruit_next_nn(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = F,
                          options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec$list_to))
  expect_true(length(res$rec$vec_from) == options$rec_options$num_rec)
  
  # works with enforce_matched
  res <- .recruit_next_nn(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = T,
                          options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec$list_to))
  expect_true(length(res$rec$vec_from) == options$rec_options$num_rec)
})

#########################

## .recruit_next_distant_cor is correct

test_that(".recruit_next_distant_cor works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 20
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)
  
  vec_start <- which(dat$df_info$time <= 0.1)
  list_end <- list(which(dat$df_info$time >= 0.9))
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$obs_x, dat$obs_y, df$df_x, df$df_y,
                                          vec_start, list_end, rec_method = "distant_cor",
                                          options = list())
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res <- .recruit_next_distant_cor(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                                   dim_reduc_obj, nn_mat, nn_obj, enforce_matched = F,
                                   options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec$list_to))
  expect_true(length(res$rec$vec_from) == length(res_cand$vec_cand))
  
  # works with enforce_matched
  res <- .recruit_next_distant_cor(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                                   dim_reduc_obj, nn_mat, nn_obj, enforce_matched = T,
                                   options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec$list_to))
  expect_true(length(res$rec$vec_from) == length(res_cand$vec_cand))
})

#################

## .recruit_next is correct

test_that(".recruit_next works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 20
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)
  
  vec_start <- which(dat$df_info$time <= 0.1)
  list_end <- list(which(dat$df_info$time >= 0.9))
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$obs_x, dat$obs_y, df$df_x, df$df_y,
                                          vec_start, list_end)
  
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
    
  res <- .recruit_next(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, dim_reduc_obj, 
                       nn_mat, nn_obj, enforce_matched = F, options$rec_options)
  
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true("postprocess" %in% names(res$diagnostic))
})
