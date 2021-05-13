context("Test forming method")

## .init_est_matrices is correct

test_that(".init_est_matrices works", {
  set.seed(10)
  n <- 100; p1 <- 10; p2 <- 5
  mat_x <- matrix(runif(n*p1), n, p1)
  mat_y <- matrix(runif(n*p2), n, p2)
  vec_start <- 1:10
  list_end <- list(80:90, 91:100)
  cell_name <- paste0("n", 1:n)
  df_res <- .init_chrom_df(n, vec_start, list_end, cell_name)
  
  res <- .init_est_matrices(mat_x, mat_y, df_res)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2"))))
  len <- length(c(vec_start, unlist(list_end)))
  expect_true(all(dim(res$mat_x1) == c(len, p1)))
  expect_true(all(dim(res$mat_y2) == c(len, p2)))
})

#######

## .update_estimation_literal is correct

test_that(".update_estimation_literal works", {
  set.seed(10)
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list())
  n <- 100; p1 <- 10; p2 <- 5
  mat_x <- matrix(runif(n*p1), n, p1)
  mat_y <- matrix(runif(n*p2), n, p2)
  mat_x1 <- mat_x[c(11:20, 1:10),]
  mat_y2 <- mat_y[c(1:10, 1:10),]
  vec_matched <- 11:20
  idx1 <- 11:20
  rec <- list(vec_from = 21:30, list_to = lapply(11:20,function(x){x}))
  
  res <- .update_estimation_literal(mat_x, mat_y, mat_x1, mat_y2, 
                                    rec, options$form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2"))))
  expect_true(all(dim(res$mat_x1) == c(nrow(mat_x1) + length(rec$vec_from), ncol(mat_x1))))
  expect_true(all(dim(res$mat_y2) == c(nrow(mat_y2) + length(rec$vec_from), ncol(mat_y2))))
})

## .update_estimation_average is correct

test_that(".update_estimation_average works", {
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
                                          form_method = "average",
                                          options = list(rec_num_rec = 5))
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res_rec <- .recruit_next_nn(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = F,
                          options$rec_options)
  
  mat_x1 <- mat_x[c(11:20, 1:10),]
  mat_y2 <- mat_y[c(1:10, 1:10),]
  res <- .update_estimation_average(mat_x, mat_y, mat_x1, mat_y2, 
                                    res_rec$rec, options$form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2"))))
  expect_true(all(dim(res$mat_x1) == c(nrow(mat_x1) + length(res_rec$rec$vec_from), ncol(mat_x1))))
  expect_true(all(dim(res$mat_y2) == c(nrow(mat_y2) + length(res_rec$rec$vec_from), ncol(mat_y2))))
})
