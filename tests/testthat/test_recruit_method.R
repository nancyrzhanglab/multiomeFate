context("Test recruiting method")

generate_data_test <- function(seed = 10){
  set.seed(seed)
  g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                   directed = F)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
  idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
  combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                            num_per_wave = num_per_wave, 
                                            distinct_waves = distinct_waves)
  
  res <- simulate_data_input(combn_wave_mat)
  df_x <- res$df_x; df_y <- res$df_y
  list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise
  df_cell <- simulate_df_cell(100, time_max = max(df_y$time_end_scaffold, na.rm = T),
                              num_branch = 3)
  
  set.seed(10)
  dat <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
  mat_x <- dat$obs_x; mat_y <- dat$obs_y
  
  vec_start <- which(df_cell$time <= 10)
  list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
    intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
  })
  
  list(mat_x = mat_x, mat_y = mat_y, 
       df_x = df_x, df_y = df_y,
       vec_start = vec_start, list_end = list_end)
}

## .predict_yfromx is correct

test_that(".predict_yfromx works", {
  set.seed(12)
  dat <- generate_data_test()
  prep_obj <- chromatin_potential_prepare(dat$mat_x, dat$mat_y, dat$df_x, dat$df_y,
                                          dat$vec_start, dat$list_end, rec_method = "nn",
                                          options = list(rec_num_rec = 5))
  est_options <- prep_obj$options$est_options
  
  res_g <- .estimate_g_glmnet(dat$mat_x, dat$mat_y, weights = rep(1, nrow(dat$mat_x)),
                              est_options)
  
  res <- .predict_yfromx(dat$mat_x, res_g, family = est_options$family)
  
  expect_true(is.matrix(res))
  expect_true(all(!is.na(res)))
  expect_true(all(dim(res) == dim(dat$mat_y)))
})

#############

## .recruit_next_nn is correct

test_that(".recruit_next_nn works", {
  set.seed(10)
  dat <- generate_data_test()
  p1 <- ncol(dat$mat_x)
  p2 <- ncol(dat$mat_y)
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$mat_x, dat$mat_y, dat$df_x, dat$df_y,
                                     dat$vec_start, dat$list_end, rec_method = "nn",
                                     options = list(rec_num_rec = 5))
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  mat_g <- matrix(runif(p1*p2), nrow = p1, ncol = p2)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res <- .recruit_next_nn(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = F,
                          options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec))
  expect_true(length(res$rec) == options$rec_options$num_rec)
  
  # works with enforce_matched
  res <- .recruit_next_nn(dat$mat_x, dat$mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = T,
                          options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec))
  expect_true(length(res$rec) == options$rec_options$num_rec)
  for(i in 1:length(res$rec)){
    expect_true(any(!is.na(df_res$order_rec[res$rec[[i]]$to])))
  }
})

#########################

## .recruit_next_distant_cor is correct

test_that(".recruit_next_distant_cor works", {
  set.seed(10)
  dat <- generate_data_test()
  p1 <- ncol(dat$mat_x)
  p2 <- ncol(dat$mat_y)
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$mat_x, dat$mat_y, dat$df_x, dat$df_y,
                                          dat$vec_start, dat$list_end, rec_method = "distant_cor",
                                          options = list())
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_g <- prep_obj$nn_g
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  mat_g <- matrix(runif(p1*p2), nrow = p1, ncol = p2)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res <- .recruit_next_distant_cor(dat$mat_x, dat$mat_y, res_cand$vec_cand, res_g, df_res, 
                                   dim_reduc_obj, nn_g, nn_mat, nn_obj, enforce_matched = F,
                                   prep_obj$options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec))
  expect_true(length(res$rec) == length(res_cand$vec_cand))
  
  # works with enforce_matched
  res <- .recruit_next_distant_cor(dat$mat_x, dat$mat_y, res_cand$vec_cand, res_g, df_res, 
                                   dim_reduc_obj, nn_g, nn_mat, nn_obj, enforce_matched = T,
                                   options$rec_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true(all(sort(names(res$rec)) == sort(c("vec_from", "list_to"))))
  expect_true(length(res$rec$vec_from) == length(res$rec$list_to))
  expect_true(is.list(res$rec))
  expect_true(length(res$rec) == length(res_cand$vec_cand))
  for(i in 1:length(res$rec)){
    expect_true(any(!is.na(df_res$order_rec[res$rec[[i]]$to])))
  }
})

#################

## .recruit_next is correct

test_that(".recruit_next works", {
  set.seed(10)
  dat <- generate_data_test()
  p1 <- ncol(dat$mat_x)
  p2 <- ncol(dat$mat_y)
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(dat$mat_x, dat$mat_y, dat$df_x, dat$df_y,
                                          dat$vec_start, dat$list_end, rec_method = "distant_cor",
                                          options = list())
 
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_g <- prep_obj$nn_g
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  mat_g <- matrix(runif(p1*p2), nrow = p1, ncol = p2)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
    
  res <- .recruit_next(dat$mat_x, dat$mat_y, res_cand$vec_cand, res_g, df_res, dim_reduc_obj, 
                       nn_g, nn_mat, nn_obj, enforce_matched = F, df_cell = NA,
                       prep_obj$options$rec_options)
  
  expect_true(all(sort(names(res)) == sort(c("rec", "diagnostic"))))
  expect_true("postprocess" %in% names(res$diagnostic))
})

##################################

## .find_to_list is correct

test_that(".find_to_list works", {
  tmp <- generate_data_test()
  mat_x <- tmp$mat_x; mat_y <- tmp$mat_y
  df_x <- tmp$df_x; df_y <- tmp$df_y
  vec_start <- tmp$vec_start; list_end <- tmp$list_end
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                          vec_start, list_end)
  
  vec_cand <- c(15:20); i <- 1; cell <- vec_cand[i]
  nn_g <- prep_obj$nn_g; nn_mat <- prep_obj$nn_mat
  list_nn <- .extract_nn_list(vec_cand, nn_mat)
  rec_options <- prep_obj$options$rec_options
  
  res <- .find_to_list(cell, include_idx = NA, exclude_idx = list_nn[[i]],
                       nn_g, nn_mat, rec_options)
  
  expect_true(is.list(res))
})
