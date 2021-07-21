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
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(form_bool_include_start = F))
  form_options <- options$form_options
  
  res <- .init_est_matrices(mat_x, mat_y, df_res, form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2"))))
  len <- length(unlist(list_end))
  expect_true(all(dim(res$mat_x1) == c(len, p1)))
  expect_true(all(dim(res$mat_y2) == c(len, p2)))
  
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(form_bool_include_start = T))
  form_options <- options$form_options
  
  res <- .init_est_matrices(mat_x, mat_y, df_res, form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2"))))
  len <- length(c(vec_start, unlist(list_end)))
  expect_true(all(dim(res$mat_x1) == c(len, p1)))
  expect_true(all(dim(res$mat_y2) == c(len, p2)))
})

#######

## .update_estimation_average is correct

test_that(".update_estimation_average works", {
  set.seed(10)
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
  p1 <- ncol(mat_x); p2 <- ncol(mat_y)
  
  vec_start <- which(df_cell$time <= 10)
  list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
    intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
  })
  
  set.seed(10)
  prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                          vec_start, list_end,
                                          rec_method = "nn")
  df_res <- prep_obj$df_res; nn_mat <- prep_obj$nn_mat
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_obj <- prep_obj$nn_obj; options <- prep_obj$options
  mat_g <- matrix(runif(p1*p2), nrow = p1, p2)
  res_g <- list(mat_g = mat_g, vec_g = rep(0, p2))
  
  res_cand <- .candidate_set_nn_any(prep_obj$df_res, prep_obj$nn_mat, options$cand_options)
  
  res_rec <- .recruit_next_nn(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                          dim_reduc_obj, nn_obj, enforce_matched = F,
                          options$rec_options)
  
  mat_x1 <- mat_x[c(11:20, 1:10),]
  mat_y2 <- mat_y[c(1:10, 1:10),]
  res <- .update_estimation_average(mat_x, mat_y, mat_x1, mat_y2, 
                                    res_rec$rec, options$form_options)
  
  expect_true(all(sort(names(res)) == sort(c("mat_x1", "mat_y2", "weights"))))
  expect_true(all(dim(res$mat_x1) == c(nrow(mat_x1) + sum(sapply(res_rec$rec, function(x){length(x$from)})), ncol(mat_x1))))
  expect_true(all(dim(res$mat_y2) == c(nrow(mat_y2) + sum(sapply(res_rec$rec, function(x){length(x$from)})), ncol(mat_y2))))
})
