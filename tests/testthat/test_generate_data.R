context("Test generating data")

## .generate_ygivenx is correct

test_that(".generate_ygivenx works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 5
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj))
  
  x_true <- obj_next$vec_startx
  x_obs <- stats::rbinom(length(x_true), size = 1, prob = x_true)
  res <- .generate_ygivenx(obj_next, x_true = x_true, x_obs = x_obs)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("y_true", "y_obs", "idx_time"))))
  expect_true(length(res$y_true) == p2)
  expect_true(all(res$y_true >= 0))
  expect_true(length(res$y_obs) == p2)
  expect_true(all(res$y_obs >= 0))
  expect_true(length(res$y_obs) == length(res$y_true))
})

#########################

## .generate_xgiveny is correct

test_that(".generate_xgiveny works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 5
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj))
  
  y_true <- obj_next$vec_starty
  y_obs <- stats::rpois(length(y_true), lambda = y_true)
  res <- .generate_xgiveny(obj_next, y_true = y_true, y_obs = y_obs, idx_time = NA)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_true", "x_obs", "time"))))
  expect_true(length(res$x_true) == p1)
  expect_true(all(res$x_true >= 0))
  expect_true(all(res$x_true <= 1))
  expect_true(length(res$x_obs) == p1)
  expect_true(all(res$x_obs >= 0))
  expect_true(all(res$x_obs <= 1))
  expect_true(length(res$x_obs) == length(res$x_true))
  expect_true(res$time >= 0)
  expect_true(res$time <= 1)
})

##################################

## .generate_data_single is correct

test_that(".generate_data_single works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 100
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  
  res <- .generate_data_single(obj_next, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("obs_x", "obs_y", "true_x", "true_y", "df_info"))))
  idx <- which(names(res) %in% c("obs_x", "obs_y", "true_x", "true_y", "df_info"))
  for(i in 1:length(idx)){
    if(names(res)[idx[i]] == "df_info"){
      expect_true(is.data.frame(res[[idx[i]]]))
    } else {
      expect_true(is.matrix(res[[idx[i]]]))
    }
  }
  for(i in 2:length(idx)){
    expect_true(nrow(res[[idx[i]]]) == nrow(res[[idx[1]]]))
  }
  expect_true(ncol(res$obs_x) == p1)
  expect_true(ncol(res$true_x) == p1)
  expect_true(ncol(res$obs_y) == p2)
  expect_true(ncol(res$true_y) == p2)
})

#########################

## generate_data is correct

test_that("generate_data works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 100
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  
  res <- generate_data(obj_next, number_runs = 5, sample_perc = 0.9, verbose = F)
  
  expect_true(class(res) == "mf_simul")
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("obs_x", "obs_y", "true_x", "true_y", "df_info", "df_x", "df_y"))))
  idx <- which(names(res) %in% c("obs_x", "obs_y", "true_x", "true_y", "df_info"))
  for(i in 1:length(idx)){
    if(names(res)[idx[i]] == "df_info"){
      expect_true(is.data.frame(res[[idx[i]]]))
    } else {
      expect_true(is.matrix(res[[idx[i]]]))
    }
  }
  for(i in 2:length(idx)){
    expect_true(nrow(res[[idx[i]]]) == nrow(res[[idx[1]]]))
  }
  expect_true(ncol(res$obs_x) == p1)
  expect_true(ncol(res$true_x) == p1)
  expect_true(ncol(res$obs_y) == p2)
  expect_true(ncol(res$true_y) == p2)
  expect_true(nrow(res$df_x) == p1)
  expect_true(nrow(res$df_y) == p2)
})
