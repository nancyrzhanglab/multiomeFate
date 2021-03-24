context("Test generating data")

## .generate_ygivenx is correct

test_that(".generate_ygivenx works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 5
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj))
  
  res <- .generate_ygivenx(obj_next, obj_next$vec_startx)
  
  expect_true(length(res) == p2)
  expect_true(all(res >= 0))
})

#########################

## .generate_xgiveny is correct

test_that(".generate_xgiveny works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 5
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj))
  
  res <- .generate_xgiveny(obj_next, obj_next$vec_starty)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "time"))))
  expect_true(length(res$x) == p1)
  expect_true(all(res$x >= 0))
  expect_true(all(res$x <= 1))
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
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)*3
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
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)*3
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
