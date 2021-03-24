context("Test generating data")

## .generate_ygivenx is correct

test_that(".generate_ygivenx works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 5
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), bool_traj_y = T)
  
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
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), bool_traj_y = T)
  
  res <- .generate_xgiveny(obj_next, obj_next$vec_starty)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "time"))))
  expect_true(length(res$x) == p1)
  expect_true(all(res$x >= 0))
  expect_true(all(res$x <= 1))
  expect_true(res$time >= 0)
  expect_true(res$time <= 1)
})
