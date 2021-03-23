context("Test preparing object containing next-cell-generation")

## .compute_xfromy_starting is correct

test_that(".compute_xfromy_starting works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000; window = 10
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  vec_y <- seq(1, 0, length.out = p2)
  
  res <- .compute_xfromy_starting(vec_y, mat_g)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == p1)
  expect_true(all(abs(res %*% mat_g - vec_y) <= 1e-4))
})

#################################

## .compute_xfromy_next is correct

test_that(".compute_xfromy_next works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000; window = 10
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  
  timepoints <- 20
  traj_mat <- generate_traj_cascading(df$df_y, timepoints = timepoints)*2
  vec_x <- .compute_xfromy_starting(traj_mat[1,], mat_g)
  
  res <- .compute_xfromy_next(vec_x, traj_mat[2,], mat_g)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == p1)
  expect_true(all(abs(res %*% mat_g - traj_mat[2,]) <= 1e-4))
})

############################

## .compute_xfromy is correct

test_that(".compute_xfromy_next works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000; window = 10
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  
  timepoints <- 30
  traj_mat <- generate_traj_cascading(df$df_y, timepoints = timepoints)*2
  res <- .compute_xfromy(list(traj_mat), mat_g)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(timepoints, p1)))
  
  bool_vec <- sapply(1:timepoints, function(i){
    all(abs(res[i,] %*% mat_g - traj_mat[i,]) <= 1e-4)
  })
  expect_true(all(bool_vec))
})

