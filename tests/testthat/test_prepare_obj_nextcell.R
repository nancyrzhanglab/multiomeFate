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
  traj_mat <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(1), min_val = 1)
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
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(1), min_val = 1)
  res <- .compute_xfromy(list(mat_traj), mat_g)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(timepoints, p1)))
  
  bool_vec <- sapply(1:timepoints, function(i){
    all(abs(res[i,] %*% mat_g - log(mat_traj[i,])) <= 1e-4)
  })
  expect_true(all(bool_vec))
})

#########################################

## prepare_obj_nextcell is correct

test_that("prepare_obj_nextcell works", {
  set.seed(10)
  p1 <- 20; p2 <- 5; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 30
  
  mat_traj <- generate_traj_cascading(df$df_x, timepoints = timepoints)
  res <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj), verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "mf_obj_next")
  expect_true(all(sort(names(res)) == sort(c("df_x", "df_y", "mat_g", "ht",
                                             "obj_blueprint", "vec_startx", "vec_starty"))))
  
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  res <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj), verbose = F)
})

##########################

## .possion_ygivenx is correct

test_that(".possion_ygivenx works", {
  set.seed(10)
  p1 <- 100; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  x <- runif(p1)
  
  res <- .possion_ygivenx(x, mat_g)
  expect_true(length(res) == p2)
  expect_true(all(res >= 0))
})

test_that(".possion_ygivenx roughly gives the correct mean", {
  set.seed(10)
  p1 <- 100; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  x1 <- seq(0, 1, length.out = p1)
  x2 <- rep(0.5, length = p1)
  
  trials <- 100
  res1 <- sapply(1:trials, function(i){
    .possion_ygivenx(x1, mat_g)
  })
  row1 <- rowMeans(res1)
  res2 <- sapply(1:trials, function(i){
    .possion_ygivenx(x2, mat_g)
  })
  row2 <- rowMeans(res2)
  
  row_target <- exp(x1 %*% mat_g)
  
  expect_true(sum(abs(row1 - row_target)) <= sum(abs(row2 - row_target)))
})

###########################

## .glmnet_logistic is correct

test_that(".glmnet_logistic works", {
  set.seed(10)
  n <- 50
  covariate <- rbind(matrix(rnorm(3*n), n, 3), matrix(rnorm(3*n, mean = 1), n, 3))
  response_prob <- cbind(c(runif(n, max = 0.1), runif(n, min = 0.9)), runif(2*n))
  
  res <- .glmnet_logistic(covariate, response_prob)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_coef", "vec_intercept"))))
  expect_true(all(dim(res$mat_coef) == c(3,2)))
  expect_true(length(res$vec_intercept) == 2)
})

test_that(".glmnet_logistic doesn't crash when all the values are the same", {
  set.seed(10)
  n <- 50
  covariate <- rbind(matrix(rnorm(3*n), n, 3), matrix(rnorm(3*n, mean = 1), n, 3))
  response_prob <- cbind(rep(0, 2*n), rep(0.9, 2*n))
  
  res <- .glmnet_logistic(covariate, response_prob)
  
  expect_true(all(dim(res$mat_coef) == c(3,2)))
  expect_true(all(is.na(res$mat_coef)))
  expect_true(all(res$vec_intercept == c(0, 0.9)))
})

#########################

## .bernoulli_xgiveny is correct

test_that(".bernoulli_xgiveny gives meaningful answers", {
  set.seed(10)
  n <- 50
  covariate <- rbind(matrix(rnorm(3*n), n, 3), matrix(rnorm(3*n, mean = 1), n, 3))
  response_prob <- cbind(c(runif(n, max = 0.1), runif(n, min = 0.9)), runif(2*n))
  obj <- .glmnet_logistic(covariate, response_prob)
  
  res <- sapply(1:nrow(covariate), function(i){
    .bernoulli_xgiveny(covariate[i,], obj$mat_coef, obj$vec_intercept)
  })
  
  expect_true(nrow(res) == 2)
  expect_true(sd(res[1,]) >= sd(res[2,]))
  expect_true(median(res[1,1:n]) <= median(res[1,(n+1):(2*n)]))
})

test_that(".bernoulli_xgiveny gives meaningful answers under constant vectors", {
  set.seed(10)
  n <- 50
  covariate <- rbind(matrix(rnorm(3*n), n, 3), matrix(rnorm(3*n, mean = 1), n, 3))
  response_prob <- cbind(rep(0, 2*n), rep(0.9, 2*n))
  obj <- .glmnet_logistic(covariate, response_prob)
  
  res <- sapply(1:nrow(covariate), function(i){
    .bernoulli_xgiveny(covariate[i,], obj$mat_coef, obj$vec_intercept)
  })
  
  expect_true(all(res[1,] == 0))
  expect_true(all(res[2,] == 0.9))
})

#########################



