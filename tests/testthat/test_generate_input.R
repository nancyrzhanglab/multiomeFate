context("Test generating inputs")

## generate_df_simple is correct

test_that("generate_df_simple works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000
  res <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  
  expect_true(is.list(res))
  expect_true(all(sapply(res, is.data.frame)))
  expect_true(all(sort(names(res)) == sort(c("df_x", "df_y"))))
  expect_true(nrow(res$df_x) == p1)
  expect_true(nrow(res$df_y) == p2)
  expect_true(all(res$df_x$location > 0))
  expect_true(all(res$df_x$location <= genome_length))
  expect_true(all(res$df_y$location > 0))
  expect_true(all(res$df_y$location <= genome_length))
})

#####################

## generate_gmat_simple is correct

test_that("generate_gmat_simple works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000; window = 10
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = window)
  res <- generate_gmat_simple(df$df_x, df$df_y, window = window)
  
  expect_true(all(dim(res) == c(nrow(df$df_x), nrow(df$df_y))))
  bool_vec <- sapply(1:nrow(df$df_y), function(i){
    idx <- which(abs(df$df_x$location - df$df_y$location[i]) > window)
    if(length(idx) == 0) return(TRUE)
    
    all(res[idx,i] == 0)
  })
  expect_true(all(bool_vec))
  expect_true(all(rownames(res) == df$df_x$name))
  expect_true(all(colnames(res) == df$df_y$name))
})

#####################

## generate_traj_cascading is correct

test_that("generate_traj_cascading works", {
  p1 <- 20; p2 <- 5; genome_length <- 1000; window = 10
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = window)
  
  resolution <- 20
  res <- generate_traj_cascading(df$df_x, resolution = resolution)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2"))))
  expect_true(all(sapply(res, ncol) == p1))
  expect_true(all(sapply(res, nrow) == resolution-1))
})




