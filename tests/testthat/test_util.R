context("Testing utils")

## .log_sum_exp is correct

test_that(".log_sum_exp is correct", {
  trials <- 100
  n <- 20
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    vec <- runif(n)
    res1 <- .log_sum_exp(vec)
    res2 <- log(sum(exp(vec)))
    
    abs(res1 - res2) <= 1e-6
  })
  
  expect_true(all(bool_vec))
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i+1)
    vec <- runif(n)
    res1 <- .log_sum_exp(log(vec))
    res2 <- log(sum(vec))
    
    abs(res1 - res2) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that(".log_sum_exp handles NA and Infs", {
  set.seed(10)
  vec <- c(runif(10), -Inf, NA, runif(10), NA, -Inf)
  res <- .log_sum_exp(vec)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(!is.na(res))
  expect_true(!is.infinite(res))
  expect_true(!is.null(res))
})

############################

## .exp_ratio is correct

test_that(".exp_ratio is correct", {
  trials <- 100
  n <- 20
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    vec <- runif(n)
    res1 <- .exp_ratio(vec)
    tmp <- exp(vec)
    res2 <- tmp/sum(tmp)
    
    abs(sum(res1 - res2)) <= 1e-6 & abs(sum(res1) - 1) <= 1e-6
  })
  
  expect_true(all(bool_vec))
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i+1)
    vec <- runif(n)
    res1 <- .exp_ratio(log(vec))
    res2 <- vec/sum(vec)
    
    abs(sum(res1 - res2)) <= 1e-6 & abs(sum(res1) - 1) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that(".exp_ratio handles NA and Infs", {
  set.seed(10)
  vec <- c(runif(10), -Inf, NA, runif(10), NA, -Inf)
  res <- .exp_ratio(vec)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == length(vec))
  expect_true(all(is.na(res[c(11,12,23,24)])))
  expect_true(all(!is.na(res[-c(11,12,23,24)])))
  expect_true(abs(sum(res[-c(11,12,23,24)]) - 1) <= 1e-6)
  
  vec <- c(rep(NA, 5), 0.5, rep(NA, 5))
  res <- .exp_ratio(vec)
  expect_true(all(is.na(res[-6])))
  expect_true(res[6] == 1)
})
