context("Testing weighted grenander")

## .weighted_cdf is correct

test_that(".weighted_cdf works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  res <- .weighted_cdf(values = values, weights = weights)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("cdf", "x"))))
  expect_true(sum(abs(range(res$cdf) - c(0,1))) <= 1e-6)
  expect_true(sum(abs(range(res$x) - c(0, max(values)))) <= 1e-6)
  expect_true(all(diff(res$x) >= 0))
  expect_true(all(diff(res$cdf) >= 0))
})

######################

## .compute_decreasing_density is correct
test_that(".compute_decreasing_density works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  cdf_res <- .weighted_cdf(values = values, weights = weights)
  res <- .compute_decreasing_density(cdf = cdf_res$cdf,
                                     x = cdf_res$x)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "pdf"))))
  expect_true(all(diff(res$x) >= 0))
  expect_true(all(diff(res$pdf) <= 0))
  
  area <- sum(diff(res$x)*res$pdf[-length(res$pdf)])
  expect_true(abs(area - 1) <= 1e-6)
})

#############################

## estimate_grenander is correct

test_that("estimate_grenander works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  
  res <- estimate_grenander(values = values,
                            weights = weights)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "grenander"))
  expect_true(all(sort(names(res)) == sort(c("x", "pdf", "log_pdf", "scaling_factor"))))
  expect_true(all(diff(res$pdf) <= 1e-6))
  expect_true(all(diff(res$x) >= 0))
})

test_that("estimate_grenander works for a non-decreasing density", {
  set.seed(10)
  n <- 1000
  values <- sapply(1:n, function(i){
    mixture <- stats::rbinom(1, 1, 0.5)
    if(mixture == 0) stats::rexp(1) else stats::rnorm(1, mean = 5, sd = 0.5)
  })
  weights <- rep(1, n)
  
  res <- estimate_grenander(values = values,
                            weights = weights)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "grenander"))
  expect_true(all(diff(res$pdf) <= 1e-6))
  expect_true(all(diff(res$x) >= 0))
})

test_that("estimate_grenander yields a density that integrates to 1", {
  set.seed(10)
  n <- 1000
  values <- sapply(1:n, function(i){
    mixture <- stats::rbinom(1, 1, 0.5)
    if(mixture == 0) stats::rexp(1) else stats::rnorm(1, mean = 5, sd = 0.5)
  })
  weights <- rep(1, n)
  
  res <- estimate_grenander(values = values,
                            weights = weights)
  
  area <- sum(diff(res$x)*res$pdf[-length(res$pdf)])
  expect_true(abs(area - 1) <= 1e-6)
})

test_that("estimate_grenander is the MLE for the unweighted samples", {
  trials <- 10
  nboot <- 100
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    
    n <- 100
    values <- sapply(1:n, function(i){
      mixture <- stats::rbinom(1, 1, 0.5)
      if(mixture == 0) stats::rexp(1) else stats::rnorm(1, mean = 5, sd = 0.5)
    })
    weights <- rep(1, n)
    
    res <- estimate_grenander(values = values,
                              weights = weights)
    
    ll <- sum(log(sapply(values, function(x){
      evaluate_grenander(obj = res, x = x)
    })))
    
    ll_others <- sapply(1:nboot, function(j){
      set.seed(j)
      weights <- stats::runif(n)
      res_alt <- estimate_grenander(values = values,
                                    weights = weights)
      sum(log(sapply(values, function(x){
        evaluate_grenander(obj = res_alt, x = x)
      })))
    })
    
    ll >= max(ll_others)
  })
  
  expect_true(all(bool_vec))
})

test_that("estimate_grenander is the MLE for the weighted samples", {
  trials <- 10
  nboot <- 1000
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    
    n <- 10
    values <- sapply(1:n, function(i){
      mixture <- stats::rbinom(1, 1, 0.5)
      if(mixture == 0) stats::rexp(1) else stats::rnorm(1, mean = 5, sd = 0.5)
    })
    weights <- stats::runif(n)
    
    res <- estimate_grenander(values = values,
                              weights = weights)
    
    ll <- sum(weights*log(sapply(values, function(x){
      evaluate_grenander(obj = res, x = x)
    })))
    
    ll_others <- sapply(1:nboot, function(j){
      set.seed(j)
      weights2 <- stats::runif(n)
      res_alt <- estimate_grenander(values = values,
                                    weights = weights2)
      sum(weights*log(sapply(values, function(x){
        evaluate_grenander(obj = res_alt, x = x)
      })))
    })
    
    ll >= max(ll_others)
  })
  
  expect_true(all(bool_vec))
})

#######################

## evaluate_grenander is correct
test_that("evaluate_grenander works", {
  set.seed(10)
  n <- 1000
  values <- sapply(1:n, function(i){
    mixture <- stats::rbinom(1, 1, 0.5)
    if(mixture == 0) stats::rexp(1) else stats::rnorm(1, mean = 5, sd = 0.5)
  })
  weights <- runif(n)
  
  obj1 <- estimate_grenander(values = values*10,
                             weights = weights,
                             scaling_factor = 10)
  stopifnot(abs(obj1$scaling_factor - 10) <= 1e-6)
  
  obj2 <- estimate_grenander(values = values,
                             weights = weights,
                             scaling_factor = 1)
  stopifnot(abs(obj2$scaling_factor - 1) <= 1e-6)
  
  obj3 <- estimate_grenander(values = values*10,
                             weights = weights,
                             scaling_factor = 1)
  stopifnot(abs(obj3$scaling_factor - 1) <= 1e-6)
  
  test_values <- runif(100)
  bool_vec <- sapply(test_values, function(x){
    res1 <- evaluate_grenander(obj = obj1, x = 10*x)
    res2 <- evaluate_grenander(obj = obj2, x = x)
    res3 <- evaluate_grenander(obj = obj3, x = x)
    
    abs(res1 - res2) <= 1e-6 & abs(res1 - res3) >= 1e-6
  })
  
  expect_true(all(bool_vec))
})


