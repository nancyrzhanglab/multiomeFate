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
test_that(".weighted_cdf works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  cdf_res <- .weighted_cdf(values = values, weights = weights)
  res <- .compute_decreasing_density(cdf = cdf_res$cdf,
                                     x = cdf_res$x)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x.knots", "y.knots", "slope.knots"))))
  expect_true(all(diff(res$x.knots) >= 0))
  expect_true(all(diff(res$y.knots) >= 0))
  expect_true(sum(abs(range(res$y.knots) - c(0,1))) <= 1e-6)
  expect_true(all(diff(res$slope.knots) <= 0))
})

#####################

## .epanechnikov_weights is correct

test_that(".epanechnikov_weights works", {
  res <- .epanechnikov_weights(c(-1, -.5, .2, .7, 1))
  expect_true(abs(sum(res)-1) <= 1e-6)
  expect_true(all(abs(res[c(1,5)]) <= 1e-6))
})

#####################

## .epanechnikov_kernel_all is 
test_that(".epanechnikov_kernel_all works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  cdf_res <- .weighted_cdf(values = values, weights = weights)
  stepdensity_res <- .compute_decreasing_density(cdf = cdf_res$cdf,
                                             x = cdf_res$x)
  x <- stepdensity_res$x.knots 
  pdf <- stepdensity_res$slope.knots
  pdf <- c(pdf[1], pdf)
  res <- .epanechnikov_kernel_all(bandwidth = 0.5,
                                  x = x,
                                  pdf = pdf)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "grenander"))
  expect_true(all(sort(names(res)) == sort(c("x", "pdf"))))
  expect_true(sum(abs(res$x - x)) <= 1e-6)
  expect_true(all(diff(res$pdf) <= 0))
})

########################

## .smooth_stepdensity is correct
test_that(".smooth_stepdensity works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  cdf_res <- .weighted_cdf(values = values, weights = weights)
  stepdensity_res <- .compute_decreasing_density(cdf = cdf_res$cdf,
                                                 x = cdf_res$x)
  res <- .smooth_stepdensity(bandwidth = 0.5,
                             stepdensity_res = stepdensity_res)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "grenander"))
  expect_true(all(sort(names(res)) == sort(c("x", "pdf"))))
  expect_true(all(diff(res$pdf) <= 0))
})

test_that(".smooth_stepdensity roughly integrates to 1", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  cdf_res <- .weighted_cdf(values = values, weights = weights)
  stepdensity_res <- .compute_decreasing_density(cdf = cdf_res$cdf,
                                                 x = cdf_res$x)
  res <- .smooth_stepdensity(bandwidth = 0.005,
                             stepdensity_res = stepdensity_res)
  
  right_area <- sum(res$pdf[-length(res$pdf)] * diff(res$x))
  left_area <- sum(res$pdf[-1] * diff(res$x))
  
  expect_true(abs(right_area - 1) <= .001)
  expect_true(abs(left_area - 1) <= .001)
  expect_true(right_area >= left_area)
})
