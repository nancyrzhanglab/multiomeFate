context("Test vertex hunting")

## .simplex_dist is correct

test_that(".simplex_dist works", {
  vertices <- cbind(rbind(0, diag(3)), 0)
  coef <- c(0,rep(1/3, 3))
  theta <- as.numeric(coef%*%vertices)
  theta[4] <- 0.5
  
  res <- .simplex_dist(theta, vertices)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("combination", "value"))))
  expect_true(length(res$value) == 1)
  expect_true(is.numeric(res$value))
  expect_true(res$value > 0)
  
  expect_true(sum(abs(res$combination - coef)) <= 1e-6)
  expect_true(abs(res$value - 0.5) <= 1e-6)
})

test_that(".simplex_dist returns 0 correctly", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 5; K <- 4
    vertices <- matrix(runif(n*K), nrow = n, ncol = K)
    coef <- runif(n); coef <- coef/sum(coef)
    theta <- as.numeric(coef%*%vertices)
    
    res <- .simplex_dist(theta, vertices)
    
    bool1 <- sum(abs(res$combination - coef)) <= 1e-6
    bool2 <- abs(res$value) <= 1e-6
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that(".simplex_dist returns at most the added distance", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 5; K <- 4
    vertices <- matrix(runif(n*K), nrow = n, ncol = K)
    coef <- runif(n); coef <- coef/sum(coef)
    theta <- as.numeric(coef%*%vertices)
    dist_add <- rnorm(K)
    theta <- theta + dist_add
    
    res <- .simplex_dist(theta, vertices)
    abs(res$value) <= .l2norm(dist_add)
  })
  
  expect_true(all(bool_vec))
})

test_that(".simplex_dist works when the perturbation is in an extra dimension", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 5; K <- 4
    vertices <- matrix(runif(n*K), nrow = n, ncol = K)
    vertices <- cbind(vertices, 0)
    coef <- runif(n); coef <- coef/sum(coef)
    theta <- as.numeric(coef%*%vertices)
    dist_add <- runif(1)
    theta[K+1] <- theta[K+1] + dist_add
    
    res <- .simplex_dist(theta, vertices)
    abs(res$value - dist_add) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

#############

## .vertex_hunting is correct

test_that(".vertex_hunting works", {
  set.seed(10)
  n <- 500
  vertices <- cbind(rbind(0, diag(2)))
  mat <- t(sapply(1:n, function(x){
    coef <- runif(nrow(vertices)); coef <- coef/sum(coef)
    as.numeric(coef %*% vertices)
  }))
  mat <- mat + matrix(rnorm(prod(dim(mat)), sd = 0.01), nrow = nrow(mat), ncol = ncol(mat))
  
  res <- .vertex_hunting(mat, fixed_clustering = list(),
                         m = 20, K0 = 8, num_restart = 10, max_tries = 500)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(ncol(mat)+1, ncol(mat))))
})



