context("Test markov cluster")

generate_markov <- function(){
  set.seed(10)
  n <- 20; K <- 3
  
  # generate U and V
  U <- do.call(rbind, lapply(1:K, function(k){
    t(sapply(1:n, function(x){
      tmp <- runif(K); tmp[k] <- runif(1, max = 3)
      tmp <- tmp/sum(tmp)
      tmp
    }))
  }))
  
  V <- do.call(rbind, lapply(1:K, function(k){
    t(sapply(1:n, function(x){
      tmp <- runif(K); tmp[k] <- runif(1, max = 3)
      tmp <- tmp/sum(tmp)
      tmp
    }))
  }))

  # set anchor points
  V[n,] <- c(1,0,0)
  V[2*n,] <- c(0,1,0)
  V[3*n,] <- c(0,0,1)
  
  # normalize V properly
  for(j in 1:3){
    V[,j] <- V[,j]/sum(V[,j])
  }
    
  # compute the transition matrix
  P <- tcrossprod(U, V)
  
  # sample with multiple restarts
  tmp <- 10000*abs(Re(eigen(t(P))$vectors[,1]))
  N <- round(.mult_mat_vec(P, tmp))
  N
}

##########################

## .optim_fn is correct

test_that(".optim_fn works", {
  vertices <- cbind(rbind(0, diag(3)), 0)
  theta <- c(0, 1/3, 1/3, 1/3)
  x <- as.numeric(theta %*% vertices)
  vec <- c(0.1, 0, 0.4, 0.6)
  
  res <- .optim_fn(x, vec, vertices)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})

#################################

## .optim_gn is correct

test_that(".optim_gn works", {
  vertices <- cbind(rbind(0, diag(3)), 0)
  theta <- c(0, 1/3, 1/3, 1/3)
  vec <- as.numeric(theta %*% vertices)
  x <- c(0.1, 0, 0.4, 0.6)
  
  res <- .optim_gn(x, vec, vertices)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 4)
})

test_that(".optim_gn provides reasonable gradients inequalities", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    vertices <- matrix(rnorm(4*3), nrow = 4, ncol = 3)
    vec <- rnorm(3)
    x1 <- runif(4); x1 <- x1/sum(x1)
    x2 <- runif(4); x2 <- x2/sum(x2)
    
    res <- .optim_gn(x1, vec, vertices)
    val1 <- .optim_fn(x2, vec, vertices)
    val2 <-  .optim_fn(x1, vec, vertices) + as.numeric(res %*% (x2-x1))
    
    val1 > val2
  })
  
  expect_true(all(bool_vec))
})

test_that(".optim_gn provides reasonable gradients with numDeriv::grad", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    vertices <- matrix(rnorm(4*3), nrow = 4, ncol = 3)
    vec <- rnorm(3)
    x <- runif(4); x <- x/sum(x)
    
    res1 <- .optim_gn(x, vec, vertices)
    res2 <- numDeriv::grad(
      .optim_fn, x, 
      side = NULL, vec = vec, vertices = vertices
    )
    
    sum(abs(res1 - res2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

########################

## .weighting_optimization is correct

test_that(".weighting_optimization works", {
  set.seed(10)
  vertices <- matrix(rnorm(4*3), nrow = 4, ncol = 3)
  vec <- rnorm(3)
  res <- .weighting_optimization(vec, vertices)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 4)
})

test_that(".weighting_optimization is correct", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(10)
    vertices <- rbind(0, diag(3))
    theta <- runif(4); theta <- theta/sum(theta)
    vec <- as.numeric(theta %*% vertices)
    
    res <- .weighting_optimization(vec, vertices)
    
    sum(abs(res - theta)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})


