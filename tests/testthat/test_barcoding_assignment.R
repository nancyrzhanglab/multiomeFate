context("Testing barcoding assignments")

## .multinomial_posterior is correct

test_that(".multinomial_posterior works", {
  set.seed(10)
  nlineages <- 20
  n <- 100
  gamma <- seq(0.1,10,length.out=nlineages)
  lin_mat <- matrix(sample(0:5, n*nlineages, replace = T, prob = 6:1),
                    nrow = nlineages, ncol = n)
  res <- .multinomial_posterior(
    bool_force_rebase = T,
    gamma = gamma,
    lin_mat = lin_mat
  )
  
  expect_true(is.matrix(res))
  expect_true(sum(abs(colSums(res) - 1)) <= 1e-3)
  expect_true(all(dim(res) == dim(lin_mat)))
})

##############################

## .multinomial_posterior_vector is correct

test_that(".multinomial_posterior_vector is correct", {
  trials <- 2000
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    nlineages <- 20
    gamma <- seq(0.1,10,length.out=nlineages)
    lgamma <- log(gamma)
    lin_count <- sample(0:20, nlineages, replace = T, prob = 21:1)
    res <- .multinomial_posterior_vector(
      bool_force_rebase = F,
      lgamma = lgamma,
      lin_count = lin_count
    )
    res2 <- .multinomial_posterior_vector(
      bool_force_rebase = T,
      lgamma = lgamma,
      lin_count = lin_count
    )
    
    res3 <- gamma^(lin_count)
    res3 <- res3/sum(res3)
    
    bstarc <-  which.max(lin_count)
    ldeltabc <- lgamma - lgamma[bstarc]
    diff_vec <- lin_count - lin_count[bstarc]
    tmp <- exp(lin_count*ldeltabc + diff_vec*lgamma[bstarc])
    res4 <- tmp/sum(tmp)
    
    sum(abs(res - res2)) <= 1e-5 & sum(abs(res - res3)) <= 1e-5 & sum(abs(res - res4)) <= 1e-5
  })
  
  expect_true(all(bool_vec))
})
