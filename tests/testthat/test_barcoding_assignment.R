context("Testing barcoding assignments")

## .assignment_multinomial_posterior is correct

test_that(".assignment_multinomial_posterior works", {
  set.seed(10)
  nlineages <- 20
  n <- 100
  gamma <- seq(0.1,10,length.out=nlineages)
  lin_mat <- matrix(sample(0:5, n*nlineages, replace = T, prob = 6:1),
                    nrow = nlineages, ncol = n)
  res <- .assignment_multinomial_posterior(
    bool_force_rebase = T,
    lgamma = ,
    lin_mat
  )
})