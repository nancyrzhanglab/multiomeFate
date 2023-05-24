context("Test supervised PCA")

test_that("supervised_pca yields uncorrelated variables", {
  set.seed(10)
  n <- 100
  x <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = matrix(c(1,0.9,0.1,0.1,
                                                         0.9,1,0.1,0.1,
                                                         0.1,0.1,1,0.9,
                                                         0.1,0.1,0.9,1), nrow = 4, ncol = 4))
  y <- form_onehot_classification_mat(sample(1:3,n,replace = T))
  res <- supervised_pca(x,y)
  
  cor_mat <- crossprod(res$dimred)
})