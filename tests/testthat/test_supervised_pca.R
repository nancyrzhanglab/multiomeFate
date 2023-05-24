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
  diag(cor_mat) <- 0
  expect_true(sum(abs(cor_mat)) <= 1e-6)
})

test_that("supervised_pca finds meaningful correlation with the response", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    n <- 100
    x <- MASS::mvrnorm(n, mu = rep(0, 3), Sigma = diag(c(5,2,1)))
    y <- sapply(1:n, function(i){ifelse(x[i,3] >= stats::median(x[,3]), 0, 1)})
    y_mat <- form_onehot_classification_mat(y)
    res <- supervised_pca(x,y_mat,k=1)
    spca_res <- res$dimred
    pca_res <- x %*% svd(x)$v[,1]
    
    # par(mfrow = c(1,2)); plot(spca_res[,1], col = y+1); plot(pca_res[,1], col = y+1)
    pos_idx <- which(y == 1)
    neg_idx <- which(y == 0)
    spca_separation <- abs(mean(spca_res[pos_idx]) - mean(spca_res[neg_idx]))/(sd(spca_res[pos_idx]) + sd(spca_res[neg_idx]))
    pca_separation <- abs(mean(pca_res[pos_idx]) - mean(pca_res[neg_idx]))/(sd(pca_res[pos_idx]) + sd(pca_res[neg_idx]))
    
    spca_separation > pca_separation
  })
  
  expect_true(all(bool_vec))
})
