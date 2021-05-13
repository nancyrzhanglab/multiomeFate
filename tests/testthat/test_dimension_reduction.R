context("Test dimension reduction")

## dimension_reduction is correct

test_that("dimension_reduction works and has correct calculations", {
  set.seed(10)
  k <- 2
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(dim_nlatent_x = k))
  dim_options <- options$dim_options
  
  n <- 100; p <- 50
  cov_mat <- matrix(0, p, p)
  cov_mat[1:(p/2), 1:(p/2)] <- 0.5
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- 1
  mat <- MASS::mvrnorm(n, mu = rep(1,p), Sigma = cov_mat)
  mode = "x"
  
  res <- dimension_reduction(mat, mode = "x", dim_options)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dimred", "vec_mean", "vec_sd", "mat_proj"))))
  expect_true(is.matrix(res$dimred))
  expect_true(all(dim(res$dimred) == c(n, k)))
  expect_true(!is.matrix(res$vec_mean))
  expect_true(length(res$vec_mean) == p)
  expect_true(length(res$vec_sd) == p)
  expect_true(is.matrix(res$mat_proj))
  expect_true(all(dim(res$mat_proj) == c(p, k)))
  
  expect_true(all(abs(colMeans(res$dimred)) <= 1e-4))
  expect_true(sum(abs(res$dimred - scale(mat, center = T, scale = T)%*%res$mat_proj)) <= 1e-4)
  
  tmp <- stats::prcomp(mat, center = T, scale. = T)$x[,1:k]
  tmp <- tmp/svd(tmp)$d[1]
  sign_mat <- matrix(c(1,1, 1,-1, -1,1, -1,-1), nrow = 4, ncol = 2, byrow = T)
  bool_vec <- sapply(1:4, function(i){
    sum(abs(res$dimred - tmp%*%diag(sign_mat[i,]))) <= 1e-4
  })
  expect_true(any(bool_vec))
})

##########################

## .dimension_reduction_pca is correct

test_that(".dimension_reduction_pca works on sparse matrices", {
  set.seed(10)
  k <- 2
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(dim_nlatent_x = k))
  dim_options <- options$dim_options
  
  n <- 100; p <- 50
  cov_mat <- matrix(0, p, p)
  cov_mat[1:(p/2), 1:(p/2)] <- 0.5
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- 1
  mat <- MASS::mvrnorm(n, mu = rep(1,p), Sigma = cov_mat)
  mat[sample(1:prod(dim(mat)), size = round(0.9*prod(dim(mat))))] <- 0
  mat2 <- Matrix::Matrix(mat, sparse = TRUE)
  mode = "x"
  
  res <- dimension_reduction(mat2, mode = "x", dim_options)
  
  expect_true(class(mat2) == "dgCMatrix")
  expect_true(all(sort(names(res)) == sort(c("dimred", "vec_mean", "vec_sd", "mat_proj"))))
  
  expect_true(all(abs(colMeans(res$dimred)) <= 1e-4))
  expect_true(sum(abs(res$dimred - scale(mat, center = T, scale = T)%*%res$mat_proj)) <= 1e-4)
  
  tmp <- stats::prcomp(mat, center = T, scale. = T)$x[,1:k]
  tmp <- tmp/svd(tmp)$d[1]
  sign_mat <- matrix(c(1,1, 1,-1, -1,1, -1,-1), nrow = 4, ncol = 2, byrow = T)
  bool_vec <- sapply(1:4, function(i){
    sum(abs(res$dimred - tmp%*%diag(sign_mat[i,]))) <= 1e-4
  })
  expect_true(any(bool_vec))
})

####################

## .apply_dimred is correct

test_that(".apply_dimred works", {
  set.seed(10)
  k <- 2
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(dim_nlatent_x = k))
  dim_options <- options$dim_options
  
  n <- 100; p <- 50
  cov_mat <- matrix(0, p, p)
  cov_mat[1:(p/2), 1:(p/2)] <- 0.5
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- 1
  mat <- MASS::mvrnorm(n, mu = rep(1,p), Sigma = cov_mat)
  mode = "x"
  
  dim_reduc_obj <- vector("list", 0)
  tmp <- dimension_reduction(mat, mode = "x", dim_options)
  x_dimred <- tmp$dimred
  dim_reduc_obj$x_mean <- tmp$vec_mean; dim_reduc_obj$x_sd <- tmp$vec_sd
  dim_reduc_obj$x_proj <- tmp$mat_proj
  
  res <- t(sapply(1:n, function(i){
    .apply_dimred(mat[i,], mode = "x", dim_reduc_obj)
  }))
  
  expect_true(sum(abs(res - x_dimred)) <= 1e-4)
})

