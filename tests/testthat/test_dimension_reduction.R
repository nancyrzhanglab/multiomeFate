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
  expect_true(all(sort(names(res)) == sort(c("scores", "dim_reduc_obj"))))
  expect_true(is.matrix(res$scores))
  expect_true(all(dim(res$scores) == c(n, k)))
  expect_true(!is.matrix(res$dim_reduc_obj$vec_mean))
  expect_true(length(res$dim_reduc_obj$vec_mean) == p)
  expect_true(length(res$dim_reduc_obj$vec_sd) == p)
  expect_true(is.matrix(res$dim_reduc_obj$mat_proj))
  expect_true(all(dim(res$dim_reduc_obj$mat_proj) == c(p, k)))
  
  expect_true(all(abs(colMeans(res$scores)) <= 1e-4))
  expect_true(sum(abs(res$scores - scale(mat, center = T, scale = T)%*%res$dim_reduc_obj$mat_proj)) <= 1e-4)
  
  tmp <- stats::prcomp(mat, center = T, scale. = T)$x[,1:k]
  tmp <- tmp/svd(tmp)$d[1]
  sign_mat <- matrix(c(1,1, 1,-1, -1,1, -1,-1), nrow = 4, ncol = 2, byrow = T)
  bool_vec <- sapply(1:4, function(i){
    sum(abs(res$scores - tmp%*%diag(sign_mat[i,]))) <= 1e-4
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
  expect_true(all(sort(names(res)) == sort(c("scores", "dim_reduc_obj"))))
  
  expect_true(all(abs(colMeans(res$scores)) <= 1e-4))
  expect_true(sum(abs(res$scores - scale(mat, center = T, scale = T)%*%res$dim_reduc_obj$mat_proj)) <= 1e-4)
  
  tmp <- stats::prcomp(mat, center = T, scale. = T)$x[,1:k]
  tmp <- tmp/svd(tmp)$d[1]
  sign_mat <- matrix(c(1,1, 1,-1, -1,1, -1,-1), nrow = 4, ncol = 2, byrow = T)
  bool_vec <- sapply(1:4, function(i){
    sum(abs(res$scores - tmp%*%diag(sign_mat[i,]))) <= 1e-4
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
  x_scores <- tmp$scores; dim_reduc_obj <- tmp$dim_reduc_obj
  
  res <- t(sapply(1:n, function(i){
    .apply_dimred(mat[i,], dim_reduc_obj)
  }))
  
  expect_true(sum(abs(res - x_scores)) <= 1e-4)
})

