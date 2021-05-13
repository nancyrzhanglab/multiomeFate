context("Test nearest neighbor")

## nearest_neighbor is correct

test_that("nearest_neighbor works", {
  set.seed(10)
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list())
  n <- 100; p <- 20
  mat <- MASS::mvrnorm(n, rep(0,p), diag(p))
  
  res <- nearest_neighbor(mat, options$nn_options)
  
  expect_true(class(res) == "Rcpp_AnnoyEuclidean")
  expect_true(is.environment(res))
  
  vec <- rnorm(p)
  tmp1 <- res$getNNsByVectorList(vec, 5, search_k = -1, include_distances = T)
  expect_true(is.list(tmp1))
  expect_true(all(sort(names(tmp1)) == sort(c("item", "distance"))))
  expect_true(all(order(tmp1$distance, decreasing = F) == 1:5))
  
  tmp2 <- res$getNNsByVector(vec, 5)
  expect_true(length(tmp2) == 5)
  expect_true(all(tmp2 == tmp1$item))
})

####################

## .query_nn is correct

test_that(".query_nn works", {
  set.seed(10)
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list(nn_nn = 5))
  nn_options <- options$nn_options
  
  n <- 100; p <- 20
  mat <- t(sapply(1:n, function(i){
    stats::rnorm(p, mean = i, sd = 0.1)
  }))
  
  nn_obj <- nearest_neighbor(mat, options$nn_options)
  
  res <- .query_nn(nn_obj, options$nn_options)
  
  expect_true(all(dim(res) == c(n,5)))
  bool_vec <- sapply(1:n, function(i){
    all(abs(res[i,]-i) <= 5) & !i %in% res[i,] & length(unique(res[i,])) == length(res[i,])
  })
  expect_true(all(bool_vec))
})
