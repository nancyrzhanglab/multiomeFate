# see Algorithm 1 in "Supervised principal component analysis: Visualization, classification andregression on subspaces and submanifolds"
supervised_pca <- function(x, y){
  stopifnot(is.matrix(x), is.matrix(y),
            nrow(x) == nrow(y))
  
  x <- t(x); y <- t(y)
  
  n <- ncol(x)
  L <- crossprod(y)
  H <- diag(n) - tcrossprod(rep(1,n))/n
  Q <- tcrossprod(x %*% H %*% L %*% H, x)
  eigen_res <- eigen(Q)
  U <- eigen_res$vectors[,1:2]
  
  res <- t(crossprod(U,x))
  colnames(res) <- c("SPCA_1", "SPCA_2")
  colnames(U) <- c("SPCA_1", "SPCA_2")
  rownames(U) <- rownames(x)
  
  list(dimred = res, U = U)
}

form_onehot_classification_mat <- function(y){
  uniq_val <- sort(unique(y))
  k <- length(uniq_val)
  n <- length(y)
  
  mat <- matrix(0, nrow = n, ncol = k)
  for(j in 1:k){
    mat[which(y == uniq_val[j]),j] <- 1
  }
  colnames(mat) <- uniq_val
  if(length(names(y)) > 0) rownames(mat) <- names(y)
  
  mat
}