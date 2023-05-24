# see Algorithm 1 in "Supervised principal component analysis: Visualization, classification andregression on subspaces and submanifolds"
supervised_pca <- function(x, y, 
                           k = 2,
                           orthogonalize = T){
  stopifnot(is.matrix(x), is.matrix(y),
            nrow(x) == nrow(y))
  
  n <- nrow(x)
  H <- diag(n) - matrix(1/n, nrow = n, ncol = n)
  half_mat <- crossprod(y, H %*% x)
  Q <- crossprod(half_mat)
  eigen_res <- eigen(Q)
  U <- Re(eigen_res$vectors[,1:k])
  
  res <- x %*% U
  
  if(orthogonalize){
    svd_res <- svd(res)
    rotation_mat <- svd_res$v
    U <- U %*% rotation_mat
    res <- x %*% U
  }
  
  colnames(res) <- paste0("SPCA_", 1:k)
  colnames(U) <- paste0("SPCA_", 1:k)
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