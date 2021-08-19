.l2norm <- function(x){sqrt(sum(x^2))}

# return idx such that vec1[idx] == vec2
.matching_idx <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  
  ord_1 <- order(vec1, decreasing = F)
  rank_2 <- rank(vec2)
  
  ord_1[rank_2]
}

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == ncol(mat))
  mat * rep(vec, rep(nrow(mat), length(vec)))
}

# [[note to self: clean this function up]]
.svd_truncated <- function(mat, K, K_full_rank, vec_mean, vec_sd){
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  if(K == min(dim(mat))) K_full_rank <- T
  
  if(min(dim(mat)) > 2*(K+2)){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      irlba::irlba(mat, nv = ifelse(K_full_rank, K, K+2), center = vec_mean,
                   scale = vec_sd)
    }, warning = function(e){
      if(all(!is.null(vec_mean))) mat <- mat - vec_mean[col(mat)]
      if(all(!is.null(vec_sd))) mat <- .mult_mat_vec(mat, 1/vec_sd)
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
      
    }, error = function(e){
      if(all(!is.null(vec_mean))) mat <- mat - vec_mean[col(mat)]
      if(all(!is.null(vec_sd))) mat <- .mult_mat_vec(mat, 1/vec_sd)
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    })
  } else {
    if(all(!is.null(vec_mean))) mat <- mat - vec_mean[col(mat)]
    if(all(!is.null(vec_sd))) mat <- .mult_mat_vec(mat, 1/vec_sd)
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  
  # pass row-names and column-names
  if(length(rownames(mat)) != 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) != 0) rownames(res$v) <- colnames(mat)
  
  res
}