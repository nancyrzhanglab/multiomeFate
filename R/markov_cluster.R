markov_cluster <- function(count_mat, K, 
                           m = round(min(c(10*K, nrow(count_mat)/5))), 
                           K0 = min(c(2*K, m-1)), 
                           delta = 0.9,
                           num_restart = 10, max_tries = 500){
  stopifnot(nrow(count_mat) == ncol(count_mat), all(count_mat >= 0),
            inherits(count_mat, c("dgCMatrix", "dgTMatrix", "matrix")),
            K >= 2, K %% 1 == 0)
  n <- nrow(count_mat)
  
  # step 1: SVD
  if(inherits(count_mat, c("dgCMatrix", "dgTMatrix"))){
    col_sum <- sparseMatrixStats::colSums2(count_mat)
  } else {
    col_sum <- matrixStats::colSums2(count_mat)
  }
  
  stopifnot(all(col_sum > 0))
  tmp <- .mult_mat_vec(count_mat, 1/col_sum^(0.5))
  svd_res <- .svd_truncated(tmp, K = K, K_full_rank = F, 
                            vec_mean = NULL, vec_sd = NULL)
  d_mat <- .mult_vec_mat(1/svd_res$v[,1], svd_res$v[,-1])
  
  # step 2: vertex hunting
  vertices <- .vertex_hunting(d_mat, m = m, K0 = K0, 
                              num_restart = num_restart, max_tries = max_tries)
  
  # step 3: optimization
  w_mat <- t(sapply(1:n, function(i){
    .weighting_optimization(d_mat[i,], vertices)
  }))
  w_mat[w_mat < 0] <- 0
  w_mat <- t(apply(w_mat, 1, function(x){x/sum(x)}))
  v_est <- .mult_vec_mat(svd_res$v[,1]*col_sum^(1/2), w_mat)
  v_est <- apply(v_est, 2, function(x){x/sum(x)})
  
  # step 4: completion
  if(inherits(count_mat, c("dgCMatrix", "dgTMatrix"))){
    row_sum <- sparseMatrixStats::rowSums2(count_mat)
  } else {
    row_sum <- matrixStats::rowSums2(count_mat)
  }
  p_mat <- .mult_vec_mat(1/row_sum, count_mat)
  u_est <- p_mat %*% v_est %*% solve(crossprod(v_est))
  
  anchor_idx <- which(apply(w_mat, 1, max) > delta)
  
  list(u = u_est, v = v_est, anchor_idx = anchor_idx)
}

#########################################

.optim_fn <- function(x, vec, vertices){
  sum((vec - as.numeric(x %*% vertices))^2) + (1-sum(x))^2
}

.optim_gn <- function(x, vec, vertices){
  vec1 <- vec - as.numeric(x %*% vertices)
  val2 <- 1-sum(x)
  
  -2*as.numeric(vertices %*% vec1) - 2*val2*rep(1, length(x))
}

.weighting_optimization <- function(vec, vertices){
  stopifnot(nrow(vertices) == ncol(vertices)+1)
  n <- nrow(vertices)
  
  res <- stats::optim(rep(1/n, n), fn = .optim_fn, gr = .optim_gn, method = "BFGS",
                      vec = vec, vertices = vertices)
  res$par
}