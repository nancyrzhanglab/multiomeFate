dimension_reduction <- function(mat, mode, dim_options){
  stopifnot(mode %in% c("x", "y"))
  
  if(dim_options[["method"]] == "pca"){
    res <- .dimension_reduction_pca(mat, mode, dim_options)
  } else {
    stop("Dimension reduction method not found")
  }
  
  res
}

###########

.dimension_reduction_pca <- function(mat, mode, dim_options, tol = 1e-4){
  ## [[note to self: only this sparse class is supported currently. Shouldn't be
  ## too hard to extend to other sparse classes though]]
  bool_sparse <- any(class(mat) == "dgCMatrix")
  
  if(dim_options$mean){
    if(bool_sparse){
      vec_mean <- sparseMatrixStats::colMeans2(mat)
    } else {
      vec_mean <- matrixStats::colMeans2(mat)
    }
  } else {
    vec_mean <- rep(0, ncol(mat))
  }
  
  if(dim_options$sd){
    if(bool_sparse){
      vec_sd <- sparseMatrixStats::colSds(mat)
    } else {
      vec_sd <- matrixStats::colSds(mat)
      vec_sd[which(abs(vec_sd) <= tol)] <- 1
    }
  } else {
    vec_sd <- rep(1, ncol(mat))
  }
    
  K <- ifelse(mode == "x", dim_options$nlatent_x, dim_options$nlatent_y)
  svd_res <- .svd_truncated(mat, K = K, K_full_rank = F, vec_mean = vec_mean, vec_sd = vec_sd)
  dimred <- .mult_mat_vec(svd_res$u, svd_res$d/svd_res$d[1])
  
  list(dimred = dimred, vec_mean = vec_mean, vec_sd = vec_sd, 
       mat_proj = .mult_mat_vec(svd_res$v, rep(1/svd_res$d[1], ncol(svd_res$v))))
}

####################

.apply_dimred <- function(vec, mode, dim_reduc_obj){
  if(mode == "x"){
    stopifnot(length(vec) == length(dim_reduc_obj$x_mean))
    
    vec <- (vec - dim_reduc_obj$x_mean)/dim_reduc_obj$x_sd
    as.numeric(vec %*% dim_reduc_obj$x_proj)
  } else {
    stopifnot(length(vec) == length(dim_reduc_obj$y_mean))
    
    vec <- (vec - dim_reduc_obj$y_mean)/dim_reduc_obj$y_sd
    as.numeric(vec %*% dim_reduc_obj$y_proj)
  }
}

.apply_dimred_mat <- function(mat, mode, dim_reduc_obj){
  n <- nrow(mat)
  t(sapply(1:n, function(i){
    .apply_dimred(mat[i,], mode, dim_reduc_obj)
  }))
}