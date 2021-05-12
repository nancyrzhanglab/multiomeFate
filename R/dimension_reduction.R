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

.dimension_reduction_pca <- function(mat, mode, dim_options){
  ## [[note to self: only this sparse class is supported currently. Shouldn't be
  ## too hard to extend to other sparse classes though]]
  bool_sparse <- class(mat) == "dgCMatrix"
  
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
    }
  } else {
    vec_sd <- rep(1, ncol(mat))
  }
    
  svd_res <- irlba::irlba(mat, nv = ifelse(mode == "x", dim_options$nlatent_x, dim_options$nlatent_y), 
                          center = vec_mean, scale = vec_sd)
  vec_sing <- svd_res$d; vec_sing <- vec_sing/vec_sing[1]
  dimred <- .mult_mat_vec(svd_res$u, vec_sing)
  
  
  list(dimred = dimred, vec_mean = vec_mean, vec_sd = vec_sd, 
       mat_proj = .mult_mat_vec(svd_res$v, vec_sing))
}

####################

.apply_dimred <- function(vec, mode, dim_reduc_obj){
  if(mode == "x"){
    stopifnot(length(vec) == length(dim_reduc_obj$x_mean))
    
    vec <- (vec - dim_reduc_obj$x_mean)/dim_reduc_obj$x_sd
    as.numeric(dim_reduc_obj$x_proj %*% vec)
  } else {
    stopifnot(length(vec) == length(dim_reduc_obj$y_mean))
    
    vec <- (vec - dim_reduc_obj$y_mean)/dim_reduc_obj$y_sd
    as.numeric(dim_reduc_obj$y_proj %*% vec)
  }
}