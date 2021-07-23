#' Compute the dimension reduction
#'
#' This function returns two objects: \code{scores} (the 
#' matrix with \code{nrow(mat)} rows and \code{dim_options$dims_x} or 
#' \code{dim_options$dims_y} columns, depending on \code{mode})
#' representing the projection of \code{mat} into the lower-dimensional
#' space, and \code{dim_reduc_obj} (a list containing \code{vec_mean}
#' and \code{vec_sd} both as 
#' vectors of length \code{ncol(mat)} and a matrix \code{mat_proj}
#' with \code{ncol(mat)} rows and \code{dim_options$dims_x} or 
#' \code{dim_options$dims_y} columns, all three of which are needed
#' to project a new vector into said lower-dimensional space)
#'
#' @param mat full data for a modality, where each row is a cell and each column is a variable
#' @param mode string, either \code{"x"} or \code{"y"}, which dictates
#' if  \code{dim_options$dims_x} or 
#' \code{dim_options$dims_y} is used as the number of latent dimensions
#' @param dim_options one of the outputs from \code{.chrom_options}
#'
#' @return a list
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
  
  if(mode == "x"){
    dims_vec <- dim_options$dims_x
  } else {
    dims_vec <- dim_options$dims_y
  }

  svd_res <- .svd_truncated(mat, K = max(dims_vec), K_full_rank = F, vec_mean = vec_mean, vec_sd = vec_sd)
  scores <- .mult_mat_vec(svd_res$u[,dims_vec,drop = F], svd_res$d[dims_vec]/svd_res$d[min(dims_vec)])
  
  list(scores = scores, dim_reduc_obj = list(vec_mean = vec_mean, vec_sd = vec_sd, 
       mat_proj = .mult_mat_vec(svd_res$v[,dims_vec,drop = F], rep(1/svd_res$d[min(dims_vec)], length(dims_vec)))))
}

####################

.apply_dimred <- function(vec, dim_reduc_obj){
  stopifnot(length(vec) == length(dim_reduc_obj$vec_mean))
  
  vec <- (vec - dim_reduc_obj$vec_mean)/dim_reduc_obj$vec_sd
  as.numeric(vec %*% dim_reduc_obj$mat_proj)
}

# [[note to self: this could probably be done more efficiently]]
.apply_dimred_mat <- function(mat, dim_reduc_obj){
  n <- nrow(mat)
  t(sapply(1:n, function(i){
    .apply_dimred(mat[i,], dim_reduc_obj)
  }))
}