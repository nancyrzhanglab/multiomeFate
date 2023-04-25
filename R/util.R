.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

##########################

# computes log(exp(x_1) + exp(x_2) + exp(x_3) + ...)
# via  log(exp(x_1+C) + exp(x_2+C) + exp(x_3+C) + ...) - C
.log_sum_exp <- function(vec, max_val = 25){
  vec <- vec[!is.na(vec)]
  vec <- vec[!is.infinite(vec)]
  if(length(vec) == 0) return(NA)
  C <- max_val - max(vec)
  log(sum(exp(vec+C))) - C
}

# computes exp(x_1)/(exp(x_1) + exp(x_2) + ...)
# via exp(x_1+C)/(exp(x_1+C) + exp(x_2+C) + ...)
.exp_ratio <- function(vec, max_val = 25){
  vec <- vec[!is.na(vec)]
  vec <- vec[!is.infinite(vec)]
  if(length(vec) == 0) return(NA)
  C <- max_val - max(vec)
  tmp <- exp(vec+C)
  tmp/sum(tmp)
}