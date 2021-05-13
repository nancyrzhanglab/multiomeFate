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
