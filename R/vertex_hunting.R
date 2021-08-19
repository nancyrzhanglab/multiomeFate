.vertex_hunting <- function(mat, m, K0, num_restart, max_tries){
  stopifnot(m > K0, K0 > ncol(mat))
  
  # initialize centers
  K <- ncol(mat)+1
  res <- stats::kmeans(mat, m, iter.max = 100, nstart = num_restart)
  theta <- res$centers
  dist_mat <- as.matrix(stats::dist(theta, method = "euclidean"))
  
  # greedy selection
  idx <- unique(as.numeric(which(dist_mat == max(dist_mat), arr.ind = T)))
  if(K0 > 2){
    while(length(idx) <= K0) {
      remaining_idx <- c(1:m)[-idx]
      tmp <- remaining_idx[which.max(apply(dist_mat[remaining_idx,idx], 1, mean))]
      idx <- c(idx, tmp)
    }
  }
  
  # find the best subset
  comb_mat <- utils::combn(1:K0, K)
  if(max_tries < ncol(combn_mat)){
    combn_mat <- combn_mat[,sample(1:ncol(combn_mat), max_tries)]
  }
  max_values <- rep(0, ncol(comb_mat))
  for (i in 1:ncol(comb_mat)){
    for (j in setdiff(1:m, comb_mat[,i])){
      max_values[i] <- max(.simplex_dist(theta[j,], theta[comb_mat[,i],])$value, max_values[i])
    }
  }
  
  min_idx <- which.min(max_values)
  theta[comb_mat[,min_idx],]
}

#' The l_2 distance between a point and a simplex
#' 
#' This function computes the l_2 distance between a point and a simplex, 
#' that is the shortest l_2 distance between the given point and any point in the simplex.
#' @param theta A  vector, representing a point.
#' @param V The matrix, with each row being a vertex.
#' @return The l_2 distance.
#' @references Modified from \url{https://github.com/cran/TopicScore/blob/master/R/topic_score.R#L145}
.simplex_dist <- function(theta, V){
  stopifnot(length(theta) == ncol(V), nrow(V) >= 2, !is.matrix(theta),
            nrow(V) <= ncol(V)+1)
  n <- nrow(V); K <- ncol(V)
  
  idx <- n
  theta2 <- theta - V[idx,]; V2 <- V[-idx,,drop = F] - rep(V[idx,], each = n-1)
  
  D <- tcrossprod(V2)
  d <- V2 %*% theta2
  b0 <- c(rep(0,n-1), -1)
  A <- cbind(diag(n-1), rep(-1, n-1))

  obj <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b0)
  value <- sqrt(max(2*obj$value + .l2norm(theta2)^2, 0))

  comb <- rep(0, n)
  comb[-idx] <- obj$solution; comb[idx] <- 1-sum(obj$solution)
  
  list(combination = comb, value = value)
}