.init_est_matrices <- function(mat_x, mat_y, vec_start, list_end){
  # initialize
  vec_onlyend <- unlist(list_end)
  vec <- c(vec_start, vec_onlyend)
  n0 <- length(vec); p1 <- ncol(mat_x); p2 <- ncol(mat_y)
  mat_x1 <- matrix(NA, nrow = n0, ncol = p1)
  mat_y1 <- matrix(NA, nrow = n0, ncol = p2)
  mat_y2 <- matrix(NA, nrow = n0, ncol = p2)
  
  # fill in steady-states
  mat_x1 <- mat_x[vec,,drop = F]; mat_y2 <- mat_y[vec,,drop = F]
  mat_y1 <- mat_y[vec_onlyend,,drop = F]
  
  list(mat_x1 = mat_x1, mat_y1 = mat_y1, mat_y2 = mat_y2, 
       idx1 = vec_onlyend)
}

# [[note to self: there must be a better way to write this...]]
.update_estimation_matrices <- function(mat_x, mat_y,
                                        mat_x1, mat_y1, mat_y2, idx1,
                                        rec, form_options){
  stopifnot(form_options[["method"]] == "literal")
  
  p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
  stopifnot(ncol(mat_x1) == p1, ncol(mat_y1) == p2, ncol(mat_y2) == p2,
            nrow(mat_y) == n, nrow(mat_x1) == nrow(mat_y2))
  stopifnot(all(idx1 <= n), length(idx1) == nrow(mat_y1), 
            length(idx1) == length(unique(idx1)), all(idx1 > 0), 
            all(idx1 %% 1 == 0))
  
  n_org1 <- nrow(mat_x1); len <- sum(sapply(rec$list_to, length))
  mat_x1 <- rbind(mat_x1, matrix(NA, nrow = len, ncol = p1))
  mat_y2 <- rbind(mat_y2, matrix(NA, nrow = len, ncol = p2))
  
  counter <- 1
  for(i in 1:length(rec$list_to)){
    for(j in 1:length(rec$list_to[[i]])){
      mat_x1[n_org1+counter,] <- mat_x[rec$vec_from[i],]
      mat_y2[n_org1+counter,] <- mat_y[rec$list_to[[i]][j],]
      
      counter <- counter+1
    }
  }
  
  ###
  
  n_org2 <- nrow(mat_y1)
  mat_y1 <- rbind(mat_y1, matrix(NA, nrow = length(rec$vec_from), ncol = p2))
  counter <- 1
  for(i in 1:length(rec$vec_from)){
    mat_y1[n_org2+counter,] <- mat_y[rec$vec_from[i],]
  }
  idx1 <- c(idx1, rec$vec_from)
  
  stopifnot(all(idx1 <= n), length(idx1) == nrow(mat_y1), 
            length(idx1) == length(unique(idx1)), all(idx1 > 0), 
            all(idx1 %% 1 == 0))
  
  list(mat_x1 = mat_x1, mat_y1 = mat_y1, mat_y2 = mat_y2, 
       idx1 = idx1)
}