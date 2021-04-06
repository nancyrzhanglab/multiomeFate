.init_est_matrices <- function(mat_x, mat_y, vec_start, list_end){
  stopifnot(nrow(mat_x) == nrow(mat_y))
  n <- nrow(mat_x)
  stopifnot(all(vec_start > 0), all(vec_start %% 1 == 0), all(vec_start <= n))
  for(i in 1:length(list_end)){
    stopifnot(all(list_end[[i]] > 0), all(list_end[[i]] %% 1 == 0), all(list_end[[i]] <= n))
  }
  
  # initialize
  vec_onlyend <- unlist(list_end)
  vec <- c(vec_start, vec_onlyend)
  stopifnot(length(vec) == length(unique(vec)))
  
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

# [[note to self: design a lot of tests for this]]
.update_estimation_matrices <- function(mat_x, mat_y,
                                        mat_x1, mat_y1, mat_y2, idx1,
                                        rec, form_options){
  p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
  stopifnot(c("vec_from", "list_to") %in% names(rec))
  stopifnot(ncol(mat_x1) == p1, ncol(mat_y1) == p2, ncol(mat_y2) == p2,
            nrow(mat_y) == n, nrow(mat_x1) == nrow(mat_y2))
  stopifnot(all(idx1 <= n), length(idx1) == nrow(mat_y1), 
            length(idx1) == length(unique(idx1)), all(idx1 > 0), 
            all(idx1 %% 1 == 0))
  stopifnot(!any(rec$vec_from %in% idx1), all(unlist(rec$list_to) %in% idx1))
  
  if(form_options[["method"]] == "literal"){
    res <- .update_estimation_literal(mat_x, mat_y,
                                      mat_x1, mat_y1, mat_y2, idx1,
                                      rec, form_options)
  } else if(form_options[["method"]] == "average"){
    res <- .update_estimation_average(mat_x, mat_y,
                                      mat_x1, mat_y1, mat_y2, idx1,
                                      rec, form_options)
  } else {
    stop("Forming method not found")
  }
  
  stopifnot(all(res$idx1 <= n), length(res$idx1) == nrow(res$mat_y1), 
            length(res$idx1) == length(unique(res$idx1)), all(res$idx1 > 0), 
            all(res$idx1 %% 1 == 0))
  
  res
}

########################

.update_estimation_literal <- function(mat_x, mat_y,
                                       mat_x1, mat_y1, mat_y2, idx1,
                                       rec, form_options){
  p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
  
  # for mat_x1 and mat_y2
  idx_from <- unlist(lapply(1:length(rec$list_to), function(i){
    rep(rec$vec_from, length = length(rec$list_to[[i]]))
  }))
  idx_to <- unlist(rec$list_to)
  stopifnot(length(idx_from) == length(idx_to))
  mat_x1 <- rbind(mat_x1, mat_x[idx_from,,drop = F])
  mat_y2 <- rbind(mat_y2, mat_y[idx_to,,drop = F])
  
  # for mat_y1
  n_org2 <- nrow(mat_y1)
  mat_y1 <- rbind(mat_y1, mat_y[rec$vec_from,,drop = F])
  idx1 <- c(idx1, rec$vec_from)
  
  list(mat_x1 = mat_x1, mat_y1 = mat_y1, mat_y2 = mat_y2, 
       idx1 = idx1)
}

.update_estimation_average <- function(mat_x, mat_y,
                                       mat_x1, mat_y1, mat_y2, idx1,
                                       rec, form_options){
  # p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
  # 
  # # for mat_x1 and mat_y2
  # idx_from <- unlist(lapply(1:length(rec$list_to), function(i){
  #   rep(rec$vec_from, length = length(rec$list_to[[i]]))
  # }))
  # idx_to <- unlist(rec$list_to)
  # stopifnot(length(idx_from) == length(idx_to))
  # mat_x1 <- rbind(mat_x1, mat_x[idx_from,,drop = F])
  # mat_y2 <- rbind(mat_y2, mat_y[idx_to,,drop = F])
  # 
  # # for mat_y1
  # n_org2 <- nrow(mat_y1)
  # mat_y1 <- rbind(mat_y1, mat_y[rec$vec_from,,drop = F])
  # idx1 <- c(idx1, rec$vec_from)
  # 
  # list(mat_x1 = mat_x1, mat_y1 = mat_y1, mat_y2 = mat_y2, 
  #      idx1 = idx1)
}