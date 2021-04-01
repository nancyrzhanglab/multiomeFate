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

.update_estimation_matrices <- function(mat_x1, mat_y1, mat_y2, idx1, idx2, 
                                        res, form_options){
  stopifnot(form_options[["method"]] == "literal")
  
}