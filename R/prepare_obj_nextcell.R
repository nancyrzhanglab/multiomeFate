# some inputs: 
# - df_x (data frame of inputs for modal 1: name, location)
# - df_y (name, location, baseline expression [the intercept])
# - matrix of regression coef for g, 
# - list of mat_1 [can be proportions between [0,1]], 
# - list of mat_2 [can be proportions between [0,1]],
# - igraph for branching structure
#
# some outputs: essentially all the ingredients for the simulation:
# - df_x, df_y
# - the matrix of coefficients formed by g
# - a huge hash table that stores each unique element of "g(mat_1)"
# and contains 1) probability of which traj to use, and 2) that traj's logistic function coefficients,
# and 3) the psuedotime
# - a huge matrix of all the elements in "g(mat_1)", corresponding to the hash table 
# [note: in the future, replace this with an exposed C++ obj from RANN: https://github.com/jefferislab/RANN/blob/master/R/nn.R]
# WARNING: We'll code as if there's no branching for now. But in the future, it'll prob require putting information in the nodes of 
#  \code{branching_graph}, and we'll need fancy functions to grab the correct rows, etc.
prepare_obj_nextcell <- function(df_x, df_y, mat_g, list_traj_mat, bool_traj_y = T, 
                                 branching_graph = NA, coarseness = 0.1, max_y = 1e5, verbose = T){
  
  if(is.na(branching_graph)) stopifnot(length(list_traj_mat) == 1) else {
    stopifnot(class(branching_graph) == "igraph", igraph::vcount(branching_graph) == length(list_traj_mat), igraph::components(branching_graph)$no == 1)
  }
  stopifnot(all(mat_g >= 0)) # [note to self: needed for now for simplicity -- should be removed in the future]
  stopifnot(all(matrix(1, nrow = nrow(list_traj_mat[[1]]), ncol = nrow(mat_g)) %*% mat_g >= list_traj_mat[[1]])) # [note to self : needed for now for simplicity -- assumes only one trajectory]
    
  # [note to self: need a check to make sure resolution is not too large wrt number of rows in list_x1's matrices]
  
  # [note to self: put a function here to check branching_graph wrt list_traj_mat.
  # for example, the first row of list_traj_mat when a branch occurs are the same]
  
  # enumerate all unique df_y's needed for the hash table
  ## [note to self: hard-code the fact it's a linear trajectory. we'll need to use paths from start to end later on -- how would this capture cycles?]
  if(!bool_traj_y){
    # if trajectory is for df_x
    mat_startx <- do.call(rbind, list_traj_mat) #x1
    mat_starty <- t(sapply(1:nrow(list_traj_mat[[1]]), function(i){
      .possion_ygivenx(list_traj_mat[[1]][i,], mat_g, max_val = max_y)
    })) # y2
  } else {
    # if trajectory is for df_y
    mat_startx <- .compute_xfromy(list_traj_mat, mat_g) #x1
    mat_starty <- do.call(rbind, list_traj_mat) #y2
  }
  n_total <- nrow(mat_starty) # count how many unique rows there are
  list_time <- list(seq(0, 1, length.out = n_total)) # [note to self: currently hard-coded for linear trajectory]
   
  # initialize hash table
  ht <- hash::hash()
  counter <- 1
  for(k in 1:length(list_traj_mat)){
    for(i in 1:nrow(list_traj_mat[[k]])){
      if(verbose && i %% floor(nrow(list_traj_mat[[k]])/10) == 0) cat('*')
      
      ## grab the relevant rows
      ## [note to self: hard-code the fact it's a linear trajectory]
      idx <- c(max(round(i-coarseness*n_total), 1):min(round(i+coarseness*n_total), n_total-1))
      mat_y2 <- mat_starty[idx,,drop = F] 
      mat_x2 <- mat_startx[idx+1,,drop = F] 
     
      ## perform logistic regression
      list_coef <- .glmnet_logistic(mat_y2, mat_x2)
      
      ## store the values 
      ht[[as.character(counter)]] <- list(list_coef = list_coef, time = list_time[[k]][i])
      
      counter <- counter+1
    }
  }
  stopifnot(counter-1 == n_total)
 
  # prepare outputs
  start_x <- list_x1[[1]][1,]
  start_y <- .possion_ygivenx(start_x, mat_g, max_val = max_y)
  structure(list(df_x = df_x, df_y = df_y, mat_g = mat_g, ht = ht, mat_starty = mat_starty,
                 start_x = start_x, start_y = start_y),
            class = "obj_next")
}

####################################

# [note to self: this entire function is currently coded assuming one trajectory]
.compute_xfromy <- function(list_traj_mat, mat_g){
  p1 <- nrow(mat_g)
  
  # compute starting x
  mat_startx <- matrix(0, nrow = nrow(list_traj_mat[[1]]), ncol = nrow(mat_g))
  mat_startx[1,] <- .compute_xfromy_starting(list_traj_mat[[1]][1,], mat_g)
  
  # compute the next x's
  for(i in 2:p1){
    mat_startx[i,] <- .compute_xfromy_next(mat_startx[i-1,], list_traj_mat[[1]][i,], mat_g)
  }
  
  mat_startx
}

.compute_xfromy_starting <- function(vec_y, mat_g, tol = 1e-6){
  stopifnot(ncol(mat_g) == length(vec_y))
  p1 <- nrow(mat_g); p2 <- ncol(mat_g)
  
  obj <- rep(1, p1)
  A <- t(mat_g)
  constr_lb <- vec_y - tol
  constr_ub <- vec_y + tol
  var_lb <- rep(0, p1)
  var_ub <- rep(1, p1)
  
  res <- suppressWarnings(clplite::clp_solve(obj, A, constr_lb, constr_ub, 
                                             var_lb, var_ub, max = FALSE))
  stopifnot(res$status == 0)
  vec <- res$solution
  vec[vec <= tol] <- 0
  vec
}

# min |x2 - x1| : G*x2 = y2, and x2 in [0,1]
#  where x1 (vector of length p1) is the previous ATAC, x2 (vector of length p1)
#    is the current ATAC we wish to solve for
#  and G is the linear function that maps the current ATAC to RNA,
#  and y (vector of length p2) is the current RNA we're given
#
# reformulate to:
# min e_1 + ... + e_p1 : e_i >= x2_i - x1_i and  -e_i <= -(x2_i - x1_i) for all i in {1,...,p1}
#  and G*x2 = y2, and x2 in [0,1]
# 
# reformulate to:
# min e_1 + ... + e_p1 : x1_i >= x2_i - e_i and -x1_i <= -x2_i + e_i for all i in {1,...,p1}
#  and G*x2 = y2, and x2 in [0,1]
# put the x2's first, then e's
.compute_xfromy_next <- function(vec_x, vec_y, mat_g, tol = 1e-6){
  stopifnot(all(vec_x >= 0), all(vec_x <= 1), length(vec_x) == nrow(mat_g), length(vec_y) == ncol(mat_g))
  p1 <- nrow(mat_g); p2 <- ncol(mat_g)
  
  obj <- c(rep(0, p1), rep(1, p1))
  
  A <- matrix(0, nrow = 2*p1+p2, ncol = 2*p1)
  for(i in 1:p1){
    A[i, c(i, i+p1)] <- c(1,-1)
  }
  for(i in 1:p1){
    A[i+p1, c(i, i+p1)] <- c(-1,1)
  }
  A[((2*p1+1):(2*p1+p2)),1:p1] <- t(mat_g)
  
  constr_lb <- c(rep(-Inf, p1), -vec_x, vec_y - tol)
  constr_ub <- c(vec_x, rep(Inf, p1), vec_y + tol)
  var_lb <- rep(0, 2*p1)
  var_ub <- c(rep(1, p1), rep(Inf, p1))
  
  res <- suppressWarnings(clplite::clp_solve(obj, A, constr_lb, constr_ub, 
                                             var_lb, var_ub, max = FALSE))
  stopifnot(res$status == 0)
  vec <- res$solution[1:p1]
  vec[vec <= tol] <- 0
  vec
}

####################

.possion_ygivenx <- function(x, mat_g, max_val = 1e5){
  pmin(as.numeric(exp(x %*% mat_g)), max_val)
}

.bernoulli_xgiveny <- function(y, mat_coef, vec_intercept){
  p1 <- ncol(mat_coef)
  
  sapply(1:p1, function(i){
    val <- as.numeric(y%*%mat_coef[,i]) + vec_intercept[i]
    .sigmoid(val)
  })
}

.glmnet_logistic <- function(covariate, response_prob){
  x <- covariate
  p1 <- ncol(response_prob); p2 <- ncol(covariate)
  mat_coef <- matrix(NA, nrow = p2, ncol = p1)
  vec_intercept <- rep(NA, length = ncol(response_prob))
  
  for(i in 1:p1){
    y_mat <- cbind(1-response_prob[,i], response_prob[,i])
    fit <- glmnet::glmnet(x = x, y = y_mat, family = "binomial",
                          standardize = FALSE, intercept = TRUE)
    
    len <- length(fit$lambda)
    mat_coef[,i] <- fit$beta[,len]
    vec_intercept[i] <- fit$a0[len]
  }
  
  list(mat_coef = mat_coef, vec_intercept = vec_intercept)
}