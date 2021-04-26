#' Generate object to prepare for simulating data
#' 
#' If \code{list_traj_mat} is used to describe the evoluation of \code{df_x}, then the values in \code{list_traj_mat}
#' currently denote the expected probability of a Bernoulli.
#' If \code{list_traj_mat} is used  to describe the evoluation of \code{df_y}, then the values in \code{list_traj_mat}
#' currently denote the expected probability of a Poisson. Importantly, this means
#' \code{log} of these values denote the natural parameters. It is important then for
#' \code{list_traj_mat} to have values 1 or more in this case.
#' 
#' Additional requirements: \code{df_x} and \code{df_y} are required to be data frames
#' with a column named \code{name}. The output of this function is meant as the input for
#' \code{generate_data}. 
#' 
#' Notes to self: 1) All entries in \code{mat_g} currently need to be non-negative,
#' 2) Only linear trajectories are allowed currently, 3) Modality 1 (represented by
#' \code{df_x}) is assumed be Bernoulli and Modality 2 (represented by \code{df_y})
#' is assumed by Poisson. 4) Currently, the order of names in \code{mat_g} and 
#' \code{list_traj_mat} must match the exact order of names in \code{df_x} and \code{df_y}
#'
#' @param df_x data frame where each row contains information about each variable in Modality 1
#' @param df_y data frame where each row contains information about each variable in Modality 2
#' @param mat_g matrix with \code{nrow(df_x)} rows and \code{nrow(df_y)} columns containing the linear coefficients on how
#' Modality 1's value affects the mean of Modality 2's values (i.e., the deterministic
#' relation from Modality 1 to Modality 2)
#' @param list_traj_mat list of matrices, where each matrix contains either \code{nrow(df_x)}
#' or \code{nrow(df_y)} columns. Each element of the list represents a different
#' segment in the trajectories. Each of the matrices within the list must have
#' column names matching either \code{df_x$name} or \code{df_y$name}.
#' @param branching_graph \code{igraph} object describing the branching structure
#' @param coarseness number between 0 and 1, describing how many rows in \code{list_traj_mat}
#' are used to fit the probabilistic relation from Modality 2 to Modality 1. The larger the number,
#' the more coarse (i.e., more approximate) the probabilistic relation is, but potentially
#' computationally more stable.
#' @param max_y positive integer, denoting the maximum value in Modality 2
#' @param verbose boolean
#'
#' @return object of class \code{mf_obj_next} containing the following components:
#' \code{df_x} and \code{df_y} and \code{mat_g} (all same as input), 
#' \code{ht} (a \code{hash} object storing all the pseudotimes and logistic regression coefficient
#' comprising the probabilistic relation from Modality 2 to Modality 1),
#' \code{mat_y2all} (a matrix with number of rows equal to \code{length(ht)},
#' where row \code{i} in this matrix is the Modality 2 vector corresponding to
#' the information in \code{ht[[as.character(i)]]}, and number of columns equal to
#' \code{nrow(df_y)})
#' \code{vec_startx} (the initial value for the cell in Modality 1),
#' \code{vec_starty} (the initial value for the cell in Modality 2)
#' @export
prepare_obj_nextcell <- function(df_x, df_y, mat_g, list_traj_mat, 
                                 branching_graph = NA, coarseness = 0.1, 
                                 max_y = 1e5, verbose = T){
  # run a lot of checks
  stopifnot(coarseness > 0, coarseness <= 1, is.data.frame(df_x), is.data.frame(df_y),
            "name" %in% names(df_x), "name" %in% names(df_y), 
            length(unique(sapply(list_traj_mat, ncol))) == 1)
  stopifnot(length(c(df_x$name, df_y$name)) == length(unique(c(df_x$name, df_y$name))),
            all(colnames(list_traj_mat[[1]]) == df_x$name) || 
              all(colnames(list_traj_mat[[1]]) == df_y$name))
  if(length(list_traj_mat) > 1){
    for(i in 2:length(list_traj_mat)){
      all(colnames(list_traj_mat[[1]]) == colnames(list_traj_mat[[i]]))
    }
  }
  stopifnot(all(rownames(mat_g) == df_x$name), all(colnames(mat_g) == df_y$name))
  
  if(is.na(branching_graph)) stopifnot(length(list_traj_mat) == 1) else {
    stopifnot(class(branching_graph) == "igraph", igraph::vcount(branching_graph) == length(list_traj_mat), igraph::components(branching_graph)$no == 1)
  }
  
  # [note to self vv: needed for now for simplicity -- should be removed in the future]
  stopifnot(all(mat_g >= 0)) 
  if(all(sort(colnames(list_traj_mat[[1]])) == sort(df_y$name))) {
    bool_traj_y <- T
    # [note to self vv: needed for now for simplicity -- assumes only one trajectory, and poisson data]
    tmp <- matrix(1, nrow = nrow(list_traj_mat[[1]]), ncol = nrow(mat_g)) 
    stopifnot(all(tmp %*% mat_g >= log(list_traj_mat[[1]])))
  } else {
    bool_traj_y <- F
  }
    
  # [note to self: need a check to make sure resolution is not too large wrt number of rows in list_x1's matrices]
  
  # [note to self: put a function here to check branching_graph wrt list_traj_mat.
  # for example, the first row of list_traj_mat when a branch occurs are the same]
  
  ####
  # some preparation
  
  # enumerate all unique df_y's needed for the hash table
  ## [note to self: hard-code the fact it's a linear trajectory. we'll need to use paths from start to end later on -- how would this capture cycles?]
  if(!bool_traj_y){
    # if trajectory is for df_x
    mat_x1all <- do.call(rbind, list_traj_mat) 
    mat_y2all <- t(sapply(1:nrow(list_traj_mat[[1]]), function(i){
      .possion_ygivenx(list_traj_mat[[1]][i,], mat_g, max_val = max_y)
    })) 
  } else {
    # if trajectory is for df_y
    mat_x1all <- .compute_xfromy(list_traj_mat, mat_g) #x1
    mat_y2all <- do.call(rbind, list_traj_mat) #y2
  }
  n_total <- nrow(mat_y2all) # count how many unique rows there are
  list_time <- list(seq(0, 1, length.out = n_total)) # [note to self: currently hard-coded for linear trajectory]
   
  ################
  
  # initialize hash table
  ht <- hash::hash()
  counter <- 1
  for(k in 1:length(list_traj_mat)){
    for(i in 1:nrow(list_traj_mat[[k]])){
      if(verbose && nrow(list_traj_mat[[k]]) > 10 && i %% floor(nrow(list_traj_mat[[k]])/10) == 0) cat('*')
      
      ## grab the relevant rows
      ## [note to self: hard-code the fact it's a linear trajectory]
      idx <- c(max(round(i-coarseness*n_total), 1):min(round(i+coarseness*n_total), n_total-1))
      mat_y2 <- mat_y2all[idx,,drop = F] 
      mat_x2 <- mat_x1all[idx+1,,drop = F] 
     
      ## perform logistic regression
      list_coef <- .glmnet_logistic(mat_y2, mat_x2)
      
      ## store the values 
      ht[[as.character(counter)]] <- list(list_coef = list_coef, time = list_time[[k]][i])
      
      counter <- counter+1
    }
  }
  stopifnot(counter-1 == n_total)
 
  # prepare outputs
  # [note to self: make this more flexible]
  vec_startx <- mat_x1all[2,]
  vec_starty <- mat_y2all[2,]
  structure(list(df_x = df_x, df_y = df_y, mat_g = mat_g, ht = ht, mat_y2all = mat_y2all,
                 vec_startx = vec_startx, vec_starty = vec_starty),
            class = "mf_obj_next")
} 

####################################

# [note to self: this entire function is currently coded assuming one trajectory,
#  and also assumes poisson data]
.compute_xfromy <- function(list_traj_mat, mat_g){
  p1 <- nrow(mat_g); n <- nrow(list_traj_mat[[1]])
  
  # compute starting x
  mat_startx <- matrix(0, nrow = n, ncol = nrow(mat_g))
  mat_startx[1,] <- .compute_xfromy_starting(log(list_traj_mat[[1]][1,]), mat_g)
  
  # compute the next x's
  for(i in 2:n){
    mat_startx[i,] <- .compute_xfromy_next(mat_startx[i-1,], log(list_traj_mat[[1]][i,]), mat_g)
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

# [note to self: not sure what sensible tests there are for this. maybe tests set up for stationary points?]
.bernoulli_xgiveny <- function(y, mat_coef, vec_intercept){
  p1 <- ncol(mat_coef)
  
  sapply(1:p1, function(i){
    if(all(is.na(mat_coef[,i]))){
      vec_intercept[i]
    } else {
      val <- as.numeric(y%*%mat_coef[,i]) + vec_intercept[i]
      .sigmoid(val)
    }
  })
}

.glmnet_logistic <- function(covariate, response_prob, tol = 1e-3){
  x <- covariate
  p1 <- ncol(response_prob); p2 <- ncol(covariate)
  mat_coef <- matrix(NA, nrow = p2, ncol = p1)
  vec_intercept <- rep(NA, length = ncol(response_prob))
  
  for(i in 1:p1){
    if (diff(range(response_prob[,i])) <= tol){
      mat_coef[,i] <- NA; vec_intercept[i] <- stats::median(response_prob[,i])
      
    } else {
      y_mat <- cbind(1-response_prob[,i], response_prob[,i])
      fit <- glmnet::glmnet(x = x, y = y_mat, family = "binomial",
                            standardize = FALSE, intercept = TRUE, alpha = 0)
      
      len <- length(fit$lambda)
      mat_coef[,i] <- fit$beta[,len]
      vec_intercept[i] <- fit$a0[len]
    }
  }
  
  list(mat_coef = mat_coef, vec_intercept = vec_intercept)
}