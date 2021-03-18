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
prepare_obj_nextcell <- function(df_x, df_y, mat_g, list_x1, list_x2, 
                                 branching_graph = NA, coarseness = 0.1, max_y = 1e5, verbose = T){
  
  stopifnot(length(list_x1) == length(list_x2), all(sapply(length(list_x1), function(i){all(dim(list_x1[[i]]) == dim(list_x2[[i]]))})))
  if(is.na(branching_graph)) stopifnot(length(list_x1) == 1) else {
    stopifnot(class(branching_graph) == "igraph", igraph::vcount(branching_graph) == length(list_x1), igraph::components(branching_graph)$no == 1)
  }
  # [need a check to make sure resolution is not too large wrt number of rows in list_x1's matrices]
  
  # [in the future: put a function here to check branching_graph wrt list_x1 and list_x2.
  # for example, the first row of list_x1 when a branch occurs are the same]
  
  # compute pseudotime of each rows in list_x1
  ## [for now: hard-code the fact it's a linear trajectory. we'll need to use paths from start to end later on -- how would this capture cycles?]
  ## count how many unique rows there are
  mat_starty <- t(sapply(1:nrow(list_x1[[1]]), function(i){
    .possion_ygivenx(x = list_x1[[1]][i,], mat_g, max_val = max_y)
  }))
  n_total <- nrow(mat_starty)
  list_time <- list(seq(0, 1, length.out = n_total)) # !!
  
  # initialize hash table
  ht <- hash::hash()
  counter <- 1
  for(k in 1:length(list_x1)){
    for(i in 1:nrow(list_x1[[k]])){
      if(verbose && i %% floor(nrow(list_x1[[k]])/10) == 0) cat('*')
      
      ## grab the relevant rows 
      ## [for now: hard-code the fact it's a linear trajectory]
      idx <- c(max(round(i-coarseness*n_total), 1):min(round(i+coarseness*n_total), n_total))
      mat_x1 <- list_x1[[k]][idx,,drop = F]
      mat_y2 <- t(apply(mat_x1, 1, function(x){.possion_ygivenx(x, mat_g, max_val = max_y)}))
      mat_x2 <- list_x2[[k]][idx,,drop = F]
      
      ## perform logistic regression
      list_coef <- .glmnet_logistic(mat_y2, mat_x2)
      
      ## store the values 
      ht[[as.character(counter)]] <- list(list_coef = list_coef, time = list_time[[k]][i])
      
      counter <- counter+1
    }
  }
  stopifnot(counter-1 == n_total)
 
  # prepare outputs
  structure(list(df_x = df_x, df_y = df_y, mat_g = mat_g, ht = ht, mat_starty = mat_starty),
            class = "obj_next")
}

generate_ygivenx <- function(obj_next, x){
  stopifnot(class(obj_next) == "obj_next")
  
  # generate y from x
  y <- .possion_ygivenx(x, obj_next$mat_g)
  
  # add the intercepts given by dat_y
  y <- y + obj_next$df_y$baseline
  
  # generate poisson
  stats::rpois(length(y), lambda = y)
}

generate_xgiveny <- function(obj_next, y){
  stopifnot(class(obj_next) == "obj_next")
  
  # find the nearest neighbor 
  tmp <- matrix(y, nrow = 1, ncol = length(y))
  idx <- RANN::nn2(obj_next$mat_starty, query = tmp, k = 1)$nn.idx[1,1]
  
  # grab the information from the hash table
  info <- obj_next$ht[[as.character(idx)]]
  
  # [to come: determine which traj to use]
  
  # use logistic regression
  x <- .bernoulli_xgiveny(y, info$list_coef$mat_coef, info$list_coef$vec_intercept)
  x <- stats::rbinom(length(x), size = 1, prob = x)
  
  # prepare output (include the pseudotime from the hashtable)
  list(x = x, time = info$time)
}

####################################

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