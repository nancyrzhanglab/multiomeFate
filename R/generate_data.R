# some inputs: obj_nextcell, the initial point, max_cells
# to add in future: non-uniform density of cells along time
#
# some outputs: data frame of cells [uncorrupted atac and rna] as well as its branch and psuedotime,
# true datapoint mother
# and another data of simply cells with the noisy data
generate_data <- function(obj_next, max_n = 2*length(obj_next$ht), number_runs = 1,
                          time_tol = 0.01, verbose = T){
  # [!! handle number_runs in the future]
  # initialize the noiseless matrix
  p1 <- nrow(obj_next$df_x); p2 <- nrow(obj_next$df_y)
  init_row <- 10
  mat_x <- matrix(NA, nrow = init_row, ncol = p1); colnames(mat_x) <- obj_next$df_x$name
  mat_y <- matrix(NA, nrow = init_row, ncol = p2); colnames(mat_y) <- obj_next$df_y$name
  df_info <- data.frame(time = rep(NA, length = init_row), counter = rep(NA, length = init_row))
  
  n <- 1
  tmp <- obj_next$start_x
  mat_x[n,] <- stats::rbinom(length(tmp), size = 1, prob = tmp)
  tmp <- obj_next$start_y
  mat_y[n,] <- stats::rpois(length(tmp), lambda = tmp)
  df_info[n,"time"] <- 0; df_info[n,"counter"] <- n
  
  # while loop
  while(n < max_n){
    if(verbose) print(paste0("n: ", n, ", time:", round(df_info[n,"time"], 2)))
    if(df_info[n,"time"] > 1-time_tol) break()
      
    ## generate next y, based on x
    mat_y[n+1,] <- generate_ygivenx(obj_next, mat_x[n,])
    
    ## based on the next y, generate that corresponding x
    tmp <- generate_xgiveny(obj_next, y = mat_y[n+1,])
    mat_x[n+1,] <- tmp$x
    df_info[n+1, "time"] <- tmp$time; df_info[n+1, "counter"] <- n+1
    
    ## expand the matrix if necessary
    if(sum(is.na(mat_x[,1])) < 2){
      mat_x <- .expand_matrix(mat_x); mat_y <- .expand_matrix(mat_y)
      df_info <- .expand_df(df_info)
    }
    
    n <- n+1
  }
 
  # [to do: add technical noise to matrices. dropout?]
  # [in the future: df_info should contain who is the mother, what traj?]
  idx <- which(is.na(mat_x[,1]))
  if(length(idx) > 0){
    mat_x <- mat_x[-idx,]; mat_y <- mat_y[-idx,]; df_info <- df_info[-idx,]
  }
  obs_x <- mat_x; obs_y <- mat_y
  
  # prepare output
  list(df_x = obj_next$df_x, df_y = obj_next$df_y,
       obs_x = obs_x, obs_y = obs_y, true_x = mat_x, true_y = mat_y, 
       df_info = df_info)
}

################

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
  # [note: in the future, replace this with an exposed C++ obj from RANN: https://github.com/jefferislab/RANN/blob/master/R/nn.R]
  tmp <- matrix(y, nrow = 1, ncol = length(y))
  idx <- RANN::nn2(obj_next$mat_starty, query = tmp, k = 1)$nn.idx[1,1]
  
  # grab the information from the hash table
  info <- obj_next$ht[[as.character(idx)]]
  
  # [to come: determine which traj to use]
  
  # use logistic regression
  x <- .bernoulli_xgiveny(y, info$list_coef$mat_coef, info$list_coef$vec_intercept)
  x <- stats::rbinom(length(x), size = 1, prob = x)
  
  # prepare output (include the pseudotime from the hashtable)
  # [should record trajectory in the future also]
  list(x = x, time = info$time)
}

#############################

.expand_matrix <- function(mat, scaling = .5){
  n <- nrow(mat); p <- ncol(mat)
  rbind(mat, matrix(NA, nrow = ceiling(scaling*n), ncol = p))
}

.expand_df <- function(df, scaling = .5){
  n <- nrow(df); p <- ncol(df)
  df2 <- matrix(NA, nrow = ceiling(scaling*n), ncol = p)
  df2 <- as.data.frame(df2)
  colnames(df2) <- colnames(df)
  
  rbind(df, df2)
}
