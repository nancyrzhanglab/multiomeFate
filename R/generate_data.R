#' Generating data
#' 
#' This is under development
#' 
#' Notes to self: 1) Currently, the function does not add technical noise.
#' 2) It currently only generates data for Modality 1 as Bernoullis and for Modality 2 as Poissons.
#'
#' @param obj_next returned by \code{prepare_obj_nextcell}
#' @param max_n positive integer, for maximum number of cells for each of the \code{number_runs} runs
#' @param number_runs positive integer, for number of times we start at \code{obj_next$vec_startx} and follow that cell's trajectory
#' @param sample_perc number between \code{0} and \code{1}, denoting what percentage of cells we've generated in this simulation
#' are returned in the final output. Here, \code{1} means all the cells generated are returned.
#' @param time_tol positive integer between  \code{0} and \code{1}, typically close to \code{0},
#' denoting the stopping criteria for each of the \code{number_runs} runs. A run
#' stops when it generates a cell whose pseudotime (as dictated in \code{obj_next})
#' is larger than \code{1-time_tol}.
#' @param verbose boolean
#'
#' @return object of class \code{mf_simul} with the following components: 
#' \code{df_x} (the data frame from \code{obj_next} containing information of Modality 1),
#' \code{df_y} (the data frame from \code{obj_next} containing information of Modality 2),
#' \code{obs_x} (the observed data matrix where each row is a cell and each column is one of \code{nrow(df_x)} variables),
#' \code{obs_y} (the observed data matrix where each row is a cell and each column is one of \code{nrow(df_y)} variables),
#' \code{true_x} (the true data matrix, which is \code{obs_x} without the technical noise),
#' \code{true_y} (the true data matrix, which is \code{obs_y} without the technical noise), and
#' \code{df_info} (the data frame with the same number of rows as \code{obs_x} and \code{obs_y}, where each row
#' contains meta-information about the cell)
#' @export
generate_data <- function(obj_next, max_n = 2*length(obj_next$ht), number_runs = 5,
                          sample_perc = 1,
                          time_tol = 0.01, verbose = T){
  stopifnot(sample_perc > 0, sample_perc <= 1, class(obj_next) == "mf_obj_next")
  list_out <- lapply(1:number_runs, function(i){
    .generate_data_single(obj_next, max_n, time_tol, verbose)
  })
  
  # merge all outputs
  n_tot <- sum(sapply(list_out, function(x){nrow(x$df_info)}))
  idx <- sort(sample(1:n_tot, replace = F, size = ceiling(n_tot*sample_perc)))
  res <- .merge_run_outputs(list_out, idx)
  
  # prepare output
  structure(list(df_x = obj_next$df_x, df_y = obj_next$df_y,
       obs_x = res$obs_x, obs_y = res$obs_y, 
       true_x = res$true_x, true_y = res$true_y, 
       df_info = res$df_info), class = "mf_simul")
}

################

.generate_data_single <- function(obj_next, max_n = 2*length(obj_next$ht), time_tol = 0.01,
                                  verbose = T){
  stopifnot(time_tol > 0, time_tol <= 1)
  
  # initialize
  p1 <- nrow(obj_next$df_x); p2 <- nrow(obj_next$df_y)
  init_row <- 10
  mat_x <- matrix(NA, nrow = init_row, ncol = p1); colnames(mat_x) <- obj_next$df_x$name
  mat_y <- matrix(NA, nrow = init_row, ncol = p2); colnames(mat_y) <- obj_next$df_y$name
  df_info <- data.frame(time = rep(NA, length = init_row), counter = rep(NA, length = init_row))
  
  n <- 1
  tmp <- obj_next$vec_startx
  mat_x[n,] <- stats::rbinom(length(tmp), size = 1, prob = tmp)
  tmp <- obj_next$vec_starty
  mat_y[n,] <- stats::rpois(length(tmp), lambda = tmp)
  df_info[n,"time"] <- 0; df_info[n,"counter"] <- n
  
  # while loop
  while(n < max_n){
    # [note to self: This print statement could be improved]
    if(verbose) print(paste0("n: ", n, ", time:", round(df_info[n,"time"], 2)))
    if(df_info[n,"time"] > 1-time_tol) break()
    
    ## generate next y, based on x
    mat_y[n+1,] <- .generate_ygivenx(obj_next, mat_x[n,])
    
    ## based on the next y, generate that corresponding x
    tmp <- .generate_xgiveny(obj_next, y = mat_y[n+1,])
    mat_x[n+1,] <- tmp$x
    df_info[n+1, "time"] <- tmp$time; df_info[n+1, "counter"] <- n+1
    
    ## expand the matrix if necessary
    if(sum(is.na(mat_x[,1])) < 2){
      mat_x <- .expand_matrix(mat_x); mat_y <- .expand_matrix(mat_y)
      df_info <- .expand_df(df_info)
    }
    
    n <- n+1
  }
  
  # [note to self: df_info should contain who is the mother, what traj?]
  idx <- which(is.na(mat_x[,1]))
  if(length(idx) > 0){
    mat_x <- mat_x[-idx,]; mat_y <- mat_y[-idx,]; df_info <- df_info[-idx,]
  }
  # [note to self: need to add technical noise]
  obs_x <- mat_x; obs_y <- mat_y
  
  list(obs_x = obs_x, obs_y = obs_y, true_x = mat_x, true_y = mat_y, 
        df_info = df_info)
}

.generate_ygivenx <- function(obj_next, x){
  stopifnot(class(obj_next) == "mf_obj_next")
  
  # generate y from x
  y <- .possion_ygivenx(x, obj_next$mat_g)
  
  # add the intercepts given by dat_y
  y <- y + obj_next$df_y$baseline
  
  # generate poisson
  stats::rpois(length(y), lambda = y)
}

.generate_xgiveny <- function(obj_next, y){
  stopifnot(class(obj_next) == "mf_obj_next")
  
  # find the nearest neighbor 
  # [note: in the future, replace this with an exposed C++ obj from RANN: https://github.com/jefferislab/RANN/blob/master/R/nn.R]
  tmp <- matrix(y, nrow = 1, ncol = length(y))
  idx <- RANN::nn2(obj_next$mat_y2all, query = tmp, k = 1)$nn.idx[1,1]
  
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

.merge_run_outputs <- function(list_out, idx = NA){
  obs_x <- do.call(rbind, lapply(list_out, function(x){x$obs_x}))
  obs_y <- do.call(rbind, lapply(list_out, function(x){x$obs_y}))
  true_x <- do.call(rbind, lapply(list_out, function(x){x$true_x}))
  true_y <- do.call(rbind, lapply(list_out, function(x){x$true_y}))
  df_info <- do.call(rbind, lapply(list_out, function(x){x$df_info}))
  
  if(!any(is.na(idx))){
    obs_x <- obs_x[idx,,drop = F]; obs_y <- obs_y[idx,,drop = F]
    true_x <- true_x[idx,,drop = F]; true_y <- true_y[idx,,drop = F]
    df_info <- df_info[idx,,drop = F]
  }
  
  list(obs_x = obs_x, obs_y = obs_y, true_x = true_x, true_y = true_y,
       df_info = df_info)
}

.expand_matrix <- function(mat, scaling = .5){
  stopifnot(scaling > 0)
  
  n <- nrow(mat); p <- ncol(mat)
  rbind(mat, matrix(NA, nrow = ceiling(scaling*n), ncol = p))
}

.expand_df <- function(df, scaling = .5){
  stopifnot(scaling > 0)
  
  n <- nrow(df); p <- ncol(df)
  df2 <- matrix(NA, nrow = ceiling(scaling*n), ncol = p)
  df2 <- as.data.frame(df2)
  colnames(df2) <- colnames(df)
  
  rbind(df, df2)
}
