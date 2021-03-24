#' Generate a simple genome
#'
#' @param p1 positive integer, number of variables in Modality 1
#' @param p2 positive integer, number of variables in Modality 2
#' @param genome_length positive integer, for length of the genome
#' @param window positive integer, for the length of each the \code{p2} variables "shoulder" 
#' where the \code{round(p2/p1)} variables are located
#' @param intercept2 numeric, for baseline signal for each of the \code{p2} variables
#' @param prefix1 string, for the prefix for the \code{p1} variable names
#' @param prefix2 string, for the prefix for the \code{p2} variable names
#'
#' @return two data frames, \code{df_x} and \code{df_y}, having \code{p1} and \code{p2}
#' rows for information of their respective variables
#' @export
generate_df_simple <- function(p1, p2, genome_length = 10*(p1+p2), window = max(floor(genome_length/(4*p2)),1), 
                               intercept2 = 0, prefix1 = "peak", prefix2 = "gene"){
  stopifnot(all(c(p1,p2,genome_length,window) > 0), all(c(p1,p2,genome_length,window) %% 1 == 0), 
            window <= genome_length/2, length(intercept2) %in% c(1, p2), genome_length >= 2*(p1+p2),
            ceiling(p2/p1) < 2*window)
  
  # generate locations for p2
  y_loc <- round(seq(0, genome_length, length.out = p2+2)); y_loc <- y_loc[-c(1,length(y_loc))]
  stopifnot(all(diff(sort(y_loc)) > 0), all(y_loc > 0))
  
  # generate locations for p1
  num_p1inp2 <- rep(floor(p1/p2), length = p2)
  if(sum(num_p1inp2) < p1) num_p1inp2[1:(p1-sum(num_p1inp2))] <- num_p1inp2[1:(p1-sum(num_p1inp2))]+1
  x_loc <- unlist(lapply(1:p2, function(i){
     round(seq(max(y_loc[i]-window, 1), min(y_loc[i]+window,genome_length), 
               length.out = num_p1inp2[i]))
  }))
  stopifnot(all(diff(sort(x_loc)) > 0), all(x_loc > 0))
  
  # generate dfs
  df_x <- data.frame(name = paste0(prefix1, 1:p1), location = x_loc)
  df_y <- data.frame(name = paste0(prefix2, 1:p2), location = y_loc, baseline = intercept2)
  
  list(df_x = df_x, df_y = df_y)
}


#' Generate coefficients forming the deterministic
#' relation from Modality 1 to Modality 2
#' 
#' The data frames \code{df_x} and \code{df_y} can be from \code{generate_df_simple}.
#' The input \code{window} below should ideally that used in \code{generate_df_simple}.
#'
#' @param df_x data frame where each row contains information about each variable in Modality 1
#' @param df_y data frame where each row contains information about each variable in Modality 2
#' @param window positive integer, for the genomic window that variables in Modality 1
#' can affect variables in Modality 2.
#' @param sparsity positive value between 0 and 1, denoting the percentage of 
#' variables in Modality 2 in within \code{window} of the variable in Modality 1 are "turned off"
#' @param signal_mean mean value of each non-zero value
#' @param signal_sd standard deviation of each non-zero value
#'
#' @return matrix with \code{nrow(df_x)} rows and \code{nrow(df_y)} columns 
#' @export
generate_gcoef_simple <- function(df_x, df_y, window = max(floor(0.5*(nrow(df_y)/nrow(df_x))),1), 
                                 sparsity = 0, signal_mean = 3*nrow(df_y)/nrow(df_x), signal_sd = 0){
  stopifnot(length(signal_mean) == 1, length(signal_sd) == 1, signal_sd >= 0,
            is.data.frame(df_x), is.data.frame(df_y), 
            c("name", "location") %in% names(df_x), c("name", "location") %in% names(df_y))
  
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  
  mat_g <- sapply(1:p2, function(i){
    signal_vec <- rep(0, p1)
    
    # find all p1's within window of a particular p2[i]
    idx <- which(abs(df_x$location - df_y$location[i]) <= window)
    if(length(idx) == 0) {return(signal_vec)}
    
    if(sparsity > 0){
      chosen <- stats::rbinom(length(idx), size = 1, prob = 1-sparsity)
      idx <- idx[which(chosen == 1)]
      if(length(idx) == 0) {return(signal_vec)}
    }
    
    signal_vec[idx] <- stats::rnorm(length(idx), mean = signal_mean, sd = signal_sd)
    signal_vec
  })
  
  colnames(mat_g) <- df_y$name
  rownames(mat_g) <- df_x$name
  
  mat_g
  
}

# cascading sigmoids with slope controled by \code{steepness}.
# output mat_1 and mat_2 with \code{nrow(mat_1) = nrow(mat_2) = timepoints}.

#' Generate cascading-from-the-left pattern
#'
#' @param df_mod data frame, for example, either \code{df_x} or \code{df_y} used as
#' input in \code{generate_gcoef_simple}.
#' @param steepness positive numeric denoting the steepness of the logistic function
#' (when the domain is from 0 to 1)
#' @param start_midpoint numeric denoting where the logistic function's mid starts off
#' @param end_midpoint  numeric denoting where the logistic function's mid ends up
#' @param timepoints discretization of how many rows are outputted by this function.
#' The larger this number is, the higher the resolution and the more rows there are 
#' @param max_val maximum value of the output matrix. (The values are rescaled by multiplying
#' by this value)
#'
#' @return matrix with \code{timepoints} rows and \code{nrow(df_mod)} columns
#' @export
generate_traj_cascading <- function(df_mod, steepness = 10, 
                                    start_midpoint = 0, end_midpoint = 1, timepoints = 10*nrow(df_mod),
                                    max_val = 1){
  stopifnot(end_midpoint > start_midpoint, timepoints > 3, timepoints %% 1 == 0)
  
  start <- min(df_mod$location); end <- max(df_mod$location); len <- end-start
  p1 <- nrow(df_mod)
  
  # set up the values on the [0,1] scale
  mid_vec <- seq(start_midpoint, end_midpoint, length.out = timepoints)
  eval_vec <- (df_mod$location-start)/len
  
  mat <- t(sapply(1:timepoints, function(i){
    1-.sigmoid(eval_vec, x0 = mid_vec[i], k = steepness)
  }))
  mat <- mat*max_val
  colnames(mat) <- df_mod$name
  
  mat
}

#####################

# https://en.wikipedia.org/wiki/Logistic_function
.sigmoid <- function(x, x0 = 0, k = 1){
  1/(1+exp(-k*(x-x0)))
}

