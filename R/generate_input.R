# generate a fake genome of length \code{genome_length} where we put 
# \code{p2} genes uniformly spaced out, and \code{round(p2/p1)} peaks around each gene
generate_df_simple <- function(p1, p2, genome_length = 10*(p1+p2), window = max(floor(genome_length/(4*p2)),1), 
                                intercept2 = 0, prefix1 = "peak", prefix2 = "gene"){
  stopifnot(all(c(p1,p2,genome_length,window) > 0), all(c(p1,p2,genome_length,window) %% 1 == 0), 
            window <= genome_length, length(intercept2) %in% c(1, p2), genome_length >= 2*(p1+p2))
  
  # generate locations for p2
  y_loc <- round(seq(0, genome_length, length.out = p2+2)); y_loc <- y_loc[-c(1,length(y_loc))]
  stopifnot(all(diff(sort(y_loc)) > 0))
  
  # generate locations for p1
  num_p1inp2 <- rep(floor(p1/p2), length = p2)
  if(sum(num_p1inp2) < p1) num_p1inp2[1:(p1-sum(num_p1inp2))] <- num_p1inp2[1:(p1-sum(num_p1inp2))]+1
  x_loc <- unlist(lapply(1:p2, function(i){
    y_loc[i] + round(seq(-window, window, length.out = num_p1inp2[i]))
  }))
  stopifnot(all(diff(sort(x_loc)) > 0))
  
  # generate dfs
  df_x <- data.frame(name = paste0(prefix1, 1:p1), location = x_loc)
  df_y <- data.frame(name = paste0(prefix1, 1:p2), location = y_loc, baseline = intercept2)
  
  list(df_x = df_x, df_y = df_y)
}

# generate a simple matrix of regression coefficients
generate_gmat_simple <- function(df_x, df_y, window = max(floor(0.5*(nrow(df_y)/nrow(df_x))),1), 
                                 sparsity = 0, signal_mean = 5*nrow(df_y)/nrow(df_x), signal_sd = 0){
  stopifnot(length(signal_mean) == 1, length(signal_sd) == 1, signal_sd >= 0)
  
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  
  g_mat <- sapply(1:p2, function(i){
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
  
}

# cascading sigmoids with slope controled by \code{steepness}.
# output mat_1 and mat_2 with \code{nrow(mat_1) = nrow(mat_2) = resolution}.
generate_traj_cascading <- function(df_x, steepness = 1, 
                                    start_midpoint = 0, end_midpoint = 1, resolution = 10){
  
}