# generate a fake genome of length \code{genome_length} where we put 
# \code{p2} genes uniformly spaced out, and \code{round(p2/p1)} peaks around each gene
generaate_df_simple <- function(p1, p2, window = max(floor(0.5*(p2/p1)),1), genome_length = 100){
  
}

# generate a simple matrix of regression coefficients
generate_gmat_simple <- function(df_x, df_y, window = max(floor(0.5*(nrow(df_y)/nrow(df_x))),1), 
                                 sparsity = 0, signal_mean = 1, signal_sd = 0){
  
}

# cascading sigmoids with slope controled by \code{steepness}.
# output mat_1 and mat_2 with \code{nrow(mat_1) = nrow(mat_2) = resolution}.
generate_traj_cascading <- function(df_x, steepness = 1, 
                                    start_midpoint = 0, end_midpoint = 1, resolution = 10){
  
}