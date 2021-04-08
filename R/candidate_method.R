.candidate_set <- function(mat_x, df_res, cand_options){
  if(cand_options[["method"]] == "nn_xonly"){
    res <- .candidate_set_nn_xonly(mat_x, df_res, cand_options)
  } else {
    stop("Candidate method not found")
  }
  
  res
}

#######################

.candidate_set_nn_xonly <- function(mat_x, df_res, cand_options){
  # extract the indices already recruited
  n <- nrow(df_res)
  idx_free <- which(is.na(df_res$order_rec))
  idx_rec <- which(!is.na(df_res$order_rec))
  
  if(length(idx_free) == 0) return(numeric(0))
  if(length(idx_free) <= cand_options$nn) return(idx_free)
  
  # find the free points that are nearest neighbors to any of the recruited points
  res <- RANN::nn2(mat_x[idx_free,,drop = F], query = mat_x[idx_rec,,drop = F], 
                   k = cand_options$nn)
  
  idx_free[sort(unique(as.numeric(res$nn.idx)))]
}