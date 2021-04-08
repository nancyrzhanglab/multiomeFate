#' Find candidate cells
#' 
#' Based on the cells previously recruited (recorded in \code{df_res}),
#' find other cells based on \code{mat_x} (the data for Modality 1) and
#' \code{df_res} to potentially be recruited in this current
#' iteration.  
#' 
#' The options are:
#' \itemize{
#' \item \code{nn_xonly_any}: Find the candidates based only on which
#' unrecruited cells are one of the \code{cand_options$nn} nearest neighbors
#' to any previously-recruited cells in Modality 1
#' \item \code{nn_xonly_avg}: Find the \code{cand_options$num_cand} candidates 
#' which are \code{cand_options$num_cand} unrecruited cells which have the closest
#' \code{cand_options$nn} nearest neighbors of previously-recruited cells in Modality 1,
#' where we compute the average distance to said nearest neighbors using
#' \code{cand_options$average}
#' \item \code{all}: All the unrecruited cells are considered as candidates
#' }
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param df_res data frame recording the current results, generated within \code{chromatin_potential}
#' @param cand_options one of the outputs from \code{.chrom_options}
#'
#' @return a vector of integers between 1 and \code{nrow(mat_x)}
.candidate_set <- function(mat_x, df_res, cand_options){
  if(cand_options[["method"]] == "nn_xonly_any"){
    res <- .candidate_set_nn_xonly_any(mat_x, df_res, cand_options)
  } else if(cand_options[["method"]] == "nn_xonly_avg"){
    res <- .candidate_set_nn_xonly_avg(mat_x, df_res, cand_options)
  } else if(cand_options[["method"]] == "all"){
    res <- .candidate_set_all(df_res)
  } else {
    stop("Candidate method not found")
  }
  
  res
}

#######################

.candidate_set_nn_xonly_any <- function(mat_x, df_res, cand_options){
  # extract the indices already recruited
  n <- nrow(df_res)
  idx_free <- which(is.na(df_res$order_rec))
  idx_rec <- which(!is.na(df_res$order_rec))
  
  if(length(idx_free) == 0) return(numeric(0))
  if(length(idx_free) <= cand_options$nn) return(idx_free)
  
  # find the free points that are nearest neighbors to any of the recruited points
  # [[note to self: check if it's worthwhile to expose the ANN KD-tree here]]
  res <- RANN::nn2(mat_x[idx_free,,drop = F], query = mat_x[idx_rec,,drop = F], 
                   k = cand_options$nn)
  
  idx_free[sort(unique(as.numeric(res$nn.idx)))]
}

.candidate_set_nn_xonly_avg <- function(mat_x, df_res, cand_options){
  # extract the indices already recruited
  n <- nrow(df_res)
  idx_free <- which(is.na(df_res$order_rec))
  idx_rec <- which(!is.na(df_res$order_rec))
  
  if(length(idx_free) == 0) return(numeric(0))
  if(length(idx_free) <= cand_options$num_cand) return(idx_free)
  
  # find the free points that are nearest neighbors to any of the recruited points
  # [[note to self: check if it's worthwhile to expose the ANN KD-tree here]]
  res <- RANN::nn2(mat_x[idx_rec,,drop = F], query = mat_x[idx_free,,drop = F], 
                   k = cand_options$nn)
  
  if(cand_options$average == "mean"){
    func <- mean
  } else if(cand_options$average == "median"){
    func <- stats::median
  }else {
    stop("Candidate method (option: 'average') not found")
  }
  
  idx <- order(apply(res$nn.dist, 1, func), decreasing = F)[1:cand_options$num_cand]
  
  idx_free[idx]
}

.candidate_set_all <- function(df_res){
  sort(which(is.na(df_res$order_rec)))
}