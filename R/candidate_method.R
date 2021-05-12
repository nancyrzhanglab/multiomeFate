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
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param res_g output of \code{.estimate_g}
#' @param df_res data frame recording the current results, generated within \code{chromatin_potential}
#' @param cand_options one of the outputs from \code{.chrom_options}
#'
#' @return a list containing \code{vec_cand} (vector of integers between 1 and \code{nrow(mat_x)})
#' and a list \code{diagnostic} for possible diagnostics
.candidate_set <- function(mat_x, mat_y, df_res, nn_mat, cand_options){
  if(cand_options[["method"]] == "nn_any"){
    res <- .candidate_set_nn_any(df_res, nn_mat, cand_options)
  } else if(cand_options[["method"]] == "nn_freq"){
    res <- .candidate_set_nn_freq(df_res, nn_mat, cand_options)
  } else if(cand_options[["method"]] == "all"){
    res <- .candidate_set_all(df_res, cand_options)
  } else {
    stop("Candidate method not found")
  }
  
  if(cand_options$run_diagnostic){
    # currently nothing here
    res$diagnostic$postprocess <- NA
  }
  
  res
}

#######################

.candidate_set_nn_any <- function(df_res, nn_mat, cand_options){
  # extract the indices of cells recruited in the previous iteration
  idx_free <- which(is.na(df_res$order_rec))
  
  if(cand_options$only_latest){
    val <- max(df_res$order, na.rm = T)
    idx_rec <- which(df_res$order_rec == val)
  } else {
    idx_rec <- which(!is.na(df_res$order_rec))
  }
 
  if(length(idx_free) == 0) return(list(vec_cand = numeric(0), diagnostic = list()))
  if(length(idx_free) <= cand_options$num_cand) return(list(vec_cand = idx_free, diagnostic = list()))
  
  # find the free points that are nearest neighbors to any of the recruited points
  idx_nn <- unique(as.vector(nn_mat[idx_rec,]))
  
  vec_cand <- intersect(idx_nn, idx_free)
  stopifnot(length(vec_cand) > 0)
  
  # run the diagnostic
  list_diagnos <- list() 
  if(cand_options$run_diagnostic){
    # nothing currently here
  }
  
  list(vec_cand = vec_cand, diagnostic = list_diagnos)
}

.candidate_set_nn_freq <- function(df_res, nn_mat, cand_options){
  # extract the indices already recruited
  idx_free <- which(is.na(df_res$order_rec))
  idx_rec <- which(!is.na(df_res$order_rec))
  
  if(length(idx_free) == 0) return(list(vec_cand = numeric(0), diagnostic = list()))
  if(length(idx_free) <= cand_options$num_cand) return(list(vec_cand = idx_free, diagnostic = list()))
  
  # find the free points that are nearest neighbors to any of the recruited points
  res <- as.vector(nn_mat[idx_rec,])
  res <- res[res %in% idx_free]
  tab <- table(res); idx_order <- order(tab, decreasing = T)
  vec_cand <- as.numeric(names(tab[idx_order][1:cand_options$num_cand]))
  
  # run the diagnostic
  list_diagnos <- list() 
  if(cand_options$run_diagnostic){
    # nothing currently here
  }
  
  list(vec_cand = vec_cand, diagnostic = list_diagnos)
}

.candidate_set_all <- function(df_res, cand_options){
  vec_cand <- sort(which(is.na(df_res$order_rec)))
  
  # run the diagnostic
  list_diagnos <- list() 
  if(cand_options$run_diagnostic){
    # nothing currently here
  }
  
  list(vec_cand = vec_cand, diagnostic = list_diagnos)
}