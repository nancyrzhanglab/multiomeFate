#' Based on the candidate cells, recruit a specific cell
#' 
#' Here, \code{mat_x} are all the cells in our dataset, among which
#' we are saying that \code{idx1} are cells that have previously been
#' recruited already, and \code{vec_cand} are the currently-selected
#' candidates that we are hoping to recruit and match to one of the
#' previously-recruited cells in \code{idx1}. This 
#' matching is done by using our estimation function \code{res_g}
#' on the candidate cell's Modality 1 expression in \code{mat_x}
#' and seeing if it's close to any of the Modality 2 expressions
#' in \code{mat_y1}.
#' 
#' The options are:
#' \itemize{
#' \item \code{nn_yonly}: Recruit \code{rec_options$num_rec} cells
#' whose predicted Modality 2 expression has the smallest average
#' (for example, mean, or in general, depending on what \code{rec_options$average}
#' is set to) distance to its \code{rec_options$nn} previously-recruited
#' nearest neighbors based on \code{mat_y1}
#' }
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param vec_cand output of \code{.candidate_set}
#' @param mat_y1 output of \code{.init_est_matrices} or \code{.update_estimation_matrices}, representing
#' the data for Modality 1
#' @param idx1 output of \code{.init_est_matrices} or \code{.update_estimation_matrices}, representing
#' the data for Modality 1
#' @param res_g output of \code{.estimate_g}
#' @param df_res data frame recording the current results, generated within \code{chromatin_potential}
#' @param rec_options one of the outputs from \code{.chrom_options}
#'
#' @return a list of vector of integers \code{vec_from} and a list
#' of integers \code{list_to}
.recruit_next <- function(mat_x, vec_cand, mat_y1, idx1, res_g, df_res, 
                          rec_options){
  stopifnot(all(idx1 <= nrow(mat_x)), length(idx1) == nrow(mat_y1), length(idx1) == length(unique(idx1)),
            all(idx1 %% 1 == 0), all(idx1 > 0), all(idx1 <= nrow(mat_x)))
  stopifnot(all(vec_cand %% 1 == 0), all(vec_cand > 0), all(vec_cand <= nrow(mat_x)),
            length(vec_cand) == length(unique(vec_cand)))
  stopifnot(!any(vec_cand %in% idx1))
  stopifnot(all(is.na(df_res$order_rec[vec_cand])), !any(is.na(df_res$order_rec[idx1])))
  
  if(rec_options[["method"]] == "nn_yonly"){
    res <- .recruit_next_nn_yonly(mat_x, vec_cand, mat_y1, idx1, res_g, rec_options)
  } else {
    stop("Recruit method not found")
  }

  res
}

###################

.recruit_next_nn_yonly <- function(mat_x, vec_cand, mat_y1, idx1, res_g, 
                                    rec_options){
  num_rec <- min(rec_options$num_rec, length(vec_cand))
  nn <- min(c(rec_options$nn, ceiling(nrow(mat_y1)/2)))
  
  # apply mat_g to mat_x
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g)
  
  # see which prediction is closest to mat_y1
  # [[note to self: check if it's worthwhile to expose the ANN KD-tree here]]
  res <- RANN::nn2(mat_y1, query = pred_y, k = nn)
  if(rec_options$average == "mean"){
    func <- mean
  } else if(rec_options$average == "median"){
    func <- stats::median
  }else {
    stop("Recruiting method (option: 'average') not found")
  }
  
  idx <- order(apply(res$nn.dist, 1, func), decreasing = F)[1:num_rec]
  
  vec_from <- vec_cand[idx]
  list_to <- lapply(idx, function(i){idx1[res$nn.idx[i,]]})
  list(vec_from = vec_from, list_to = list_to)
}

#########################

# [[note to self: I'm not sure about this function name, also, Poisson hard-coded right now]]
.predict_yfromx <- function(mat_x, res_g){
  stopifnot(c("vec_g", "mat_g") %in% names(res_g))
  
  p2 <- ncol(res_g$mat_g)
  nat_param <- mat_x %*% res_g$mat_g
  
  #[[note to self: There's gotta be a cleaner way to do this]]
  for(j in 1:p2){
    nat_param[,j] <- nat_param[,j] + res_g$vec_g[j]
  }
  
  exp(nat_param)
}