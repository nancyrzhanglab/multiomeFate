#' Based on the candidate cells, recruit a specific cell
#' 
#' Here, \code{mat_x} are all the cells in our dataset, among which
#' we are saying that \code{vec_matched} are cells that have previously been
#' recruited already, and \code{vec_cand} are the currently-selected
#' candidates that we are hoping to recruit and match to one of the
#' previously-recruited cells acccording to \code{df_res}. This 
#' matching is done by using our estimation function \code{res_g}
#' on the candidate cell's Modality 1 expression in \code{mat_x}
#' and seeing if it's close to any prreviously-recruited cells's Modality 2 expressions
#' in \code{mat_y}.
#' 
#' The options are:
#' \itemize{
#' \item \code{nn_yonly}: Recruit \code{rec_options$num_rec} cells
#' whose predicted Modality 2 expression has the smallest average
#' (for example, mean, or in general, depending on what \code{rec_options$average}
#' is set to) distance to its \code{rec_options$nn} previously-recruited
#' nearest neighbors (dictated by the information in \code{df_res})
#' }
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param vec_cand output of \code{.candidate_set}
#' @param res_g output of \code{.estimate_g}
#' @param df_res data frame recording the current results, generated within \code{chromatin_potential}
#' @param rec_options one of the outputs from \code{.chrom_options}
#'
#' @return a list of two things: a list called \code{rec} that contains
#' a  vector of integers \code{vec_from} and a list
#' of integers \code{list_to}, and a list called \code{diagnostic} that 
#' contains optionally-computed diagnostics to better-understand the recruitment
.recruit_next <- function(mat_x, mat_y, vec_cand, res_g, df_res, 
                          rec_options, dir_back){
  stopifnot(all(vec_cand %% 1 == 0), all(vec_cand > 0), all(vec_cand <= nrow(mat_x)),
            length(vec_cand) == length(unique(vec_cand)))
  vec_matched <- which(!is.na(df_res$order_rec))
  stopifnot(!any(vec_cand %in% vec_matched))
  stopifnot(all(is.na(df_res$order_rec[vec_cand])), !any(is.na(df_res$order_rec[vec_matched])))
  
  if(rec_options[["method"]] == "nn_yonly"){
    res <- multiomeFate:::.recruit_next_nn_yonly(mat_x, mat_y, vec_cand, res_g, df_res, rec_options, dir_back)
  } else {
    stop("Recruit method not found")
  }

  res
}

###################

.recruit_next_nn_yonly <- function(mat_x, mat_y, vec_cand, res_g, df_res,
                                    rec_options, dir_back){
  if(dir_back==TRUE){
  vec_matched <- which(!is.na(df_res$order_rec))###
  }else if(dir_back==FALSE){
  vec_matched <- .find_vector_matched(rec_options = rec_options, vec_cand = vec_cand, df_res = df_res, mat_x=mat_x)
  }else{
    stop("Please specify \"dir_back\" as TRUE/FALSE.")
  }
  num_rec <- min(rec_options$num_rec, length(vec_cand))
  nn <- min(c(rec_options$nn, ceiling(length(vec_matched)/2)))
  
  

  # apply mat_g to mat_x
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g)
  
  # see which prediction is closest to mat_y[vec_matched,]
  # [[note to self: check if it's worthwhile to expose the ANN KD-tree here]]
  res <- RANN::nn2(mat_y[vec_matched,], query = pred_y, k = nn)
  if(rec_options$average == "mean"){
    func <- mean
  } else if(rec_options$average == "median"){
    func <- stats::median
  }else {
    stop("Recruiting method (option: 'average') not found")
  }
  
  idx <- order(apply(res$nn.dist, 1, func), decreasing = F)[1:num_rec]
  
  # run the diagnostic
  list_diagnos <- list()
  if(rec_options$run_diagnostic){
    # [[note to self: put diagnostics here]]
  }
  
  vec_from <- vec_cand[idx]
  list_to <- lapply(idx, function(i){vec_matched[res$nn.idx[i,]]})
  list(rec = list(vec_from = vec_from, list_to = list_to),
       diagnostic = list_diagnos)
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
  
  #exp(nat_param)
  nat_param
}

#########################

.find_vector_matched=function(rec_options=NULL, vec_cand=NULL, df_res=NULL, mat_x=NULL){
  n <- nrow(df_res)
  idx_free <- which(! ((1:n) %in% vec_cand))
  idx_rec <- vec_cand
  if(length(idx_free) == 0) return(numeric(0))
  if(length(idx_free) <= rec_options$num_rec) return(idx_free)
  nn <- length(vec_cand)

  res <- RANN::nn2(mat_x[idx_rec,,drop = F], query = mat_x[idx_free,,drop = F],
                   k = nn)

  if(rec_options$average == "mean"){
    func <- mean
  } else if(rec_options$average == "median"){
    func <- stats::median
  }else {
    stop("Candidate method (option: 'average') not found")
  }

  idx <- order(apply(res$nn.dist, 1, func), decreasing = F)[1:rec_options$search_num]

  idx_free[idx]
}

