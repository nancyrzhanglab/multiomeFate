#' Initialize the matrices for estimation later
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param df_res data frame formed during the \code{chromatin_potential_prepare} function
#' 
#' @return list of 2 matrices, \code{mat_x1} and \code{mat_y2}, as
#' well as vector of indices \code{vec_matched}
.init_est_matrices <- function(mat_x, mat_y, df_res){
  stopifnot(nrow(mat_x) == nrow(mat_y))
  
  # initialize
  vec_start <- which(df_res$init_state == -1)
  vec_onlyend <- which(df_res$init_state > 0)
  vec <- c(vec_start, vec_onlyend)
  stopifnot(length(vec) == length(unique(vec)))

  # fill in steady-states
  mat_x1 <- mat_x[vec,,drop = F]; mat_y2 <- mat_y[vec,,drop = F]
  
  list(mat_x1 = mat_x1, mat_y2 = mat_y2)
}

#' Update the matrices for estimation later
#' 
#' Using the output of \code{.recruit_next} (in \code{rec}), update
#' \code{mat_x1}, \code{mat_y2} by grabbing
#' the appropriate cells in \code{mat_x} and \code{mat_y} based on the
#' options in \code{form_options}
#' 
#' The options are:
#' \itemize{
#' \item \code{literal}: If \code{rec$list_to[[i]]} matches many cells from a single cell in \code{rec$vec_from[i]},
#' form one row in \code{mat_y2} for each matched cell in \code{rec$list_to[[i]]} (and 
#' duplicate the cell in \code{rec$vec_from[i]} in \code{mat_x1} as many times as needed)
#' \item \code{average}: If \code{rec$list_to[[i]]} matches many cells from a single cell in \code{rec$vec_from[i]},
#' average all the variables among the matched cell in \code{rec$list_to[[i]]} based on
#' \code{form_options$average} to form one row in \code{mat_y2} (and one row in 
#' \code{mat_x1} for \code{rec$vec_from[i]})
#' }
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param mat_x1 previous \code{mat_x1}, say, from an earlier call to \code{.init_est_matrices} or \code{.update_estimation_matrices}
#' @param mat_y2 previous \code{mat_y2}, say, from an earlier call to \code{.init_est_matrices} or \code{.update_estimation_matrices}
#' @param rec output \code{rec} from \code{.recruit_next}
#' @param form_options one of the outputs from \code{.chrom_options}
#'
#' @return list of 2 matrices, \code{mat_x1} and \code{mat_y2}
.update_estimation_matrices <- function(mat_x, mat_y,
                                        mat_x1, mat_y2,
                                        rec, form_options){
  p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
  stopifnot(ncol(mat_x1) == p1, ncol(mat_y2) == p2,
            nrow(mat_y) == n, nrow(mat_x1) == nrow(mat_y2))
  
  if(form_options[["method"]] == "average"){
    res <- .update_estimation_average(mat_x, mat_y, mat_x1, mat_y2, 
                                      rec, form_options)
  } else {
    stop("Forming method not found")
  }
  
  res
}

########################

.update_estimation_average <- function(mat_x, mat_y, mat_x1, mat_y2,
                                       rec, form_options){
  p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)

  # for mat_x1 and mat_y2
  tmp_x <- t(sapply(1:length(rec), function(i){
    if(form_options$average == "mean"){
      if("dgcMatrix" %in% class(mat_x)){
        sparseMatrixStats::colMeans2(mat_x[rec[[i]]$from,,drop=F])
      } else {
        matrixStats::colMeans2(mat_x[rec[[i]]$from,,drop=F])
      }
    } else if(form_options$average == "median") {
      if("dgcMatrix" %in% class(mat_x)){
        sparseMatrixStats::colMedians(mat_x[rec[[i]]$from,,drop=F])
      } else {
        matrixStats::colMedians(mat_x[rec[[i]]$from,,drop=F])
      }
    } else {
      stop("Forming method (option: 'average') not found")
    }
  }))
  
  tmp_y <- t(sapply(1:length(rec), function(i){
    if(form_options$average == "mean"){
      if("dgcMatrix" %in% class(mat_y)){
        sparseMatrixStats::colMeans2(mat_y[rec[[i]]$to,,drop=F])
      } else {
        matrixStats::colMeans2(mat_y[rec[[i]]$to,,drop=F])
      }
    } else if(form_options$average == "median") {
      if("dgcMatrix" %in% class(mat_y)){
        sparseMatrixStats::colMedians(mat_y[rec[[i]]$to,,drop=F])
      } else {
        matrixStats::colMedians(mat_y[rec[[i]]$to,,drop=F])
      }
    } else {
      stop("Forming method (option: 'average') not found")
    }
  }))
  
  mat_x1 <- rbind(mat_x1, tmp_x)
  mat_y2 <- rbind(mat_y2, tmp_y)

  list(mat_x1 = mat_x1, mat_y2 = mat_y2)
}