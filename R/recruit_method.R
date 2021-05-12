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
                          nn_mat, nn_obj, rec_options){
  stopifnot(all(vec_cand %% 1 == 0), all(vec_cand > 0), all(vec_cand <= nrow(mat_x)),
            length(vec_cand) == length(unique(vec_cand)))
  vec_matched <- which(!is.na(df_res$order_rec))
  stopifnot(!any(vec_cand %in% vec_matched))
  stopifnot(all(is.na(df_res$order_rec[vec_cand])), !any(is.na(df_res$order_rec[vec_matched])))
  
  if(rec_options[["method"]] == "nn"){
    res <- .recruit_next_nn(mat_x, vec_cand, res_g, df_res, nn_obj, rec_options)
  } else if(rec_options[["method"]] == "distant_cor"){
    res <- .recruit_next_distant_cor(mat_x, mat_y, vec_cand, res_g, df_res, 
                                         nn_mat, nn_obj, rec_options)
  } else {
    stop("Recruit method not found")
  }
  
  if(rec_options$run_diagnostic){
    res_diag <- .recruit_diagnostic_global(mat_x, mat_y, vec_cand, res_g, 
                                          df_res, res, rec_options)
    res$diagnostic$postprocess <- res_diag
  }

  res
}

###################

.recruit_next_nn <- function(mat_x, vec_cand, res_g, df_res, nn_obj, rec_options){
  num_rec <- min(rec_options$num_rec, length(vec_cand))
  nn <- min(c(rec_options$nn, ceiling(length(vec_target)/2)))
  
  # apply mat_g to mat_x
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
  
  # see which cells are closest to the prediction 
  nn_res <- lapply(1:nrow(pred_y), function(i){
    res <- nn_obj$getNNsByItemList(pred_y[i,], nn, search_k = -1, include_distances = T)
    res$item <- res$item+1
    
    if(i %in% res$item){
      idx <- which(res$item == i)
      res$item <- res$item[-i]; res$distance <- res$distance[-i]
    }
    
    res
  })
  
  if(rec_options$average == "mean"){
    func <- mean
  } else if(rec_options$average == "median"){
    func <- stats::median
  }else {
    stop("Recruiting method (option: 'average') not found")
  }
  
  idx <- order(sapply(nn_res, function(tmp){func(tmp$distance)}), decreasing = F)[1:num_rec]
  
  # run the diagnostic
  list_diagnos <- list() 
  if(rec_options$run_diagnostic){
    # nothing currently here
  }
  
  vec_from <- vec_cand[idx] 
  list_to <- lapply(idx, function(i){nn_res[[i]]$item})
  list(rec = list(vec_from = vec_from, list_to = list_to),
       diagnostic = list_diagnos)
}

.recruit_next_distant_cor <- function(mat_x, mat_y, vec_cand, res_g, df_res, 
                                          nn_mat, nn_obj, rec_options){
  # apply mat_g to mat_x
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
  
  list_to <- lapply(1:length(vec_cand), function(i){
    cell <- vec_cand[i]
    nn_size <- ncol(nn_mat)
    nn_cand <- nn_mat[cell, ]
    nn_pred <- nn_obj$getNNsByVector(pred_y[i,], round(rec_options$inflation*nn_size)) + 1
    
    if(length(intersect(nn_pred[1:nn_size], nn_cand)) > 0) {
      nn_pred <- nn_pred[1:nn_size]
    }
    
    # find all nn's that aren't too close to cell itself
    nn_pred <- intersect(nn_pred[1:nn_size], nn_cand)
    stopifnot(length(nn_pred) > 0, !vec_cand[i] %in% nn_pred)
    
    # from this set of cells, find the ones with highest pearson
    pred_diff <- pred_y[i,] - mat_y[vec_cand[i],]
    pearson_vec <- sapply(nn_pred, function(j){
      matched_diff <- mat_y[j,] - mat_y[vec_cand[i],]
      stats::cor(pred_diff, matched_diff, method = rec_options$method)
    })
    
    idx <- nn_pred[which.max(pearson_vec)]
    tmp <- intersect(nn_mat[idx, ], nn_cand)
    if(length(tmp) >= rec_options$nn){
      vec_to <- c(idx, tmp[1:rec_options$nn])
    } else {
      vec_to <- c(idx, tmp)
    }
   
    vec_to
  })
  
  # run the diagnostic
  list_diagnos <- list() 
  if(rec_options$run_diagnostic){
    # nothing currently here
  }
  
  list(rec = list(vec_from = vec_cand, list_to = list_to),
       diagnostic = list_diagnos)
}

##############################3

.recruit_diagnostic_global <- function(mat_x, mat_y, vec_cand, res_g, df_res, 
                                       res_rec, rec_options){
  vec_matched <- which(!is.na(df_res$order_rec))
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g)
  
  nn <- min(c(rec_options$nn, nrow(mat_y)-1))
  nn_res <- RANN::nn2(mat_y, query = pred_y, k = nn+1)
  
  mat_diag <- sapply(1:length(vec_cand), function(i){
    tmp <- nn_res$nn.idx[i,]
    tmp <- tmp[tmp != vec_cand[i]] # ignore loops to itself
    stopifnot(length(tmp) > 0)
    forward_num <- length(which(tmp %in% vec_matched))
    current_num <- length(which(tmp %in% vec_cand))
    backward_num <- length(tmp) - forward_num - current_num
    
    c(forward_num = forward_num, current_num = current_num, backward_num = backward_num)
  })
  
  lis_nn <- lapply(1:length(vec_cand), function(i){
    tmp <- nn_res$nn.idx[i,]
    tmp[tmp != vec_cand[i]] # ignore loops to itself
  })
  names(lis_nn) <- as.character(vec_cand)
  
  selected <- vec_cand %in% res_rec$rec$vec_from
  df_diag <- data.frame(idx = vec_cand, forward_num = mat_diag["forward_num",],
             current_num = mat_diag["current_num",], 
             backward_num = mat_diag["backward_num",],
             selected = selected)
  
  list(df_diag = df_diag, lis_nn = lis_nn)
  
}

#########################

# [[note to self: I'm not sure about this function name]]
.predict_yfromx <- function(mat_x, res_g, family){
  stopifnot(c("vec_g", "mat_g") %in% names(res_g))
  
  p2 <- ncol(res_g$mat_g)
  nat_param <- mat_x %*% res_g$mat_g
  
  #[[note to self: There's gotta be a cleaner way to do this]]
  for(j in 1:p2){
    nat_param[,j] <- nat_param[,j] + res_g$vec_g[j]
  }
  
  if(family == "gaussian"){
    nat_param
  } else if(family == "poisson"){
    exp(nat_param)
  } else {
    stop("family not found")
  }
  
}