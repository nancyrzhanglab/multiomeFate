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
#' \item \code{nn}: Recruit \code{rec_options$num_rec} cells
#' whose combined Modality 1 and predicted Modality 2 expression
#' vector has the smallest average 
#' (for example, mean, or in general, depending on what \code{rec_options$average}
#' is set to) distance to its \code{rec_options$nn} 
#' nearest neighbors (based on \code{nn_obj})
#' \item \code{distant_cor}: Recruit all the cells in \code{vec_cand}
#' and match cell \code{i} to a supposed cell \code{j}
#' (and its \code{rec_options$nn} nearest-neighbors)
#' whose Modality 2 difference between cell \code{j} and \code{i} 
#' has the highest correlation (via \code{rec_optiosn$cor_method})
#' with the cell \code{i}'s change in Modality 2 (computed via \code{res_g}),
#' searching among the \code{rec_otions$inflation*ncol(nn_mat)}
#' nearest-neighbors of the combined Modality 1 and predicted Modality 2 expression
#' vector but excluding the nearest neighbors of cell \code{i}
#' (based on \code{nn_mat})
#' }
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param vec_cand output of \code{.candidate_set}
#' @param res_g output of \code{.estimate_g}
#' @param df_res data frame recording the current results, generated within \code{chromatin_potential}
#' @param dim_reduc_obj object to compute the dimension reduction for a given vector,
#' computed by \code{chromatin_potential_prepare}
#' @param nn_g nearest neighbor graph, represented as an \code{igraph} object
#' @param nn_mat the nearest-neighbor matrix, output from \code{chromatin_potential_prepare}
#' @param nn_obj the exposed C++ \code{RcppAnnoy} that encodes the nearest neighbor
#' information for the \code{n} cells
#' @param enforce_matched boolean, where if \code{TRUE}, recruited cells are matched
#' to only cells that previously-recruited
#' @param df_cell (optional) data frame containing information of each cell, primarily used
#' for simulations where the true cell types are known (for \code{.recruit_next})
#' @param rec_options one of the outputs from \code{.chrom_options}
#'
#' @return a list of two things: a list called \code{rec} that contains
#' a  vector of integers \code{vec_from} and a list
#' of integers \code{list_to}, and a list called \code{diagnostic} that 
#' contains optionally-computed diagnostics to better-understand the recruitment
.recruit_next <- function(mat_x, mat_y, vec_cand, res_g, df_res, dim_reduc_obj, 
                          nn_g, nn_mat, nn_obj, enforce_matched, df_cell,
                          rec_options){
  stopifnot(all(vec_cand %% 1 == 0), all(vec_cand > 0), all(vec_cand <= nrow(mat_x)),
            length(vec_cand) == length(unique(vec_cand)))
  vec_matched <- which(!is.na(df_res$order_rec))
  stopifnot(!any(vec_cand %in% vec_matched))
  stopifnot(all(is.na(df_res$order_rec[vec_cand])), !any(is.na(df_res$order_rec[vec_matched])))
  
  if(rec_options[["method"]] == "nn"){
    res <- .recruit_next_nn(mat_x, mat_y, vec_cand, res_g, df_res, dim_reduc_obj, 
                            nn_obj, enforce_matched, rec_options)
  } else if(rec_options[["method"]] == "distant_cor"){
    res <- .recruit_next_distant_cor(mat_x, mat_y, vec_cand, res_g, df_res, 
                                     dim_reduc_obj, nn_g, nn_mat, nn_obj, 
                                     enforce_matched, rec_options)
  } else if(rec_options[["method"]] == "distant_cor_oracle"){
    res <- .recruit_next_distant_cor_oracle(mat_x, mat_y, vec_cand, res_g, df_res, 
                                     dim_reduc_obj, nn_g, nn_mat, nn_obj, enforce_matched, 
                                     df_cell, rec_options)
  } else {
    stop("Recruit method not found")
  }
  
  if(rec_options$run_diagnostic){
    res$diagnostic$postprocess <- res$rec
  }

  res
}

###################
# [[note to self: a lot of these routines should be refactored...]]

.recruit_next_nn <- function(mat_x, mat_y, vec_cand, res_g, df_res, dim_reduc_obj, nn_obj, 
                             enforce_matched, rec_options){
  num_rec <- min(rec_options$num_rec, length(vec_cand))
  
  # apply mat_g to mat_x
  len <- length(vec_cand)
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
  
  # initialize variables for the loop
  if(!rec_options$parallel && future::nbrOfWorkers() == 1){
    my_lapply <- pbapply::pblapply
    if(rec_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
  } else {
    my_lapply <- future.apply::future_lapply
  }
  
  if(enforce_matched){
    matched_idx <- which(!is.na(df_res$order_rec))
    matched_x <- .apply_dimred_mat(mat_x[matched_idx,,drop = F], dim_reduc_obj$x)
    matched_y <- .apply_dimred_mat(mat_y[matched_idx,,drop = F], dim_reduc_obj$y)
    matched_mat <- cbind(matched_x, matched_y)
    nn <- min(c(rec_options$nn, sum(!is.na(df_res$order_rec))))
  } else {
    nn <- min(c(rec_options$nn, ceiling(nrow(mat_x)/2)))
  }
  
  # see which cells are closest to the prediction 
  nn_res <- my_lapply(1:len, function(i){
    vec <- c(.apply_dimred(mat_x[vec_cand[i],], dim_reduc_obj$x),
             .apply_dimred(pred_y[i,], dim_reduc_obj$y))
    
    # allow cell to be matched to any other cell
    if(!enforce_matched){
      res <- nn_obj$getNNsByVectorList(vec, nn, search_k = -1, include_distances = T)
      res$item <- res$item+1
     
    # only match a cell to another previously-matched cell
    } else {
      tmp <- RANN::nn2(matched_mat, query = matrix(vec, nrow = 1), k = nn)
      res <- list(item = matched_idx[tmp$nn.idx[1,]], distance = tmp$nn.dist[1,])
    }
    
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
  
  rec_list <- lapply(1:length(idx), function(i){
    list(from = vec_cand[idx[i]], to = nn_res[[idx[i]]]$item)
  })
  list(rec = rec_list, diagnostic = list_diagnos)
}

.recruit_next_distant_cor <- function(mat_x, mat_y, vec_cand, res_g, df_res, 
                                      dim_reduc_obj, nn_g, nn_mat, nn_obj, 
                                      enforce_matched, rec_options){
  nn_size <- ncol(nn_mat)
  
  if(enforce_matched) matched_idx <- which(!is.na(df_res$order_rec))
  
  # apply mat_g to mat_x
  if(rec_options$bool_avg_from){
    list_nn <- .extract_nn_list(vec_cand, nn_mat)
    mat_avg_x <- .construct_avg_mat(mat_x, list_nn)
    mat_avg_y <- .construct_avg_mat(mat_y, list_nn)
    pred_y <- .predict_yfromx(mat_avg_x, res_g, rec_options$family)
  } else {
    list_nn <- lapply(vec_cand, function(x){x})
    mat_avg_x <- mat_x[vec_cand,]
    mat_avg_y <- mat_y[vec_cand,]
    pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
  }
  
  # initialize variables for the loop
  if(!rec_options$parallel && future::nbrOfWorkers() == 1){
    my_lapply <- pbapply::pblapply
    if(rec_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
  } else {
    my_lapply <- future.apply::future_lapply
  }

  list_to <- my_lapply(1:length(vec_cand), function(i){
    cell <- vec_cand[i]
    vec <- c(.apply_dimred(mat_avg_x[i,], dim_reduc_obj$x),
             .apply_dimred(pred_y[i,], dim_reduc_obj$y))
  
    # allow cell to be matched to potentially any other cell
    if(!rec_options$bool_pred_nn){
      if(!enforce_matched){
        # consider the nearest neighbors
        nn_pred <- nn_obj$getNNsByVector(vec, nn_size) + 1
      } else {
        # consider any cells, as long as it's matched
        nn_pred <- sample(matched_idx, size = ceiling(rec_options$matched_sampling_rate * length(matched_idx)))
      }
      list_nn_to <- lapply(nn_pred, function(x){x})
    
    # match to only cells near target cell
    } else {
      # [[note to self: the following could definitely be better implemented,
      #  as we recompute a lot of things unnecessarily]]
      if(!enforce_matched){
        # consider nearest neighbors that are not the immediate nn's (i.e., exclude_idx)
        list_nn_to <- .find_to_list(cell, include_idx = NA,
                                    exclude_idx = list_nn[[i]],
                                    nn_g, nn_mat, rec_options)
      } else{
        # consider nearest neighbors that are not the immediate nn's (i.e., exclude_idx)
        #  but must only include cells that have been previously matched (i.e., include_idx)
        list_nn_to <- .find_to_list_matched(mat_x, mat_y, cell, dim_reduc_obj, 
                              include_idx = setdiff(matched_idx, list_nn[[i]]), 
                              exclude_idx = list_nn[[i]], nn_mat,
                              rec_options)
      }
    }
    
    # from this set of cells, find the ones with highest pearson
    pred_diff <- pred_y[i,] - mat_avg_y[i,]
    mat_nn_to <- .construct_avg_mat(mat_y, list_nn_to)
    cor_vec <- sapply(1:nrow(mat_nn_to), function(j){
      matched_diff <- mat_nn_to[j,] - mat_avg_y[i,]
      stats::cor(pred_diff, matched_diff, method = rec_options$cor_method)
    })
    
    list_nn_to[[which.max(cor_vec)]]
  })
  
  # run the diagnostic
  list_diagnos <- list() 
  if(rec_options$run_diagnostic){
    list_diagnos[["pred_y"]] <- .predict_yfromx(mat_x, res_g, family = "gaussian")
  }
  
  rec_list <- lapply(1:length(list_nn), function(i){
    list(from = list_nn[[i]], to = list_to[[i]])
  })
  list(rec = rec_list, diagnostic = list_diagnos)
}

.recruit_next_distant_cor_oracle <- function(mat_x, mat_y, vec_cand, res_g, df_res, 
                                      dim_reduc_obj, nn_g, nn_mat, nn_obj, 
                                      enforce_matched, df_cell, rec_options){
  nn_size <- ncol(nn_mat)

  # apply mat_g to mat_x
  if(rec_options$bool_avg_from){
    list_nn <- .extract_nn_list(vec_cand, nn_mat)
    mat_avg_x <- .construct_avg_mat(mat_x, list_nn)
    mat_avg_y <- .construct_avg_mat(mat_y, list_nn)
    pred_y <- .predict_yfromx(mat_avg_x, res_g, rec_options$family)
  } else {
    list_nn <- lapply(vec_cand, function(x){x})
    mat_avg_X <- mat_x[vec_cand,]
    mat_avg_y <- mat_y[vec_cand,]
    pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
  }
  
  # initialize variables for the loop
  if(!rec_options$parallel && future::nbrOfWorkers() == 1){
    my_lapply <- pbapply::pblapply
    if(rec_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
  } else {
    my_lapply <- future.apply::future_lapply
  }
  
  list_to <- my_lapply(1:length(vec_cand), function(i){
    cell <- vec_cand[i]
    vec <- c(.apply_dimred(mat_avg_x[i,], dim_reduc_obj$x),
             .apply_dimred(pred_y[i,], dim_reduc_obj$y))
    
    # allow cell to be matched to any other cell
    if(!rec_options$bool_pred_nn){
      nn_pred <- nn_obj$getNNsByVector(vec, nn_size) + 1
      list_nn_to <- lapply(nn_pred, function(x){x})
      
      # match to only cells near target cell
    } else {
      include_idx <- intersect(which(df_cell$branch == df_cell$branch[cell]),
                              which(df_cell$time >= df_cell$time[cell]))
      list_nn_to <- .find_to_list(cell, include_idx = setdiff(include_idx, list_nn[[i]]),
                                  exclude_idx = list_nn[[i]], 
                                  nn_g, nn_mat, rec_options)
    }
    
    # from this set of cells, find the ones with highest pearson
    pred_diff <- pred_y[i,] - mat_avg_y[i,]
    mat_nn_to <- .construct_avg_mat(mat_y, list_nn_to)
    cor_vec <- sapply(1:nrow(mat_nn_to), function(j){
      matched_diff <- mat_nn_to[j,] - mat_avg_y[i,]
      stats::cor(pred_diff, matched_diff, method = rec_options$cor_method)
    })
    
    list_nn_to[[which.max(cor_vec)]]
  })
  
  # run the diagnostic
  list_diagnos <- list() 
  if(rec_options$run_diagnostic){
    # nothing currently here
  }
  
  rec_list <- lapply(1:length(list_nn), function(i){
    list(from = list_nn[[i]], to = list_to[[i]])
  })
  list(rec = rec_list, diagnostic = list_diagnos)
}

##############################3

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
    res <- nat_param
  } else if(family == "poisson"){
    res <- exp(nat_param)
  } else {
    stop("family not found")
  }
  
  if("vec_threshold" %in% names(res_g)){
    for(j in 1:p2){
      val <- res_g$vec_threshold[j]
      res[res[,j] <= val, j] <- val
    }
  }
  
  pmax(res, 0)
}

.extract_nn_list <- function(vec_cand, nn_mat){
  lapply(vec_cand, function(cell){
    c(cell, nn_mat[cell,])
  })
}

.construct_avg_mat <- function(mat, list_idx){
  t(sapply(list_idx, function(vec_idx){
    if("dgcMatrix" %in% class(mat)){
      sparseMatrixStats::colMeans2(mat[vec_idx,,drop=F])
    } else {
      matrixStats::colMeans2(mat[vec_idx,,drop=F])
    }
  }))
}

# include_idx are indices we intersect
# exclude_idx are indices we setdiff out
.find_to_list <- function(cell, include_idx, 
                          exclude_idx,
                          nn_g, nn_mat, rec_options){
  if(!all(is.na(include_idx)) & !all(is.na(exclude_idx))){
    stopifnot(length(intersect(include_idx, exclude_idx)) == 0)
  }
    
  if(rec_options$bool_avg_from){
    nn_pred <- as.numeric(igraph::ego(nn_g, order = 4, nodes = cell, mindist = 4)[[1]])
    if(length(nn_pred) == 0){ # happens if the graph doesn't have large enough radius
      # [[note to self: refactor this better...]]
      nn_pred <- as.numeric(igraph::ego(nn_g, order = 2, nodes = cell, mindist = 2)[[1]])
    }
  } else {
    nn_pred <- as.numeric(igraph::ego(nn_g, order = 2, nodes = cell, mindist = 2)[[1]])
  }
  
  if(length(nn_pred) > 0){
    nn_pred <- lapply(nn_pred, function(x){
      c(x, nn_mat[x,])
    })
  }
  
  # deal with includes
  if(length(nn_pred) > 0 & all(!is.na(include_idx))){
    nn_pred <- lapply(nn_pred, function(vec){
      intersect(vec, include_idx)
    })
  }
  
  # deal with excludes
  if(length(nn_pred) > 0 & all(!is.na(exclude_idx))){
    nn_pred <- lapply(nn_pred, function(vec){
      setdiff(vec, exclude_idx)
    })
  }
  
  # remove empty sets
  nn_pred <- nn_pred[which(sapply(nn_pred, length) > 0)]
    
  # if all else fails, match to a cell previously matched
  if(length(nn_pred) == 0){
    stopifnot(!all(is.na(include_idx)))
    nn_pred <- sample(include_idx, size = ceiling(rec_options$matched_sampling_rate * length(include_idx)),
                      replace = F)
    nn_pred <- lapply(nn_pred, function(x){
      c(x, nn_mat[x,])
    })
  }
  
  nn_pred
}

.find_to_list_matched <- function(mat_x, mat_y, cell, dim_reduc_obj, 
                                  include_idx, exclude_idx, nn_mat,
                                  rec_options){
  stopifnot(length(intersect(include_idx, exclude_idx)) == 0)
  
  matched_x <- .apply_dimred_mat(mat_x[include_idx,,drop = F], dim_reduc_obj$x)
  matched_y <- .apply_dimred_mat(mat_y[include_idx,,drop = F], dim_reduc_obj$y)
  matched_mat <- cbind(matched_x, matched_y)
  
  cell_x <- .apply_dimred_mat(mat_x[cell,,drop = F], dim_reduc_obj$x)
  cell_y <- .apply_dimred_mat(mat_y[cell,,drop = F], dim_reduc_obj$y)
  cell_vec <- c(cell_x, cell_y)
  
  tmp <- RANN::nn2(matched_mat, query = matrix(cell_vec, nrow = 1), k = rec_options$enforced_nn)$nn.idx[1,]
  nn_pred <- include_idx[tmp]
  nn_pred <- lapply(nn_pred, function(x){
    setdiff(c(x, intersect(nn_mat[x,], include_idx)), exclude_idx)
  })
  
  nn_pred
}
