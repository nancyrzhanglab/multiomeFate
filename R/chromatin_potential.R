#' Main function to estimate chromatin potential
#' 
#' \code{.init_est_matrices} forms 2 different matrices
#' that this function updates every iteration (via \code{.update_estimation_matrices}):
#' \itemize{
#' \item \code{mat_x1}: This is the matrix about Modality 1
#' (i.e., has \code{ncol(mat_x)} variables)
#' where each row is a cell that has been previously-recruited, grabbed
#' from \code{mat_x}
#' \item \code{mat_y2}: This is the matrix about Modality 2
#' (i.e., has \code{ncol(mat_y)} variables)
#' where each row is the "future" cell matched
#' to the corresponding cell in \code{mat_x1}, grabbed
#' from \code{mat_y}. That is, this matrix has the same number of 
#' rows as \code{mat_x1}, but the \code{i}th row in
#' \code{mat_x1} might represent a different cell than the
#' \code{i}th row in \code{mat_y2}
#' }
#' In short, \code{mat_x1} and \code{mat_y2} always have the same number
#' of rows (but might represent different cells), and these two matrices
#' are used for \code{.estimate_g} to estimate the link from Modality 1 
#' to Modality 2. 
#'
#' @param prep_obj object of class \code{chromatin_potential_prep}, created using the
#' \code{chromatin_potential_prepare} function
#' @param mat_g_init (optional) initial estimate of the matrix of coefficients linking Modality 1 to Modality 2
#' @param vec_g_init (optional) initial estimate of the vector of intercepts linking Modality 1 to Modality 2
#' @param verbose boolean
#'
#' @return object of class \code{chromatin_potential}
#' @export
chromatin_potential <- function(prep_obj, mat_g_init = NA, vec_g_init = rep(0, ncol(mat_y)),
                                vec_threshold_init = rep(0, ncol(mat_y)),
                                df_cell = NA, bool_oracle = F, verbose = T){
  # pull the appropriate objects for convenience
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
  df_res <- prep_obj$df_res; dim_reduc_obj <- prep_obj$dim_reduc_obj
  nn_g <- prep_obj$nn_g; nn_mat <- prep_obj$nn_mat; nn_obj <- prep_obj$nn_obj
  list_diagnos <- prep_obj$list_diagnos; options <- prep_obj$options
  
  dim_options <- options$dim_options; nn_options <- options$nn_options
  form_options <- options$form_options; est_options <- options$est_options
  cand_options <- options$cand_options; rec_options <- options$rec_options
  
  # initialize
  n <- nrow(mat_x)
  
  ht_neighbor <- .init_chrom_ht(which(df_res$order_rec == 0))
  tmp <- .init_est_matrices(mat_x, mat_y, df_res, form_options)
  mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
  list_diagnos <- list()
  iter <- 1
  weights <- .init_weights(nrow(mat_x1), form_options)
  
  # while:
  while(length(ht_neighbor) < n){
    # [[note to self: put a better statement here]]
    if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                             round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), ")"))
    if(verbose & !any(is.na(weights))) print(paste0("Weights range from ", round(min(weights),2), " to ", round(max(weights),2)))
    ## estimate res_g
    if((iter == 1 | est_options$hold_initial) && !any(is.na(mat_g_init)) && !any(is.na(vec_g_init))){
      res_g <- list(mat_g = mat_g_init, vec_g = vec_g_init, vec_threshold = vec_threshold_init)
    } else {
      res_g <- .estimate_g(mat_x1, mat_y2, weights, est_options)
    }
   
    ## construct candidate set
    res_cand <- .candidate_set(mat_x, mat_y, df_res, nn_mat, cand_options)
    df_res <- .update_chrom_df_cand(df_res, res_cand$vec_cand)
    stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
    list_diagnos[[as.character(iter)]]$candidate <- res_cand$diagnostic
    
    ## recruit an element from the candidate se
    enforce_matched <- length(which(df_res$order_rec == 0)) > length(which(df_res$order_rec > 0)) & !bool_oracle
    res_rec <- .recruit_next(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                             dim_reduc_obj, nn_g, nn_mat, nn_obj, enforce_matched,
                             df_cell, rec_options)
    list_diagnos[[as.character(iter)]]$recruit <- res_rec$diagnostic
    
    ## update
    tmp <- .update_estimation_matrices(mat_x, mat_y, mat_x1, mat_y2, weights,
                                       res_rec$rec, form_options)
    mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2; weights <- tmp$weights
    ht_neighbor <- .update_chrom_ht(ht_neighbor, res_rec$rec, enforce_matched)
    df_res <- .update_chrom_df_rec(df_res, res_rec$rec, iter)
    
    iter <- iter+1
  }

  # output
  structure(list(res_g = res_g, mat_x = mat_x, mat_y = mat_y, 
                 df_x = df_x, df_y = df_y, df_res = df_res, 
                 df_cell = df_cell,
                 dim_reduc_obj = dim_reduc_obj,
                 ht_neighbor = ht_neighbor, 
                 nn_mat = nn_mat, nn_obj = nn_obj,
                 list_diagnos = list_diagnos, options = options),
       class = "chromatin_potential")
}

#########################


.init_chrom_ht <- function(vec_end){
  ht_neighbor <- hash::hash()
  for(i in vec_end){
    ht_neighbor[[as.character(i)]] <- c(neighbor = i)
  }
  
  ht_neighbor
}


.update_chrom_df_cand <- function(df_res, vec_cand){
  stopifnot(all(is.na(df_res$order_rec[vec_cand])))
  tmp <- df_res$init_state[vec_cand]; tmp <- tmp[!is.na(tmp)]
  if(length(tmp) > 0) stopifnot(all(tmp < 0))
  stopifnot(all(vec_cand <= nrow(df_res)), all(vec_cand %% 1 == 0), all(vec_cand > 0),
            length(vec_cand) == length(unique(vec_cand)))
  
  df_res$num_cand[vec_cand] <-  df_res$num_cand[vec_cand]+1
  df_res
}

.update_chrom_ht <- function(ht_neighbor, rec, enforce_matched){
  if(enforce_matched) {
    for(i in 1:length(rec)){
      tmp <- as.character(rec[[i]]$to)
      stopifnot(all(tmp %in% hash::keys(ht_neighbor)))
    }
  }
  
  for(i in 1:length(rec)){
    for(idx in rec[[i]]$from){
      tmp <- ht_neighbor[[as.character(idx)]]
      ht_neighbor[[as.character(idx)]] <- c(tmp, rec[[i]]$to)
    }
  }
  
  ht_neighbor
}

# [[note to self: currently does not check that iter is the latest/largest iteration in df_res]]
.update_chrom_df_rec <- function(df_res, rec, iter){
  for(i in 1:length(rec)){
    df_res$order_rec[rec[[i]]$from] <- iter
  }
  
  df_res
}