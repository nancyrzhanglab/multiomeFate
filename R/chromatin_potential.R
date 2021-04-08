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
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param df_x the data frame containing information of Modality 1
#' @param df_y the data frame containing information of Modality 2
#' @param vec_start integers between 1 and \code{nrow(mat_x)} to denote the cells at the start state
#' @param list_end list of integers between 1 and \code{nrow(mat_x)} to denote the cells any of the end states
#' @param form_method string
#' @param est_method string
#' @param cand_method string
#' @param rec_method string
#' @param options list
#' @param verbose boolean
#'
#' @return object of class \code{chromatin_potential}
#' @export
chromatin_potential <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                                form_method = "average", est_method = "glmnet",
                                cand_method = "nn_xonly_avg", rec_method = "nn_yonly", 
                                options = list(),
                                verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), ncol(mat_y) == nrow(df_y),
            is.list(options))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # check all the options
  full_options <- .chrom_options(form_method, est_method, cand_method, rec_method, options)
  form_options <- full_options$form_options; est_options <- full_options$est_options
  cand_options <- full_options$cand_options; rec_options <- full_options$rec_options
  
  # initialize
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
  df_res <- .init_chrom_df(n, vec_start, list_end, cell_name)
  ht_neighbor <- .init_chrom_ht(list_end)
  list_diagnos <- list()
  iter <- 1
  if(est_options$enforce_cis){
    est_options <- .gene_peak_map(df_x, df_y, est_options)
  }
  
  # while:
  while(length(ht_neighbor) < n){
    # [[note to self: put a better statement here]]
    if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                             round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), ")"))
    ## estimate res_g
    res_g <- .estimate_g(mat_x1, mat_y2, est_options)
    
    ## construct candidate set
    vec_cand <- .candidate_set(mat_x, mat_y, df_res, cand_options)
    df_res <- .update_chrom_df_cand(df_res, vec_cand)
    stopifnot(all(is.na(df_res$order_rec[vec_cand])))
    
    ## recruit an element from the candidate set
    res <- .recruit_next(mat_x, mat_y, vec_cand, 
                         res_g, df_res, rec_options)
    stopifnot(all(is.na(df_res$order_rec[res$rec$vec_from])))
    
    
    ## update
    tmp <- .update_estimation_matrices(mat_x, mat_y, mat_x1, mat_y2, 
                                       res$rec, form_options)
    mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
    ht_neighbor <- .update_chrom_ht(ht_neighbor, res$rec$vec_from, res$rec$list_to)
    df_res <- .update_chrom_df_rec(df_res, res$rec$vec_from, iter)
    list_diagnos[[as.character(iter)]] <- res$diagnostic
    
    iter <- iter+1
  }

  # output
  structure(list(res_g = res_g, df_res = df_res, ht_neighbor = ht_neighbor, 
                 mat_x = mat_x, mat_y = mat_y, df_x = df_x, df_y = df_y,
                 list_diagnos = list_diagnos,
                 options = full_options),
       class = "chromatin_potential")
}

#########################

.init_chrom_df <- function(n, vec_start, list_end, cell_name){
  stopifnot(all(vec_start %% 1 == 0), all(vec_start > 0), all(vec_start <= n))
  stopifnot(all(sapply(list_end, function(vec){all(vec %% 1 == 0) & all(vec > 0) & all(vec <= n)})))
  tmp <- c(vec_start, unlist(list_end))
  stopifnot(length(tmp) == length(unique(tmp)))
  
  df_res <- data.frame(idx = 1:n, init_state = rep(NA, n), num_cand = rep(0, n),
                       order_rec = rep(NA, n))
  if(length(cell_name) == n) rownames(df_res) <- cell_name
  
  df_res$init_state[vec_start] <- -1
  for(i in 1:length(list_end)){
    df_res$init_state[list_end[[i]]] <- i
    df_res$order_rec[list_end[[i]]] <- 0
  }
  
  df_res
}

.init_chrom_ht <- function(list_end){
  ht_neighbor <- hash::hash()
  vec <- unlist(list_end)
  for(i in vec){
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

.update_chrom_ht <- function(ht_neighbor, vec_from, list_to){
  tmp <- as.character(unlist(list_to))
  stopifnot(all(tmp %in% hash::keys(ht_neighbor)))
  
  for(i in 1:length(vec_from)){
    ht_neighbor[[as.character(vec_from[i])]] <- list_to[[i]]
  }
  
  ht_neighbor
}

# [[note to self: currently does not check that iter is the latest/largest iteration in df_res]]
.update_chrom_df_rec <- function(df_res, vec_from, iter){
  stopifnot(all(is.na(df_res$order_rec[vec_from])))
  
  df_res$order_rec[vec_from] <- iter
  df_res
}