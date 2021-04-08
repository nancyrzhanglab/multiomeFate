# output: mat_g, dataframe of when each row got recruited, # of times it was a candidate, order of recruitment, and 
# hash table of who its nearest neighbors are
chromatin_potential <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                                form_method = "average", est_method = "glmnet",
                                cand_method = "nn_xonly", rec_method = "nn_yonly", 
                                options = list(),
                                verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), ncol(mat_y) == nrow(df_y),
            is.list(options))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # check all the options
  tmp <- .chrom_options(form_method, est_method, cand_method, rec_method, options)
  form_options <- tmp$form_options; est_options <- tmp$est_options
  cand_options <- tmp$cand_options; rec_options <- tmp$rec_options
  
  # initialize
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
  idx1 <- tmp$idx1
  df_res <- .init_chrom_df(n, vec_start, list_end, cell_name)
  ht_neighbor <- .init_chrom_ht(list_end)
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
    res_g <- .estimate_g(mat_x1, mat_y2, df_y, est_options)
    
    ## construct candidate set
    vec_cand <- .candidate_set(mat_x, df_res, cand_options)
    df_res <- .update_chrom_df_cand(df_res, vec_cand)
    stopifnot(all(is.na(df_res$order_rec[vec_cand])), !any(vec_cand %in% idx1))
    
    ## recruit an element from the candidate set
    rec <- .recruit_next(mat_x, vec_cand, mat_y1, idx1,
                         res_g, df_res, rec_options)
    stopifnot(all(is.na(df_res$order_rec[rec$vec_from])), !any(rec$vec_from %in% idx1))
    
    
    ## update
    tmp <- .update_estimation_matrices(mat_x, mat_y,
                                       mat_x1, mat_y1, mat_y2, idx1,
                                       rec, form_options)
    mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
    idx1 <- tmp$idx1
    ht_neighbor <- .update_chrom_ht(ht_neighbor, rec$vec_from, rec$list_to)
    df_res <- .update_chrom_df_rec(df_res, rec$vec_from, iter)
    
    iter <- iter+1
  }

  # output
  list(res_g = res_g, df_res = df_res, ht_neighbor = ht_neighbor, 
       options = list(form_options = form_options, est_options = est_options,
                      cand_options = cand_options, rec_options = rec_options))
}

#########################

# columns: steady-state (neg for initial, pos for end, or NA)
# num times was a candidate, and order of recruitment
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
  expect_true(all(tmp %in% hash::keys(ht_neighbor)))
  
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