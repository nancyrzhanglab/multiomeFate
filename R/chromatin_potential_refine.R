chromatin_potential_refine <- function(chrom_obj, iter_max = 10, tol = 1e-4,
                                       verbose = T){
  # set up
  mat_x <- chrom_obj$mat_x; mat_y <- chrom_obj$mat_y
  df_x <- chrom_obj$df_x; df_y <- chrom_obj$df_y
  df_res <- prep_obj$df_res; dim_reduc_obj <- chrom_obj$dim_reduc_obj
  ht_neighbor <- chrom_obj$ht_neighbor
  nn_mat <- chrom_obj$nn_mat; nn_obj <- chrom_obj$nn_obj
  res_g <- res$res_g
  list_diagnos <- chrom_obj$list_diagnos; options <- chrom_obj$options
  
  dim_options <- options$dim_options; nn_options <- options$nn_options
  form_options <- options$form_options; est_options <- options$est_options
  cand_options <- options$cand_options; rec_options <- options$rec_options
  
  n <- nrow(mat_x); iter <- 1
  df_res$order_rec <- NA
  
  while(iter <= iter_max){
    if(verbose) print(paste0("Refining: On iteration ", iter))
    
    # rematch all the cells based on g
    res_rec <- .recruit_next(mat_x, mat_y, vec_cand = 1:n, res_g, df_res, 
                             dim_reduc_obj, nn_mat, nn_obj, enforce_matched = F,
                             rec_options)
    ht_neighbor <- .refine_chrom_ht(res_rec$rec)
    
    # reestimate g
    tmp <- .update_estimation_matrices(mat_x, mat_y, 
                                       mat_x1 = matrix(NA, nrow = 0, ncol = ncol(mat_x)), 
                                       mat_y2 = matrix(NA, nrow = 0, ncol = ncol(mat_y)), 
                                       res_rec$rec, form_options)
    mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
    next_res_g <- .estimate_g(mat_x1, mat_y2, est_options)
    
    # check if terminate
    if(all(abs(next_res_g$mat_g - res_g$mat_g) <= tol)) break()
    res_g <- next_res_g
    
    iter <- iter + 1
  }
  
  structure(list(res_g = res_g, mat_x = mat_x, mat_y = mat_y, 
                 df_x = df_x, df_y = df_y, df_res = df_res, 
                 dim_reduc_obj = dim_reduc_obj,
                 ht_neighbor = ht_neighbor, 
                 nn_mat = nn_mat, nn_obj = nn_obj,
                 list_diagnos = list_diagnos, options = options),
            class = "chromatin_potential")
}

#####################

.refine_chrom_ht <- function(list_res){
  n <- length(list_res$vec_from)
  stopifnot(length(unique(list_res$vec_from)) == n, max(list_res$vec_from) == n)
  
  ht_neighbor <- hash::hash()
  for(i in 1:length(list_res$vec_from)){
    ht_neighbor[[as.character(list_res$vec_from[i])]] <- list_res$list_to[[i]]
  }
  
  ht_neighbor
}