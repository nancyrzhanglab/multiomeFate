# assumes that mat_x and mat_y are normalized
preparation <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                        nn_method = "pca",
                        form_method = "average", est_method = "glmnet",
                        cand_method = "nn_xonly_avg", rec_method = "nn_yonly", 
                        options = list(), verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), ncol(mat_y) == nrow(df_y),
            is.list(options))
  stopifnot(all(mat_x >= 0), all(mat_y >= 0))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # check all the options
  full_options <- .chrom_options(form_method, est_method, cand_method, rec_method, options)

  # compute the dimension reduction
  tmp <- dimension_reduction(mat_x, nn_method)
  x_dimred <- tmp$x_dimred; x_mean <- tmp$x_mean; x_sd <- tmp$x_sd
  tmp <- dimension_reduction(mat_y, nn_method)
  y_dimred <- tmp$y_dimred; y_mean <- tmp$y_mean; y_sd <- tmp$y_sd
  
  # form the nn
  n <- nrow(mat_x)
  latent_dim <- length(x_mean) + length(y_mean)
  nn_obj <- new(RcppAnnoy::AnnoyEuclidean, latent_dim)
  for(i in 1:n){
    
  }
  
  structure(list(mat_x = mat_x, mat_y = mat_y, df_x = df_x, df_y = df_y,
                 list_diagnos = list_diagnos),
            class = "chromatin_potential_prep")
}