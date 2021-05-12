# assumes that mat_x and mat_y are normalized
chromatin_potential_prepare <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                        dim_method = "pca", nn_method = "annoy",
                        form_method = "literal", est_method = "glmnet",
                        cand_method = "nn_any", rec_method = "distant_pearson", 
                        options = list(), verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), ncol(mat_y) == nrow(df_y),
            is.list(options))
  stopifnot(all(mat_x >= 0), all(mat_y >= 0))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # check all the options
  full_options <- .chrom_options(form_method, est_method, 
                                 cand_method, rec_method, 
                                 options)
  dim_options <- full_options$dim_options; nn_options <- full_options$nn_options

  # compute the dimension reduction
  dim_reduc_obj <- vector("list", 0)
  tmp <- dimension_reduction(mat_x, mode = "x", dim_options)
  x_dimred <- tmp$dimred
  dim_reduc_obj$x_mean <- tmp$vec_mean; dim_reduc_obj$x_sd <- tmp$vec_sd
  dim_reduc_obj$x_proj <- tmp$mat_proj
  tmp <- dimension_reduction(mat_y, mode = "y", dim_options)
  y_dimred <- tmp$dimred
  dim_reduc_obj$y_mean <- tmp$vec_mean; dim_reduc_obj$y_sd <- tmp$vec_sd
  dim_reduc_obj$y_proj <- tmp$mat_proj
  
  # form the nn
  n <- nrow(mat_x)
  all_dimred <- cbind(x_dimred, y_dimred)
  nn_obj <- nearest_neighbor(all_dimred, nn_options)
  
  # query each point's nn's
  nn_mat <- .query_nn(nn_obj, nn_options)
  
  structure(list(mat_x = mat_x, mat_y = mat_y, df_x = df_x, df_y = df_y,
                 dim_reduc_obj = dim_reduc_obj, 
                 nn_obj = nn_obj, nn_mat = nn_mat,
                 list_diagnos = list_diagnos),
            class = "chromatin_potential_prep")
}