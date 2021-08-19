#' Prepare chromatin potential
#' 
#' To use this function properly, we assume that assumes that \code{mat_x} and 
#' \code{mat_y} are normalized
#'
#' @param mat_x full data for Modality 1, where each row is a cell and each column is a variable
#' @param mat_y full data for Modality 2, where each row is a cell and each column is a variable
#' @param df_x the data frame containing information of Modality 1
#' @param df_y the data frame containing information of Modality 2
#' @param vec_start integers between 1 and \code{nrow(mat_x)} to denote the cells at the start state
#' @param list_end list of integers between 1 and \code{nrow(mat_x)} to denote the cells any of the end states
#' @param dim_method string
#' @param nn_method string
#' @param form_method string
#' @param est_method string
#' @param cand_method string
#' @param rec_method string
#' @param ht_map \code{NA} or a \code{hash} object denoting which variables
#' in \code{mat_x} correspond to each variable (the keys) in \code{mat_y}
#' @param options list
#' @param verbose boolean
#'
#' @return object of class \code{chromatin_potential_prep}
#' @export
chromatin_potential_prepare <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                        dim_method = "pca", nn_method = "annoy",
                        form_method = "average", est_method = "glmnet",
                        cand_method = "nn_any", rec_method = "distant_cor", 
                        ht_map = NA, options = list(), verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), 
            ncol(mat_y) == nrow(df_y), is.list(options),
            length(unlist(list_end)) > 0)
  stopifnot(all(mat_x >= 0), all(mat_y >= 0))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # check all the options
  full_options <- .chrom_options(dim_method, nn_method,
                                 form_method, est_method, 
                                 cand_method, rec_method, 
                                 options)
  dim_options <- full_options$dim_options; nn_options <- full_options$nn_options

  # compute the dimension reduction
  dim_reduc_obj <- vector("list", 0)
  tmp <- dimension_reduction(mat_x, mode = "x", dim_options)
  x_scores <- tmp$scores; dim_reduc_obj$x <- tmp$dim_reduc_obj
  tmp <- dimension_reduction(mat_y, mode = "y", dim_options)
  y_scores <- tmp$scores; dim_reduc_obj$y <- tmp$dim_reduc_obj
  
  # form the nn
  n <- nrow(mat_x)
  all_scores <- cbind(x_scores, y_scores)
  nn_obj <- nearest_neighbor(all_scores, nn_options)
  
  # query each point's nn's
  nn_mat <- .query_nn(nn_obj, nn_options)
  nn_g <- .nn_to_igraph(nn_mat, name_vec = rownames(mat_x))
  
  # initialize
  list_diagnos <- list()
  df_res <- .init_chrom_df(n, vec_start, list_end, cell_name)
  if(full_options$est_options$enforce_cis){
    if(class(ht_map) != "hash" && all(is.na(ht_map))){
      full_options$est_options <- .gene_peak_map(df_x, df_y, full_options$est_options)
    } else {
      full_options$est_options$ht_map <- ht_map
    }
  }
  
  structure(list(mat_x = mat_x, mat_y = mat_y, df_x = df_x, df_y = df_y,
                 df_res = df_res, dim_reduc_obj = dim_reduc_obj, 
                 nn_g = nn_g, nn_mat = nn_mat, nn_obj = nn_obj, 
                 list_diagnos = list_diagnos, options = full_options),
            class = "chromatin_potential_prep")
}

#################################

.init_chrom_df <- function(n, vec_start, list_end, cell_name){
  if(!all(is.na(vec_start))) stopifnot(all(vec_start %% 1 == 0), all(vec_start > 0), all(vec_start <= n))
  stopifnot(all(sapply(list_end, function(vec){all(vec %% 1 == 0) & all(vec > 0) & all(vec <= n)})))
  if(!all(is.na(vec_start))) {
    tmp <- c(vec_start, unlist(list_end))
  } else {
    tmp <- unlist(list_end)
  }
  stopifnot(length(tmp) == length(unique(tmp)))
  
  df_res <- data.frame(idx = 1:n, init_state = rep(NA, n), num_cand = rep(0, n),
                       order_rec = rep(NA, n))
  if(length(cell_name) == n) rownames(df_res) <- cell_name
  
  if(!all(is.na(vec_start))) df_res$init_state[vec_start] <- -1
  for(i in 1:length(list_end)){
    df_res$init_state[list_end[[i]]] <- i
    df_res$order_rec[list_end[[i]]] <- 0
  }
  
  df_res
}

.nn_to_igraph <- function(nn_mat, name_vec){
  # generate the edge-mat
  edge_mat <- do.call(rbind, lapply(1:nrow(nn_mat), function(i){
    matrix(c(rep(i, ncol(nn_mat)), nn_mat[i,]), ncol = 2, nrow = ncol(nn_mat))
  }))
  
  nn_g <- igraph::graph_from_edgelist(edge_mat, directed = F)
  nn_g <- igraph::simplify(nn_g)
  
  if(!all(is.na(name_vec))){
    stopifnot(length(name_vec) == nrow(nn_mat))
    nn_g <- igraph::set_vertex_attr(nn_g, "name", value = name_vec)
  }
  
  nn_g
}