#' Plot UMAP embedding (for \code{mf_simul})
#'
#' @param obj object of class \code{mf_simul}, from the output of \code{generate_data}
#' @param mode_x boolean, where \code{TRUE} means the data from Modality 1 is involved in the UMAP
#' @param mode_y boolean, where \code{TRUE} means the data from Modality 2 is involved in the UMAP
#' @param noiseless boolean, where \code{TRUE} means information from 
#' \code{dat$true_x} and/or \code{dat$true_y} is plotted, and \code{FALSE}
#' means information from \code{dat$obs_x} and/or \code{dat$obs_y} is plotted
#' @param k positive numeric for the number of principal scores extracted
#' @param reorder boolean. If \code{TRUE}, color the cells by their pseudotime.
#' If \code{FALSE}, ccolor the cells by the order they appear in simulation
#' @param num_col positive numeric, for the number of distinct colors used
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_umap.mf_simul <- function(obj, mode_x = T, mode_y = T, noiseless = F, k = 10,
                               reorder = T, num_col = 10, ...){
  stopifnot(class(obj) == "mf_simul", mode_x | mode_y)
  
  if(noiseless){
    embedding <- .umap_embedding(obj$true_x, obj$true_y, mode_x, mode_y, k)
  } else {
    embedding <- .umap_embedding(obj$obs_x, obj$obs_y, mode_x, mode_y, k)
  }
  
  # extract color
  if(reorder){
    vec <- obj$df_info$time
  } else {
    vec <- obj$df_info$counter
  }
  col_palette <- grDevices::colorRampPalette(c("red", "blue"))(num_col)
  vec_val <- seq(min(vec), max(vec), length.out = 10)
  col_vec <- sapply(vec, function(x){
    col_palette[which.min(abs(x - vec_val))]
  })
  
  graphics::plot(embedding[,1], embedding[,2], pch = 16, col = col_vec, ...)
  invisible()
}

#' Plot UMAP embedding (for \code{chromatin_potential})
#'
#' @param obj object of class \code{chromatin_potential}, from the output of \code{chromatin_potential}
#' @param mode_x boolean, where \code{TRUE} means the data from Modality 1 is involved in the UMAP
#' @param mode_y boolean, where \code{TRUE} means the data from Modality 2 is involved in the UMAP
#' @param k positive numeric for the number of principal scores extracted
#' @param percent_arrows numeric
#' @param arrow_length numeric
#' @param col_vec \code{NA} or vector of colors
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_umap.chromatin_potential <- function(obj, mode_x = T, mode_y = T, k = 10, 
                                          percent_arrows = 0.3,
                                          arrow_length = 0.05,
                                          col_vec = NA, ...){
  stopifnot(class(obj) == "chromatin_potential", mode_x | mode_y, percent_arrows >= 0,
            percent_arrows <= 1)
  if(all(is.na(col_vec))){
    col_vec <- grDevices::rgb(0.5, 0.5, 0.5, 0.5)
  } else {
    stopifnot(length(col_vec) == nrow(obj$mat_x))
  }
  
  embedding <- .umap_embedding(obj$mat_x, obj$mat_y, mode_x, mode_y, k)
  graphics::plot(embedding[,1], embedding[,2], pch = 16, col = col_vec, ...)
  
  # add arrows
  if(percent_arrows > 0){
    n <- nrow(obj$mat_x)
    idx <- sample(1:n, size = round(n * percent_arrows))
    
    for(i in idx){
      vec_from <- embedding[i,]
      neigh <- obj$ht_neighbor[[as.character(i)]]
      if(length(neigh) == 1){
        vec_to <- embedding[neigh[1],]
      } else {
        vec_to <- colMeans(embedding[neigh,])
      }
      
      suppressWarnings(graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2],
                       length = arrow_length))
    }
  }
  
  invisible()
}

############################

.umap_embedding <- function(mat_x, mat_y, mode_x, mode_y, k){
  mat <- numeric(0)
  if(mode_x){
    if(k > ncol(mat_x)){
      tmp <- scale(mat_x)
    } else {
      tmp <- stats::prcomp(mat_x)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  if(mode_y){
    if(k > ncol(mat_y)){
      tmp <- scale(mat_y)
    } else {
      tmp <- stats::prcomp(mat_y)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  # extract embedding
  embedding <- suppressWarnings(Seurat::RunUMAP(mat, verbose = F)@cell.embeddings)
}