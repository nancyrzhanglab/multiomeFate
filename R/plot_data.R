#' Plot activation 
#'
#' @param dat object of class \code{mf_simul}, from the output of \code{generate_data}
#' @param mode \code{1} or \code{2}, denoting whether or not the plot is of Modality 1 or 2 respecitvely
#' @param cutoff numeric, for the value that needs to be exceeded for a point to be plotted
#' @param reorder boolean. If \code{TRUE}, order the cells by their pseudotime.
#' If \code{FALSE}, order the cells by the order they appear in simulation
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_activation <- function(dat, mode = 1, cutoff = 0.5, reorder = T, ...){
  stopifnot(class(dat) == "mf_simul", mode %in% c(1,2), length(mode) == 1)
  
  if(mode == 1){
    mat <- dat$obs_x; loc <- dat$df_x$location
  } else {
    mat <- dat$obs_y; loc <- dat$df_y$location
  }
  n <- nrow(mat)
 
  if(reorder){
    mat <- mat[order(dat$df_info$time, decreasing = F), ]
  } else {
    mat <- mat[order(dat$df_info$counter, decreasing = F), ]
  }
  
  graphics::plot(NA, xlim = range(loc), ylim = c(0,n+1), xlab = "Position", ylab = "Cell")
  for(i in 1:n){
    vec <- mat[i,]
    idx <- which(vec >= cutoff)
    if(length(idx) > 0){
      graphics::points(x = loc[idx], y = rep(n-i+1, length(idx)), pch = 16, ...)
    }
  }
  
  invisible()
}

#' Plot heatmap
#'
#' @param dat object of class \code{mf_simul}, from the output of \code{generate_data}
#' @param mode \code{1} or \code{2}, denoting whether or not the plot is of Modality 1 or 2 respecitvely
#' @param reorder boolean. If \code{TRUE}, order the cells by their pseudotime.
#' If \code{FALSE}, order the cells by the order they appear in simulation
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_heatmap <- function(dat, mode = 1, reorder = T, ...){
  stopifnot(class(dat) == "mf_simul")
  
  if(mode == 1){
    mat <- dat$obs_x
  } else {
    mat <- dat$obs_y
  }
  
  if(reorder){
    mat <- mat[order(dat$df_info$time, decreasing = F), ]
  } else {
    mat <- mat[order(dat$df_info$counter, decreasing = F), ]
  }
  
  graphics::image(.rotate(mat), ...)
  
  invisible()
}

#' Plot UMAP embedding
#'
#' @param dat object of class \code{mf_simul}, from the output of \code{generate_data}
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
plot_umap <- function(dat, mode_x = T, mode_y = T, noiseless = T, k = 10,
                      reorder = T, num_col = 10, ...){
  stopifnot(class(dat) == "mf_simul", mode_x | mode_y)
  
  # [note to self: currently noiseless does nothing]
  mat <- numeric(0)
  if(mode_x){
    if(k > ncol(dat$obs_x)){
      tmp <- scale(dat$obs_x)
    } else {
      tmp <- stats::prcomp(dat$obs_x)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  if(mode_y){
    if(k > ncol(dat$obs_y)){
      tmp <- scale(dat$obs_y)
    } else {
      tmp <- stats::prcomp(dat$obs_y)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  # extract embedding
  embedding <- Seurat::RunUMAP(mat, verbose = F)@cell.embeddings
  
  # extract color
  if(reorder){
    vec <- dat$df_info$time
  } else {
    vec <- dat$df_info$counter
  }
  col_palette <- grDevices::colorRampPalette(c("red", "blue"))(num_col)
  vec_val <- seq(min(vec), max(vec), length.out = 10)
  col_vec <- sapply(vec, function(x){
    col_palette[which.min(abs(x - vec_val))]
  })

  graphics::plot(embedding[,1], embedding[,2], pch = 16, col = col_vec, ...)
  invisible()
}

##############

.rotate <- function(mat){t(mat)[,nrow(mat):1]}