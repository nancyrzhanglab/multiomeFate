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

##############

.rotate <- function(mat){t(mat)[,nrow(mat):1]}