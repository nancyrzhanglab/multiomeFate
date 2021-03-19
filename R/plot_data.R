plot_activation <- function(dat, ...){
  n <- nrow(dat$obs_x)
  graphics::plot(NA, xlim = range(dat$df_x$location), ylim = c(0,n+1), xlab = "Position", ylab = "Cell")
  for(i in 1:n){
    vec <- dat$obs_x[i,]
    idx <- which(vec > 0)
    if(length(idx) > 0){
      graphics::points(x = dat$df_x$location[idx], y = rep(n-i+1, length(idx)), pch = 16, ...)
    }
  }
  
  invisible()
}

plot_heatmap <- function(dat, ...){
  image(.rotate(dat$obs_y), ...)
}

plot_umap <- function(dat, mode_x = T, mode_y = T, noiseless = T){
  
}

##############

.rotate <- function(mat){t(mat)[,nrow(mat):1]}