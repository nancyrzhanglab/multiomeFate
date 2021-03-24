plot_activation <- function(dat, reorder = F, ...){
  n <- nrow(dat$obs_x)
  mat <- dat$obs_x
  if(reorder){
    mat <- mat[order(dat$df_info$time, decreasing = F), ]
  }
  graphics::plot(NA, xlim = range(dat$df_x$location), ylim = c(0,n+1), xlab = "Position", ylab = "Cell")
  for(i in 1:n){
    vec <- mat[i,]
    idx <- which(vec > 0)
    if(length(idx) > 0){
      graphics::points(x = dat$df_x$location[idx], y = rep(n-i+1, length(idx)), pch = 16, ...)
    }
  }
  
  invisible()
}

plot_heatmap <- function(dat, reorder = F, ...){
  mat <- dat$obs_y
  if(reorder){
    mat <- mat[order(dat$df_info$time, decreasing = F), ]
  }
  
  graphics::image(.rotate(mat), ...)
  
  invisible()
}

plot_umap <- function(dat, mode_x = T, mode_y = T, noiseless = T, k = 10,
                      color_by = "time", num_col = 10, ...){
  # [[currently noiseless does nothing]]
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
  if(color_by == "time"){
    vec <- dat$df_info$time
  } else {
    vec <- dat$df_info$counter
  }
  col_palette <- grDevices::colorRampPalette(c("red", "blue"))(num_col)
  vec_val <- seq(min(vec), max(vec), length.out = 10)
  col_vec <- sapply(vec, function(x){
    col_palette[which.min(abs(x - vec_val))]
  })

  plot(embedding[,1], embedding[,2], pch = 16, col = col_vec, ...)
  invisible()
}

##############

.rotate <- function(mat){t(mat)[,nrow(mat):1]}