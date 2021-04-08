#' Plot pseudotime direction over recruitment index
#'
#' @param res object of class \code{chromatin_potential}, from the output of \code{chromatin_potential}
#' @param vec_time numeric vector of pseudotimes, with length \code{nrow(res$mat_x)}
#' @param jitter boolean
#' @param arrow_length numeric
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_arrow_iteration <- function(res, vec_time, jitter = T,
                                 arrow_length = 0.1, ...){
  n <- nrow(res$mat_x)
  stopifnot(length(vec_time) == n)
  
  if(jitter){vec_time <- vec_time + stats::rnorm(n, sd = diff(range(vec_time))/100)}
  vec_from <- rep(NA, n); vec_to <- rep(NA, n)
  key_vec <- as.character(sort(as.numeric(hash::keys(res$ht_neighbor))))
  key_vec <- key_vec[order(res$df_res$order_rec, decreasing = F)]
  
  for(i in 1:length(key_vec)){
    vec_from[i] <- vec_time[as.numeric(key_vec[i])]
    vec_to[i] <- vec_time[res$ht_neighbor[[key_vec[i]]][1]]
  }
  
  graphics::plot(NA, xlim = c(1,n), ylim = range(vec_time), ...)
  for(i in 1:n){
    if(res$df_res$order_rec[as.numeric(key_vec[i])] == 0) {
      graphics::points(i, vec_from[i], pch = 16, col = "red")
    } else {
      col <- ifelse(vec_from[i] < vec_to[i], "black", "green")
      suppressWarnings(graphics::arrows(x0 = i, y0 = vec_from[i], x1 = i, y1 = vec_to[i], 
                       length = arrow_length, col = col))
    }
  }
  
  invisible()
}