# return a smoothed density returned at a lot of different points
estimate_grenander <- function(values,
                               weights,
                               scaling_factor = 1){
  cdf_empirical <- .weighted_cdf(
    values = values/scaling_factor,
    weights = weights
  )
  
  obj <- .compute_decreasing_density(
    cdf = cdf_empirical$cdf,
    x = cdf_empirical$x
  )
  
  obj$log_pdf <- log(obj$pdf)
  obj$scaling_factor <- scaling_factor
  class(obj) <- "grenander"
  obj
}

evaluate_grenander <- function(obj,
                               x,
                               bool_log = F){
  stopifnot(inherits(obj, "grenander"), x>=0)
  
  idx <- max(which(obj$x < x/obj$scaling_factor))
  if(bool_log) obj$log_pdf[idx] else obj$pdf[idx]
}

##########################################

# construct the weighted CDF function
.weighted_cdf <- function(values,
                          weights,
                          tol = 1e-6){
  stopifnot(all(values >= 0), length(values) == length(weights))
  
  mat <- cbind(c(0, values), c(0, weights))
  colnames(mat) <- c("x", "w")
  
  # compute cdf
  mat <- mat[order(mat[,"x"]),]
  cdf_vec <- cumsum(mat[,"w"])
  cdf_vec <- cdf_vec/max(cdf_vec)
  
  # clean up duplicated values
  tmp <- .remove_duplicates(associated_vec = cdf_vec,
                            target_vec = mat[,"x"])
  
  list(cdf = tmp$associated_vec,
       x = tmp$target_vec)
}

# run LCM from https://search.r-project.org/CRAN/refmans/fdrtool/html/gcmlcm.html
# outputs are: x.knots, y.knots, slope.knots
.compute_decreasing_density <- function(cdf, x){
  res <- fdrtool::gcmlcm(x = x,
                         y = cdf,
                         type = "lcm")
  
  list(x = res$x,
       pdf = c(res$slope.knots,0))
}

#######################################

# assumes target_vec is sorted already, and associated_vec is either NULL or of equal length to target_vec
.remove_duplicates <- function(associated_vec,
                               target_vec,
                               tol = 1e-6){
  stopifnot(length(associated_vec) == length(target_vec))
  
  diff_vec <- diff(target_vec)
  idx <- which(diff_vec <= tol)
  if(length(idx) > 0){
    target_vec <- target_vec[-idx]
    associated_vec <- associated_vec[-idx]
  }
  
  list(associated_vec = associated_vec,
       target_vec = target_vec)
}
