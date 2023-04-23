# return a smoothed density returned at a lot of different points
estimate_grenander <- function(values,
                               weights,
                               bandwidth = diff(stats::quantile(values, probs = c(0.25,0.75)))/10,
                               discretization_stepsize = bandwidth/5,
                               scaling_factor = 1){
  cdf_empirical <- .weighted_cdf(
    values = values/scaling_factor,
    weights = weights
  )
  
  stepdensity_res <- .compute_decreasing_density(
    cdf = cdf_empirical$cdf,
    x = cdf_empirical$x
  )
  
  smooth_pdf_est <- .smooth_stepdensity(
    bandwidth = bandwidth/scaling_factor,
    discretization_stepsize = discretization_stepsize/scaling_factor,
    stepdensity_res = stepdensity_res
  )
  
  res <- smooth_pdf_est
  right_area <- sum(res$pdf[-length(res$pdf)] * diff(res$x))
  left_area <- sum(res$pdf[-1] * diff(res$x))
  res$param <- c(res$param, 
                 left_area = left_area,
                 right_area = right_area,
                 scaling_factor = scaling_factor)

  res
}

evaluate_grenander <- function(obj,
                               x){
  stopifnot(inherits(obj, "grenander"), x>=0)
  idx <- which.min(abs(obj$x - x/obj$param$scaling_factor))
  obj$pdf[idx]
}

##########################################

# construct the weighted CDF function
.weighted_cdf <- function(values,
                          weights,
                          tol = 1e-6){
  stopifnot(all(values >= 0))
  
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
  fdrtool::gcmlcm(x = x,
                  y = cdf,
                  type = "lcm")
}

# see https://en.wikipedia.org/wiki/Kernel_(statistics)
# inspired by https://arxiv.org/pdf/1512.07445.pdf
# assumes x is already sorted
.epanechnikov_kernel_all <- function(bandwidth,
                                     x,
                                     pdf,
                                     tol = 1e-6){
  stopifnot(length(x) == length(pdf))
  n <- length(x)
  pdf_smooth <- rep(NA, n)
  x_trans <- x/bandwidth
  for(i in 1:n){
    idx <- which(abs(x_trans - x_trans[i]) <= 1)
    x_segment <- x_trans[idx]
    weight_vec <- .epanechnikov_weights(x_segment - x_trans[i])
    pdf_smooth[i] <- sum(weight_vec*pdf[idx])
  }
  
  structure(list(x = x,
                 pdf = pdf_smooth,
                 param = list(bandwidth = bandwidth)),
            class = "grenander")
}

# assumes all values in x_vec are within +/- 1
.epanechnikov_weights <- function(x_vec){
  weight_vec <- 3/4*(1-x_vec^2)
  weight_vec/sum(weight_vec)
}

.smooth_stepdensity <- function(bandwidth,
                                discretization_stepsize,
                                stepdensity_res){
  x <- stepdensity_res$x.knots # assumed to be ordered and starts at 0
  pdf <- stepdensity_res$slope.knots # has a length 1 shorter than x
  pdf <- c(pdf[1], pdf)
  
  # create a lot of discretizations
  x1 <- seq(0, max(x)+2*bandwidth, by = discretization_stepsize)
  x1 <- c(x1, x)
  
  # remove duplicated x1's
  x1 <- sort(unique(x1))
  tmp <- .remove_duplicates(associated_vec = NULL,
                            target_vec = x1)
  x1 <- tmp$target_vec
  
  # now compute all the y's
  pdf1 <- sapply(1:length(x1), function(i){
    idx <- which(x >= x1[i]) 
    if(length(idx) > 0) {
      pdf[min(idx)]
    } else {
      # length of idx will be 0 if x1[i] is larger than all other values
      0
    }
  })
  
  # smooth
  res <- .epanechnikov_kernel_all(
    bandwidth = bandwidth,
    x = x1,
    pdf = pdf1
  )
  
  structure(list(x = res$x,
                 pdf = res$pdf,
                 param = list(
                   bandwidth = bandwidth,
                   discretization_stepsize = discretization_stepsize
                 )),
            class = "grenander")
}

#######################################

# assumes target_vec is sorted already, and associated_vec is either NULL or of equal length to target_vec
.remove_duplicates <- function(associated_vec,
                               target_vec,
                               tol = 1e-6){
  stopifnot(all(is.null(associated_vec)) | length(associated_vec) == length(target_vec))
  
  diff_vec <- diff(target_vec)
  idx <- which(diff_vec <= tol)
  if(length(idx) > 0){
    target_vec <- target_vec[-idx]
    if(!all(is.null(associated_vec))){
      associated_vec <- associated_vec[-idx]
    }
  }
  
  list(associated_vec = associated_vec,
       target_vec = target_vec)
}
