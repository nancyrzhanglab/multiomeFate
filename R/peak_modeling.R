peak_mixture_modeling <- function(bin_limits, # vector of length 2
                                  bin_midpoints, # midpoints of each bin
                                  cutmat, # rows = cells, columns = basepairs
                                  peak_locations,
                                  peak_prior,
                                  bool_freeze_prior = F, # set to T if optimization is done to too few fragments. this is the prior over peaks
                                  max_iter = 100,
                                  return_bin_mat = F, # set to T usually for only debugging purposes
                                  tol = 1e-6,
                                  verbose = 1){
  stopifnot(length(bin_midpoints) %% 2 == 1,
            inherits(cutmat, c("dgCMatrix", "matrix")))
  
  # initial assignment
  num_bins <- length(bin_midpoints)
  bin_mat <- .compute_bin_matrix(
    bin_limits = bin_limits,
    bin_midpoints = bin_midpoints,
    cutmat = cutmat,
    peak_locations = peak_locations
  )
 
  theta_vec <- .initialize_theta(bin_mat = bin_mat,
                                 num_bins = num_bins)
  if(verbose > 2) print(round(theta_vec, 2))
  
  # start iteration
  if(verbose > 0) print("Starting initialization")
  iter <- 1
  likelihood_vec <- .compute_loglikelihood(
    bin_mat = bin_mat,
    prior_vec = peak_prior,
    theta_vec = theta_vec
  )
  theta_diff <- NA
  prior_vec <- peak_prior
  while(TRUE){
    if(iter > max_iter) break()
    
    if(verbose > 1) print("E-step")
    assignment_mat <- .e_step(
      bin_mat = bin_mat,
      prior_vec = prior_vec,
      theta_vec = theta_vec
    )
    
    if(verbose > 1) print("M-step")
    theta_vec_new <- .m_step(
      assignment_mat = assignment_mat,
      bin_mat = bin_mat,
      num_bins = num_bins
    )
    if(bool_freeze_prior) { prior_vec <- colSums(assignment_mat); prior_vec <- prior_vec/sum(prior_vec) }
    if(verbose > 2) print(round(theta_vec_new, 2))
    
    if(verbose) print("Computing likelihood")
    likelihood_val <- .compute_loglikelihood(
      bin_mat = bin_mat,
      prior_vec = prior_vec,
      theta_vec = theta_vec_new
    )
    
    likelihood_vec <- c(likelihood_vec, likelihood_val)
    theta_diff <- c(theta_diff, sum(abs(theta_vec_new - theta_vec)))
    iter <- length(likelihood_vec)
    if(length(likelihood_vec) >= 2){
      if(verbose > 0) print(paste0("Iteration: ", iter, ", likelihood: ", round(likelihood_vec[iter],2)))
      if(abs(likelihood_vec[iter] - likelihood_vec[iter-1]) <= tol &
         theta_diff[iter] <= tol) break()
    }
    theta_vec <- theta_vec_new
  }
  
  if(return_bin_mat){
    bin_mat2 <- bin_mat
  } else {
    bin_mat2 <- NULL
  }
 
  structure(list(assignment_mat = assignment_mat,
                 bin_mat = bin_mat2,
                 iter = length(likelihood_vec),
                 likelihood_vec = likelihood_vec,
                 prior_vec = prior_vec,
                 theta_diff = theta_diff,
                 theta_vec = theta_vec),
            class = "peakDistribution")
}

compute_bin_midpoints <- function(peak_mat,
                                  num_bins = 7){
  width_vec <- sapply(1:nrow(peak_mat), function(i){
    peak_mat[i,"end"] - peak_mat[i,"start"] + 1
  })
  width <- stats::median(width_vec)
  
  res <- seq(-(num_bins-1)/2, (num_bins-1)/2, by=1)*width
  names(res) <- paste0("bin:", (-(num_bins-1)/2):((num_bins-1)/2))
  res
}

compute_peak_locations <- function(peak_mat){
  res <- sapply(1:nrow(peak_mat), function(i){
    round(mean(peak_mat[i,]))
  })
  names(res) <- paste0("p:", 1:nrow(peak_mat))
  res
}

compute_peak_prior <- function(mat,
                               peak_mat){
  peak_bp <- as.numeric(colnames(mat))
  count_vec <- sapply(1:nrow(peak_mat), function(i){
    idx <- intersect(which(peak_bp >= peak_mat[i,"start"]),
                     which(peak_bp <= peak_mat[i,"end"]))
    sum(sapply(idx, function(j){
      length(.nonzero_col(mat = mat,
                          col_idx = j,
                          bool_value = F))
    }))
  })
  
  count_vec <- count_vec/sum(count_vec)
  names(count_vec) <- paste0("p:", 1:nrow(peak_mat))
  count_vec
}

##################

.compute_bin_matrix <- function(bin_limits,
                                bin_midpoints,
                                cutmat,
                                peak_locations){
  stopifnot(all(bin_midpoints >= bin_limits[1]), 
            all(bin_midpoints <= bin_limits[2]))
  # a matrix with nrow = number of fragments, and ncol = number of peaks
  n <- nrow(cutmat)
  p <- length(peak_locations)
  cutmat_t <- Matrix::t(cutmat)
  num_bins <- length(bin_midpoints)
  bin_index <- (-(num_bins-1)/2):((num_bins-1)/2)
  
  # a list, where each entry is a matrix of "fragments by peaks" where its entry is which distance-bin it is
  cell_list <- lapply(1:n, function(i){
    fragment_idx <- .nonzero_col(mat = cutmat_t,
                                 col_idx = i,
                                 bool_value = F)
    if(length(fragment_idx) == 0) return(numeric(0))
    fragment_locations <- as.numeric(colnames(cutmat))[fragment_idx]
    
    tmp <- sapply(fragment_locations, function(j){
      distance_vec <- peak_locations - j
      sapply(distance_vec, function(dist){
        min_val <- min(abs(dist - bin_midpoints))
        if(min_val >= bin_limits[1] & (min_val <= bin_limits[2])){
          bin_index[which.min(abs(dist - bin_midpoints))]
        } else {
          NA
        }
      })
    })
    if(!is.matrix(tmp)) tmp <- matrix(tmp, nrow = length(tmp), ncol = p)
    if(all(dim(tmp) == c(p,length(fragment_locations)))) tmp <- t(tmp)
    rownames(tmp) <- paste0("c:", i, "_loc:", fragment_locations)
    tmp
  })
  cell_list <- cell_list[sapply(cell_list, length) > 0]
  
  res <- do.call(rbind, cell_list)
  colnames(res) <- paste0("p:", 1:p)
  res
}

.initialize_theta <- function(bin_mat,
                              num_bins){
  min_vec <- sapply(1:nrow(bin_mat), function(i){
    vec <- bin_mat[i,]
    if(all(is.na(vec))) return(NA)
    
    avec <- abs(vec)
    min_val <- min(avec, na.rm = T)
    idx <- which(avec == min_val)
    if(length(idx) == 1) {
      return(vec[idx])
    } else {
      return(vec[sample(idx,1)])
    }
  })
  
  bin_index <- (-(num_bins-1)/2):((num_bins-1)/2)
  theta_vec <- sapply(bin_index, function(i){
    length(which(min_vec == i))
  })
  theta_vec <- theta_vec/sum(theta_vec)
  names(theta_vec) <- paste0("bin:", bin_index)
  theta_vec
}

.compute_loglikelihood <- function(bin_mat,
                                   prior_vec,
                                   theta_vec,
                                   tol = 1e-6){
  m <- nrow(bin_mat) # number of fragments
  p <- ncol(bin_mat) # number of peaks
  stopifnot(length(theta_vec) >= max(bin_mat),
            abs(sum(prior_vec) - 1) <= tol)
  
  mask_idx <- which(is.na(bin_mat))
  bin_mat[mask_idx] <- 0
  num_bins <- length(theta_vec)
  reidx <- (num_bins-1)/2+1
  # remember values in bin_mat are from -(num_bins-1)/2 to (num_bins-1)/2, so we need to re-index them
  tmp1 <- matrix(prior_vec[bin_mat + reidx], nrow = m, ncol = p)
  tmp1[mask_idx] <- 0
  tmp2 <- matrix(theta_vec[bin_mat + reidx], nrow = m, ncol = p)
  tmp2[mask_idx] <- 0
  tmp3 <- tmp1 * tmp2
  sum_vec <- Matrix::rowSums(tmp3)
  sum(log(sum_vec[which(sum_vec >= tol)]))
}

# compute assignment_mat
.e_step <- function(bin_mat,
                    prior_vec,
                    theta_vec){
  m <- nrow(bin_mat) # number of fragments
  p <- ncol(bin_mat) # number of peaks
  stopifnot(length(prior_vec) == p)
  
  num_bins <- length(theta_vec)
  reidx <- (num_bins-1)/2+1
  # remember values in bin_mat are from -(num_bins-1)/2 to (num_bins-1)/2, so we need to re-index them
  tmp <- matrix(theta_vec[bin_mat + reidx], nrow = m, ncol = p)
  tmp <- tmp %*% diag(prior_vec)
  sum_vec <- rowSums(tmp)
  assignment_mat <-  diag(1/sum_vec) %*% tmp
  
  rownames(assignment_mat) <- rownames(bin_mat)
  colnames(assignment_mat) <- colnames(bin_mat)
  assignment_mat
}

.m_step <- function(assignment_mat,
                    bin_mat,
                    num_bins){
  stopifnot(all(dim(assignment_mat) == dim(bin_mat)))
  
  bin_index <- (-(num_bins-1)/2):((num_bins-1)/2)
  beta_vec <- sapply(bin_index, function(i){
    # find all the peaks that are the appropriate distance
    idx <- which(bin_mat == i)
    sum(assignment_mat[idx])
  })
  
  beta_vec <- beta_vec/sum(beta_vec)
  names(beta_vec) <- paste0("bin:", bin_index)
  beta_vec
}