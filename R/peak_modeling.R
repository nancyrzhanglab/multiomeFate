peak_mixture_modeling <- function(bandwidth,
                                  cutmat, # rows = cells, columns = basepairs
                                  peak_locations,
                                  peak_prior,
                                  peak_width,
                                  discretization_stepsize = bandwidth/5,
                                  bool_lock_within_peak = T, 
                                  bool_freeze_prior = F, # set to T if optimization is done to too few fragments. this is the prior over peaks
                                  max_iter = 100,
                                  min_prior = 0.01,
                                  num_peak_limit = 4,
                                  return_assignment_mat = F, # set to T usually for only debugging purposes
                                  return_bin_mat = F, # set to T usually for only debugging purposes
                                  tol = 1e-6,
                                  verbose = 1){
  stopifnot(length(bin_midpoints) %% 2 == 1,
            inherits(cutmat, c("dgCMatrix", "matrix")))
  
  # initial assignment
  num_bins <- length(bin_midpoints)
  bin_mat <- .compute_frag_peak_matrix(
    bool_lock_within_peak = bool_lock_within_peak,
    cutmat = cutmat,
    num_peak_limit = num_peak_limit,
    peak_locations = peak_locations,
    peak_width = peak_width
  )
  num_frags <- nrow(bin_mat)
  
  grenander_obj <- .initialize_grenander(bandwidth = bandwidth,
                                         bin_mat = bin_mat,
                                         discretization_stepsize = discretization_stepsize)
  
  # start iteration
  if(verbose > 0) print("Starting initialization")
  iter <- 1
  loglikelihood_vec <- .compute_loglikelihood(
    bin_mat = bin_mat,
    grenander_obj = grenander_obj,
    prior_vec = peak_prior
  )
  if(verbose > 1) print(paste0("Initial log-likelihood: ", round(loglikelihood_vec[1],2)))
  prior_vec <- peak_prior
  
  # TODO: Return if there are no fragments
  
  while(TRUE){
    if(iter > max_iter) break()
    
    if(verbose > 1) print("E-step")
    assignment_mat <- .e_step(
      bin_mat = bin_mat,
      grenander_obj = grenander_obj,
      prior_vec = prior_vec
    )
    
    if(verbose > 1) print("M-step")
    grenander_obj_new <- .m_step(
      assignment_mat = assignment_mat,
      bin_mat = bin_mat,
      num_bins = num_bins
    )
    if(!bool_freeze_prior) { prior_vec <- .compute_prior(assignment_mat = assignment_mat, min_prior = min_prior) }
   
    if(verbose) print("Computing likelihood")
    loglikelihood_val <- .compute_loglikelihood(
      bin_mat = bin_mat,
      grenander_obj = grenander_obj,
      prior_vec = peak_prior
    )
    
    loglikelihood_vec <- c(loglikelihood_vec, loglikelihood_val)
    iter <- length(loglikelihood_vec)
    if(length(loglikelihood_vec) >= 2){
      if(verbose > 0) print(paste0("Iteration: ", iter, ", log-likelihood: ", round(loglikelihood_vec[iter],2)))
      if(abs(loglikelihood_vec[iter] - loglikelihood_vec[iter-1]) <= tol) break()
    }
    grenander_obj <- grenander_obj_new
  }
  
  if(!return_bin_mat) bin_mat <- NULL
  if(!return_assignment_mat) assignment_mat <- NULL
  
  structure(list(assignment_mat = assignment_mat,
                 bin_mat = bin_mat,
                 grenander_obj = grenander_obj,
                 iter = length(loglikelihood_vec),
                 loglikelihood_val = loglikelihood_vec[length(loglikelihood_vec)],
                 loglikelihood_vec = loglikelihood_vec,
                 num_frags = num_frags,
                 prior_vec = prior_vec),
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
                               peak_mat,
                               min_prior = 0.01){
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
  
  prior_vec <- count_vec/sum(count_vec)
  if(any(prior_vec <= min_prior)){
    if(min_prior*length(prior_vec) > 1) {
      min_prior <- 1/(2*length(prior_vec))
    }
    prior_vec <- prior_vec + min_prior/(1-min_prior*length(prior_vec))
    prior_vec <- prior_vec/sum(prior_vec)
  }
  names(prior_vec) <- paste0("p:", 1:nrow(peak_mat))
  
  stopifnot(all(prior_vec > 0))
  prior_vec
}

##################

# compute a frag-by-peak matrix that denotes the distance from each frag to each peak, and NA otherwise
.compute_frag_peak_matrix <- function(bool_lock_within_peak,
                                      cutmat,
                                      num_peak_limit,
                                      peak_locations,
                                      peak_width){
  # a matrix with nrow = number of fragments, and ncol = number of peaks
  n <- nrow(cutmat)
  p <- length(peak_locations)
  cutmat_t <- Matrix::t(cutmat) # basepair-by-cell
  
  # a list, where each entry is a matrix of "fragments by peaks" where its entry is which distance-bin it is
  cell_list <- lapply(1:n, function(i){
    fragment_idx <- .nonzero_col(mat = cutmat_t,
                                 col_idx = i,
                                 bool_value = F)
    if(length(fragment_idx) == 0) return(numeric(0))
    fragment_locations <- as.numeric(colnames(cutmat))[fragment_idx]
    frag_len <- length(fragment_locations)
    
    # for each fragment, return a vector (length of peaks) of the distance to the peak
    tmp <- sapply(fragment_locations, function(j){
      peak_locations - j
    })
    
    if(!is.matrix(tmp)) tmp <- matrix(tmp, nrow = frag_len, ncol = p)
    stopifnot(all(dim(tmp) == c(frag_len,p)))
    rownames(tmp) <- paste0("c:", i, "_loc:", fragment_locations)
    tmp
  })
  cell_list <- cell_list[sapply(cell_list, length) > 0]
  
  res <- do.call(rbind, cell_list)
  colnames(res) <- paste0("p:", 1:p)
  
  # apply bool_lock_within_peak
  if(bool_lock_within_peak) {
    for(i in 1:nrow(res)){
      idx <- which(abs(res[i,]) <= peak_width/2)
      if(length(idx) > 0){
        if(length(idx) > 1) idx <- sample(idx, 1)
        res[i,-idx] <- NA
      }
    }
  }
  
  # apply num_peak_limit
  if(!is.na(num_peak_limit)){
    for(i in 1:nrow(res)){
      # fragments to the left of the peak
      idx <- which(res[i,] >= 0)
      if(length(idx) > num_peak_limit) {
        dist_vec <- abs(res[i,idx])
        cutoff_dist <- sort(dist_vec, decreasing = F)[num_peak_limit]
        rm_idx <- idx[which(dist_vec > cutoff_dist)]
        res[i,rm_idx] <- NA
      }
      
      # fragments to the right of the peak
      idx <- which(res[i,] <= 0)
      if(length(idx) > num_peak_limit) {
        dist_vec <- abs(res[i,idx])
        cutoff_dist <- sort(dist_vec, decreasing = F)[num_peak_limit]
        rm_idx <- idx[which(dist_vec > cutoff_dist)]
        res[i,rm_idx] <- NA
      }
    }
  }
  
  # convert all signed distances into unsigned distances
  abs(res)
}

.initialize_grenander <- function(bandwidth,
                                  bin_mat,
                                  discretization_stepsize){
  idx <- which(!is.na(bin_mat))
  values <- as.numeric(bin_mat[idx])
  weights <- rep(1, length(values))
  
  estimate_grenander(values = values,
                     weights = weights,
                     bandwidth = bandwidth,
                     discretization_stepsize = discretization_stepsize)
}

.compute_loglikelihood <- function(bin_mat,
                                   grenander_obj,
                                   prior_vec,
                                   tol = 1e-6){
  m <- nrow(bin_mat) # number of fragments
  p <- ncol(bin_mat) # number of peaks
  stopifnot(abs(sum(prior_vec) - 1) <= tol)
  
  # compute all the densities based on grenander
  idx <- which(!is.na(bin_mat))
  if(length(idx) == 0 || m == 0) return(NA) # if no fragments
  for(i in 1:idx){
    bin_mat[i] <- evaluate_grenander(
      obj = grenander_obj,
      x = bin_mat[i]
    )
  }
  if(length(idx) != prod(dim(bin_mat))) {
    bin_mat[-idx] <- 0
  }
  
  tmp <- bin_mat %*% diag(prior_vec)
  sum_vec <- Matrix::rowSums(tmp)
  sum(log(sum_vec[sum_vec >= tol]))
}

# compute assignment_mat
.e_step <- function(bin_mat,
                    grenander_obj,
                    prior_vec,
                    tol = 1e-6){
  m <- nrow(bin_mat) # number of fragments
  p <- ncol(bin_mat) # number of peaks
  stopifnot(abs(sum(prior_vec) - 1) <= tol, m > 0, p > 0)
  
  # compute all the densities based on grenander
  idx <- which(!is.na(bin_mat))
  for(i in 1:idx){
    bin_mat[i] <- evaluate_grenander(
      obj = grenander_obj,
      x = bin_mat[i]
    )
  }
  if(length(idx) != prod(dim(bin_mat))) {
    bin_mat[-idx] <- 0
  }
  
  tmp <- bin_mat %*% diag(prior_vec)
  sum_vec <- rowSums(tmp)
  idx <- which(sum_vec <= tol)
  if(length(idx) > 0) sum_vec[sum_vec <= tol] <- 1
  assignment_mat <-  diag(1/sum_vec) %*% tmp
  if(length(idx) > 0) assignment_mat[idx,] <- NA
  
  rownames(assignment_mat) <- rownames(bin_mat)
  colnames(assignment_mat) <- colnames(bin_mat)
  assignment_mat
}

.m_step <- function(assignment_mat,
                    bandwidth,
                    bin_mat,
                    discretization_stepsize){
  idx <- which(!is.na(bin_mat))
  stopifnot(all(dim(assignment_mat) == dim(bin_mat)), length(idx) > 0)
  
  values <- as.numeric(bin_mat[idx])
  weights <- as.numeric(assignment_mat[idx])
  
  estimate_grenander(values = values,
                     weights = weights,
                     bandwidth = bandwidth,
                     discretization_stepsize = discretization_stepsize)
}

.compute_prior <- function(assignment_mat,
                           min_prior){
  prior_vec <- Matrix::colSums(assignment_mat, na.rm = T)
  prior_vec <- prior_vec/sum(prior_vec)
  
  if(any(prior_vec <= min_prior)){
    if(min_prior*length(prior_vec) > 1) {
      min_prior <- 1/(2*length(prior_vec))
    }
    prior_vec <- prior_vec + min_prior/(1-min_prior*length(prior_vec))
    prior_vec <- prior_vec/sum(prior_vec)
  }
  
  stopifnot(all(prior_vec > 0))
  prior_vec
}