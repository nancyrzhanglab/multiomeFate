peak_mixture_modeling <- function(cutmat, # rows = cells, columns = basepairs
                                  peak_locations,
                                  peak_prior,
                                  peak_width,
                                  bool_lock_within_peak = T, 
                                  bool_freeze_prior = F, # set to T if optimization is done to too few fragments. this is the prior over peaks
                                  fragment_locations = NULL, # one of cutmat and fragment_locations should be NULL
                                  max_iter = 100,
                                  min_prior = 0.01,
                                  num_peak_limit = 4,
                                  return_assignment_mat = F, # set to T usually for only debugging purposes
                                  return_dist_mat = F, # set to T usually for only debugging purposes
                                  tol = 1e-6,
                                  verbose = 1){
  stopifnot(all(is.null(cutmat)) || inherits(cutmat, c("dgCMatrix", "matrix")))
  
  # initial assignment
  dist_mat <- .compute_frag_peak_matrix(
    bool_lock_within_peak = bool_lock_within_peak,
    cutmat = cutmat,
    fragment_locations = fragment_locations,
    num_peak_limit = num_peak_limit,
    peak_locations = peak_locations,
    peak_width = peak_width
  )
  num_frags <- nrow(dist_mat)
  scaling_factor <- max(dist_mat@x)/10
  
  grenander_obj <- .initialize_grenander(dist_mat = dist_mat,
                                         scaling_factor = scaling_factor)
  
  # start iteration
  if(verbose > 0) print("Starting initialization")
  iter <- 1
  loglikelihood_vec <- .compute_loglikelihood(
    dist_mat = dist_mat,
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
      dist_mat = dist_mat,
      grenander_obj = grenander_obj,
      prior_vec = prior_vec
    )
    
    if(verbose > 1) print("M-step")
    grenander_obj_new <- .m_step(
      assignment_mat = assignment_mat,
      dist_mat = dist_mat,
      scaling_factor = scaling_factor
    )
    if(!bool_freeze_prior) { prior_vec <- .compute_prior(assignment_mat = assignment_mat, min_prior = min_prior) }
    
    if(verbose) print("Computing likelihood")
    loglikelihood_val <- .compute_loglikelihood(
      dist_mat = dist_mat,
      grenander_obj = grenander_obj_new,
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
  
  if(!return_dist_mat) dist_mat <- NULL
  if(!return_assignment_mat) assignment_mat <- NULL
  
  structure(list(assignment_mat = assignment_mat,
                 dist_mat = dist_mat,
                 grenander_obj = grenander_obj,
                 iter = length(loglikelihood_vec),
                 loglikelihood_val = loglikelihood_vec[length(loglikelihood_vec)],
                 loglikelihood_vec = loglikelihood_vec,
                 num_frags = num_frags,
                 prior_vec = prior_vec),
            class = "peakDistribution")
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
                                      fragment_locations, # one of cutmat and fragment_locations should be NULL
                                      num_peak_limit,
                                      peak_locations,
                                      peak_width){
  stopifnot((!is.null(cutmat) & is.null(fragment_locations)) || (is.null(cutmat) & !is.null(fragment_locations)))
  
  # a matrix with nrow = number of fragments, and ncol = number of peaks
  n <- nrow(cutmat)
  p <- length(peak_locations)
  
  # a list, where each entry is a matrix of "fragments by peaks" where its entry is the distance
  if(is.null(fragment_locations)){
    fragment_locations <- .extract_fragment_from_cutmat(cutmat)
  }
  
  frag_len <- length(fragment_locations)
  res <- sapply(fragment_locations, function(j){
    peak_locations - j
  })
  if(!is.matrix(res)) {
    stopifnot(p == 1)
    res <- matrix(res, nrow = frag_len, ncol = 1)
  } else {
    res <- t(res)
  }
  stopifnot(all(dim(res) == c(frag_len,p)))
  rownames(res) <- sapply(1:nrow(res), function(i){
    paste0("frag:", i, "_loc:", fragment_locations[i])
  })
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
  
  # convert into sparse matrix
  # also convert all signed distances into unsigned distances
  idx_vec <- which(!is.na(res), arr.ind = F)
  idx_mat <- which(!is.na(res), arr.ind = T)
  values <- res[idx_vec]
  values[values == 0] <- 1 # prevent any distances of 0
  res2 <- Matrix::sparseMatrix(
    i = idx_mat[,1],
    j = idx_mat[,2],
    x = abs(values),
    dims = dim(res)
  )
  
  res2
}

.initialize_grenander <- function(dist_mat,
                                  scaling_factor){
  values <- dist_mat@x
  weights <- rep(1, length(values))
  
  estimate_grenander(values = values,
                     weights = weights,
                     scaling_factor = scaling_factor)
}

.compute_loglikelihood <- function(dist_mat,
                                   grenander_obj,
                                   prior_vec,
                                   tol = 1e-6){
  m <- nrow(dist_mat) # number of fragments
  p <- ncol(dist_mat) # number of peaks
  stopifnot(abs(sum(prior_vec) - 1) <= tol)
  
  # compute all the densities based on grenander
  prob_mat <- dist_mat
  vec <- prob_mat@x
  if(length(vec) == 0 || m == 0) return(NA) # if no fragments
  for(i in 1:length(vec)){
    vec[i] <- evaluate_grenander(
      obj = grenander_obj,
      x = vec[i]
    )
  }
  prob_mat@x <- vec
  
  tmp <- .mult_mat_vec(prob_mat, prior_vec)
  sum_vec <- Matrix::rowSums(tmp)
  sum(log(sum_vec[sum_vec >= tol]))
}

# compute assignment_mat
.e_step <- function(dist_mat,
                    grenander_obj,
                    prior_vec,
                    tol = 1e-6){
  m <- nrow(dist_mat) # number of fragments
  p <- ncol(dist_mat) # number of peaks
  stopifnot(abs(sum(prior_vec) - 1) <= tol, m > 0, p > 0)
  
  # compute all the densities based on grenander
  prob_mat <- dist_mat
  vec <- prob_mat@x
  if(length(vec) == 0 || m == 0) return(NA) # if no fragments
  for(i in 1:length(vec)){
    vec[i] <- evaluate_grenander(
      obj = grenander_obj,
      x = vec[i]
    )
  }
  prob_mat@x <- vec
  
  tmp <- .mult_mat_vec(prob_mat, prior_vec)
  sum_vec <- Matrix::rowSums(tmp)
  idx <- which(sum_vec <= tol)
  if(length(idx) > 0) sum_vec[sum_vec <= tol] <- 1
  assignment_mat <-  .mult_vec_mat(1/sum_vec, tmp)
  if(length(idx) > 0) assignment_mat[idx,] <- 0
  
  rownames(assignment_mat) <- rownames(dist_mat)
  colnames(assignment_mat) <- colnames(dist_mat)
  assignment_mat
}

.m_step <- function(assignment_mat,
                    dist_mat,
                    scaling_factor){
  stopifnot(all(dim(assignment_mat) == dim(dist_mat)), length(dist_mat@x) > 0)
  
  values <- dist_mat@x
  weights <- assignment_mat@x
  if(length(values) != length(weights)){
    idx <- which(as.matrix(dist_mat) != 0)
    weights <- assignment_mat[idx]
  }
  
  estimate_grenander(values = values,
                     weights = weights,
                     scaling_factor = scaling_factor)
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

.extract_fragment_from_cutmat <- function(cutmat){
  n <- nrow(cutmat)
  cutmat_t <- Matrix::t(cutmat) # basepair-by-cell
  
  fragment_idx <- unlist(lapply(1:n, function(i){
    .nonzero_col(mat = cutmat_t,
                 col_idx = i,
                 bool_value = F)
  }))
  
  as.numeric(colnames(cutmat))[fragment_idx]
}