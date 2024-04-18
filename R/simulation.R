generate_simulation<- function(embedding_mat, 
                               bool_add_randomness = TRUE, 
                               coefficient_intercept = 0, 
                               embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
                               fatefeatures_coefficient_vec = NULL,
                               fatefeatures_mat = NULL, 
                               lineage_spread = 1, 
                               lineage_prior = NA, 
                               num_lineages = 10, 
                               tol = 1e-06, 
                               verbose = 0) 
{
  if (all(is.na(lineage_prior))) 
    lineage_prior <- rep(1/num_lineages, length = num_lineages)
  lineage_prior <- lineage_prior/sum(lineage_prior)
  K <- num_lineages
  n <- nrow(embedding_mat)
  d <- ncol(embedding_mat)
  if(all(!is.null(fatefeatures_mat))){
    d2 <- ncol(fatefeatures_mat)
  } else {
    d2 <- 0
  }
  
  rho <- lineage_spread
  if (length(names(lineage_prior)) > 0) {
    warning("Overwriting names in lineage_prior")
  }
  
  names(lineage_prior) <- paste0("lineage:", 1:K)
  stopifnot(d > 1, 
            length(embedding_coefficient_vec) == d, 
            length(fatefeatures_coefficient_vec)==d2, 
            length(lineage_prior) == K, 
            all(lineage_prior >= 0), 
            abs(sum(lineage_prior) - 1) <= tol, 
            rho >= 0)
  
  if (verbose > 0) 
    print("Step 1: Selecting seed cells")
  cluster_idxs <- .kmeans_seed(embedding_mat = embedding_mat, 
                               K = K)
  if (verbose > 0) 
    print("Step 2: Computing Gaussian distributions")
  gaussian_list <- .form_gaussian_distributions(cluster_idxs = cluster_idxs, 
                                                embedding_mat = embedding_mat, rho = rho)
  if (verbose > 0) 
    print("Step 3: Computing posterior distributions")
  prob_mat <- .compute_posteriors(embedding_mat = embedding_mat, 
                                  gaussian_list = gaussian_list, lineage_prior = lineage_prior, 
                                  verbose = verbose - 1)
  if (verbose > 0) 
    print("Step 4: Sampling lineages")
  lineage_assignment <- sapply(1:n, function(i) {
    sample(1:K, size = 1, prob = prob_mat[i, ])
  })
  lineage_assignment <- factor(paste0("lineage:", lineage_assignment), 
                               levels = colnames(prob_mat))
  if (length(rownames(embedding_mat)) > 0) 
    names(lineage_assignment) <- rownames(embedding_mat)
  
  if (verbose > 0) 
    print("Step 5: Computing future lineage size")
  cell_contribution <- as.numeric(embedding_mat %*% embedding_coefficient_vec)
  if(d2 > 0){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution_truth <- exp(cell_contribution + coefficient_intercept)
  if (length(rownames(embedding_mat)) > 0) 
    names(cell_contribution_truth) <- rownames(embedding_mat)
  cell_contribution_random <- cell_contribution_truth
  
  if (bool_add_randomness) {
    if (verbose > 0) 
      print("Step 5b: (Optional) Adding randomness")
    cell_contribution_random <- stats::rpois(n = length(cell_contribution_random), 
                                             lambda = cell_contribution_random)
    if (length(rownames(embedding_mat)) > 0) 
      names(cell_contribution_random) <- rownames(embedding_mat)
  }
  
  lineage_future_size <- sapply(levels(lineage_assignment), 
                                function(lev) {
                                  idx <- which(lineage_assignment == lev)
                                  round(sum(cell_contribution_random[idx]))
                                })
  names(lineage_future_size) <- levels(lineage_assignment)
  
  if (verbose > 0) 
    print("Step 6: Outputting")
  list(cell_fate_potential = log10(cell_contribution_random + 1), 
       cell_fate_potential_truth = log10(cell_contribution_truth), 
       coefficient_intercept = coefficient_intercept,
       embedding_coefficient_vec = embedding_coefficient_vec, 
       fatefeatures_coefficient_vec = fatefeatures_coefficient_vec, 
       fatefeatures_mat = fatefeatures_mat,
       gaussian_list = gaussian_list, 
       lineage_assignment = lineage_assignment, 
       lineage_future_size = lineage_future_size, 
       prob_mat = prob_mat)
}


#############

.kmeans_seed <- function(
    embedding_mat,
    K
){
  kmeans_res <- suppressWarnings(stats::kmeans(embedding_mat, centers = 2*K))
  cluster_idxs <- sapply(1:K, function(k){
    sample(which(kmeans_res$cluster == k), 1)
  })
  names(cluster_idxs) <- paste0("lineage:", 1:K)
  cluster_idxs
}

.gaussian <- function(cov_mat,
                      mean_vec){
  stopifnot(nrow(cov_mat) == ncol(cov_mat),
            nrow(cov_mat) == length(mean_vec))
  
  structure(list(cov = cov_mat, 
                 mean = mean_vec),
            class = "gaussian")
}

.form_gaussian_distribution <- function(
    cluster_idx,
    embedding_mat,
    rho
){
  stopifnot(length(cluster_idx) == 1,
            cluster_idx <= nrow(embedding_mat),
            cluster_idx > 0,
            cluster_idx %% 1 == 0)
  
  mean_vec <- embedding_mat[cluster_idx,]
  sd_vec <- apply(embedding_mat, 2, stats::sd)
  
  cov_mat <- diag(rho*sd_vec^2)
  
  .gaussian(cov_mat = cov_mat, 
            mean_vec = mean_vec)
}

.form_gaussian_distributions <- function(
    cluster_idxs,
    embedding_mat,
    rho
){
  K <- length(cluster_idxs)
  gaussian_list <- lapply(1:K, function(k){
    .form_gaussian_distribution(
      cluster_idx = cluster_idxs[k],
      embedding_mat = embedding_mat,
      rho = rho
    )
  })
  names(gaussian_list) <- paste0("lineage:", 1:K)
  
  gaussian_list
}

.compute_posteriors <- function(
    embedding_mat,
    gaussian_list,
    lineage_prior,
    verbose = 0
){
  n <- nrow(embedding_mat)
  K <- length(lineage_prior)
  
  # all the calculations are done on the log scale
  # we use the log-sum-exp trick: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  prob_mat <- matrix(NA, nrow = n, ncol = K)
  if(length(rownames(embedding_mat)) > 0) 
    rownames(prob_mat) <- rownames(embedding_mat)
  colnames(prob_mat) <- names(lineage_prior)
  
  for(i in 1:n){
    if(verbose > 0 && n > 10 && i %% floor(n/10) == 0) cat('*')
    d_vec <- sapply(1:K, function(k){
      .dmvnorm(x = embedding_mat[i,],
               mean = gaussian_list[[k]]$mean, 
               sigma = gaussian_list[[k]]$cov, 
               log = TRUE, 
               checkSymmetry = FALSE)
    })
    stopifnot(length(lineage_prior) == length(d_vec))
    log_vec <- log(lineage_prior) + d_vec
    
    prob_mat[i,] <- .log_sum_exp_normalization(log_vec)
  }
  
  prob_mat
}

# we use the log-sum-exp trick: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
# vec is the vector of log(un-normalized probabilities), 
# and we wish to output the probabilities (non-negative, sums to 1) 
.log_sum_exp_normalization <- function(x, tol = 1e-6){
  c <- max(x)
  y <- c + log(sum(exp(x-c)))
  res <- exp(x-c)
  res <- res/sum(res)
  
  stopifnot(all(res >= -tol, abs(sum(res)-1) <= tol))
  
  res
}

# from the mvtnorm package: https://github.com/cran/mvtnorm/blob/master/R/mvnorm.R
.dmvnorm <- function (x, 
                      mean = rep(0, p), 
                      sigma = diag(p), 
                      log = FALSE, 
                      checkSymmetry = TRUE)
{
  
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  
  if(!missing(mean)) {
    if(!is.null(dim(mean))) dim(mean) <- NULL
    if (length(mean) != p)
      stop("x and mean have non-conforming size")
  }
  if(!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (checkSymmetry && !Matrix::isSymmetric(
      sigma, 
      tol = sqrt(.Machine$double.eps), 
      check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  
  ## <faster code contributed by Matteo Fasiolo mf364 at bath.ac.uk
  dec <- tryCatch(base::chol(sigma), error=function(e)e)
  if (inherits(dec, "error")) {
    ## warning("cannot compute chol(sigma)"); return(NaN)
    ## behave the same as dnorm(): return Inf or 0
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf # and all other f(.) == 0
  } else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp ^ 2)
    logretval <- - sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if(log) logretval else exp(logretval)
}
