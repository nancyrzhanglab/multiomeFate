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
  if(!is.null(fatefeatures_mat)){
    d2 <- ncol(fatefeatures_mat)
  } else {
    d2 <- 0
  }
  
  rho <- lineage_spread
  if (length(names(lineage_prior)) > 0) {
    warning("Overwriting names in lineage_prior")
  }
  
  names(lineage_prior) <- paste0("lineage:", 1:K)
  stopifnot(d > 1, length(coefficient_vec) == d, 
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
  cell_contribution <-as.numeric(embedding_mat %*% coefficient_vec)
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
       coefficient_intercept = coefficient_intercept, coefficient_vec = coefficient_vec, 
       gaussian_list = gaussian_list, lineage_assignment = lineage_assignment, 
       lineage_future_size = lineage_future_size, prob_mat = prob_mat)
}
