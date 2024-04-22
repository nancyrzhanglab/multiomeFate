generate_simulation_plastic <- function(embedding_mat, 
                                        bool_add_randomness = TRUE, 
                                        coefficient_intercept = 0, 
                                        embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
                                        fatefeatures_coefficient_vec = NULL,
                                        fatefeatures_mat = NULL, 
                                        lineage_spread = NA, # NA or a value 1 or larger. "1" means no spread
                                        num_lineages = 10, 
                                        tol = 1e-06, 
                                        verbose = 0) 
{
  K <- num_lineages
  n <- nrow(embedding_mat)
  d <- ncol(embedding_mat)
  if(all(!is.null(fatefeatures_mat))){
    d2 <- ncol(fatefeatures_mat)
    stopifnot(nrow(fatefeatures_mat) == nrow(embedding_mat),
              length(fatefeatures_coefficient_vec) == d2)
  } else {
    d2 <- 0
  }
  
  rho <- lineage_spread
  if (length(rownames(embedding_mat)) == 0){
    rownames(embedding_mat) <- paste0("cell:", 1:n)
  }
  
  stopifnot(is.na(rho) || rho >= 1, 
            d > 1, 
            length(embedding_coefficient_vec) == d)
  
  if (verbose > 0) 
    print("Step 1: Computing fate potential of all the cells")
  cell_contribution <- as.numeric(embedding_mat %*% embedding_coefficient_vec)
  if(d2 > 0){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution_truth <- exp(cell_contribution + coefficient_intercept)
  names(cell_contribution_truth) <- rownames(embedding_mat) 
  cell_contribution_random <- cell_contribution_truth
  
  if (bool_add_randomness) {
    if (verbose > 0) 
      print("Step 1b: (Optional) Adding randomness")
    cell_contribution_random <- stats::rpois(n = length(cell_contribution_random), 
                                             lambda = cell_contribution_random)
    names(cell_contribution_random) <- rownames(embedding_mat) 
  }
  
  if (verbose > 0) 
    print("Step 2: Computing probability of cells in each lineage")
  tmp <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    rho = rho,
    verbose = verbose - 1
  )
  prob_mat <- tmp$prob_mat

  if (verbose > 0) 
    print("Step 3: Assigning cells to lineages")
  lineage_assignment <- .assign_plastic_lineages(prob_mat)
  # reorder the lineage assignment
  lineage_assignment <- lineage_assignment[names(cell_contribution_truth)]
  
  lineage_future_size <- sapply(levels(lineage_assignment), 
                                function(lev) {
                                  idx <- which(lineage_assignment == lev)
                                  round(sum(cell_contribution_random[idx]))
                                })
  names(lineage_future_size) <- levels(lineage_assignment)
 
  if (verbose > 0) 
    print("Step 4: Outputting")
  return(
    structure(list(cell_fate_potential = log10(cell_contribution_random + 1), 
                   cell_fate_potential_truth = log10(cell_contribution_truth), 
                   coefficient_intercept = coefficient_intercept,
                   embedding_coefficient_vec = embedding_coefficient_vec, 
                   fatefeatures_coefficient_vec = fatefeatures_coefficient_vec, 
                   fatefeatures_mat = fatefeatures_mat,
                   lineage_assignment = lineage_assignment, 
                   lineage_future_size = lineage_future_size,
                   prob_mat = prob_mat),
              class = "multiomeFate_simulation_plastic")
  )
}

#################

.compute_plastic_probabilities <- function(
    cell_contribution_truth,
    num_lineages,
    rho,
    verbose = 0
){
  
  n <- length(cell_contribution_truth)
  K <- num_lineages
  mean_val <- mean(cell_contribution_truth)
  sd_val <- stats::sd(cell_contribution_truth)
 
  # compute rho if it's NA
  if(is.na(rho)){
    rho <- (max(abs(cell_contribution_truth - mean_val))/sd_val)/2
  }
  
  lineage_sd_vec <- seq(sd_val*rho, sd_val/rho, length.out = K)
  names(lineage_sd_vec) <- paste0("lineage:", 1:K) 
  
  # reorder the cell contribution
  idx <- .reorder_by_contribution(abs(cell_contribution_truth - mean_val))
  cell_contribution_truth <- cell_contribution_truth[idx]
 
  prob_mat <- matrix(0, nrow = n, ncol = K)
  rownames(prob_mat) <- names(cell_contribution_truth)
  colnames(prob_mat) <- names(lineage_sd_vec)
  
  for(j in 1:K){
    if(verbose > 0 && K > 10 && j %% floor(K/10) == 0) cat('*')
    lineage <- colnames(prob_mat)[j]
    
    prob_mat[,lineage] <- stats::dnorm(
      x = cell_contribution_truth,
      mean = mean_val,
      sd = lineage_sd_vec[lineage],
      log = TRUE
    )
  }
  
  for(i in 1:n){
    prob_mat[i,] <- .log_sum_exp_normalization(prob_mat[i,])
  }
  
  return(
    list(lineage_sd_vec = lineage_sd_vec,
         prob_mat = prob_mat)
  )
}

.reorder_by_contribution <- function(vec){
  n <- length(vec)
  order_dec <- order(vec, decreasing = TRUE)
  order_inc <- order(vec, decreasing = FALSE)
  
  vec <- as.numeric(rbind(order_inc, order_dec))
  vec <- vec[!duplicated(vec)]
  
  return(vec)
}

.assign_plastic_lineages <- function(prob_mat){
  n <- nrow(prob_mat)
  K <- ncol(prob_mat)
  
  maximum_lineage_size <- ceiling(n/K)
  current_size <- rep(0, K)
  names(current_size) <- colnames(prob_mat)
  
  lineage_assignment <- rep(NA, n)
  for(i in 1:n) {
    stopifnot(is.matrix(prob_mat), ncol(prob_mat) >= 1)
    
    prob_vec <- prob_mat[i,]
    if(all(prob_vec <= 1e-6)) 
      prob_vec <- rep(1/ncol(prob_mat), length = ncol(prob_mat))
    
    ## [[ruh-oh -- testing]]
    # prob_vec <- prob_vec^5
    # prob_vec <- prob_vec/sum(prob_vec)
    
    lineage <- colnames(prob_mat)[sample(1:ncol(prob_mat), size = 1, prob = prob_vec)]
    current_size[lineage] <- current_size[lineage] + 1
    lineage_assignment[i] <- lineage
    
    if(current_size[lineage] >= maximum_lineage_size){
      prob_mat <- prob_mat[,-which(colnames(prob_mat) == lineage), drop = FALSE]
    }
  }
  
  lineage_assignment <- factor(lineage_assignment, 
                               levels = names(current_size))
  names(lineage_assignment) <- rownames(prob_mat)
  
  return(lineage_assignment)
}