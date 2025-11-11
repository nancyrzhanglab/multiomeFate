generate_simulation_plastic <- function(embedding_mat, 
                                        bool_add_randomness = TRUE, 
                                        coefficient_intercept = 0, 
                                        embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
                                        fatefeatures_coefficient_vec = NULL,
                                        fatefeatures_mat = NULL, 
                                        lineage_mean_spread = 1, # NA or a value 1 or larger. "1" means no spread
                                        lineage_sd_spread = NA, # NA or a numeric. Lineage 1 is lineage_sd_spread, and Lineage num_lineage is 1/lineage_sd_spread.
                                        num_lineages = 10, 
                                        tol = 1e-06, 
                                        verbose = 0) {
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
  
  gamma <- lineage_mean_spread
  rho <- lineage_sd_spread
  if (length(rownames(embedding_mat)) == 0){
    rownames(embedding_mat) <- paste0("cell:", 1:n)
  }
  
  stopifnot(d > 1, 
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
    gamma = gamma,
    rho = rho,
    verbose = verbose - 1
  )
  prob_mat <- tmp$prob_mat; rho <- tmp$rho
  
  if (verbose > 0) 
    print("Step 3: Assigning cells to lineages")
  enforce_equal_size <- !is.na(gamma)
  lineage_assignment <- .assign_plastic_lineages(enforce_equal_size = enforce_equal_size,
                                                 prob_mat = prob_mat)
  # reorder the lineage assignment
  lineage_assignment <- lineage_assignment[names(cell_contribution_truth)]
  
  lineage_future_size <- sapply(levels(lineage_assignment), 
                                function(lev) {
                                  idx <- which(lineage_assignment == lev)
                                  round(sum(cell_contribution_random[idx]))
                                })
  names(lineage_future_size) <- levels(lineage_assignment)
  
  summary_mat <- .compute_summary_lineages(cell_fate_potential_truth = log10(cell_contribution_truth),
                                           lineage_assignment = lineage_assignment,
                                           lineage_future_size = lineage_future_size)
  
  if (verbose > 0) 
    print("Step 4: Outputting")
  return(
    structure(list(cell_fate_potential = log10(cell_contribution_random + 1), 
                   cell_fate_potential_truth = log10(cell_contribution_truth), 
                   coefficient_intercept = coefficient_intercept,
                   embedding_mat = embedding_mat,
                   fatefeatures_coefficient_vec = fatefeatures_coefficient_vec, 
                   fatefeatures_mat = fatefeatures_mat,
                   lineage_assignment = lineage_assignment, 
                   lineage_future_size = lineage_future_size,
                   lineage_sd_spread = rho,
                   prob_mat = prob_mat,
                   summary_mat = summary_mat),
              class = "multiomeFate_simulation_plastic")
  )
}

#################

.compute_plastic_probabilities <- function(
    cell_contribution_truth,
    num_lineages,
    gamma,
    rho,
    verbose = 0
){
  
  stopifnot(all(cell_contribution_truth > 0))
  
  if(is.na(rho) & is.na(gamma)){
    warning("lineage_mean_spread and lineage_sd_spread are both NA. Be aware the method will set lineages have large means AND large variances, and lineages to have small means and small variances.")
  }
  if(!is.na(gamma) && abs(gamma - 1) > 1e-4){
    warning("lineage_mean_spread can only handle NA or 1. Defaulting to 1")
  }
  
  cell_contribution_truth <- log(cell_contribution_truth)
  n <- length(cell_contribution_truth)
  K <- num_lineages
  mean_val <- mean(cell_contribution_truth)
  sd_val <- stats::sd(cell_contribution_truth)

  # compute gamma if it's NA
  if(is.na(gamma)){
    lineage_mean_vec <- stats::quantile(cell_contribution_truth, 
                                        probs = seq(1, 0, length.out = num_lineages))
    sd_val <- sd_val/4
  } else {
    lineage_mean_vec <- rep(mean_val, length = num_lineages)
  }
  
  # compute rho if it's NA
  if(is.na(rho)){
    rho <- (max(abs(cell_contribution_truth - mean_val))/sd_val)/2
  }
  lineage_sd_vec <- seq(sd_val*rho, sd_val/rho, length.out = K)
  
  
  names(lineage_mean_vec) <- paste0("lineage:", 1:K) 
  names(lineage_sd_vec) <- names(lineage_mean_vec)
  
  if(verbose > 1){
    print("Lineage means: ")
    print(round(lineage_mean_vec,2))
    print("Lineage std's: ")
    print(lineage_sd_vec)
  }

  
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
      mean = lineage_mean_vec[lineage],
      sd = lineage_sd_vec[lineage],
      log = TRUE
    )
  }
  
  for(i in 1:n){
    prob_mat[i,] <- .log_sum_exp_normalization(prob_mat[i,])
  }
  
  return(
    list(lineage_sd_vec = lineage_sd_vec,
         prob_mat = prob_mat,
         rho = rho)
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

.assign_plastic_lineages <- function(enforce_equal_size,
                                     prob_mat){
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
    
    lineage <- colnames(prob_mat)[sample(1:ncol(prob_mat), size = 1, prob = prob_vec)]
    current_size[lineage] <- current_size[lineage] + 1
    lineage_assignment[i] <- lineage
    
    if(enforce_equal_size & current_size[lineage] >= maximum_lineage_size){
      prob_mat <- prob_mat[,-which(colnames(prob_mat) == lineage), drop = FALSE]
    }
  }
  
  lineage_assignment <- factor(lineage_assignment, 
                               levels = names(current_size))
  names(lineage_assignment) <- rownames(prob_mat)
  
  return(lineage_assignment)
}

.compute_summary_lineages <- function(cell_fate_potential_truth,
                                      lineage_assignment,
                                      lineage_future_size){
  stopifnot(all(names(cell_fate_potential_truth) == names(lineage_assignment)))
  stopifnot(all(sort(unique(lineage_assignment)) == names(lineage_future_size)))
  
  mean_val <- sapply(levels(lineage_assignment), function(lineage){
    mean(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  med_val <- sapply(levels(lineage_assignment), function(lineage){
    stats::median(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  sd_val <- sapply(levels(lineage_assignment), function(lineage){
    stats::sd(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  range_val <- sapply(levels(lineage_assignment), function(lineage){
    diff(range(cell_fate_potential_truth[lineage_assignment == lineage]))
  })
  
  lineage_future_size <- lineage_future_size[levels(lineage_assignment)]
  
  mat <- rbind(mean_val, med_val, sd_val, range_val, lineage_future_size)
  rownames(mat) <- c("mean", "median", "sd", "range", "future_size")
  colnames(mat) <- levels(lineage_assignment)
  
  return(mat)
}