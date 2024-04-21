generate_simulation_attachFuture <- function(
    coefficient_intercept,
    embedding_coefficient_vec,
    future_cell_embedding_mat,
    lineage_assignment,
    previous_cell_embedding_mat,
    fatefeatures_coefficient_vec = NULL,
    fatefeatures_mat = NULL, 
    lineage_spread = 1,
    num_pushforward_training_iter = 20,
    num_subsamples = 200,
    verbose = 0
){
  stopifnot(ncol(future_cell_embedding_mat) > 1,
            is.factor(lineage_assignment))
  
  cell_contribution <- coefficient_intercept + as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) 
  if(!all(is.null(fatefeatures_coefficient_vec))){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution <- exp(cell_contribution)
  
  if(verbose > 0) print("Step 1: Adjusting all the ingredients")
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- new_coefficient_intercept + as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) 
  names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  if(!all(is.null(fatefeatures_coefficient_vec))){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution <-  exp(cell_contribution)
  cell_contribution_rounded <- round(cell_contribution)
  
  # remove all the cells that don't contribute to the future timepoint
  non_zero_idx <- which(cell_contribution_rounded > 0)
  stopifnot(length(non_zero_idx) > 1)
  cell_contribution_rounded <- cell_contribution_rounded[non_zero_idx]
  lineage_assignment <- droplevels(lineage_assignment[non_zero_idx])
  previous_cell_embedding_mat <- previous_cell_embedding_mat[non_zero_idx,]
  if(verbose > 0) print(paste0("There are ", num_future_cells, " future cells, and the sum of potentials is ", sum(cell_contribution_rounded)))
  
  if(verbose > 0) print("Step 2: Recompute all the adjusted ingredients from generate_simulation")
  cell_fate_potential <- log10(cell_contribution)
  if(length(rownames(previous_cell_embedding_mat)) > 0) 
    names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  future_lineage_size <- sapply(levels(lineage_assignment), function(lev){
    idx <- which(lineage_assignment == lev)
    round(sum(cell_contribution[idx]))
  })
  names(future_lineage_size) <- levels(lineage_assignment)
  
  if(verbose > 0) print("Step 3: Determining the push-forward function")
  pushforward_res <- .compute_pushforward(
    cell_contribution = cell_contribution_rounded,
    future_cell_embedding_mat = future_cell_embedding_mat,
    num_pushforward_training_iter = num_pushforward_training_iter,
    num_subsamples = num_subsamples,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    verbose = verbose - 1
  )
  
  if(verbose > 0) print("Step 4: Computing the mapping of each previous cell to a future cell")
  sd_vec <- apply(future_cell_embedding_mat, 2, stats::sd)
  mapping_mat <- .compute_previous_to_future_mapping(
    future_cell_embedding_mat = future_cell_embedding_mat,
    lineage_spread = lineage_spread,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    pushforward_func = pushforward_res$pushforward_func,
    sd_vec = sd_vec,
    verbose = verbose - 1
  )
  
  if(verbose > 0) print("Step 5: Assigning future cells to a lineage")
  tmp <- .assign_future_to_previous(
    mapping_mat = mapping_mat,
    previous_cell_contribution = cell_contribution_rounded,
    verbose = verbose - 1
  )
  future_cell_assignment <- tmp$future_cell_assignment
  prev_cell_num_progenitor <- tmp$prev_cell_num_progenitor
  
  stopifnot(length(prev_cell_num_progenitor) == length(lineage_assignment))
  future_lineage_size <- sapply(levels(lineage_assignment), function(lev){
    idx <- which(lineage_assignment == lev)
    sum(prev_cell_num_progenitor[idx])
  })
  names(future_lineage_size) <- levels(lineage_assignment)
  
  return(
    structure(list(cell_fate_potential = cell_fate_potential,
                   coefficient_intercept = new_coefficient_intercept,
                   future_cell_assignment = future_cell_assignment,
                   future_lineage_size = future_lineage_size,
                   mapping_mat = round(mapping_mat*1e3),
                   prev_cell_num_progenitor = prev_cell_num_progenitor),
              class = "multiomeFate_simulation_future")
  )
}

#######################################
#######################################

.adjust_coefficient_intercept <- function(
    cell_contribution,
    coefficient_intercept,
    num_future_cells,
    interval_add = 0.1,
    max_iter = 100
){
  tmp <- log(num_future_cells) - log(sum(cell_contribution))
  new_cell_contribution <- cell_contribution * exp(tmp)
  
  # now make sure when rounded, it still is larger
  iter <- 0
  interval_value <- 0
  while(sum(round(new_cell_contribution)) < sum(num_future_cells)){
    iter <- iter + 1
    interval_value <- interval_value + interval_add
    new_cell_contribution <- cell_contribution * exp(tmp+interval_value)
    
    if(iter > max_iter) stop("Error with computing intercept")
  }
  
  return(coefficient_intercept + tmp + interval_value)
}

########

.compute_pushforward <- function(
    cell_contribution,
    future_cell_embedding_mat,
    num_pushforward_training_iter,
    num_subsamples,
    previous_cell_embedding_mat,
    previous_cell_potential,
    verbose = 0
){
  n <- nrow(previous_cell_embedding_mat)
  m <- nrow(future_cell_embedding_mat)
  previous_idx <- sample(1:n, 
                         size = num_subsamples, 
                         prob = cell_contribution,
                         replace = TRUE)
  previous_mat <- previous_cell_embedding_mat[previous_idx,]
  
  if(verbose > 0) print("Compute the pushforward functions")
  pushforward_list <- lapply(num_pushforward_training_iter, function(kk){
    future_idx <- sample(1:m, size = num_subsamples, replace = TRUE)
    future_mat <- future_cell_embedding_mat[future_idx,]
    pushforward_res <- .compute_pushforward_fit(
      future_mat = future_mat,
      previous_mat = previous_mat
    )
    
    tmp <- sapply(1:num_subsamples, function(i){
      vec <- pushforward_res$pushforward_func(previous_mat[i,])
      .l2norm(vec - future_mat[i,])^2
    })
    fit <- sum(tmp)
    
    list(fit = fit,
         pushforward_res = pushforward_res)
  })
  fit_vec <- sapply(pushforward_list, function(x){x$fit})
  
  return(pushforward_list[[which.min(fit_vec)]]$pushforward_res)
}

.pushforward_func_constructor <- function(a, b){
  stopifnot(length(a) == 1)
  
  return(
    function(vec){
      stopifnot(length(vec) == length(b))
      a * vec + b
    }
  )
}

.compute_pushforward_fit <- function(
    future_mat,
    previous_mat
){
  stopifnot(all(dim(future_mat) == dim(previous_mat)))
  d <- ncol(future_mat)
  n <- nrow(future_mat)
  
  coef_mat <- sapply(1:d, function(j){
    df <- data.frame(
      x = previous_mat[,j],
      y = future_mat[,j]
    )
    lm_res <- stats::lm(y ~ ., data = df)
    stats::coef(lm_res)
  })
  coef_mat <- t(coef_mat)
  colnames(coef_mat) <- c("b", "a")
  
  a <- stats::median(coef_mat[,"a"])
  b <- sapply(1:d, function(j){
    stats::median(future_mat[,j] - a*previous_mat[,j])
  })
  
  return(
    list(a = a,
         b = b,
         pushforward_func = .pushforward_func_constructor(a = a, b = b))
  )
}

########

.compute_previous_to_future_mapping <- function(
    future_cell_embedding_mat,
    lineage_spread,
    previous_cell_embedding_mat,
    pushforward_func,
    sd_vec,
    verbose = 0
){
  stopifnot(ncol(previous_cell_embedding_mat) == ncol(future_cell_embedding_mat))
  
  n <- nrow(previous_cell_embedding_mat)
  m <- nrow(future_cell_embedding_mat)
  mapping_mat <- matrix(NA, nrow = n, ncol = m)
  if(length(rownames(future_cell_embedding_mat)) > 0){
    rownames(mapping_mat) <- rownames(previous_cell_embedding_mat)
  }
  if(length(rownames(previous_cell_embedding_mat)) > 0){
    colnames(mapping_mat) <- rownames(future_cell_embedding_mat)
  }
  
  for(i in 1:n){
    if(verbose == 1 && n > 10 && i %% floor(n/10) == 0) cat('*')
    if(verbose > 1) print(paste0("i:", i))
    mean_vec <- pushforward_func(previous_cell_embedding_mat[i,])
    mapping_mat[i,] <- .dmvnorm_log_many_samples(
      mean = mean_vec,
      sigma = lineage_spread*diag(sd_vec),
      x_mat = future_cell_embedding_mat
    )
  }
  
  # rearrange the columns
  col_sum <- colSums(mapping_mat)
  mapping_mat <- mapping_mat[,order(col_sum, decreasing = TRUE)]
  
  for(j in 1:m){
    mapping_mat[,j] <- .log_sum_exp_normalization(mapping_mat[,j])
  }
  
  return(mapping_mat)
}

# returns a vector of length nrow(x_mat), the log density for each row of x_mat
.dmvnorm_log_many_samples <- function(mean,
                                      sigma,
                                      x_mat){
  stopifnot(ncol(sigma) == nrow(sigma),
            Matrix::rankMatrix(sigma) == ncol(sigma))
  
  p <- ncol(sigma)
  n <- nrow(x_mat)
  sigma_inv <- solve(sigma)
  determinant_value <- determinant(sigma,
                                   logarithm = TRUE)$modulus
  determinant_value <- as.numeric(determinant_value)
  
  x_mat <- sweep(x_mat, 
                 MARGIN = 2,
                 STATS = mean, 
                 FUN = "-")
  lhs <- x_mat %*% sigma_inv
  rss <- sapply(1:n, function(i){
    lhs[i,] %*% x_mat[i,]
  })
  
  return(- 0.5 * determinant_value - 0.5 * p * log(2 * pi) - 0.5 * rss)
}


########

.assign_future_to_previous <- function(
    mapping_mat,
    previous_cell_contribution, 
    verbose = 0
){
  stopifnot(sum(previous_cell_contribution) >= ncol(mapping_mat),
            length(rownames(mapping_mat)) > 0,
            length(colnames(mapping_mat)) > 0)
  
  m <- ncol(mapping_mat)
  n <- nrow(mapping_mat)
  prev_cell_num_progenitor <- rep(0, length = n)
  names(prev_cell_num_progenitor) <- rownames(mapping_mat)
  future_cell_assignment <- rep(NA, length = m)
  names(future_cell_assignment) <- colnames(mapping_mat)
  
  for(j in 1:m){
    if(verbose > 0 && m > 10 && j %% floor(m/10) == 0) cat('*')
    
    prob_vec <- mapping_mat[,j]
    sample_prev_name <- sample(rownames(mapping_mat), size = 1, prob = prob_vec)
    future_cell_assignment[colnames(mapping_mat)[j]] <- sample_prev_name
    prev_cell_num_progenitor[sample_prev_name] <- prev_cell_num_progenitor[sample_prev_name] + 1
    
    if(prev_cell_num_progenitor[sample_prev_name] >= previous_cell_contribution[sample_prev_name]){
      rm_idx <- which(rownames(mapping_mat) == sample_prev_name)
      mapping_mat <- mapping_mat[-rm_idx,,drop = FALSE]
    }
  }
  
  return(
    list(future_cell_assignment = future_cell_assignment,
         prev_cell_num_progenitor = prev_cell_num_progenitor)
  )
}




