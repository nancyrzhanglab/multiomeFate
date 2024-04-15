generate_simulation_attachFuture <- function(
    coefficient_intercept,
    coefficient_vec,
    future_cell_embedding_mat,
    lineage_assignment,
    previous_lineage_assignment,
    previous_cell_embedding_mat,
    lineage_spread = 1,
    num_pushforward_training_iter = 20,
    num_subsamples = 200,
    verbose = 0
){
  stopifnot(ncol(future_cell_embedding_mat) > 1,
            length(rownames(mapping_mat)) > 0,
            length(colnames(mapping_mat)) > 0)
  
  cell_contribution <- ceiling(exp(as.numeric(embedding_mat %*% coefficient_vec) + coefficient_intercept))
  
  if(verbose > 0) print("Step 1: Adjusting all the ingredients")
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- exp(as.numeric(embedding_mat %*% coefficient_vec) + new_coefficient_intercept)
  names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  cell_contribution_rounded <- ceiling(cell_contribution)
  non_zero_idx <- which(cell_contribution_rounded > 0)
  stopifnot(length(non_zero_idx) > 1)
  cell_contribution_rounded <- cell_contribution_rounded[non_zero_idx]
  previous_cell_embedding_mat <- previous_cell_embedding_mat[non_zero_idx,]
  if(verbose > 0) print(paste0("There are ", num_future_cells, " future cells, and the sum of potentials is ", sum(cell_contribution_rounded)))
  
  if(verbose > 0) print("Step 2: Recompute all the adjusted ingredients from generate_simulation")
  cell_fate_potential <- log10(cell_contribution)
  if(length(rownames(embedding_mat)) > 0) 
    names(cell_contribution) <- rownames(embedding_mat)
  lineage_future_size <- sapply(levels(lineage_assignment), function(lev){
    idx <- which(lineage_assignment == lev)
    round(sum(cell_contribution[idx]))
  })
  names(lineage_future_size) <- levels(lineage_assignment)
  
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
    sd_vec = sd_vec
  )
  
  if(verbose > 0) print("Step 5: Assigning future cells to a lineage")
  tmp <- .assign_future_to_previous(
    mapping_mat = mapping_mat,
    previous_cell_contribution = cell_contribution_rounded
  )
  prev_lineage_assignment = tmp$prev_lineage_assignment
  
  list(cell_fate_potential = cell_fate_potential,
       coefficient_intercept = new_coefficient_intercept,
       lineage_future_size = lineage_future_size,
       prev_lineage_assignment = prev_lineage_assignment,)
}

#######################################
#######################################

.adjust_coefficient_intercept <- function(
    cell_contribution,
    coefficient_intercept,
    num_future_cells
){
  ratio <- sum(cell_contribution)/num_future_cells
  coefficient_intercept - log(ratio)
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
  
  pushforward_list[[which.min(fit_vec)]]$pushforward_res
}

.pushforward_func_constructor <- function(a, b){
  stopifnot(length(a) == 1)
  
  function(vec){
    stopifnot(length(vec) == length(b))
    a * vec + b
  }
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
  
  list(a = a,
       b = b,
       pushforward_func = .pushforward_func_constructor(a = a, b = b))
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
    if(verbose > 0 && n > 10 && i %% floor(n/10) == 0) cat('*')
    mean_vec <- pushforward_func(previous_cell_embedding_mat[i,])
    
    for(j in 1:m){
      mapping_mat[i,j] <- .dmvnorm(x = future_cell_embedding_mat[j,],
                                   mean = mean_vec,
                                   sigma = lineage_spread*diag(sd_vec),
                                   log = TRUE)
    }
  }
  
  mapping_mat
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
  prev_lineage_size <- rep(0, length = n)
  names(prev_lineage_size) <- rownames(mapping_mat)
  prev_lineage_assignment <- rep(NA, length = m)
  names(prev_lineage_assignment) <- colnames(mapping_mat)
  
  # rearrange the columns of mapping_mat
  colsum_vec <- colSums(mapping_mat)
  mapping_mat <- mapping_mat[,order(colsum_vec, decreasing = TRUE)]
  for(j in 1:m){
    sample_prev_name <- sample(rownames(mapping_mat), size = 1, prob = exp(mapping_mat[,j]))
    prev_lineage_assignment[colnames(mapping_mat)[j]] <- sample_prev_name
    prev_lineage_size[sample_prev_name] <- prev_lineage_size[sample_prev_name] + 1
    
    if(prev_lineage_size[sample_prev_name] >= previous_cell_contribution[sample_prev_name]){
      rm_idx <- which(rownames(mapping_mat) == sample_prev_name)
      mapping_mat <- mapping_mat[-rm_idx,,drop = FALSE]
    }
  }
  
  list(prev_lineage_assignment = prev_lineage_assignment,
       prev_lineage_size = prev_lineage_size)
}




