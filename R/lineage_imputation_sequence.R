lineage_imputation_sequence <- function(cell_features,
                                        cell_lineage,
                                        lineage_future_count,
                                        lambda_initial = NA,
                                        lambda_sequence_length = 50,
                                        multipler = 1e4,
                                        verbose = 1){
  res <- .compute_initial_parameters(cell_features = cell_features,
                                     cell_lineage = cell_lineage,
                                     lineage_future_count = lineage_future_count,
                                     multipler = multipler)
  coefficient_initial <- res$coefficient_initial
  if(is.na(lambda_initial)) lambda_initial <- res$lambda_initial
  
  lambda_sequence <- exp(seq(log(lambda_initial), 0, length.out = lambda_sequence_length))-1
  fit_list <- vector("list", length = lambda_sequence_length)
  
  for(i in 1:lambda_sequence_length){
    if(verbose > 0) print(paste0("Working on lambda in sequence ", i, " out of ", lambda_sequence_length))
    if(i == 1){
      coefficient_vec <- coefficient_initial
    } else {
      coefficient_vec <- fit_list[[i-1]]$coefficient_vec
    }
    
    tmp <- lineage_imputation(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              coefficient_initial_list = coefficient_vec,
                              lineage_future_count = lineage_future_count,
                              lambda = lambda_sequence[i],
                              random_initializations = 10,
                              upper_randomness = 5,
                              verbose = verbose-1)
    fit_list[[i]] <- tmp$fit
  }
  
  list(fit_list = fit_list,
       lambda_sequence = lambda_sequence)
}


# cell_features simply included for convenience
.compute_initial_parameters <- function(cell_features,
                                        cell_lineage,
                                        lineage_future_count,
                                        lambda_min = 101,
                                        multipler = 10){
  
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  num_lineages <- length(lineage_future_count)
  
  lineage_current_count <- sapply(cell_lineage_idx_list, length)
  stopifnot(all(names(lineage_current_count) == names(lineage_future_count)))
  
  future_total <- sum(lineage_future_count)
  current_total <- sum(lineage_current_count)
  
  term1 <- future_total*(1-log(future_total/current_total))
  term2 <- sum(lineage_future_count * log(lineage_current_count))
  
  lambda_initial <- -multipler*(term1 - term2)/num_lineages
  lambda_initial <- pmax(lambda_initial, lambda_min)
  
  coefficient_initial <- rep(0, ncol(cell_features))
  names(coefficient_initial) <- colnames(cell_features)
  stopifnot("Intercept" %in% names(coefficient_initial))
  coefficient_initial["Intercept"] <- log(future_total/current_total)
  
  list(coefficient_initial = coefficient_initial,
       lambda_initial = lambda_initial)
}


