lineage_imputation <- function(cell_features,
                               cell_lineage,
                               coefficient_initial_list,
                               lineage_future_count,
                               lambda = 0,
                               random_initializations = 10,
                               upper_randomness = 5,
                               verbose = 1){
  # do some preliminary formatting
  if(!is.list(coefficient_initial_list)) coefficient_initial_list <- list(coefficient_initial_list)
  list_len <- length(coefficient_initial_list)
  
  # some cleanup
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count,
                          verbose = verbose)
  cell_features <- tmp$cell_features
  cell_lineage <- tmp$cell_lineage
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  uniq_lineages <- tmp$uniq_lineages
  coefficient_initial_list <- .append_intercept_term(coefficient_initial_list)
  p <- ncol(cell_features)
  
  stopifnot(all(sort(unique(cell_lineage)) == 
                  sort(unique(names(lineage_future_count)))),
            is.matrix(cell_features), nrow(cell_features) == length(cell_lineage),
            all(sapply(coefficient_initial_list, length) == ncol(cell_features)),
            sum(is.na(cell_features)) == 0,
            sum(is.na(cell_lineage)) == 0,
            sum(is.na(lineage_future_count)) == 0)
  for(i in 1:list_len){
    if(length(names(coefficient_initial_list[[i]])) != 0){
      stopifnot(all(names(coefficient_initial_list[[i]]) == colnames(cell_features)))
    } else {
      names(coefficient_initial_list[[i]]) <- colnames(cell_features)
    }
  }
  
  # rearrange arguments
  optim_fn <- function(coefficient_vec,
                       cell_features,
                       cell_lineage,
                       cell_lineage_idx_list,
                       lambda,
                       lineage_future_count){
    .lineage_objective(cell_features = cell_features,
                       cell_lineage = cell_lineage,
                       cell_lineage_idx_list = cell_lineage_idx_list,
                       coefficient_vec = coefficient_vec,
                       lambda = lambda,
                       lineage_future_count = lineage_future_count)
  }
  
  optim_gr <- function(coefficient_vec,
                       cell_features,
                       cell_lineage,
                       cell_lineage_idx_list,
                       lambda,
                       lineage_future_count){
    .lineage_gradient(cell_features = cell_features,
                      cell_lineage = cell_lineage,
                      cell_lineage_idx_list = cell_lineage_idx_list,
                      coefficient_vec = coefficient_vec,
                      lambda = lambda,
                      lineage_future_count = lineage_future_count)
  }
  
  res_list <- vector("list", length = list_len+random_initializations)
  
  for(i in 1:list_len){
    if(verbose > 0) print(paste0("On provided initialization ", i))
    res <- stats::optim(
      par = coefficient_initial_list[[i]],
      fn = optim_fn,
      gr = optim_gr,
      method = "BFGS",
      cell_features = cell_features,
      cell_lineage = cell_lineage,
      cell_lineage_idx_list = cell_lineage_idx_list,
      lambda = lambda,
      lineage_future_count = lineage_future_count
    )
    
    res_vec <- res$par
    names(res_vec) <- colnames(cell_features)
    
    res_list[[i]] <- list(coefficient_initial = coefficient_initial_list[[i]],
                          coefficient_vec = res_vec,
                          convergence = res$convergence,
                          lambda = lambda,
                          objective_val = res$value)
  }
  
  if(random_initializations > 0){
    max_feature <- stats::quantile(abs(cell_features), probs = 0.95)
    num_cells_per_lineage <- sapply(cell_lineage_idx_list, length)
    names(num_cells_per_lineage) <- uniq_lineages
    max_count_ratio <- max(lineage_future_count[uniq_lineages]/num_cells_per_lineage[uniq_lineages])
    max_limit <- 2*log(max_count_ratio)/(p*max_feature)
    for(i in 1:random_initializations){
      if(verbose > 0) print(paste0("On random initialization ", i))
      coef_vec <- pmin(stats::runif(p, min = 0, max = max_limit), upper_randomness)
      names(coef_vec) <- colnames(cell_features)
      
      res <- stats::optim(
        par = coef_vec,
        fn = optim_fn,
        gr = optim_gr,
        method = "BFGS",
        cell_features = cell_features,
        cell_lineage = cell_lineage,
        cell_lineage_idx_list = cell_lineage_idx_list,
        lambda = lambda,
        lineage_future_count = lineage_future_count
      )
      
      res_vec <- res$par
      names(res_vec) <- colnames(cell_features)
      
      res_list[[i+list_len]] <- list(coefficient_initial = coef_vec,
                                     coefficient_vec = res_vec,
                                     convergence = res$convergence,
                                     lambda = lambda,
                                     objective_val = res$value)
    }
  }
  
  obj_vec <- sapply(res_list, function(lis){lis$objective_val})
  if(verbose > 1){
    print("Quantile of all the objective scores")
    print(stats::quantile(obj_vec))
  }
  
  list(fit =  res_list[[which.min(obj_vec)]],
       res_list = res_list)
}

evaluate_loglikelihood <- function(cell_features,
                                   cell_lineage,
                                   coefficient_vec,
                                   lineage_future_count,
                                   lambda = 0){
  
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)
  cell_features <- tmp$cell_features
  cell_lineage <- tmp$cell_lineage
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  uniq_lineages <- tmp$uniq_lineages
  
  .lineage_objective(cell_features = cell_features,
                     cell_lineage = cell_lineage,
                     cell_lineage_idx_list = cell_lineage_idx_list,
                     coefficient_vec = coefficient_vec,
                     lambda = lambda,
                     lineage_future_count = lineage_future_count)
}

#################################

.lineage_cleanup <- function(cell_features,
                             cell_lineage,
                             lineage_future_count,
                             verbose = 0){
  stopifnot(length(colnames(cell_features)) == ncol(cell_features))
  
  # some cleanup
  if(all(sort(unique(names(lineage_future_count))) != sort(unique(cell_lineage)))){
    if(verbose > 0) warning("Lineages in `lineage_future_count` are not the same as those in `cell_lineage`")
    
    uniq_lineages <- sort(intersect(unique(names(lineage_future_count)), unique(cell_lineage)))
    lineage_future_count <- lineage_future_count[names(lineage_future_count) %in% uniq_lineages]
    rm_cell_idx <- which(!cell_lineage %in% uniq_lineages)
    if(length(rm_cell_idx) > 0){
      cell_lineage <- cell_lineage[-rm_cell_idx]
      cell_features <- cell_features[-rm_cell_idx,,drop=F]
    }
  }
  
  # reorganize everything to be in the same order
  uniq_lineages <- sort(unique(names(lineage_future_count)))
  lineage_future_count <- lineage_future_count[uniq_lineages]
  cell_lineage_idx_list <- lapply(uniq_lineages, function(lineage){
    which(cell_lineage == lineage)
  })
  names(cell_lineage_idx_list) <- uniq_lineages

  
  # add intercept to cell_features
  if(!"Intercept" %in% colnames(cell_features)){
    cell_features <- cbind(1, cell_features)
    colnames(cell_features)[1] <- "Intercept"
  }
  
  list(cell_features = cell_features,
       cell_lineage = cell_lineage,
       cell_lineage_idx_list = cell_lineage_idx_list,
       lineage_future_count = lineage_future_count,
       uniq_lineages = uniq_lineages)
}

.append_intercept_term <- function(coefficient_initial_list){
  stopifnot(is.list(coefficient_initial_list))
  
  for(i in 1:length(coefficient_initial_list)){
    if(!"Intercept" %in% names(coefficient_initial_list[[i]])){
      coefficient_initial_list[[i]] <- c(0, coefficient_initial_list[[i]])
      names(coefficient_initial_list[[i]])[1] <- "Intercept"
    }
  }
  
  coefficient_initial_list
}

.lineage_objective <- function(cell_features,
                               cell_lineage,
                               cell_lineage_idx_list,
                               coefficient_vec,
                               lambda,
                               lineage_future_count){
  uniq_lineages <- names(lineage_future_count)
  num_lineages <- length(uniq_lineages)
  cell_names <- rownames(cell_features)
  
  scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec))
  names(scalar1) <- cell_names
  scalar2 <- sapply(uniq_lineages, function(lineage){
    log(sum(scalar1[cell_lineage_idx_list[[lineage]]]))
  })
  
  idx_notintercept <- which(names(coefficient_vec) != "Intercept")
  scalar3 <- .l2norm(coefficient_vec[idx_notintercept])
  
  (sum(scalar1) - sum(lineage_future_count*scalar2))/num_lineages + lambda*scalar3^2
}

.lineage_gradient <- function(cell_features,
                              cell_lineage,
                              cell_lineage_idx_list,
                              coefficient_vec,
                              lambda,
                              lineage_future_count){
  stopifnot(colnames(cell_features) == names(coefficient_vec))
  
  uniq_lineages <- names(cell_lineage_idx_list)
  num_lineages <- length(uniq_lineages)
  cell_names <- rownames(cell_features)
  lineage_future_count_full <- lineage_future_count[cell_lineage]
  
  # keep track of colnames(cell_features) that is not the intercept
  colname_vec <- colnames(cell_features)
  colname_vec <- colname_vec[colname_vec != "Intercept"]
  
  # construct scalar_vec, which is a vector with length of nrow(cell_features)
  scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec)) 
  names(scalar1) <- cell_names
  scalar2a <- lineage_future_count_full * scalar1 
  denom_vec <- sapply(uniq_lineages, function(lineage){
    sum(scalar1[cell_lineage_idx_list[[lineage]]]) 
  })
  names(denom_vec) <- uniq_lineages
  scalar2b <- denom_vec[cell_lineage]
  scalar_vec <- scalar1 - scalar2a/scalar2b
  
  # gradient of the intercept
  res1 <- sum(scalar_vec)/num_lineages
  
  # gradient of the other terms
  weighted_features <- sweep(cell_features[,colname_vec,drop=F], 
                             MARGIN = 1, 
                             STATS = scalar_vec, 
                             FUN = "*")
  res2 <- Matrix::colSums(weighted_features)/num_lineages + 2*lambda*coefficient_vec[colname_vec]
  
  res <- c(res1, res2)
  names(res) <- colnames(cell_features)
  
  res
}

.l2norm <- function(x){sqrt(sum(x^2))}