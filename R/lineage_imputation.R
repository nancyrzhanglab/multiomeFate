lineage_imputation <- function(cell_features,
                               cell_lineage,
                               coefficient_initial,
                               lineage_future_count){
  
  stopifnot(all(sort(unique(cell_lineage)) == 
                  sort(unique(names(lineage_future_count)))),
            is.matrix(cell_features), nrow(cell_features) == length(cell_lineage),
            ncol(cell_features) == length(coefficient_initial))
  
  uniq_lineages <- sort(unique(names(lineage_future_count)))
  cell_lineage_idx_list <- lapply(uniq_lineages, function(lineage){
    which(cell_lineage == lineage)
  })
  names(cell_lineage_idx_list) <- uniq_lineages
}

#################################

.lineage_objective <- function(cell_features,
                               cell_lineage,
                               cell_lineage_idx_list,
                               coefficient_vec,
                               lineage_future_count){
  uniq_lineages <- names(cell_lineage_idx_list)
  cell_names <- rownames(cell_features)
  lineage_future_count_full <- lineage_future_count[cell_lineage]
  
  scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec))
  names(scalar1) <- cell_names
  scalar2 <- sapply(uniq_lineages, function(lineage){
    log(sum(scalar1[cell_lineage_idx_list[[lineage]]]))
  })
  scalar_vec <- scalar1 - lineage_future_count_full*scalar2
  sum(scalar_vec)
}

.lineage_gradient <- function(cell_features,
                              cell_lineage,
                              cell_lineage_idx_list,
                              coefficient_vec,
                              lineage_future_count){
  uniq_lineages <- names(cell_lineage_idx_list)
  cell_names <- rownames(cell_features)
  lineage_future_count_full <- lineage_future_count[cell_lineage]
  
  scalar1 <- exp(coefficient_vec %*% cell_features)
  names(scalar1) <- cell_names
  scalar2a <- lineage_future_count_full * scalar1
  denom_vec <- sapply(uniq_lineages, function(lineage){
    sum(scalar1[cell_lineage_idx_list[[lineage]]])
  })
  names(denom_vec) <- uniq_lineages
  scalar2b <- denom_vec[cell_lineage]
  scalar_vec <- scalar1 - scalar2a/scalar2b
  
  weighted_features <- sweep(cell_features, 
                             MARGIN = 1, 
                             STATS = scalar_vec, 
                             FUN = "*")
  res <- Matrix::colSums(weighted_features)
  names(res) <- colnames(cell_features)
  res
}

.lineage_hessian <- function(cell_features,
                             cell_lineage,
                             cell_lineage_idx_list,
                             coefficient_vec,
                             lineage_future_count){
  p <- ncol(cell_features)
  n <- nrow(cell_features)
  uniq_lineages <- names(cell_lineage_idx_list)
  cell_names <- rownames(cell_features)
  lineage_future_count_full <- lineage_future_count[cell_lineage]
  
  scalar1 <- exp(coefficient_vec %*% cell_features)
  names(scalar1) <- cell_names
  
  scalar2a <- lineage_future_count_full * scalar1
  denom_vec2a <- sapply(uniq_lineages, function(lineage){
    sum(scalar1[cell_lineage_idx_list[[lineage]]])
  })
  names(denom_vec2a) <- uniq_lineages
  scalar2b <- denom_vec2a[cell_lineage]
  
  denom_vec3a <- (denom_vec2a[cell_lineage])^2
  matrix3b <- sweep(cell_features, 
                    MARGIN = 1, 
                    STATS = scalar1, 
                    FUN = "*")
  vector3b_by_lineage <- lapply(uniq_lineages, function(lineage){
    Matrix::colSums(matrix3b[cell_lineage_idx_list[[lineage]],,drop = F])
  })
  rownames(vector3b_by_lineage) <- uniq_lineages
  
  hessian_mat <- matrix(0, nrow = p, ncol = p)
  for(i in 1:n){
    base_mat <- tcrossprod(cell_features[i,,drop = F])
    hessian_mat <- hessian_mat + 
      (scalar_vec1[i] - scalar2a[i]/scalar2b[i])*base_mat + 
      (scalar2a[i]/denom_vec3a[i]) * tcrossprod(vector3b_by_lineage[cell_lineage[i],,drop = F],
                                                cell_features[i,,drop = F])
  }
  
  colnames(hessian_mat) <- ncol(cell_features)
  rownames(hessian_mat) <- ncol(cell_features)
  
  hessian_mat
}