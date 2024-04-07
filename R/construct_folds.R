construct_folds <- function(cell_lineage,
                            tab_mat,
                            future_timepoint,
                            num_folds = 10){
  stopifnot(future_timepoint %in% colnames(tab_mat))
  
  lineages_ordered <- rownames(tab_mat)[order(tab_mat[,future_timepoint], decreasing = T)]
  num_lineages <- length(lineages_ordered)
  num_per_fold <- ceiling(length(lineages_ordered)/num_folds)
  
  for(i in 1:num_per_fold){
    idx_vec <- ((i-1)*num_per_fold+1):min(i*num_per_fold, num_lineages)
    lineages_ordered[sample(idx_vec)] <- lineages_ordered[idx_vec]
  }
  fold_lineage_list <- lapply(1:num_folds, function(i){
    lineages_ordered[unique(pmin(i + (0:num_per_fold)*num_folds, num_lineages))]
  })
  names(fold_lineage_list) <- paste0("fold:", 1:num_folds)
  
  cv_cell_list <- lapply(fold_lineage_list, function(lineages){
    unlist(lapply(lineages, function(lineage){
      which(cell_lineage == lineage)
    }))
  })
  names(cv_cell_list) <- names(fold_lineage_list)
  
  list(cv_cell_list = cv_cell_list,
       fold_lineage_list = fold_lineage_list)
}