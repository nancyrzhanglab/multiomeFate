compute_entropy <- function(cell_imputation_mat,
                            later_timepoint,
                            seurat_object,
                            variable_celltype,
                            variable_lineage,
                            variable_timepoint,
                            bool_10_power = TRUE,
                            bool_jitter = TRUE,
                            entropy_bump = 0.01,
                            min_imputation = 0.01,
                            min_jitter = 0.1){
  stopifnot(
    length(rownames(cell_imputation_mat)) > 0,
    length(colnames(cell_imputation_mat)) == ncol(cell_imputation_mat)
  )
  
  k <- ncol(cell_imputation_mat)
  if(bool_10_power) {
    cell_imputation_mat <- 10^cell_imputation_mat
  }
  
  cellsize <- rowSums(cell_imputation_mat)
  n <- nrow(cell_imputation_mat)
  
  for(i in 1:n){
    tmp <- cell_imputation_mat[i,]
    if(sum(tmp) <= min_imputation){
      cell_imputation_mat[i,] <- NA
    } else {
      cell_imputation_mat[i,] <- tmp/sum(tmp)
    }
  }
  
  idx <- unique(unlist(apply(cell_imputation_mat, 2, function(x){
    which(is.na(x))
  })))
  if(length(idx) > 0) {
    cell_imputation_mat <- cell_imputation_mat[-idx,,drop = FALSE]
    cellsize <- cellsize[-idx]
  }
  
  if(bool_jitter) {
    n <- nrow(cell_imputation_mat)
    for(i in 1:n){
      if(any(is.na(cell_imputation_mat[i,]))) next()
      cell_imputation_mat[i,] <- cell_imputation_mat[i,] + stats::runif(k, min = 0, max = min_jitter)
      cell_imputation_mat[i,] <- cell_imputation_mat[i,]/sum(cell_imputation_mat[i,])
    }
  }
  
  cell_imputation_mat <- cell_imputation_mat[which(!is.na(cell_imputation_mat[,1])),]
  
  df <- as.data.frame(cell_imputation_mat)
  metadata <- seurat_object@meta.data
  df$celltype <- metadata[rownames(cell_imputation_mat), variable_celltype]
  df$celltype <- factor(df$celltype)
  df$cellsize <- cellsize[rownames(cell_imputation_mat)]
  
  # compute the dominant fate
  df$lineage <- metadata[rownames(df),variable_lineage]
  df$dominant_fate <- rep(NA, nrow(df))
  df$entropy <- rep(NA, nrow(df))
  
  for(lineage in unique(df$lineage)){
    df_idx <- which(df$lineage == lineage)
    seurat_idx <- intersect(which(metadata[,variable_timepoint] == later_timepoint),
                            which(metadata[,variable_lineage] == lineage))
    tab_vec <- table(metadata[seurat_idx, variable_celltype])
    if(length(tab_vec) == 0) next()
    df$dominant_fate[df_idx] <- names(tab_vec)[which.max(tab_vec)]
    df$entropy[df_idx] <- .shannon_entropy(tab_vec)
  }
  df$dominant_fate[which(is.na(df$dominant_fate))] <- "NA"
  df$dominant_fate <- factor(df$dominant_fate)
  
  df$entropy <- df$entropy + entropy_bump
  
  df
}

.shannon_entropy <- function(x){
  if(length(x) == 1) return(0)
  x <- x/sum(x)
  y <- log2(x)
  y[x == 0] <- 0
  -sum(x*y)
}