.construct_lineage_data <- function(){
  set.seed(10)
  L <- 10
  n_each <- 5
  variance_across_lineage <- 1
  variance_within_lineage <- 0.3
  p <- 2
  
  # construct feature matrix
  tmp <- lapply(1:L, function(lineage){
    center <- MASS::mvrnorm(n = 1, mu = rep(1,p), Sigma = variance_across_lineage*diag(p))
    mat <- MASS::mvrnorm(n = n_each, mu = center, Sigma = variance_within_lineage*diag(p))
    rownames(mat) <- paste0("c:", (lineage-1)*n_each+1:nrow(mat))
    mat
  })
  cell_features <- do.call(rbind, tmp)
  
  # construct lineage
  uniq_lineages <- paste0("lin:", 1:L)
  cell_lineage <- rep(uniq_lineages, each = n_each)
  names(cell_lineage) <- rownames(cell_features)
  
  # construct future lineage counts
  coefficient_vec <- c(1,1)
  lineage_future_count <- sapply(uniq_lineages, function(lineage){
    idx <- which(cell_lineage == lineage)
    lambda <- sum(sapply(idx, function(i){
      exp(coefficient_vec %*% cell_features[i,])
    }))
    stats::rpois(1, lambda = lambda)
  })
  
  # name things
  colnames(cell_features) <- paste0("p:", 1:p)
  names(coefficient_vec) <- colnames(cell_features)
  names(lineage_future_count) <- uniq_lineages
  
  list(cell_features = cell_features,
       cell_lineage = cell_lineage,
       coefficient_vec = coefficient_vec,
       lineage_future_count = lineage_future_count)
}
