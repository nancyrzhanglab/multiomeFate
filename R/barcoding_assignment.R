# cells as columns, lineage as rows
barcoding_assignment <- function(lin_mat,
                                 verbose = 0){
  stopifnot(is.matrix(lin_mat))
  library_size <- Matrix::colSums(lin_mat)
    
  n <- ncol(lin_mat)
  nlineages <- nrow(lin_mat)
  
  # to initialize the estimator, find the maximum count for each cell
  if(verbose) print("Starting barcoding assignment")
  cell_max_barcode <- apply(lin_mat, 2, which.max)
  cell_max_barcode_count <- apply(lin_mat, 2, max)

  beta0 <- rep(NA, nlineages)
  beta1 <- rep(NA, nlineages)
  n_won <- rep(NA, nlineages)
  n_nonzero <- rep(NA, nlineages)
  
  for(b in 1:nlineages){
    if(verbose == 1 && b %% 100==0) print(paste0(b," out of ", nlineages," done"))
    if(verbose == 2) print(paste0(b," out of ", nlineages," done"))
    
    # for barcode b, first find all the cells that have its maximum count in barcode b
    won_cells <- which(cell_max_barcode_count == b & lin_mat[b,] > 0)
    nonzero_cells <- which(lin_mat[b,] > 0)
    zero_cells <- which(lin_mat[b,] == 0)
    # cells with nonzero counts that didn't have its maximum in barcode b
    lost_cells <- which(lin_mat[b,] > 0 & cell_max_barcode_count != b) 
    
    n_won[b] <- length(won_cells)
    n_nonzero[b] <- length(nonzero_cells)
    
    # beta1 estimate among maximizing cells
    if(length(won_cells) > 0) beta1[b] <- mean(lin_mat[b,won_cells] / (library_size[won_cells]+1))
    # beta0 estimate among all non-maximizing cells
    beta0[b] <- mean(lin_mat[b,c(lost_cells,zero_cells)] / (library_size[c(lost_cells,zero_cells)]+1))
  }
  
  # global averages
  beta1_mean <- mean(beta1[!is.na(beta1)])
  
  # prevent underflow
  beta0_thresh <- pmax(pmin(beta0, quantile(beta0, 0.98)), quantile(beta0, 0.02))
  
  gamma <- beta1_mean/beta0_thresh # vector of length nlineages
  names(gamma) <- rownames(lin_mat)
  lgamma <- log(gamma)
  Bhat <- matrix(0, ncol = n, nrow = nlineages)
  colnames(Bhat) <- colnames(lin_mat)
  rownames(Bhat) <- rownames(lin_mat)
  
  # we do calculation on the log-scale and then exponentiate
  for(i in 1:n){
    if(verbose > 0 && i %% floor(n/10) == 0) cat('*')
    if(max(lin_mat[,i]) > 10){
      # high counts, avoid overflow.
      bstarc <-  which.max(lin_mat[,i])
      ldeltabc <- lgamma - lgamma[bstarc]
      diff_vec <- lin_mat[,i] - lin_mat[bstarc,i]
      temp <- exp(lin_mat[,i]*ldeltabc + diff_vec*lgamma[bstarc])
      Bhat[,i] <- temp/sum(temp)
      
    } else {
      lgammaX <- lin_mat[,i]*lgamma
      denom <- sum(exp(lgammaX))
      Bhat[,i] <- exp(lgammaX)/denom
    }
  }
  
  list(Bhat = Bhat,
       gamma = gamma)
}

barcode_clustering <- function(lin_mat,
                               cell_lower_limit = 100,
                               check_contradictions = T,
                               cor_threshold = 0.55,
                               verbose = 0){
  stopifnot(inherits(lin_mat, "dgCMatrix"))
  
  lin_mat_t <- Matrix::t(lin_mat)
  nlineages <- nrow(lin_mat)
  num_cells <- sapply(1:nlineages, function(j){
    length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
  })
  
  lin_mat_t <- lin_mat_t[,which(num_cells >= cell_lower_limit)]
  cor_mat <- stats::cor(as.matrix(lin_mat_t))
  cor_mat[lower.tri(cor_mat, diag = T)] <- NA
  # determine all the lineages to merge
  arr_idx <- which(cor_mat >= cor_threshold, arr.ind = T)
  
  lineage_list <- vector("list", length = 0)
  uniq_lineage <- sort(unique(as.numeric(arr_idx)))
  uniq_lineage <- cbind(uniq_lineage, rep(NA, length(uniq_lineage)))
  lineage_name <- rownames(cor_mat)
  
  # determine the clusters of lineages to merge
  for(i in 1:nrow(arr_idx)){
    lineage_idx <- which(uniq_lineage[,1] %in% arr_idx[i,])
    
    if(length(lineage_list) == 0) {
      lineage_list[[1]] <- sort(arr_idx[1,])
      uniq_lineage[lineage_idx,2] <- 1
    } else {
      if(all(is.na(uniq_lineage[lineage_idx,2]))) {
        lineage_list[[length(lineage_list)+1]] <- sort(arr_idx[i,])
        uniq_lineage[lineage_idx,2] <- length(lineage_list)
      } else {
        val <- uniq_lineage[lineage_idx,2]
        val <- val[!is.na(val)]
        val <- unique(val)
        if(check_contradictions & length(val) > 1) {
          stop(paste0("Contradiction among lineages", unique(sort(lineage_list[[val[1]]], lineage_list[[val[2]]], arr_idx[i,]))))
        }
        lineage_list[[val]] <- sort(unique(c(lineage_list[[val]], arr_idx[i,])))
        uniq_lineage[lineage_idx,2] <- val
      }
    }
  }
  
  # rename all the indices with their lineage name
  lineage_list_name <- lapply(lineage_list, function(vec){
    lineage_name[vec]
  })
  
  list(lineage_clusters = lineage_list_name)
}

barcode_combine <- function(lin_mat,
                            lineage_clusters){
  stopifnot(is.matrix(lin_mat))
  nlineages <- nrow(lin_mat)
  
  lineage_included_names <- sort(unique(unlist(lineage_clusters)))
  lineage_included_idx <- which(rownames(lin_mat) %in% lineage_included_names)
  lineage_excluded_idx <- setdiff(1:nlineages, lineage_included_idx)
  
  lin_untounced <- lin_mat[lineage_excluded_idx,]
  len <- length(lineage_clusters)
  lin_list <- lapply(lineage_clusters, function(vec){
    lin_idx <- which(rownames(lin_mat) %in% lineage_names)
    lin_mat[lin_idx,]
  })
  
  for(i in 1:len){
    lin_list[[i]] <- Matrix::colSums(lin_list[[i]])
  }
  
  rbind(lin_untounced, do.call(rbind, lin_list))
}

###################################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}