# cells as columns, lineage as rows
barcoding_posterior <- function(lin_mat,
                                bool_force_rebase = F,
                                tol = 1e-8,
                                verbose = 0){
  stopifnot(is.matrix(lin_mat))
  library_size <- Matrix::colSums(lin_mat)
  
  n <- ncol(lin_mat)
  nlineages <- nrow(lin_mat)
  
  # to initialize the estimator, find the maximum count for each cell
  if(verbose) print("Starting barcoding assignment")
  cell_max_barcode_count <- apply(lin_mat, 2, max)
  lineage_maximizing <- lapply(1:nlineages, function(b){
    which(lin_mat[b,] == cell_max_barcode_count & lin_mat[b,] > 0)
  })
  lin_num_winner <- sapply(lineage_maximizing, length)
  # quantile(lin_num_winner)
  
  beta0 <- rep(NA, nlineages)
  names(beta0) <- names(lin_mat)
  beta1 <- rep(NA, nlineages)
  names(beta1) <- names(lin_mat)
  
  for(b in 1:nlineages){
    if(verbose == 1 && b %% 100==0) print(paste0(b," out of ", nlineages," done"))
    if(verbose == 2) print(paste0(b," out of ", nlineages," done"))
    
    # for barcode b, first find all the cells that have its maximum count in barcode b
    won_cells <- lineage_maximizing[[b]]
    other_cells <- setdiff(1:n, lineage_maximizing[[b]])
    
    # beta1 estimate among maximizing cells
    if(length(won_cells) > 0) beta1[b] <- mean(lin_mat[b,won_cells] / (library_size[won_cells]+1))
    # beta0 estimate among all non-maximizing cells
    beta0[b] <- mean(lin_mat[b,c(other_cells)] / (library_size[c(other_cells)]+1))
  }
  
  # global averages
  beta1_mean <- mean(beta1, na.rm = T)
  
  # prevent underflow
  beta0_thresh <- pmax(pmin(beta0, quantile(beta0[beta0 > tol], 0.98)), max(quantile(beta0[beta0 > tol], 0.02)))
  
  gamma <- beta1_mean/beta0_thresh # vector of length nlineages
  names(gamma) <- rownames(lin_mat)
  
  posterior_mat <- .multinomial_posterior(bool_force_rebase = bool_force_rebase,
                                          gamma = gamma,
                                          lin_mat = lin_mat)
  
  list(beta0 = beta0,
       beta1 = beta1,
       beta1_mean = beta1_mean,
       posterior_mat = posterior_mat,
       gamma = gamma,
       lineage_num_winner = lin_num_winner)
}

barcode_clustering <- function(lin_mat,
                               cell_lower_limit = 100,
                               cor_threshold = 0.55,
                               warn_merging = TRUE,
                               verbose = 0){
  stopifnot(inherits(lin_mat, "dgCMatrix"))
  stopifnot(length(rownames(lin_mat)) > 0)
  
  lin_mat_t <- Matrix::t(lin_mat)
  nlineages <- nrow(lin_mat)
  num_cells <- sapply(1:nlineages, function(j){
    length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
  })
  
  if(verbose > 0) print("Compute correlation matrix")
  # discard any lineages in question (for the purposes of merging) that are too small
  lin_mat_t <- lin_mat_t[,which(num_cells >= cell_lower_limit)]
  # compute a correlation matrix (# rows/columsn = number of lineages)
  cor_mat <- .custom_correlation(lin_mat_t)
  cor_mat[lower.tri(cor_mat, diag = T)] <- NA
  # determine all the lineages to merge. arr_idx is a 2-column matrix
  arr_idx <- which(cor_mat >= cor_threshold, arr.ind = TRUE)
  if(verbose > 2){
    print(paste0("There are ", nrow(arr_idx), " number of highly correlated lineages to resolve."))
  }
  
  # tabulate uniq_lineage, which is going to keep track of which lineage is 
  #  assigned to which "cluster"
  lineage_list <- vector("list", length = 0)
  uniq_lineage <- sort(unique(as.numeric(arr_idx)))
  uniq_lineage <- cbind(uniq_lineage, rep(NA, length(uniq_lineage)))
  lineage_name <- rownames(cor_mat)
  
  if(verbose > 0) print("Determining how to merge lineages")
  # determine the clusters of lineages to merge
  for(i in 1:nrow(arr_idx)){
    if(verbose > 1 && nrow(arr_idx) > 10 && i %% floor(nrow(arr_idx)/10) == 0) cat('*')
    
    # find the rows in the uniq_lineage table on which we're currently working on
    # this represents a pair of lineages
    lineage_idx <- which(uniq_lineage[,1] %in% arr_idx[i,])
    
    if(length(lineage_list) == 0) {
      # if we have not yet merged any lineages, then this is straight-forward
      # create a new cluster
      lineage_list[[1]] <- sort(arr_idx[1,])
      uniq_lineage[lineage_idx,2] <- 1
      
    } else {
      # otherwise...
      
      if(all(is.na(uniq_lineage[lineage_idx,2]))) {
        # if this is a completely new cluster (i.e., all unassigned), also pretty straight-forward
        # create a new cluster
        
        lineage_list[[length(lineage_list)+1]] <- sort(arr_idx[i,])
        uniq_lineage[lineage_idx,2] <- length(lineage_list)
        if(verbose > 2) print(paste0("There are currently ", length(lineage_list), " cluster of lineages"))
        
      } else {
        # the difficulty is this step. The new lineage in question has a high
        #  correlation with an existing cluster of lineages
        
        # first find all the lineages in this existing cluster
        val <- uniq_lineage[lineage_idx,2]
        val <- val[!is.na(val)]
        val <- unique(val)
        if(length(val) > 1) {
          if(warn_merging) {stop("Merging happening")}
          ## THIS CODE ISN'T MADE YET. KL: I suspect implementing this as part of an igraph is easier
          # This code would be designed to account of the scenario that the two lineages
          #  (which are themselves correlated) are each correlated with two separate cluster of lineages
        } else {
          # just add this lineage to the existing cluster
          lineage_list[[val]] <- sort(unique(c(lineage_list[[val]], arr_idx[i,])))
          uniq_lineage[lineage_idx,2] <- val
        }
        
      }
    }
  }
  
  # rename all the indices with their lineage name
  lineage_list_name <- lapply(lineage_list, function(vec){
    lineage_name[vec]
  })
  
  list(arr_idx = arr_idx,
       lineage_clusters = lineage_list_name,
       uniq_lineage = uniq_lineage)
}

barcode_combine <- function(lin_mat,
                            lineage_clusters,
                            verbose = 0){
  stopifnot(is.matrix(lin_mat))
  nlineages <- nrow(lin_mat)
  
  print("Starting combination")
  lineage_included_names <- sort(unique(unlist(lineage_clusters)))
  lineage_included_idx <- which(rownames(lin_mat) %in% lineage_included_names)
  lineage_excluded_idx <- setdiff(1:nlineages, lineage_included_idx)
  
  print("Extracting unaffected lineages")
  lin_untounced <- lin_mat[lineage_excluded_idx,]
  len <- length(lineage_clusters)
  lin_list <- lapply(lineage_clusters, function(vec){
    lin_idx <- which(rownames(lin_mat) %in% vec)
    lin_mat[lin_idx,]
  })
  
  print("Adding lineages to be combined")
  for(i in 1:len){
    if(verbose >0 && len > 10 && i %% floor(len/10) == 0) cat('*')
    vec <- rownames(lin_list[[i]])
    colname_vec <- colnames(lin_list[[i]])
    lin_list[[i]] <- matrix(Matrix::colSums(lin_list[[i]]), ncol = ncol(lin_list[[i]]), nrow = 1)
    rownames(lin_list[[i]]) <- vec[1]
    colnames(lin_list[[i]]) <- colname_vec
  }
  
  print("Formatting final matrix")
  rbind(lin_untounced, do.call(rbind, lin_list))
}

barcoding_assignment <- function(posterior_mat,
                                 difference_val = 0.2,
                                 verbose = 0){
  n <- ncol(posterior_mat)
  
  lineage_names <- rownames(posterior_mat)
  lineage_idx <- apply(posterior_mat, 2, which.max)
  if(verbose) print("Computing difference between maximizing and second-maximizing")
  difference_vec <- apply(posterior_mat, 2, function(x){
    abs(diff(sort(x, decreasing = T)[1:2]))
  })
  
  assignment_vec <- sapply(1:n, function(i){
    if(verbose && n > 10 && i %% floor(n/10) == 0) cat('*')
    if(difference_vec[i] >= difference_val){
      return(lineage_names[lineage_idx[i]]) 
    } else {
      return(NA)
    }
  })
  
  names(assignment_vec) <- colnames(posterior_mat)
  assignment_vec
}

###################################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
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

.multinomial_posterior <- function(bool_force_rebase,
                                   gamma,
                                   lin_mat,
                                   verbose = 0){
  n <- ncol(lin_mat)
  nlineages <- nrow(lin_mat)
  
  lgamma <- log(gamma)
  Bhat <- matrix(0, ncol = n, nrow = nlineages)
  colnames(Bhat) <- colnames(lin_mat)
  rownames(Bhat) <- rownames(lin_mat)
  
  # we do calculation on the log-scale and then exponentiate
  for(i in 1:n){
    if(verbose > 0 && i %% floor(n/10) == 0) cat('*')
    Bhat[,i] <- .multinomial_posterior_vector(
      bool_force_rebase = bool_force_rebase,
      lgamma = lgamma,
      lin_count = lin_mat[,i]
    )
  }
  
  Bhat
}

.multinomial_posterior_vector <- function(bool_force_rebase,
                                          lgamma,
                                          lin_count){
  lgammaX <- lin_count*lgamma
  
  if(bool_force_rebase || max(lgammaX) > 10){
    # high counts, avoid overflow
    max_val <- max(lgammaX)
    tmp <- exp(lgammaX-max_val)
    vec <- tmp/sum(tmp)
    
  } else {
    denom <- sum(exp(lgammaX))
    vec <- exp(lgammaX)/denom
  }
  
  vec
}

# from https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
.custom_correlation <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat)) 
  cormat <- covmat/tcrossprod(sdvec)
  cormat
}