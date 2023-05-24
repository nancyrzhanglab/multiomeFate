preprocess_chromatin_peak <- function(peak_mapping_list,
                                      seurat_obj,
                                      cell_names = NULL,
                                      slot_atac = "ATAC",
                                      slot_rna = "RNA",
                                      tol = 1e-6,
                                      verbose = 0){
  if(all(is.null(cell_names))) cell_names <- colnames(seurat_obj)
  stopifnot(all(cell_names %in% colnames(seurat_obj)))

  # extract RNA matrix
  rna_mat <- Matrix::t(all_data2[[slot_rna]]@scale.data[,cell_names])
  sd_vec <- sparseMatrixStats::colSds(rna_mat)
  if(any(sd_vec <= tol)){
    rna_mat <- rna_mat[,sd_vec >= tol,drop = F]
  }
  
  # organize the set of genes 
  tmp <- sapply(peak_mapping_list, length)
  if(any(tmp == 0)) peak_mapping_list <- peak_mapping_list[-which(tmp == 0)]
  gene_vec <- sort(unique(intersect(names(peak_mapping_list), colnames(rna_mat))))
  rna_mat <- rna_mat[,gene_vec,drop = F]
  peak_mapping_list <- peak_mapping_list[gene_vec]
  len <- length(gene_vec)
  
  # construct the relevant chromatin-activity matrix by extract relevant rows from seurat_obj
  chr_peak_list <- sapply(peak_mapping_list, function(lis){
    stopifnot("overlap_idx" %in% names(lis))
    peak_idx <- lis$overlap_idx
    stopifnot(length(peak_idx) > 0)
    Matrix::t(seurat_obj[[slot_atac]]@counts[peak_idx,,drop = F])
  })
  names(chr_peak_list) <- gene_vec
  
  chract_mat <- sapply(1:len, function(i){
    sparseMatrixStats::rowSum2(chr_peak_list[[i]])
  })
  rownames(chract_mat) <- colnames(seurat_obj)
  colnames(chract_mat) <- gene_vec
  
  # remove any genes with no peaks at all
  tmp <- which(Matrix::colSums(chract_mat) <= tol)
  if(length(tmp) > 0){
    gene_vec <- gene_vec[-tmp]
    
    rna_mat <- rna_mat[,gene_vec,drop = F]
    peak_mapping_list <- peak_mapping_list[gene_vec]
    chr_peak_list <- chr_peak_list[gene_vec]
    chract_mat <- chract_mat[,gene_vec,drop = F]
    len <- length(gene_vec)
  }
  stopifnot(all(colnames(rna_mat) == names(peak_mapping_list)),
            all(colnames(rna_mat) == names(chr_peak_list)),
            all(colnames(rna_mat) == colnames(chract_mat)),
            ncol(rna_mat) == len)
  
  # now normalize the chromatin activit matrix, as it acts like a template on how to normalize the individual peaks
  lib_mat <- Matrix::rowSums(chract_mat)
  lib_mat <- pmax(lib_mat, 1)
  chract_mat_normalized <- .mult_vec_mat(1/lib_mat, chract_mat)
  chract_mat_normalized2 <- log1p(chract_mat_normalized*1e6)
  chract_mat_normalized3 <- scale(chract_mat_normalized2)
  
  set.seed(10)
  dimred_res <- irlba::irlba(chract_mat_normalized3, nv = 50)
  eigenbasis <- dimred_res$u[,2:50]

  # now actually preprocess the individual peaks associated with gene
  for(i in 1:len){
    if(verbose > 1) print(paste0("Working on gene ", i, " out of ", len))
    gene <- names(gene_vec)[i]
    mat <- chr_peak_list[[i]]
    
    # apply library size
    mat <- .mult_vec_mat(1/lib_mat, mat)
    
    # rescale based on the log-version
    # do log1p(chromatin-activity)*(fraction of how much this peak contributed to the chromatin activity)
    mat2 <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
    base <- chract_mat_normalized2[,gene]
    denom <- chract_mat_normalized[,gene]
    denom[base == 0] <- 1
    for(j in 1:ncol(mat2)){
      mat2[,j] <- base * mat[,j]/denom
    }
    
    # recenter the entire matrix based on the global mean and standard deviation
    mean_val <- mean(chract_mat_normalized2[,gene])
    sd_val <- sd(chract_mat_normalized2[,gene])
    mat3 <- (mat2 - mean_val)/sd_val
    
    # smooth the entire matrix based on the eigenbasis computed across all the chromatin activities
    # this is n-by-n projection matrix
    mat4 <- eigenbasis %*% crossprod(eigenbasis, mat3)
    rownames(mat4) <- rownames(chr_peak_list[[i]])
    colnames(mat4) <- colnames(chr_peak_list[[i]])
    
    chr_peak_list[[i]] <- mat4[cell_names,,drop = F]
  }
  
  list(chr_peak_list = chr_peak_list,
       rna_mat = rna_mat)
}

