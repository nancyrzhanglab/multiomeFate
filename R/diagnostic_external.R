diagnos_neigh <- function(chrom_fit, dat){
  stopifnot(class(chrom_fit) == "chromatin_potential", class(dat) == "mf_simul",
            !all(is.na(chrom_fit$list_diagnos[[1]]$recruit$postprocess)))

  vec_pseudotime <- dat$df_info$time
  
  res <- lapply(1:length(chrom_fit$list_diagnos), function(i){
    
    obj <- chrom_fit$list_diagnos[[i]]$recruit$postprocess$lis_nn

    mat_diag <- sapply(1:length(obj), function(j){
      idx <- as.numeric(names(obj)[j])
      cur_time <- vec_pseudotime[idx]
      
      forward_num <- length(which(vec_pseudotime[obj[[j]]] >= cur_time))
      backward_num <- length(which(vec_pseudotime[obj[[j]]] < cur_time))
      
      c(idx = idx, forward_num = forward_num, backward_num = backward_num)
    })
    
    tmp <- as.data.frame(t(mat_diag))
    tmp$selected <- chrom_fit$list_diagnos[[i]]$recruit$postprocess$df_diag$selected
    tmp
  })
  
  names(res) <- 1:length(chrom_fit$list_diagnos)
  
  structure(res, class = "diagnos_neigh")
}