compute_peak_locations <- function(peak_mat){
  res <- sapply(1:nrow(peak_mat), function(i){
    round(mean(peak_mat[i,]))
  })
  names(res) <- paste0("p:", 1:nrow(peak_mat))
  res
}

compute_peak_prior <- function(cutmat,
                               peak_mat,
                               min_prior = 0){
  peak_bp <- as.numeric(colnames(cutmat))
  count_vec <- sapply(1:nrow(peak_mat), function(i){
    idx <- intersect(which(peak_bp >= peak_mat[i,"start"]),
                     which(peak_bp <= peak_mat[i,"end"]))
    sum(sapply(idx, function(j){
      length(.nonzero_col(mat = cutmat,
                          col_idx = j,
                          bool_value = F))
    }))
  })
  if(sum(count_vec) == 0) return(rep(NA, nrow(peak_mat)))
  
  prior_vec <- count_vec/sum(count_vec)
  if(any(prior_vec <= min_prior)){
    if(min_prior*length(prior_vec) > 1) {
      min_prior <- 1/(2*length(prior_vec))
    }
    prior_vec <- prior_vec + min_prior/(1-min_prior*length(prior_vec))
    prior_vec <- prior_vec/sum(prior_vec)
  }
  names(prior_vec) <- paste0("p:", 1:nrow(peak_mat))
  
  stopifnot(all(prior_vec >= 0))
  prior_vec
}
