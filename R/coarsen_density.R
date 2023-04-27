coarsen_density <- function(obj,
                            bin_cutoff){
  stopifnot(inherits(obj, "grenander"),
            all(bin_cutoff >= min(obj$x)),
            all(bin_cutoff <= max(obj$x)),
            all(sort(bin_cutoff) == bin_cutoff))
  
  obj_new <- .add_cutoffs_to_grenander(obj = obj,
                                       bin_cutoff = bin_cutoff)
  
  area_vec_coarsen <- .compute_areas_at_cutoff(obj = obj_new, 
                                               cutoff = bin_cutoff)
  x_vec <- c(0,bin_cutoff, max(obj_new$x))
  pdf_vec <- area_vec_coarsen/diff(x_vec)
  
  .constructor_grenander(x = x_vec,
                         pdf = c(pdf_vec,0),
                         scaling_factor = obj$scaling_factor)
}

#############

.add_cutoffs_to_grenander <- function(obj,
                                      bin_cutoff){
  # add the cutoffs into obj
  x_vec <- obj$x
  pdf_vec <- obj$pdf
  
  for(val in bin_cutoff){
    idx <- which(x_vec <= val)
    stopifnot(length(idx) > 0)
    x_vec <- c(x_vec[idx], val, x_vec[-idx])
    pdf_vec <- c(pdf_vec[idx], pdf_vec[idx[length(idx)]], pdf_vec[-idx])
  }
  
  res <- .remove_duplicates(associated_vec = pdf_vec,
                            target_vec = x_vec)
  
  .constructor_grenander(x = res$target_vec,
                         pdf = res$associated_vec,
                         scaling_factor = obj$scaling_factor)
}

.compute_areas_at_cutoff <- function(obj, 
                                     cutoff,
                                     tol = 1e-6){
  area_vec <- diff(obj$x) * obj$pdf[-length(obj$pdf)]
  stopifnot(abs(sum(area_vec) - 1) <= tol)
  cutoff_idx <- sapply(cutoff, function(x){
    idx <- which(abs(obj$x - x) <= 1e-6)
    stopifnot(length(idx) == 1)
    idx
  })
  stopifnot(all(sort(cutoff_idx) == cutoff_idx), all(cutoff_idx > 1))
  cutoff_idx <- c(0, cutoff_idx, length(area_vec)+1)
  
  area_vec_coarsen <- sapply(2:length(cutoff_idx), function(i){
    sum(area_vec[max(cutoff_idx[i-1],1):(cutoff_idx[i]-1)])
  })
  stopifnot(abs(sum(area_vec_coarsen) - 1) <= tol)
  
  area_vec_coarsen
}