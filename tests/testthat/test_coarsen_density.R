context("Testing coarsen density")

## .add_cutoffs_to_grenander is correct

test_that(".add_cutoffs_to_grenander works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  obj <- estimate_grenander(values = values,
                            weights = weights)
  
  max_val <- floor(max(values))
  bin_cutoff <- seq(0,max_val,by = 1)[-1]
  res <- .add_cutoffs_to_grenander(obj,
                                   bin_cutoff = bin_cutoff)
  
  expect_true(inherits(res, "grenander"))
  expect_true(all(bin_cutoff %in% res$x))
  bool_vec <- sapply(bin_cutoff, function(x){
    val1 <- evaluate_grenander(obj = obj, x = x)
    val2 <- evaluate_grenander(obj = res, x = x)
    abs(val1 - val2) <= 1e-6
  })
  expect_true(all(bool_vec))
  area <- sum(diff(res$x) * res$pdf[-length(res$pdf)])
  expect_true(abs(area - 1) <= 1e-6)
})

######################

## .compute_areas_at_cutoff is correct

test_that(".compute_areas_at_cutoff works", {
  set.seed(10)
  values <- rexp(1000)
  weights <- runif(1000)
  obj <- estimate_grenander(values = values,
                            weights = weights)
  
  max_val <- floor(max(values))
  bin_cutoff <- seq(0,max_val,by = 1)[-1]
  obj <- .add_cutoffs_to_grenander(obj,
                                   bin_cutoff = bin_cutoff)
  
  res <- .compute_areas_at_cutoff(obj = obj,
                                  cutoff = bin_cutoff)
  expect_true(is.numeric(res))
  expect_true(length(res) == length(bin_cutoff)+1)
  expect_true(abs(sum(res)-1) <= 1e-6)
})

test_that(".compute_areas_at_cutoff is correct", {
  x_vec <- seq(0, 5, by = 0.5)
  pdf_vec <- rev(1:(length(x_vec)))
  pdf_vec <- pdf_vec/(sum(diff(x_vec)*pdf_vec[-length(pdf_vec)]))
  obj <- .constructor_grenander(x = x_vec,
                                pdf = pdf_vec,
                                scaling_factor = 1)
  cutoff <- 1:4
  
  res <- .compute_areas_at_cutoff(obj = obj,
                                  cutoff = cutoff)
  
  cutoff_augmented <- c(0, cutoff, max(obj$x)+1)
  area_vec <- diff(obj$x)*obj$pdf[-length(obj$pdf)] # this is the area started from left
  res2 <- rep(NA, 5)
  for(i in 2:length(cutoff_augmented)){
    idx <- intersect(which(obj$x >= cutoff_augmented[i-1]),
                     which(obj$x < cutoff_augmented[i]))
    res2[i-1] <- sum(area_vec[intersect(idx, 1:length(area_vec))])
  }
  
  expect_true(sum(abs(res - res2)) <= 1e-6)
})

