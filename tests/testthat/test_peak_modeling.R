context("Testing peak modeling")

## .compute_bin_matrix is correct

test_that(".compute_bin_matrix works", {
  # hypothetical genome from base-pairs 1001 to 2000
  # 120 fragments, one per cell
  # 2 peaks (at 1500 and at 1900), most fragments at 1500, 1900 or 1100
  
  set.seed(10)
  bin_midpoints <- c(-100, -50, 0, 50, 100)
  names(bin_midpoints) <- paste0("bin:", seq(-2,2,by=1))
  peak_locations <- c(1500, 1900)
  names(peak_locations) <- paste0("p:", 1:2)
  num_frags <- 120
  frag_location <- rep(c(1100,1500,1900), each = num_frags/3)
  cutmat <- sapply(frag_location, function(x){
    vec <- rep(0, length = 1000)
    names(vec) <- c(1001:2000)
    vec[as.character(x)] <- 1
    vec
  })
  cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
  rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
  
  res <- .compute_bin_matrix(bin_midpoints,
                             cutmat,
                             peak_locations)
  
  for(i in 1:40) expect_true(all(res[i,] == 5))
  for(i in 41:80) expect_true(all(res[i,] == c(3,5)))
  for(i in 81:120) expect_true(all(res[i,] == c(1,3)))
})

##########################################

## .initialize_theta is correct

test_that(".initialize_theta works", {
  # hypothetical genome from base-pairs 1001 to 2000
  # 120 fragments, one per cell
  # 2 peaks (at 1500 and at 1900), most fragments at 1500, 1900 or 1100
  
  set.seed(10)
  bin_midpoints <- c(-100, -50, 0, 50, 100)
  names(bin_midpoints) <- paste0("bin:", seq(-2,2,by=1))
  peak_locations <- c(1500, 1900)
  names(peak_locations) <- paste0("p:", 1:2)
  num_frags <- 120
  frag_location <- rep(c(1100,1500,1900), each = num_frags/3)
  cutmat <- sapply(frag_location, function(x){
    vec <- rep(0, length = 1000)
    names(vec) <- c(1001:2000)
    vec[as.character(x)] <- 1
    vec
  })
  cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
  rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
  bin_mat <- .compute_bin_matrix(bin_midpoints,
                                 cutmat,
                                 peak_locations)
  
  res <- .initialize_theta(bin_mat = bin_mat, num_bins = 5)
  expect_true(sum(abs(res - c(0,0,2/3,0,1/3))) <= 1e-6)
})
