context("Testing peak modeling")

## .compute_bin_matrix is correct

test_that(".compute_bin_matrix works", {
  # hypothetical genome from base-pairs 1001 to 2000
  # 120 fragments, one per cell
  # 2 peaks (at 1500 and at 1900), most fragments at 1500, 1900 or 1100
  
  set.seed(10)
  bin_midpoints <- c(-100, -50, 0, 50, 100)
  names(bin_midpoints) <- paste0("bin:", seq(-2,2,by=1))
  bin_limits <- c(-200,200)
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
  
  res <- .compute_bin_matrix(bin_limits = bin_limits,
                             bin_midpoints = bin_midpoints,
                             bool_lock_within_peak = T,
                             cutmat = cutmat,
                             peak_locations = peak_locations)
  
  for(i in 1:40) expect_true(res[i,1] == 2, is.na(res[i,2]))
  for(i in 41:80) expect_true(res[i,1] == 0, is.na(res[i,2]))
  for(i in 81:120) expect_true(is.na(res[i,1]), res[i,2] == 0)
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
  bin_mat <- .compute_bin_matrix(bin_limits = bin_limits,
                                 bin_midpoints = bin_midpoints,
                                 bool_lock_within_peak = T,
                                 cutmat = cutmat_dying,
                                 peak_locations = peak_locations)
  
  res <- .initialize_theta(bin_mat = bin_mat, num_bins = 5)
  expect_true(sum(abs(res - c(0,0,2/3,0,1/3))) <= 1e-6)
})

##########################################
# load("tests/assets/test.RData")
## .compute_bin_matrix is correct
test_that(".compute_bin_matrix works", {
  load("assets/test.RData")
  res <- .compute_bin_matrix(bin_limits = bin_limits,
                             bin_midpoints = bin_midpoints,
                             bool_lock_within_peak = F,
                             cutmat = cutmat_dying,
                             peak_locations = peak_locations)
  
  expect_true(all(dim(res) == c(length(cutmat_dying@x), length(peak_locations))))
})


test_that(".compute_bin_matrix works", {
  load("assets/test.RData")
  bin_mat <- .compute_bin_matrix(bin_limits = bin_limits,
                                 bin_midpoints = bin_midpoints,
                                 bool_lock_within_peak = T,
                                 cutmat = cutmat_dying,
                                 peak_locations = peak_locations)
  
  bool_vec <- sapply(1:nrow(bin_mat), function(i){
    idx <- which(bin_mat[i,] == 0)
    if(length(idx) > 0){
      if(length(idx) == 1) return(TRUE) else return(FALSE)
    } else{
      return(TRUE)
    }
  })
  
  expect_true(all(bool_vec))
})


####################################

# load("tests/assets/test.RData")
## .initialize_theta is correct
test_that(".initialize_theta works", {
  load("assets/test.RData")
  bin_mat <- .compute_bin_matrix(bin_limits = bin_limits,
                                 bin_midpoints = bin_midpoints,
                                 bool_lock_within_peak = T,
                                 cutmat = cutmat_dying,
                                 peak_locations = peak_locations)
  
  res <- .initialize_theta(bin_mat = bin_mat,
                           num_bins = length(bin_midpoints))
  
  expect_true(abs(sum(res)-1) <= 1e-5)
})

####################################

# load("tests/assets/test.RData")
## .compute_loglikelihood is correct
test_that(".compute_loglikelihood works", {
  load("assets/test.RData")
  bin_mat <- .compute_bin_matrix(bin_limits = bin_limits,
                                 bin_midpoints = bin_midpoints,
                                 bool_lock_within_peak = T,
                                 cutmat = cutmat_dying,
                                 peak_locations = peak_locations)
  theta_vec <- .initialize_theta(bin_mat = bin_mat,
                                 num_bins = length(bin_midpoints))
  res <- .compute_loglikelihood(bin_mat = bin_mat,
                                prior_vec = peak_prior,
                                theta_vec = theta_vec,
                                tol = 1e-6)
  
  expect_true(is.numeric(res))
})

####################################

# load("tests/assets/test.RData")
## .e_step is correct
test_that(".e_step works", {
  load("assets/test.RData")
  bin_mat <- .compute_bin_matrix(bin_limits = bin_limits,
                                 bin_midpoints = bin_midpoints,
                                 bool_lock_within_peak = T,
                                 cutmat = cutmat_dying,
                                 peak_locations = peak_locations)
  theta_vec <- .initialize_theta(bin_mat = bin_mat,
                                 num_bins = length(bin_midpoints))
  res <- .e_step(bin_mat = bin_mat,
                 prior_vec = peak_prior,
                 theta_vec = theta_vec)
  
  expect_true(all(abs(Matrix::rowSums(res)-1)<=1e-6))
})

####################################

# load("tests/assets/test.RData")
## peak_mixture_modeling is correct

test_that("peak_mixture_modeling works", {
  load("assets/test.RData")
  res1 <- peak_mixture_modeling(bin_limits = bin_limits,
                               bin_midpoints = bin_midpoints, 
                               cutmat = cutmat_dying, 
                               peak_locations = peak_locations,
                               peak_prior = peak_prior,
                               bool_freeze_prior = F,
                               verbose = 0)
  expect_true(inherits(res1, "peakDistribution"))
  # round(res1$theta_vec,2)
  
  res2 <- peak_mixture_modeling(bin_limits = bin_limits,
                               bin_midpoints = bin_midpoints, 
                               cutmat = cutmat_winning, 
                               peak_locations = peak_locations,
                               peak_prior = peak_prior,
                               bool_freeze_prior = T,
                               verbose = 0)
  expect_true(inherits(res2, "peakDistribution"))
  # round(res2$theta_vec,2)
  
  res3 <- peak_mixture_modeling(bin_limits = bin_limits,
                                bin_midpoints = bin_midpoints, 
                                cutmat = rbind(cutmat_dying, cutmat_winning), 
                                peak_locations = peak_locations,
                                peak_prior = peak_prior,
                                bool_freeze_prior = F,
                                verbose = 0)
  expect_true(inherits(res3, "peakDistribution"))
  # round(res3$theta_vec,2)
  
})

