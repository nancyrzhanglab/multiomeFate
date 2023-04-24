context("Testing peak modeling")

## .compute_bin_matrix is correct

test_that(".compute_frag_peak_matrix works", {
  # hypothetical genome from base-pairs 1001 to 2000
  # 120 fragments, one per cell
  # 2 peaks (at 1500 and at 1900), most fragments at 1300, 1800 or 1100
  set.seed(10)
  peak_locations <- c(1500, 1900)
  names(peak_locations) <- paste0("p:", 1:2)
  num_frags <- 120
  frag_location <- rep(c(1100,1300,1800), each = num_frags/3)
  cutmat <- sapply(frag_location, function(x){
    vec <- rep(0, length = 1000)
    names(vec) <- c(1001:2000)
    vec[as.character(x)] <- 1
    vec
  })
  cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
  rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
  
  res <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                   cutmat = cutmat,
                                   fragment_locations = NULL,
                                   num_peak_limit = 3,
                                   peak_locations = peak_locations,
                                   peak_width = 250)
  
  expect_true(inherits(res, "dgCMatrix"))
  for(i in 1:40) expect_true(all(res[i,] == c(400,800)))
  for(i in 41:80) expect_true(all(res[i,] == c(200,600)))
  for(i in 81:120) expect_true(all(res[i,] == c(0,100)))
})

test_that(".compute_frag_peak_matrix respects num_peak_limit", {
  set.seed(10)
  peak_locations <- c(1500, 1900)
  names(peak_locations) <- paste0("p:", 1:2)
  num_frags <- 120
  frag_location <- rep(c(1100,1300,1800), each = num_frags/3)
  cutmat <- sapply(frag_location, function(x){
    vec <- rep(0, length = 1000)
    names(vec) <- c(1001:2000)
    vec[as.character(x)] <- 1
    vec
  })
  cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
  rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
  
  res <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                   cutmat = cutmat,
                                   fragment_locations = NULL,
                                   num_peak_limit = 1,
                                   peak_locations = peak_locations,
                                   peak_width = 250)
  
  expect_true(inherits(res, "dgCMatrix"))
  for(i in 1:40) expect_true(all(res[i,] == c(400,0)))
  for(i in 41:80) expect_true(all(res[i,] == c(200,0)))
  for(i in 81:120) expect_true(all(res[i,] == c(0,100)))
})

# load("tests/assets/test.RData")
test_that(".compute_frag_peak_matrix works on a real cutmatrix", {
  load("../assets/test.RData")
  res <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                   cutmat = cutmat_dying,
                                   fragment_locations = NULL,
                                   num_peak_limit = 3,
                                   peak_locations = peak_locations,
                                   peak_width = peak_width)
  
  expect_true(all(dim(res) == c(length(cutmat_dying@x), length(peak_locations))))
  
  bool_vec <- sapply(1:nrow(res), function(i){
    idx <- intersect(which(res[i,] <= peak_width/2), which(res[i,] != 0))
    if(length(idx) > 0){
      if(length(idx) == 1) return(TRUE) else return(FALSE)
    } else{
      return(TRUE)
    }
  })
  
  expect_true(all(bool_vec))
})

##########################################

## .initialize_grenander is correct

test_that(".initialize_grenander works", {
  # hypothetical genome from base-pairs 1001 to 2000
  # 120 fragments, one per cell
  # 2 peaks (at 1500 and at 1900), most fragments at 1300, 1800 or 1100
  set.seed(10)
  peak_locations <- c(1500, 1900)
  names(peak_locations) <- paste0("p:", 1:2)
  num_frags <- 120
  frag_location <- rep(c(1100,1300,1800), each = num_frags/3)
  cutmat <- sapply(frag_location, function(x){
    vec <- rep(0, length = 1000)
    names(vec) <- c(1001:2000)
    vec[as.character(x)] <- 1
    vec
  })
  cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
  rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
  dist_mat <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                        cutmat = cutmat,
                                        fragment_locations = NULL,
                                        num_peak_limit = 3,
                                        peak_locations = peak_locations,
                                        peak_width = 100)
  
  res <- .initialize_grenander(dist_mat = dist_mat,
                               scaling_factor = 1)
  expect_true(inherits(res, "grenander"))
  expect_true(all(res$pdf[res$x > 900] == 0))
})

# load("tests/assets/test.RData")
test_that(".initialize_grenander works on a real cutmat", {
  load("../assets/test.RData")
  dist_mat <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                        cutmat = cutmat_dying,
                                        fragment_locations = NULL,
                                        num_peak_limit = 3,
                                        peak_locations = peak_locations,
                                        peak_width = peak_width)
  scaling_factor <- max(dist_mat@x)
  
  res <- .initialize_grenander(dist_mat = dist_mat,
                               scaling_factor = scaling_factor)
  
  expect_true(inherits(res, "grenander"))
  expect_true(all(diff(res$x) >= 0))
  expect_true(all(diff(res$pdf) <= 1e-6))
  
  area <- sum(diff(res$x)*res$pdf[-length(res$pdf)])
  expect_true(abs(area - 1) <= 1e-6)
})

####################################

## .compute_loglikelihood is correct

# load("tests/assets/test.RData")
test_that(".compute_loglikelihood is correct", {
  load("../assets/test.RData")
  dist_mat <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                        cutmat = cutmat_dying,
                                        fragment_locations = NULL,
                                        num_peak_limit = 3,
                                        peak_locations = peak_locations,
                                        peak_width = peak_width)
  scaling_factor <- max(dist_mat@x)/10
  
  grenander_obj <- .initialize_grenander(dist_mat = dist_mat,
                                         scaling_factor = scaling_factor)
  
  res <- .compute_loglikelihood(dist_mat = dist_mat,
                                grenander_obj = grenander_obj,
                                prior_vec = peak_prior)
  
  expect_true(is.numeric(res))
  expect_true(res < 0)
  
  # manual calculation
  res2 <- 0
  for(i in 1:nrow(dist_mat)){
    tmp <- 0
    peak_idxs <- which(as.numeric(dist_mat[i,]) != 0)
    stopifnot(length(peak_idxs) <= 6)
    
    for(j in 1:length(peak_idxs)){
      dist_val <- dist_mat[i,peak_idxs[j]]
      prob_dist_given_peak <- evaluate_grenander(
        obj = grenander_obj,
        x = dist_val
      )
      tmp <- tmp + peak_prior[peak_idxs[j]] * prob_dist_given_peak
    }
    
    res2 <- res2 + log(tmp)
  }
  
  expect_true(abs(res - res2) <= 1e-6)
})

####################################

# load("tests/assets/test.RData")
## .e_step is correct
test_that(".e_step works", {
  load("../assets/test.RData")
  dist_mat <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                        cutmat = cutmat_dying,
                                        fragment_locations = NULL,
                                        num_peak_limit = 3,
                                        peak_locations = peak_locations,
                                        peak_width = peak_width)
  grenander_obj <- .initialize_grenander(bandwidth = 200,
                                         dist_mat = dist_mat,
                                         discretization_stepsize = 10)
  
  res <- .e_step(dist_mat = dist_mat,
                 grenander_obj = grenander_obj,
                 prior_vec = peak_prior)
  
  expect_true(all(abs(Matrix::rowSums(res)-1)<=1e-6))
})

####################################

# load("tests/assets/test.RData")
## .m_step is correct
test_that(".m_step works", {
  load("../assets/test.RData")
  dist_mat <- .compute_frag_peak_matrix(bool_lock_within_peak = T,
                                        cutmat = cutmat_dying,
                                        fragment_locations = NULL,
                                        num_peak_limit = 3,
                                        peak_locations = peak_locations,
                                        peak_width = peak_width)
  grenander_obj <- .initialize_grenander(bandwidth = 200,
                                         dist_mat = dist_mat,
                                         discretization_stepsize = 10)
  assignment_mat <- .e_step(dist_mat = dist_mat,
                            grenander_obj = grenander_obj,
                            prior_vec = peak_prior)
  
  res <- .m_step(
    assignment_mat = assignment_mat,
    bandwidth = 200,
    dist_mat = dist_mat,
    discretization_stepsize = 10
  )
  
  expect_true(inherits(res, "grenander"))
  expect_true(all(diff(res$x) >= 0))
  expect_true(all(diff(res$pdf) <= 1e-6))
  expect_true(abs(res$param$left_area - 1) <= 0.1)
  expect_true(abs(res$param$right_area - 1) <= 0.1)
})

####################################

## peak_mixture_modeling is correct
# load("tests/assets/test.RData")
test_that("peak_mixture_modeling works", {
  load("../assets/test.RData")
  res1 <- peak_mixture_modeling(bandwidth = 200,
                                cutmat = cutmat_dying,
                                peak_locations = peak_locations,
                                peak_prior = peak_prior,
                                peak_width = peak_width,
                                bool_freeze_prior = F,
                                max_iter = 5,
                                verbose = 0)
  expect_true(inherits(res1, "peakDistribution"))
  # plot(res1$grenander_obj$x, res1$grenander_obj$pdf, main = "Dying")
  # res1$grenander_obj$param
  
  res2 <- peak_mixture_modeling(bandwidth = 200,
                                cutmat = cutmat_winning,
                                peak_locations = peak_locations,
                                peak_prior = peak_prior,
                                peak_width = peak_width,
                                bool_freeze_prior = T,
                                max_iter = 5,
                                verbose = 0)
  expect_true(inherits(res2, "peakDistribution"))
  # plot(res2$grenander_obj$x, res2$grenander_obj$pdf, main = "Winning")
  # res2$grenander_obj$param
  
  res3 <- peak_mixture_modeling(bandwidth = 200,
                                cutmat = rbind(cutmat_dying, cutmat_winning), 
                                peak_locations = peak_locations,
                                peak_prior = peak_prior,
                                peak_width = peak_width,
                                bool_freeze_prior = F,
                                max_iter = 5,
                                verbose = 0)
  expect_true(inherits(res3, "peakDistribution"))
  # plot(res3$grenander_obj$x, res3$grenander_obj$pdf, main = "Both")
  # res3$grenander_obj$param
  
  # there's at least some non-trivial optimization going on
  expect_true(max(c(length(res1$loglikelihood_vec),
                    length(res2$loglikelihood_vec),
                    length(res3$loglikelihood_vec))) > 2)
})

