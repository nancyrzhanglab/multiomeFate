context("Testing peak modeling")

## .lrt_onefold is correct

# load("tests/assets/test.RData")
test_that(".lrt_onefold works on a real cutmatrix", {
  load("../assets/test.RData")
  
  frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
  frag_die <- .extract_fragment_from_cutmat(cutmat_dying)
  
  set.seed(10)
  idx_win <- sample(1:length(frag_win), size = round(length(frag_win)/2))
  idx_die <- sample(1:length(frag_die), size = round(length(frag_die)/2))
  
  bandwidth <- 200
  discretization_stepsize <- 10
  
  res <- .lrt_onefold(
    bandwidth = bandwidth,
    frag_die = frag_die, 
    frag_win = frag_win,
    idx_die_split1 = idx_die,
    idx_win_split1 = idx_win,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = T, 
    max_iter = 100,
    min_prior = 0.01,
    num_peak_limit = 4,
    tol = 1e-6,
    verbose = 0
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("grenander_both",
                                             "grenander_die",
                                             "grenander_win",
                                             "loglikelihood_denom",
                                             "loglikelihood_num",
                                             "teststat"))))
})

test_that(".lrt_onefold fails to reject for synthetic sharp null", {
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
  peak_width <- 100
  bandwidth <- 20
  discretization_stepsize <- 2
  peak_prior <- c(0.5,0.5)
  
  frag_all <- .extract_fragment_from_cutmat(cutmat)
  
  trials <- 50
  pvalue_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    print(i)
    
    idx <- sample(1:length(frag_all), size = round(length(frag_all)/2))
    idx_die <- sample(1:length(idx), size = round(length(idx)/2))
    idx_win <- sample(1:(length(frag_all) - length(idx)),
                      size = round((length(frag_all) - length(idx))/2))
    
    res <- .lrt_onefold(
      bandwidth = bandwidth,
      frag_die = frag_all[idx], 
      frag_win = frag_all[-idx],
      idx_die_split1 = idx_die,
      idx_win_split1 = idx_win,
      peak_locations = peak_locations,
      peak_prior = peak_prior,
      peak_width = peak_width,
      discretization_stepsize = discretization_stepsize, 
      bool_lock_within_peak = T, 
      max_iter = 100,
      min_prior = 0.01,
      num_peak_limit = 4,
      tol = 1e-6,
      verbose = 0
    )
    
    min(1/res$teststat,1)
  })
  
  expect_true(all(quantile(pvalue_vec) - seq(0,1,length.out=5) >= -.1))
})

####################

# .compute_crossfit_teststat is correct

# load("tests/assets/test.RData")
test_that(".compute_crossfit_teststat works on a real cutmatrix", {
  load("../assets/test.RData")
  
  frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
  frag_die <- .extract_fragment_from_cutmat(cutmat_dying)
  
  set.seed(10)
  idx_win <- sample(1:length(frag_win), size = round(length(frag_win)/2))
  idx_die <- sample(1:length(frag_die), size = round(length(frag_die)/2))
  
  bandwidth <- 200
  discretization_stepsize <- 10
  
  res <- .compute_crossfit_teststat(
    bandwidth = bandwidth,
    frag_die = frag_die, 
    frag_win = frag_win,
    idx_die = idx_die,
    idx_win = idx_win,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = T, 
    max_iter = 100,
    min_fragments = 6,
    min_prior = 0.01,
    num_peak_limit = 4,
    tol = 1e-6,
    verbose = 0
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("grenander_both",
                                             "grenander_die",
                                             "grenander_win",
                                             "loglikelihood_denom",
                                             "loglikelihood_num",
                                             "pvalue",
                                             "teststat"))))
})

# load("tests/assets/test.RData")
test_that(".compute_crossfit_teststat fails to reject under the sharp null", {
  load("../assets/test.RData")
  bandwidth <- 200
  discretization_stepsize <- 10
  
  frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
  frag_die <- .extract_fragment_from_cutmat(cutmat_dying)
  frag_all <- c(frag_win, frag_die)
  
  trials <- 10
  
  pvalue_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    idx <- sample(1:length(frag_all), size = round(length(frag_all)/2))
    idx_die <- sample(1:length(idx), size = round(length(idx)/2))
    idx_win <- sample(1:(length(frag_all) - length(idx)),
                      size = round((length(frag_all) - length(idx))/2))
    
    res <- .compute_crossfit_teststat(
      bandwidth = bandwidth,
      frag_die = frag_all[idx], 
      frag_win = frag_all[-idx],
      idx_die = idx_die,
      idx_win = idx_win,
      peak_locations = peak_locations,
      peak_prior = peak_prior,
      peak_width = peak_width,
      discretization_stepsize = discretization_stepsize, 
      bool_lock_within_peak = T, 
      max_iter = 100,
      min_fragments = 6,
      min_prior = 0.01,
      num_peak_limit = 4,
      tol = 1e-6,
      verbose = 0
    )
    print(paste0("Trial ", i, ": p-value of ", round(min(1/res$teststat,1),5)))
    
    min(1/res$teststat,1)
  })
})

#######################

## peak_testing is correct

# load("tests/assets/test.RData")
test_that("peak_testing works on a real cutmatrix", {
  load("../assets/test.RData")
  bandwidth <- 200
  discretization_stepsize <- 10
  
  res <- peak_testing(
    bandwidth = bandwidth,
    cutmat_dying = cutmat_dying, 
    cutmat_winning = cutmat_winning,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize,
    verbose = 0
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("grenander_both1",
                                             "grenander_both2",
                                             "grenander_die1",
                                             "grenander_die2",
                                             "grenander_win1",
                                             "grenander_win2",
                                             "loglikelihood_denom1",
                                             "loglikelihood_denom2",
                                             "loglikelihood_num1",
                                             "loglikelihood_num2",
                                             "teststat"))))
})


