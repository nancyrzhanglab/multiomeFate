peak_testing <- function(bandwidth,
                         cutmat_dying, 
                         cutmat_winning,
                         peak_locations,
                         peak_prior,
                         peak_width,
                         discretization_stepsize = NA, #if NA, it defaults to max(round(max(dist_mat@x)/50),5)
                         bool_lock_within_peak = T, 
                         max_iter = 100,
                         min_fragments = 6,
                         min_prior = 0.01,
                         num_peak_limit = 4,
                         tol = 1e-6,
                         verbose = 1){
  # extract all the relevant fragments
  frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
  frag_die <- .extract_fragment_from_cutmat(cutmat_dying)
  
  if(length(frag_win) < min_fragments | length(frag_die) < min_fragments) return(NA)
  
  # do a cross-fit
  # sample split
  idx_win <- sample(1:length(frag_win), size = round(length(frag_win)/2))
  idx_die <- sample(1:length(frag_die), size = round(length(frag_die)/2))
  
  # for one cross
  .compute_crossfit_teststat(
    bandwidth = bandwidth,
    frag_die = frag_die, 
    frag_win = frag_win,
    idx_die = idx_die,
    idx_win = idx_win,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    max_iter = max_iter,
    min_fragments = min_fragments,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    tol = tol,
    verbose = verbose
  )
}

#######################################

# see https://arxiv.org/pdf/1912.11436.pdf, equation 10
.compute_crossfit_teststat <- function(bandwidth,
                                       frag_die, 
                                       frag_win,
                                       idx_die,
                                       idx_win,
                                       peak_locations,
                                       peak_prior,
                                       peak_width,
                                       discretization_stepsize, 
                                       bool_lock_within_peak, 
                                       max_iter,
                                       min_fragments,
                                       min_prior,
                                       num_peak_limit,
                                       tol,
                                       verbose){
  len_die <- length(frag_die); len_win <- length(frag_win)
  
  stopifnot(length(idx_die) < len_die - floor(min_fragments/2),
            length(idx_win) < len_win - floor(min_fragments/2),
            length(idx_die) > floor(min_fragments/2),
            length(idx_win) > floor(min_fragments/2),
            floor(min_fragments/2) > 1)
  
  fit1 <- .lrt_onefold(
    bandwidth = bandwidth,
    frag_die = frag_die, 
    frag_win = frag_win,
    idx_die_split1 = idx_die,
    idx_win_split1 = idx_win,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    max_iter = max_iter,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    tol = tol,
    verbose = verbose
  )
  
  idx_die2 <- c(1:length(frag_die))[-idx_die]
  idx_win2 <- c(1:length(frag_win))[-idx_win]
  fit2 <- .lrt_onefold(
    bandwidth = bandwidth,
    frag_die = frag_die, 
    frag_win = frag_win,
    idx_die_split1 = idx_die2,
    idx_win_split1 = idx_win2,
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    max_iter = max_iter,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    tol = tol,
    verbose = verbose
  )
  
  teststat <- (fit1$teststat + fit2$teststat)/2
  
  list(
    grenander_both1 = fit1$grenander_both,
    grenander_both2 = fit2$grenander_both,
    grenander_die1 = fit1$grenander_die,
    grenander_die2 = fit2$grenander_die,
    grenander_win1 = fit1$grenander_win,
    grenander_win2 = fit2$grenander_win,
    loglikelihood_denom1 = fit1$loglikelihood_denom,
    loglikelihood_denom2 = fit2$loglikelihood_denom,
    loglikelihood_num1 = fit1$loglikelihood_num,
    loglikelihood_num2 = fit2$loglikelihood_num,
    pvalue = min(1/teststat, 1),
    teststat = teststat
  )
}

.lrt_onefold <- function(bandwidth,
                         frag_die, 
                         frag_win,
                         idx_die_split1,
                         idx_win_split1,
                         peak_locations,
                         peak_prior,
                         peak_width,
                         discretization_stepsize, 
                         bool_lock_within_peak, 
                         max_iter,
                         min_prior,
                         num_peak_limit,
                         tol,
                         verbose){
  len_die <- length(frag_die); len_win <- length(frag_win)
 
  # p1 is for the alternative (unconstrained)
  # compute win!=die on the first fold
  fit_win <- peak_mixture_modeling(
    bandwidth = bandwidth,
    cutmat = NULL, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    bool_freeze_prior = T, # assumed to not have a lot of fragments, so we need to freeze
    fragment_locations = frag_win[idx_win_split1], 
    max_iter = max_iter,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    return_dist_mat = F,
    tol = tol,
    verbose = verbose
  )
  
  fit_die <- peak_mixture_modeling(
    bandwidth = bandwidth,
    cutmat = NULL, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    bool_freeze_prior = T,
    fragment_locations = frag_die[idx_die_split1], 
    max_iter = max_iter,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    return_dist_mat = F,
    tol = tol,
    verbose = verbose
  )
  
  # p0 is for the null
  # compute the win=die on the second fold
  fit_both <- peak_mixture_modeling(
    bandwidth = bandwidth,
    cutmat = NULL, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    peak_width = peak_width,
    discretization_stepsize = discretization_stepsize, 
    bool_lock_within_peak = bool_lock_within_peak, 
    bool_freeze_prior = T,
    fragment_locations = c(frag_win[-idx_win_split1], frag_die[-idx_die_split1]), 
    max_iter = max_iter,
    min_prior = min_prior,
    num_peak_limit = num_peak_limit,
    return_dist_mat = T, # needed for the numerator likelihood later
    tol = tol,
    verbose = verbose
  )
  stopifnot(nrow(fit_both$dist_mat) == len_win + len_die - length(idx_win_split1) - length(idx_die_split1))
  loglikelihood_denom <- fit_both$loglikelihood_val
  
  # compute the likelihood ratio of L(win!=die)/L(win=die) on the second fold
  # that is, p1(Y_second)/p0(Y_second)
  loglikelihood_outofsample_win <- .compute_loglikelihood(
    dist_mat = fit_both$dist_mat[1:(len_win - length(idx_win_split1)),],
    grenander_obj = fit_win$grenander_obj,
    prior_vec = fit_win$prior_vec
  )
  loglikelihood_outofsample_die <- .compute_loglikelihood(
    dist_mat = fit_both$dist_mat[(len_win-length(idx_win_split1)+1):nrow(fit_both$dist_mat),],
    grenander_obj = fit_die$grenander_obj,
    prior_vec = fit_die$prior_vec
  )
  loglikelihood_num <- loglikelihood_outofsample_win + loglikelihood_outofsample_die
  
  # return test statistic
  teststat <- exp(loglikelihood_num - loglikelihood_denom)
  
  list(
    grenander_both = fit_both$grenander_obj,
    grenander_die = fit_die$grenander_obj,
    grenander_win = fit_win$grenander_obj,
    loglikelihood_denom = loglikelihood_denom,
    loglikelihood_num = loglikelihood_num,
    teststat = teststat
  )
}
