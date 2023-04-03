peak_testing <- function(res_all, # outputs of peak_mixture_modeling
                         res_winning,
                         res_dying){
  stopifnot(inherits(res_all, "peakDistribution"),
            inherits(res_winning, "peakDistribution"),
            inherits(res_dying, "peakDistribution"))
  
  lrt <- -2 * (res_all$loglikelihood_val - (res_winning$loglikelihood_val + res_dying$loglikelihood_val))
  pval <- 1-stats::pchisq(lrt, df = length(res_dying$theta_vec)-1)
  
  list(lrt = lrt, pval = pval)
}