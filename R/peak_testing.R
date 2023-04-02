peak_testing <- function(res_all, # outputs of peak_mixture_modeling
                         res_winning,
                         res_dying){
  stopifnot(inherits(res_all, "peakDistribution"),
            inherits(res_winning, "peakDistribution"),
            inherits(res_dying, "peakDistribution"))
  
  lrt <- -2 * (res_all$likelihood_vec[res_all$iter] - (res_winning$likelihood_vec[res_winning$iter] + res_dying$likelihood_vec[res_dying$iter]))
  pval <- 1-stats::pchisq(lrt, df = length(res_dying$theta_vec)-1)
  
  list(lrt = lrt, pval = pval)
}