.estimate_g <- function(mat_x1, mat_y2, df_x, df_y, est_options){
  stopifnot(est_options[["method"]] == "glmnet")
  
  # check to see if enforce_cis 
}

.glmnet_fancy <- function(x, y, family, 
                          switch, switch_cutoff,
                          alpha, standardize, intercept,
                          cv, nfolds, cv_choice){
  n <- length(y); p <- ncol(x)
  
  if(switch & n > p*switch_cutoff){
    # use glm
    if(intercept) x <- cbind(x,1)
    res <- stats::glm.fit(x, y, intercept = F)
    vec_coef <- res[1:p]
    val_int <- ifelse(intercept, vec_coef[p+1], 0)
    
  } else {
    if(cv & n > 10*nfolds){
      # use cv.glmnet
      res <- glmnet::cv.glmnet(x, y, nfolds = nfolds, alpha = alpha,
                               standardize = standardize, intercept = intercept)
      lambda <- res[[cv_choice]]
      res <- glmnet::glmnet(x, y, alpha = alpha,
                            standardize = standardize, intercept = intercept)
      tmp <- as.numeric(stats::coef(res, s = lambda))
      vec_coef <- tmp[-1]; val_int <- tmp[1]
      
    } else {
      # use glmnet and pick the densest solution
      res <- glmnet::glmnet(x, y, alpha = alpha,
                            standardize = standardize, intercept = intercept)
      tmp <- res$lambda; len <- length(tmp); stopifnot(len > 1)
      lambda <- mean(tmp[c(len-1,len)])
      tmp <- as.numeric(stats::coef(res, s = lambda))
      vec_coef <- tmp[-1]; val_int <- tmp[1]
    }
  }
  
  list(val_int = val_int, vec_coef = vec_coef)
}