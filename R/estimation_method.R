.estimate_g <- function(mat_x1, mat_y2, df_y, est_options){
  stopifnot(est_options[["method"]] == "glmnet")
  if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")
  
  p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)
  res_g <- lapply(1:p2, function(j){
    
    if(est_options$enforce_cis){
      ## find the region around each peak
      idx_x <- est_options$ht_map[[as.character(j)]]
    } else {
      #[[note to self: I think there's a cleaner way to write this]]
      idx_x <- 1:p1
    }
    if(length(idx_x) == 0) next() #[[note to self: do I need to handle this better?]]
    
    idx_y <- j
  
    ## apply glmnet
    lis <- .glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
                         family = est_options$family, 
                         switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
                         alpha = est_options$alpha, standardize = est_options$standardize, intercept = est_options$intercept,
                         cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice)
    
    lis
  })

  if(length(colnames(mat_y2)) != 0) names(res_g) <- colnames(mat_y2)
  
  res_g
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