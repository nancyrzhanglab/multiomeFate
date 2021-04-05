.estimate_g <- function(mat_x1, mat_y2, df_y, est_options){
  if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")
  
  if(est_options[["method"]] == "glmnet"){
    res_g <- .estimate_g_glmnet(mat_x1, mat_y2, df_y, est_options)
  } else {
    stop("Estimation method not found")
  }
  
  ## [[note to self: res_g$mat_g should really be a sparse matrix class]]
  stopifnot(all(!is.na(res_g$mat_g)))

  if(length(colnames(mat_y2)) != 0) {
    colnames(res_g$mat_g) <- colnames(mat_y2)
  }
  if(length(colnames(mat_x1)) != 0) {
    rownames(res_g$mat_g) <- colnames(mat_x1)
  }
  
  res_g
}

##############################

.estimate_g_glmnet <- function(mat_x1, mat_y2, df_y, est_options){
  if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")
  
  p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)
  list_res <- lapply(1:p2, function(j){
    
    if(est_options$enforce_cis){
      ## find the region around each peak
      idx_x <- est_options$ht_map[[as.character(j)]]
    } else {
      #[[note to self: I think there's a cleaner way to write this]]
      idx_x <- 1:p1
    }
    if(length(idx_x) == 0) return(val_int = 0, vec_coef = rep(0, p1))
    
    idx_y <- j
    
    ## apply glmnet
    .glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
                  family = est_options$family, 
                  switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
                  alpha = est_options$alpha, standardize = est_options$standardize, intercept = est_options$intercept,
                  cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice)
  })
  
  .transform_est_matrix(list_res, est_options, p1)
}

.glmnet_fancy <- function(x, y, family, 
                          switch, switch_cutoff,
                          alpha, standardize, intercept,
                          cv, nfolds, cv_choice){
  n <- length(y); p <- ncol(x)
  
  if(switch & n > p*switch_cutoff){
    # use glm
    if(intercept) x <- cbind(x,1)
    if(family == "poisson") {
      family_func <- stats::poisson()
    } else {
      stop("family not found")
    }
    res <- stats::glm.fit(x, y, family = family_func, intercept = F)
    vec_coef <- res$coefficients[1:p]
    vec_coef[is.na(vec_coef)] <- 0
    val_int <- ifelse(intercept, vec_coef[p+1], 0)
    
  } else {
    if(cv & n > 10*nfolds){
      # use cv.glmnet
      res <- glmnet::cv.glmnet(x, y, family = family, nfolds = nfolds, alpha = alpha,
                               standardize = standardize, intercept = intercept)
      lambda <- res[[cv_choice]]
      res <- glmnet::glmnet(x, y, family = family, alpha = alpha,
                            standardize = standardize, intercept = intercept)
    } else {
      # use glmnet and pick the densest solution
      res <- glmnet::glmnet(x, y, family = family, alpha = alpha,
                            standardize = standardize, intercept = intercept)
      tmp <- res$lambda; len <- length(tmp); stopifnot(len > 1)
      lambda <- mean(tmp[c(len-1,len)])
    }
    
    tmp <- as.numeric(stats::coef(res, s = lambda))
    vec_coef <- tmp[-1]; val_int <- tmp[1]
  }
  
  list(val_int = val_int, vec_coef = vec_coef)
}

.transform_est_matrix <- function(list_res, est_options, p1){
  p2 <- length(list_res)
  
  mat_g <- matrix(0, nrow = p1, ncol = p2)
  vec_g <- rep(0, length = p2)
  
  for(j in 1:p2){
    vec_g[j] <- list_res[[j]]$val_int
  
    if(est_options$enforce_cis){
      idx_x <- est_options$ht_map[[as.character(j)]]
      mat_g[idx_x,j] <- list_res[[j]]$vec_coef
    } else {
      mat_g[,j] <- list_res[[j]]$vec_coef
    }
  }
  
  list(mat_g = mat_g, vec_g = vec_g)
}