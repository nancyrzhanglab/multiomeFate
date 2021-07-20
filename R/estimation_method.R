#' Estimate the GLM that links the two modalities
#' 
#' Estimate the GLM for each variable of \code{mat_y2} onto
#' all the variables in \code{mat_x1}. 
#' 
#' If \code{est_options$enforce_cis} is \code{TRUE}, this function
#' will assume \code{est_options$ht_map} is a code{hash} object
#' already populated with all the variables in Modality 1 that 
#' are associated with Modality 2
#' 
#' The options are:
#' \itemize{
#' \item \code{glmnet}: Estimate the link from \code{mat_x1} to \code{mat_y2}
#' using \code{glmnet::glmnet}. The options \code{est_options$family}, \code{est_options$alpha}, 
#' \code{est_options$standardize}, and \code{est_options$intercept} are the arguments
#' for \code{glmnet::glmnet}. If \code{est_options$cv} is \code{TRUE}, then
#' \code{glmnet::cv.glmnet} is used to select the regularization parameter
#' using \code{est_options$nfolds} folds via \code{est_options$cv_choice}.
#' If \code{est_options$cv} is \code{FALSE}, use the densest solution.
#' \code{est_options$switch} is a boolean that automatically switches 
#' from \code{glmnet::glmnet} to \code{stats::glm} if there are 
#' \code{est_options$switch_cutoff} times more cells in \code{mat_x1} than variables
#' in \code{mat_x1}. 
#' }
#' 
#' @param mat_x1 Output of \code{.init_est_matrices} or \code{.update_estimation_matrices}, representing
#' the data for Modality 1
#' @param mat_y2 Output of \code{.init_est_matrices} or \code{.update_estimation_matrices}, representing
#' the data for Modality 2
#' @param est_options one of the outputs from \code{.chrom_options}
#'
#' @return a list of a matrix \code{mat_g} (of dimensions \code{ncol(mat_x1)} by \code{ncol(mat_y2)})
#' and a vector \code{vec_g} (of length \code{ncol(mat_y2)})
.estimate_g <- function(mat_x1, mat_y2, weights, est_options){
  stopifnot(all(is.na(weights)) || length(weights) == nrow(mat_x1))
  if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")
  
  if(est_options[["method"]] %in% c("glmnet", "threshold_glmnet")){
    res_g <- .estimate_g_glmnet(mat_x1, mat_y2, weights, est_options)
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

.estimate_g_glmnet <- function(mat_x1, mat_y2, weights, est_options){
  stopifnot(nrow(mat_x1) == nrow(mat_y2))
  if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")
  
  p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)

  # initialize variables for the loop
  if(!est_options$parallel && future::nbrOfWorkers() == 1){
    my_lapply <- pbapply::pblapply
    if(est_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
  } else {
    my_lapply <- future.apply::future_lapply
  }
  
  list_res <- my_lapply(1:p2, function(j){
    if(est_options$enforce_cis){
      ## find the region around each peak
      idx_x <- est_options$ht_map[[as.character(j)]]
    } else {
      #[[note to self: I think there's a cleaner way to write this]]
      #[[note to self: write a test for this scenario]]
      idx_x <- 1:p1
    }
    if(length(idx_x) == 0) return(val_int = mean(mat_y2[,j]), vec_coef = rep(0, p1))
    
    idx_y <- j
    
    ## apply glmnet
    if(est_options[["method"]] == "glmnet"){
      .glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
                    weights = weights,
                    family = est_options$family, 
                    switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
                    alpha = est_options$alpha, standardize = est_options$standardize, intercept = est_options$intercept,
                    cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice,
                    bool_round = est_options$bool_round)
    } else{
      .threshold_glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
                              weights = weights,
                              family = est_options$family, 
                              switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
                              alpha = est_options$alpha, standardize = est_options$standardize, intercept = est_options$intercept,
                              cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice,
                              bool_round = est_options$bool_round,
                              num_iterations = est_options$num_iterations, 
                              initial_quantile = est_options$initial_quantile)
    }
  })
  
  .transform_est_matrix(list_res, est_options, p1)
}

#####################################################

.glmnet_fancy <- function(x, y, weights, family, 
                          switch, switch_cutoff,
                          alpha, standardize, intercept,
                          cv, nfolds, cv_choice, bool_round){
  n <- length(y); p <- ncol(x)
  if(bool_round) y <- round(y)
  
  if(switch & n > p*switch_cutoff){
    if(length(weights) <= 1) weights <- NULL
    # use glm
    if(intercept) x <- cbind(x,1)
    if(family == "poisson") {
      family_func <- stats::poisson()
    } else if(family == "gaussian") {
      family_func <- stats::gaussian()
    } else {
      stop("family not found")
    }
    res <- stats::glm.fit(x, y, weights = weights, family = family_func, intercept = F)
    vec_coef <- res$coefficients[1:p]
    vec_coef[is.na(vec_coef)] <- 0
    val_int <- ifelse(intercept, vec_coef[p+1], 0)
    
  } else {
    if(length(weights) <= 1) weights <- rep(1, length(y))
    
    if(cv & n > 10*nfolds){
      # use cv.glmnet
      res <- glmnet::cv.glmnet(x, y, weights = weights, family = family, nfolds = nfolds, alpha = alpha,
                               standardize = standardize, intercept = intercept)
      lambda <- res[[cv_choice]]
      res <- glmnet::glmnet(x, y, weights = weights, family = family, alpha = alpha,
                            standardize = standardize, intercept = intercept)
    } else {
      # use glmnet and pick the densest solution
      res <- glmnet::glmnet(x, y, weights = weights, family = family, alpha = alpha,
                            standardize = standardize, intercept = intercept)
      tmp <- res$lambda; len <- length(tmp); stopifnot(len > 1)
      lambda <- mean(tmp[c(len-1,len)])
    }
    
    tmp <- as.numeric(stats::coef(res, s = lambda))
    vec_coef <- tmp[-1]; val_int <- tmp[1]
  }
  
  list(val_int = val_int, vec_coef = vec_coef)
}

.threshold_glmnet_fancy <- function(x, y, weights, family, 
                                    switch, switch_cutoff,
                                    alpha, standardize, intercept,
                                    cv, nfolds, cv_choice, bool_round,
                                    num_iterations, initial_quantile,
                                    tol = 1e-4){
  prev_threshold <- stats::quantile(y, probs = initial_quantile)
  iter <- 1
  
  while(iter <= num_iterations){
    # update the regression
    idx <- which(y >= prev_threshold)
    if(length(weights) == 1) weight_vec <- NA else weight_vec <- weights[idx]
    res_glm <- .glmnet_fancy(x[idx,,drop = F], y[idx], weight_vec,
                             family, switch, switch_cutoff,
                             alpha, standardize, intercept,
                             cv, nfolds, cv_choice, bool_round)
    
    # update the threshold
    # [[note to self: assumes family is Gaussian essentially]]
    next_threshold <- .update_threshold_glmnet(x, y, weights, res_glm)
    
    if(abs(next_threshold - prev_threshold) <= tol) break()
    iter <- iter+1
  }
  
  list(val_int = res_glm$val_int, vec_coef = res_glm$vec_coef,
       val_threshold = next_threshold)
}

.update_threshold_glmnet <- function(x, y, weights, res_glm){
  pred_y <- as.numeric(x %*% res_glm$vec_coef) + res_glm$val_int
  
  if(length(weights) <= 1) {
    f <- function(val_threshold, y, pred_y){
      sum((y - pmax(pred_y, val_threshold))^2)
    }
  } else {
    f <- function(val_threshold, y, pred_y){
      sum(weights*(y - pmax(pred_y, val_threshold))^2)
    }
  }
  
  stats::optimize(f, interval = c(min(y), max(y)), maximum = F,
                  y = y, pred_y = pred_y)$minimum
}

#################################

.transform_est_matrix <- function(list_res, est_options, p1){
  p2 <- length(list_res)
  bool_thres <- "val_threshold" %in% names(list_res[[1]])
  
  mat_g <- matrix(0, nrow = p1, ncol = p2)
  vec_g <- rep(0, length = p2)
  if(bool_thres) {
    vec_threshold <- rep(0, length = p2)
  }
  
  for(j in 1:p2){
    vec_g[j] <- list_res[[j]]$val_int
    if(bool_thres){ vec_threshold[j] <- list_res[[j]]$val_threshold }
  
    if(est_options$enforce_cis){
      idx_x <- est_options$ht_map[[as.character(j)]]
      mat_g[idx_x,j] <- list_res[[j]]$vec_coef
    } else {
      mat_g[,j] <- list_res[[j]]$vec_coef
    }
  }
  
  if(bool_thres){
    list(mat_g = mat_g, vec_g = vec_g, vec_threshold = vec_threshold)
  } else {
    list(mat_g = mat_g, vec_g = vec_g)
  }
  
}