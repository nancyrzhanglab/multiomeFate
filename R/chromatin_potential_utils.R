# function to check all the options are there, and reformulate
.chromatin_options <- function(forming_method, estimation_method, 
                                    candidate_method, recruit_method, options){
  stopifnot(forming_method %in% c("literal"))
  stopifnot(estimation_method %in% c("glmnet"))
  stopifnot(candidate_method %in% c("nn"))
  stopifnot(recruit_method %in% c("singleton"))
  
  form_options <- .forming_options(forming_method, options)
  est_options <- .estimation_options(estimation_method, options)
  cand_options <- .candidate_options(candidate_method, options)
  rec_options <- .recruit_options(recruit_method, options)
  
  list(form_options = form_options, est_options = est_options,
       cand_options = cand_options, rec_options = rec_options)
}

.forming_options <- function(forming_method, options){
  form_options <- vector("list", 0)
  form_options$method <- forming_method
  form_options
}

.estimation_options <- function(estimation_method, options){
  est_options <- vector("list", 0)
  est_options$method <- estimation_method
  
  if(est_options$method == "glmnet"){
    # family
    if(!is.null(options$est_family)){
      est_options$family <- options$est_family
    } else {
      est_options$family <- "poisson"
    }
    
    # standardize T/F
    if(!is.null(options$est_standardize)){
      est_options$standardize <- options$est_standardize
    } else {
      est_options$standardize <- F
    }
    
    # intercept T/F
    if(!is.null(options$est_intercept)){
      est_options$intercept <- options$est_intercept
    } else {
      est_options$intercept <- T
    }
    
    # alpha
    if(!is.null(options$est_alpha)){
      est_options$alpha <- options$est_alpha
    } else {
      est_options$alpha <- 1
    }
    
    # cv T/F
    if(!is.null(options$est_cv)){
      est_options$cv <- options$est_cv
    } else {
      est_options$cv <- T
    }
    
    # cv num.folds
    if(!is.null(options$est_nfolds)){
      est_options$nfolds <- options$est_nfolds
    } else {
      est_options$nfolds <- 5
    }
  }
  
  est_options
}

.candidate_options <- function(candidate_method, options){
  cand_options <- vector("list", 0)
  cand_options$method <- candidate_method
  
  if(cand_options$method == "nn"){
    # number of nn
    if(!is.null(options$candidate_nn)){
      cand_options$nn <- options$candidate_nn
    } else {
      cand_options$nn <- 10
    }
    
    # distance metric
    if(!is.null(options$candidate_metric)){
      cand_options$metric <- options$candidate_metric
      ## [note to self: seems like there's no metric in RANN]
      stopifnot(cand_options$metric == "euclidean")
    } else {
      cand_options$metric <- "euclidean"
    }
  }
  
  cand_options
}

.recruit_options <- function(recruit_method, options){
  rec_options <- vector("list", 0)
  rec_options$method <- recruit_method
  
  if(rec_options$method == "singleton"){
    # distance metric
    if(!is.null(options$rec_metric)){
      rec_options$metric <- options$rec_metric
      ## [note to self: seems like there's no metric in RANN]
      stopifnot(rec_options$metric == "euclidean")
    } else {
      rec_options$metric <- "euclidean"
    }
  }
  
  rec_options
}
