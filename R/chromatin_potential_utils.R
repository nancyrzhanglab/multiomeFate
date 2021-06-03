#' Check all the options for chromatin potential
#' 
#' This function checks all that options in \code{options} are appropriate
#' for the method specified in \code{dim_method}, \code{nn_method}, 
#' \code{form_method}, \code{est_method}, 
#' \code{cand_method}, and \code{rec_method},
#' and fills in the default values for any options not present in \code{options}
#'
#' @param dim_method string
#' @param nn_method string
#' @param form_method string
#' @param est_method string
#' @param cand_method string
#' @param rec_method string
#' @param options list
#'
#' @return list of 6 lists
.chrom_options <- function(dim_method, nn_method, form_method, est_method, 
                           cand_method, rec_method, 
                           options){
  stopifnot(dim_method %in% c("pca"))
  stopifnot(nn_method %in% c("annoy"))
  stopifnot(form_method %in% c("literal", "average"))
  stopifnot(est_method %in% c("glmnet", "threshold_glmnet"))
  stopifnot(cand_method %in% c("nn_any", "nn_freq", "all"))
  stopifnot(rec_method %in% c("nn", "distant_cor", "distant_cor_oracle"))
  stopifnot(is.list(options))
  
  idx <- grep("^dim_*|^nn_*|^form_*|^est_*|^cand_*|^rec_*", names(options))
  if(length(idx) != length(options)){
    if(length(idx) == 0){
      warning("Option ", paste0(names(options), collapse = " and "), " not used.")
    } else {
      warning("Option ", paste0(names(options)[-idx], collapse = " and "), " not used.")
    }
  }
  
  ## [note to self: I should check the type (num,char,bool) of all the options]
  dim_options <- .dim_options(dim_method, options)
  nn_options <- .nn_options(nn_method, options)
  form_options <- .forming_options(form_method, options)
  est_options <- .estimation_options(est_method, options)
  cand_options <- .candidate_options(cand_method, options)
  rec_options <- .recruit_options(rec_method, est_options, options)
  
  list(dim_options = dim_options, nn_options = nn_options,
       form_options = form_options, est_options = est_options,
       cand_options = cand_options, rec_options = rec_options)
}

##############

.dim_options <- function(dim_method, options){
  prefix <- "dim"
  
  if(dim_method == "pca"){
    list_default <- list(mean = T, sd = T, nlatent_x = 10, nlatent_y = 10)
    dim_options <- .fill_options(options, list_default, prefix)
  } 
  
  dim_options$method <- dim_method
  dim_options
}

.nn_options <- function(nn_method, options){
  prefix <- "nn"
  
  if(nn_method == "annoy"){
    list_default <- list(nn = 20, parallel = F, ntrees = 50, metric = "euclidean",
                         verbose = F)
    nn_options <- .fill_options(options, list_default, prefix)
  } 
  
  nn_options$method <- nn_method
  nn_options
}

##############

.forming_options <- function(form_method, options){
  prefix <- "form"
  
  if(form_method == "literal"){
    form_options <- list()
  } else if(form_method == "average"){
    list_default <- list(average = "median")
    form_options <- .fill_options(options, list_default, prefix)
    
    stopifnot(form_options$average %in% c("mean", "median"))
  }
  
  form_options$method <- form_method
  form_options
}

.estimation_options <- function(est_method, options){
  prefix <- "est"

  if(est_method == "glmnet"){
    
    list_default <- list(family = "gaussian", 
                         enforce_cis = T, cis_window = 200,
                         switch = F, switch_cutoff = 10,
                         alpha = 1, standardize = F, intercept = T,
                         cv = T, nfolds = 5, cv_choice = "lambda.1se",
                         bool_round = F, run_diagnostic = T,
                         hold_initial = F, parallel = F, verbose = F)
    est_options <- .fill_options(options, list_default, prefix)
  } else if(est_method == "threshold_glmnet"){
    
    list_default <- list(family = "gaussian", 
                         enforce_cis = T, cis_window = 200,
                         switch = F, switch_cutoff = 10,
                         alpha = 1, standardize = F, intercept = T,
                         cv = T, nfolds = 5, cv_choice = "lambda.1se",
                         bool_round = F, 
                         num_iterations = 10, initial_quantile = 0.25,
                         run_diagnostic = T, hold_initial = F, 
                         parallel = F, verbose = F)
    est_options <- .fill_options(options, list_default, prefix)
    
    stopifnot(est_options$family == "gaussian")
  }
  
  if(est_options$family == "poisson") stopifnot(est_options$bool_round) # [[note to self: glmnet seems to be not handle non-integers for poisson glm...]]
  stopifnot(!est_options$standardize) # [[note to self: I should eventually code this up although it'll be quite a hassle...]]
  stopifnot(est_options$nfolds >= 3, est_options$cv_choice %in% c("lambda.1se", "lambda.min")) # requirement by glmnet
  stopifnot(est_options$cis_window > 0)
  
  est_options$method <- est_method
  est_options
}

#########

.candidate_options <- function(cand_method, options){
  prefix <- "cand"
  
  if(cand_method == "nn_any"){
    list_default <- list(num_cand = 10, only_latest = T, run_diagnostic = T)
    cand_options <- .fill_options(options, list_default, prefix)
    
  } else if(cand_method == "nn_freq"){
    list_default <- list(num_cand = 30, run_diagnostic = T)
    cand_options <- .fill_options(options, list_default, prefix)
    
  } else if(cand_method == "all"){
    list_default <- list(run_diagnostic = T)
    cand_options <- .fill_options(options, list_default, prefix)
  }
  
  cand_options$method <- cand_method
  cand_options
}

.recruit_options <- function(rec_method, est_options, options){
  prefix <- "rec"
  
  if(rec_method == "nn"){
    list_default <- list(nn = 10, num_rec = 10, average = "mean", parallel = F, 
                         run_diagnostic = T, verbose = F)
    rec_options <- .fill_options(options, list_default, prefix)

    stopifnot(rec_options$average %in% c("mean", "median"))
    
  } else if(rec_method %in% c("distant_cor", "distant_cor_oracle")){
    list_default <- list(cor_method = "pearson", parallel = F, 
                         bool_avg_from = T, bool_pred_nn = T, 
                         run_diagnostic = T, verbose = F,
                         matched_sampling_rate = 1)
    rec_options <- .fill_options(options, list_default, prefix)
    
    stopifnot(rec_options$matched_sampling_rate > 0, rec_options$matched_sampling_rate <= 1)
    stopifnot(rec_options$method %in% c("pearson", "spearman", "kendall"))
  } 
  
  rec_options$method <- rec_method
  rec_options$family <- est_options$family
  rec_options
}

########################

.fill_options <- function(options, list_default, prefix){
  org_name <- names(list_default)
  target_name <- paste0(prefix, "_", org_name)
  
  # do a warning check
  idx <- grep(paste0("^", prefix, "_*"), names(options))
  if(length(idx) > 0){
    idx2 <- which(!names(options)[idx] %in% target_name)
    if(length(idx2) > 0){
      warning("Option ", paste0(names(options)[idx[idx2]], collapse = " and "), " not used.")
    }
  }
  
  res <- vector("list", 0)
  for(i in 1:length(list_default)){
    if(!is.null(options[[target_name[i]]])){
      res[[org_name[i]]] <- options[[target_name[i]]]
    } else {
      res[[org_name[i]]] <- list_default[[org_name[i]]]
    }
  }
  
  res
}

# assumes x is peak, y is gene
.gene_peak_map <- function(df_x, df_y, est_options){
  stopifnot(est_options$cis_window > 0)
  n <- nrow(df_y)
  
  ht_map <- hash::hash()
  for(i in 1:n){
    loc <- df_y$location[i]
    idx <- which(abs(df_x$location - loc) <= est_options$cis_window)
    ht_map[[as.character(i)]] <- idx
  }
  
  est_options$ht_map <- ht_map
  est_options
}
