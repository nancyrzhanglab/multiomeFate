# function to check all the options are there, and reformulate
.chrom_options <- function(form_method, est_method, cand_method, rec_method, 
                           options){
  stopifnot(form_method %in% c("literal", "average"))
  stopifnot(est_method %in% c("glmnet_yonly"))
  stopifnot(cand_method %in% c("nn"))
  stopifnot(rec_method %in% c("nn"))
  stopifnot(is.list(options))
  
  idx <- grep("^form_*|^est_*|^cand_*|^rec_*", names(options))
  if(length(idx) != length(options)){
    warning("Option ", paste0(names(options)[-idx], collapse = " and "), " not used.")
  }
  
  ## [note to self: I should check the type (num,char,bool) of all the options]
  form_options <- .forming_options(form_method, options)
  est_options <- .estimation_options(est_method, options)
  cand_options <- .candidate_options(cand_method, options)
  rec_options <- .recruit_options(rec_method, options)
  
  list(form_options = form_options, est_options = est_options,
       cand_options = cand_options, rec_options = rec_options)
}

.forming_options <- function(form_method, options){
  form_options <- vector("list", 0)
  form_options$method <- form_method
  form_options
}

.estimation_options <- function(est_method, options){
  prefix <- "est"

  if(est_method == "glmnet_yonly"){
    
    list_default <- list(family = "poisson", 
                         enforce_cis = T, cis_window = 200,
                         switch = T, switch_cutoff = 10,
                         alpha = 1, standardize = F, intercept = F,
                         cv = T, nfolds = 5, cv_choice = "lambda.1se")
    est_options <- .fill_options(options, list_default, prefix)
  }
  
  stopifnot(!est_options$standardize) # [[note to self: I should eventually code this up although it'll be quite a hassle...]]
  stopifnot(est_options$nfolds >= 3, est_options$cv_choice %in% c("lambda.1se", "lambda.min")) # requirement by glmnet
  stopifnot(est_options$cis_window > 0)
  
  
  est_options$method <- est_method
  est_options
}

.candidate_options <- function(cand_method, options){
  prefix <- "cand"
  
  if(cand_method == "nn"){
    list_default <- list(nn = 10, metric = "euclidean")
    cand_options <- .fill_options(options, list_default, prefix)
  }
  
  ## [note to self: seems like there's no metric in RANN]
  stopifnot(cand_options$metric == "euclidean")
  
  cand_options$method <- cand_method
  cand_options
}

.recruit_options <- function(rec_method, options){
  prefix <- "rec"
  
  if(rec_method == "nn"){
    list_default <- list(nn = 10, num_rec = 10, metric = "euclidean")
    rec_options <- .fill_options(options, list_default, prefix)
  }
  
  ## [note to self: seems like there's no metric in RANN]
  stopifnot(rec_options$metric == "euclidean")
  
  rec_options$method <- rec_method
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
