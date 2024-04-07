lineage_cv <- function(cell_features,
                       cell_lineage,
                       future_timepoint,
                       lineage_future_count,
                       lambda_initial,
                       lambda_sequence_length,
                       tab_mat,
                       num_folds = 10,
                       savefile_tmp = NULL,
                       seed_number = 10,
                       verbose = 0
){
  tmp <- construct_folds(
    cell_lineage = cell_lineage,
    tab_mat = tab_mat,
    future_timepoint = future_timepoint,
    num_folds = num_folds
  )
  cv_cell_list <- tmp$cv_cell_list
  fold_lineage_list <- tmp$fold_lineage_list
  
  cv_fit_list <- vector("list", length = num_folds)
  names(cv_fit_list) <- names(cv_cell_list)
  
  if(!is.null(seed_number)) set.seed(seed_number)
  for(i in 1:num_folds){
    fold <- names(cv_cell_list)[i]
    if(verbose > 0) print(paste0("Dropping fold #", i, " out of ", num_folds))
    
    #################
    
    # training 
    cell_features_train <- cell_features[-cv_cell_list[[fold]],,drop = F]
    cell_lineage_train <- cell_lineage[-cv_cell_list[[fold]]]
    lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]
    
    if(!is.null(seed_number)) set.seed(seed_number)
    train_fit <- multiomeFate:::lineage_imputation_sequence(
      cell_features = cell_features_train,
      cell_lineage = cell_lineage_train,
      lambda_initial = lambda_initial,
      lambda_sequence_length = lambda_sequence_length,
      lineage_future_count = lineage_future_count_train,
      verbose = verbose-1
    )
    
    #################
    
    # training evaluation
    lambda_sequence <- train_fit$lambda_sequence
    
    train_loglik <- sapply(1:length(lambda_sequence), function(kk){
      multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_train,
                                            cell_lineage = cell_lineage_train,
                                            coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                            lineage_future_count = lineage_future_count_train,
                                            lambda = 0)
    })
    
    #################
    
    # testing
    cell_features_test <- cell_features[cv_cell_list[[fold]],,drop = F]
    cell_lineage_test <- cell_lineage[cv_cell_list[[fold]]]
    lineage_future_count_test <- lineage_future_count[which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]
    
    test_loglik <- sapply(1:length(lambda_sequence), function(kk){
      multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                            cell_lineage = cell_lineage_test,
                                            coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                            lineage_future_count = lineage_future_count_test,
                                            lambda = 0)
    })
    
    #################
    
    cv_fit_list[[fold]] <- list(test_loglik = test_loglik,
                                train_loglik = train_loglik,
                                train_fit = train_fit)
    
    if(all(!is.null(savefile_tmp))){
      date_of_run <- Sys.time()
      
      save(cv_fit_list, date_of_run, 
           file = savefile_tmp)
    }
  }
  
  cv_fit_list
}