lineage_cv_finalize <- function(cell_features,
                                cell_lineage,
                                fit_res,
                                lineage_future_count){
  test_vec <- sapply(fit_res, function(x){x$test_loglik})
  test_quantile <- apply(test_vec, 1, function(vec){stats::median(vec)})
  lambda_sequence <- fit_res[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(test_quantile)]
  
  final_fit <- multiomeFate:::lineage_imputation(
    cell_features = cell_features,
    cell_lineage = cell_lineage,
    coefficient_initial_list = fit_res[[1]]$train_fit$fit_list[[which.min(test_quantile)]]$coefficient_vec,
    lambda = lambda,
    lineage_future_count = lineage_future_count,
    verbose = 0
  )
  
  ########
  
  cell_features <- cbind(1, cell_features)
  colnames(cell_features)[1] <- "Intercept"
  stopifnot(all(colnames(cell_features) == names(final_fit$fit$coefficient_vec)))
  
  cell_imputed_score <- as.numeric(cell_features %*% final_fit$fit$coefficient_vec)
  names(cell_imputed_score) <- rownames(cell_features)
  
  cell_imputed_count <- exp(cell_imputed_score)
  uniq_lineage <- sort(unique(cell_lineage))
  lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
    sum(cell_imputed_count[which(cell_lineage == lineage)])
  })
  cell_imputed_score2 <- log10(exp(cell_imputed_score)) # this one is on the log10 scale
  
  list(cell_imputed_score = cell_imputed_score2,
       coefficient_vec = final_fit$fit$coefficient_vec,
       lambda = lambda,
       lineage_imputed_count = lineage_imputed_count)
}