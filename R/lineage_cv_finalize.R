#' Finalize lineage imputation after cross validation
#' 
#' This function is used after running \code{multiomeFate::lineage_cv()}.
#' Chooses \code{lambda} by minimizing the median held-out objective across folds,
#' then refits once on all cells at the chosen \code{lambda}.
#'
#' @inheritParams lineage_cv
#' @param fit_res This is the output of \code{lineage_cv()}. 
#' 
#' @returns A list with the following elements: \code{cell_imputed_score} (a
#' vector of length \code{nrow(cell_features)}) that denotes the predicted
#' progenies spawning from each particular cell, \code{coefficient_vec} (the
#' coefficient vector of lenght \code{ncol(cell_features)+1}) that denotes
#' the coefficients in the GLM, \code{lambda} (the chosen parameter after
#' cross-validation), and \code{lineage_imputed_count} (the vector of length
#' \code{lineage_future_count} that denotes the number of predicted cells 
#' at the future timepoint in each lineage).
#' @export
lineage_cv_finalize <- function(cell_features,
                                cell_lineage,
                                fit_res,
                                lineage_future_count){
  test_vec <- sapply(fit_res, function(x){x$test_loglik})
  test_quantile <- apply(test_vec, 1, function(vec){stats::median(vec)})
  lambda_sequence <- fit_res[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(test_quantile)]
  
  final_fit <- lineage_imputation(
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
  names(lineage_imputed_count) <- uniq_lineage
  cell_imputed_score2 <- log10(exp(cell_imputed_score)) # this one is on the log10 scale
  
  list(cell_imputed_score = cell_imputed_score2,
       coefficient_vec = final_fit$fit$coefficient_vec,
       lambda = lambda,
       lineage_imputed_count = lineage_imputed_count)
}