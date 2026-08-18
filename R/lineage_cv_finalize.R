#' Finalize CYFER after cross-validation
#'
#' This function is used after running \code{multiomeFate::cyfer()}.
#' Chooses \code{lambda} by minimizing the median held-out objective across folds,
#' then refits once on all cells at the chosen \code{lambda}.
#'
#' @inheritParams cyfer
#' @param fit_res This is the output of \code{cyfer()}.
#' @param seed_number Seed value for reproducibility reasons. Default is \code{10}.
#' Governs the optimizer's random restarts in the final refit. Set to \code{NULL}
#' to leave the random number stream untouched.
#'
#' @returns A list with the following elements: \code{cell_imputed_score} (a
#' named vector of length \code{nrow(cell_features)}) that denotes the predicted
#' progenies spawning from each particular cell, \code{coefficient_vec} (the
#' coefficient vector of length \code{ncol(cell_features)+1}) that denotes
#' the coefficients in the GLM, \code{lambda} (the chosen parameter after
#' cross-validation), and \code{lineage_imputed_count} (the vector of length
#' \code{lineage_future_count} that denotes the number of predicted cells
#' at the future time point in each lineage).
#'
#' Every cell in \code{cell_features} is scored and named by its row name, even
#' cells whose lineage is absent from \code{lineage_future_count}. Such cells are
#' dropped from the refit (they carry no future count to fit against) but the
#' fitted coefficients still apply to them, so the returned score covers every
#' row that was passed in.
#'
#' @section Scales --- exp() versus 10^():
#'
#' Three scales are in play here and they are easy to confuse. The GLM's linear
#' predictor is on the \bold{natural-log} scale,
#'
#' \deqn{Z_i = \beta_0 + X_{i,\cdot}^\top \beta,}{Z_i = beta_0 + x_i' beta,}
#'
#' so \code{coefficient_vec} is on that scale, and cell \code{i}'s expected
#' number of progeny at the future time point is \code{exp(Z_i)} --- with
#' \code{exp()}, not \code{10^}.
#'
#' The two vectors this function returns then sit on \emph{different} scales,
#' deliberately:
#'
#' \describe{
#'   \item{\code{cell_imputed_score}}{\bold{log10} of the expected progeny count,
#'     that is \code{log10(exp(Z_i))}. This is the per-cell fate potential the
#'     paper reports, and what the plotting helpers and the downstream
#'     selection/adaptation indices expect.}
#'   \item{\code{lineage_imputed_count}}{a count on the \bold{natural} scale,
#'     \code{sum(exp(Z_i))} over the cells of the lineage, directly comparable to
#'     \code{lineage_future_count}.}
#' }
#'
#' Because the cell scores have already been converted to base 10, the identity
#' linking the two uses \code{10^} and not \code{exp()}:
#'
#' \preformatted{lineage_imputed_count[l] == sum(10^cell_imputed_score[cells in l])}
#'
#' Applying \code{exp()} to \code{cell_imputed_score} is the mistake this note
#' exists to prevent: it computes \code{exp(log10(exp(Z_i)))}, which is a count on
#' no scale at all and is silently plausible-looking. To go back:
#' \code{10^cell_imputed_score} recovers the expected progeny count, and
#' \code{log(10)*cell_imputed_score} recovers the linear predictor \code{Z_i}.
#' @export
cyfer_finalize <- function(cell_features,
                           cell_lineage,
                           fit_res,
                           lineage_future_count,
                           seed_number = 10){
  if (is.null(rownames(cell_features))) stop("cell_features must have row names (cell IDs)")
  if (is.null(colnames(cell_features))) stop("cell_features must have column names (feature names)")

  # The intercept is added below, so a constant column supplied by the caller
  # would make the design collinear (and one named "Intercept" would duplicate
  # the column name outright).
  constant_bool_vec <- apply(cell_features, 2, function(x){length(unique(x)) == 1})
  if(any(constant_bool_vec)){
    stop("`cell_features` has a constant column (",
         paste0(colnames(cell_features)[constant_bool_vec], collapse = ", "),
         "). The intercept is added internally, so remove it.")
  }

  # Names are kept because `as.character()` drops them and `.lineage_cleanup()`
  # uses them to check row-alignment against `cell_features` in the final refit.
  cell_lineage <- stats::setNames(as.character(cell_lineage), names(cell_lineage))
  if(!is.list(fit_res) || length(fit_res) == 0){
    stop("`fit_res` must be a non-empty list, the output of `cyfer()`")
  }

  # The held-out curves are medianed position-wise and the winning position is
  # read off fold 1's path, so position kk has to mean the same lambda in every
  # fold. 
  path_list <- lapply(fit_res, function(x){x$train_fit$lambda_sequence})
  missing_vec <- which(sapply(path_list, is.null))
  if(length(missing_vec) > 0){
    stop("fold(s) ", paste0(missing_vec, collapse = ", "),
         " carry no `train_fit$lambda_sequence`: `fit_res` is not the output of `cyfer()`")
  }
  mismatch_vec <- which(!sapply(path_list, function(path){
    isTRUE(all.equal(path, path_list[[1]]))
  }))
  if(length(mismatch_vec) > 0){
    stop("fold(s) ", paste0(mismatch_vec, collapse = ", "),
         " were fit on a different lambda sequence than fold 1, so the held-out ",
         "curves are not comparable and their median is meaningless. Refit with a ",
         "single numeric `lambda_initial` passed to `cyfer()`.")
  }
  
  # `sapply()` returns a vector rather than a matrix when the lambda path has
  # length 1, so name the fold margin explicitly.
  test_mat <- sapply(fit_res, function(x){x$test_loglik})
  if(!is.matrix(test_mat)) test_mat <- matrix(test_mat, nrow = 1)
  test_quantile <- apply(test_mat, 1, function(vec){stats::median(vec)})
  lambda_sequence <- fit_res[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(test_quantile)]

  if(!is.null(seed_number)) set.seed(seed_number)
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