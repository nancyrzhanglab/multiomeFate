#' CYFER: Cell Fate via Exponential Regression (cross-validation)
#'
#' Runs K-fold CV over a decreasing sequence of \code{lambda} values produced by
#' \code{lineage_imputation_sequence()}, selecting \code{lambda} by held-out objective.
#' This implements the CYFER method described in Chen, Lin et al.
#'
#' @param cell_features A numeric matrix where each row represents a cell, and each
#' column represents a feature (for instance, the fastTopics scores).
#' Let \code{n} denote the number of cells (rows).
#' Row names (cell IDs) and column names (feature names) are required.
#' Do \bold{not} supply an intercept column, and more generally no constant
#' column: the intercept is added internally (by \code{.lineage_cleanup()} when
#' fitting, and by \code{cyfer_finalize()} when scoring). A constant column you
#' supply yourself is therefore duplicated, making the design collinear --- and a
#' column already named \code{Intercept} produces two columns of that name.
#' \code{cyfer_finalize()} errors on any constant column for this reason;
#' \code{cyfer()} tolerates one, so the error may not appear until the finalize
#' step. It is also conventional to \code{scale()} the features before fitting,
#' both to keep \code{exp()} away from overflow and to make the ridge penalty
#' comparable across features.
#' @param cell_lineage A character or factor vector of length \code{n} where
#' element \code{i} of \code{cell_lineage} denotes which lineage cell \code{i}
#' belongs to. Factors are coerced to character internally, so unused factor
#' levels are harmless.
#' @param lineage_future_count A named numeric vector (where the names are the
#' lineage names that appeared in \code{cell_lineage}) that denotes the
#' number of cells at the future time point for each lineage.
#' @param lambda_initial The initial value of lambda to perform cross-validation on.
#' @param lambda_sequence_length The number of lambdas to perform cross-validation on.
#' The search starts with \code{lambda_initial} and then decays exponentially to 0.
#' @param num_folds Number of folds to do cross-validation on. Default is \code{10}.
#' Must be at least 2 and at most the number of distinct lineages.
#' @param savefile_tmp Filepath to save files to. Default is \code{NULL} (no
#' temporary save files).
#' @param seed_number Seed value for reproducibility reasons. Default is \code{10}.
#' Governs both the fold assignment and the optimizer's random restarts.
#' @param verbose A numeric, where numbers larger than 1 successively request more
#' information to be printed out as the algorithm proceeds.
#'
#' @returns A list with the number of elements corresponding to \code{num_folds}.
#' Each element of this list contains: \code{test_loglik} (the negative log-likelihood
#' on the held-out lineages), \code{train_loglik} (the negative log-likelihood on the trained
#' lineages), and \code{train_fit} (the actual fit, after using the \code{lineage_imputation_sequence()}).
#'
#' \code{test_loglik} and \code{train_loglik} are the \emph{unpenalized} objective
#' evaluated at each lambda along \code{train_fit$lambda_sequence} (lower is
#' better), not literal log-likelihoods. Pass the result to
#' \code{\link{cyfer_finalize}} to select lambda and score the cells; note that
#' the per-cell scores it returns are on the log10 scale while the coefficients
#' here are on the natural-log scale.
#'
#' Lineages named in \code{lineage_future_count} that have no cells in
#' \code{cell_lineage} are dropped before the folds are built, so they neither
#' occupy a fold nor count towards \code{num_folds}.
#' @export
cyfer <- function(cell_features,
                  cell_lineage,
                  lineage_future_count,
                  lambda_initial,
                  lambda_sequence_length,
                  num_folds = 10,
                  savefile_tmp = NULL,
                  seed_number = 10,
                  verbose = 0
){
  if (is.null(rownames(cell_features))) stop("cell_features must have row names (cell IDs)")
  if (is.null(colnames(cell_features))) stop("cell_features must have column names (feature names)")

  # Subsetting a factor retains its unused levels, which would misalign it with
  # the per-fold `lineage_future_count`. Work in character throughout.
  cell_lineage <- as.character(cell_lineage)

  # Drop lineages with no cells at the first time point before the folds are
  # built. `construct_folds()` would otherwise deal such a lineage into a fold
  # where it contributes nothing, and a fold made up entirely of them yields
  # `cv_cell_list[[fold]] == NULL` -- so `cell_features[-NULL,,drop=F]` trains
  # on zero rows.
  lineage_future_count <- lineage_future_count[names(lineage_future_count) %in% unique(cell_lineage)]

  # `construct_folds()` calls sample(), so the seed must be set before it runs
  # for `seed_number` to make the fold assignment reproducible.
  if(!is.null(seed_number)) set.seed(seed_number)

  tmp <- construct_folds(
    cell_lineage = cell_lineage,
    lineage_future_count = lineage_future_count,
    num_folds = num_folds
  )
  cv_cell_list <- tmp$cv_cell_list
  fold_lineage_list <- tmp$fold_lineage_list

  cv_fit_list <- vector("list", length = num_folds)
  names(cv_fit_list) <- names(cv_cell_list)

  for(i in 1:num_folds){
    fold <- names(cv_cell_list)[i]
    if(verbose > 0) print(paste0("Dropping fold #", i, " out of ", num_folds))

    #################

    # training
    cell_features_train <- cell_features[-cv_cell_list[[fold]],,drop = F]
    cell_lineage_train <- cell_lineage[-cv_cell_list[[fold]]]
    lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]

    if(!is.null(seed_number)) set.seed(seed_number)
    train_fit <- lineage_imputation_sequence(
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
      evaluate_loglikelihood(cell_features = cell_features_train,
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
      evaluate_loglikelihood(cell_features = cell_features_test,
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

  structure(cv_fit_list,
            class = "cyfer")
}