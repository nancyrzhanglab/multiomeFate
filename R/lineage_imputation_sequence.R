#' Fit lineage imputation along a lambda path
#' 
#' This function calls \code{lineage_imputation()} for a sequence of lambdas.
#'
#' @inheritParams cyfer
#' @param lambda_min,lambda_max Floor and cap applied to the internal
#'   data-driven \code{lambda_initial} heuristic. Only used when
#'   \code{lambda_initial} is \code{NA}.
#' @param multipler Scaling factor for the internal \code{lambda_initial} heuristic.
#' 
#' @return A list with \code{fit_list} (solution estimated by \code{lineage_imputation()} per lambda) and \code{lambda_sequence}.
#' \code{lambda_sequence} starts at \code{lambda_initial} and decays
#' exponentially to \code{0}, so it is strictly \emph{decreasing}; \code{fit_list[[i]]}
#' is the fit at \code{lambda_sequence[i]}, warm-started from \code{fit_list[[i-1]]}.
#' Each \code{coefficient_vec} is on the \bold{natural-log} scale --- see the
#' "Scales" section of \code{\link{cyfer_finalize}}.
#' @examples
#' \donttest{
#' data(priming_simulation)
#' set.seed(10)
#' path <- lineage_imputation_sequence(
#'   cell_features = priming_simulation$cell_features,
#'   cell_lineage = priming_simulation$cell_lineage,
#'   lineage_future_count = priming_simulation$lineage_future_count,
#'   lambda_initial = 3,
#'   lambda_sequence_length = 4,
#'   verbose = 0)
#' path$lambda_sequence                       # decreasing, ends at 0
#' sapply(path$fit_list, function(fit){fit$objective_val})
#' }
#' @export
lineage_imputation_sequence <- function(cell_features,
                                        cell_lineage,
                                        lineage_future_count,
                                        lambda_initial = NA,
                                        lambda_max = 101, # only for controlling the initial lambda
                                        lambda_min = 0.01, # only for controlling the initial lambda
                                        lambda_sequence_length = 50,
                                        multipler = 1e4,
                                        verbose = 1){
  res <- .compute_initial_parameters(cell_features = cell_features,
                                     cell_lineage = cell_lineage,
                                     lineage_future_count = lineage_future_count,
                                     lambda_max = lambda_max,
                                     lambda_min = lambda_min,
                                     multipler = multipler)
  coefficient_initial <- res$coefficient_initial
  if(is.na(lambda_initial)) lambda_initial <- res$lambda_initial
  
  lambda_sequence <- exp(seq(log(lambda_initial+1), 0, 
                             length.out = lambda_sequence_length))-1
  fit_list <- vector("list", length = lambda_sequence_length)
  
  for(i in 1:lambda_sequence_length){
    if(verbose > 1) {
      if(verbose > 0 && lambda_sequence_length > 10 && i %% floor(lambda_sequence_length/10) == 0){
        cat('*')
      } else {
        print(paste0("Working on lambda in sequence ", i, " out of ", 
                     lambda_sequence_length))
      }
    }
    if(i == 1){
      coefficient_vec <- coefficient_initial
    } else {
      coefficient_vec <- fit_list[[i-1]]$coefficient_vec
    }
    
    tmp <- lineage_imputation(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              coefficient_initial_list = coefficient_vec,
                              lineage_future_count = lineage_future_count,
                              lambda = lambda_sequence[i],
                              random_initializations = 10,
                              upper_randomness = 5,
                              verbose = verbose-2)
    fit_list[[i]] <- tmp$fit
  }
  
  list(fit_list = fit_list,
       lambda_sequence = lambda_sequence)
}


#' Data-driven starting point for the lambda path and the coefficients
#'
#' Produces the two things \code{lineage_imputation_sequence()} needs before it
#' can start: where to begin the lambda path, and what coefficient vector to
#' start the first fit from.
#'
#' The coefficient start is the interpretable half. Everything but the intercept
#' starts at zero, and the intercept starts at
#' \code{log((future_total + 1) / (current_total + 1))} --- the log of the
#' overall growth ratio. That is the exact maximizer if no feature carried any
#' information, so the path begins at the null model and the features have to
#' earn their coefficients.
#'
#' \code{lambda_initial} is a heuristic, not an estimate: it scales the gap
#' between the null-model objective and a per-lineage term by \code{multipler},
#' then clamps to \code{[lambda_min, lambda_max]}. \bold{On typical data it
#' saturates at the \code{lambda_max} cap}, so any test of its \emph{value} must
#' widen the clamp first (\code{lambda_min = 0}, \code{lambda_max = 1e12},
#' \code{multipler = 1}) or it will pass against any implementation whatsoever.
#'
#' The \code{+1} smoothing in both \code{log_growth_ratio} and
#' \code{log1p(lineage_current_count)} guards the all-zero
#' \code{lineage_future_count} case, which would otherwise give \code{log(0)}
#' and carry \code{-Inf} / \code{NaN} silently into \code{optim()}. It shifts the
#' result on \emph{every} input, not only the degenerate one.
#'
#' @inheritParams cyfer
#' @param lambda_max Cap on the returned \code{lambda_initial}. Default
#'   \code{101}.
#' @param lambda_min Floor on the returned \code{lambda_initial}. Default
#'   \code{0.01}. Must not exceed \code{lambda_max}; asserted.
#' @param multipler Scaling factor applied before the clamp. Default
#'   \code{1e4}, matching its only caller. (Spelling is the existing one.)
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{coefficient_initial}}{named numeric vector of length
#'       \code{ncol(cell_features) + 1}, zero everywhere except
#'       \code{Intercept}. On the natural-log scale.}
#'     \item{\code{lambda_initial}}{a single numeric in
#'       \code{[lambda_min, lambda_max]}.}
#'   }
#'
#' @noRd
# cell_features simply included for convenience
.compute_initial_parameters <- function(cell_features,
                                        cell_lineage,
                                        lineage_future_count,
                                        lambda_max = 101,
                                        lambda_min = 0.01,
                                        multipler = 1e4){
  
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)
  cell_features <- tmp$cell_features
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  num_lineages <- length(lineage_future_count)
  
  lineage_current_count <- sapply(cell_lineage_idx_list, length)
  stopifnot(all(names(lineage_current_count) == names(lineage_future_count)))
  
  future_total <- sum(lineage_future_count)
  current_total <- sum(lineage_current_count)
  
  # +1 smoothing throughout: an all-zero `lineage_future_count` would otherwise
  # give log(0) = -Inf here and 0*Inf = NaN in `lambda_initial`, both of which
  # travel silently into `optim`
  log_growth_ratio <- log((future_total+1)/(current_total+1))
  term1 <- future_total*(1-log_growth_ratio)
  term2 <- sum(lineage_future_count * log1p(lineage_current_count))

  lambda_initial <- -multipler*(term1 - term2)/num_lineages
  # floor at lambda_min, cap at lambda_max (these were swapped, and both
  # defaulted to 101, which collapsed the clamp to the constant 101)
  stopifnot(lambda_min <= lambda_max)
  lambda_initial <- min(max(lambda_initial, lambda_min), lambda_max)
  
  coefficient_initial <- rep(0, ncol(cell_features))
  names(coefficient_initial) <- colnames(cell_features)
  stopifnot("Intercept" %in% names(coefficient_initial))
  coefficient_initial["Intercept"] <- log_growth_ratio
  
  list(coefficient_initial = coefficient_initial,
       lambda_initial = lambda_initial)
}


