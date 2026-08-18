#' Fit lineage imputation at a single lambda
#'
#' @inheritParams cyfer
#' @param coefficient_initial_list A numeric vector or list of numeric vectors
#'   of starting coefficients (names should match feature names; \code{Intercept} added if missing).
#' @param lambda Ridge penalty weight on non-intercept coefficients.
#' @param random_initializations Number of additional random starts.
#' @param upper_randomness Upper cap for random initial coefficients.
#' @param maxit Iteration cap passed to \code{optim()}. \code{NA} (the default)
#'   sets it to \code{max(100, 10*p)}, where \code{p} counts the intercept:
#'   BFGS searches a \code{p}-dimensional space, so the budget has to grow with
#'   it. \code{optim()}'s own default of 100 was silently truncating fits at
#'   realistic feature counts. This is a cap, not a target --- a fit that meets
#'   \code{reltol} stops earlier and costs nothing extra.
#' 
#' @return An object of class \code{"lineage_imputation"} with \code{fit} and \code{res_list}.
#' \code{fit} is the element of \code{res_list} with the smallest
#' \code{objective_val}; \code{res_list} holds one entry per initialization
#' (supplied first, then random).
#'
#' \code{fit$coefficient_vec} is on the \bold{natural-log} scale, so
#' \code{exp(cbind(Intercept = 1, cell_features) \%*\% coefficient_vec)} is the
#' expected number of progeny per cell. Note that \code{cyfer_finalize()} returns
#' its \code{cell_imputed_score} on the \bold{log10} scale instead --- see the
#' "Scales" section of \code{\link{cyfer_finalize}}.
#' @export
lineage_imputation <- function(cell_features,
                               cell_lineage,
                               coefficient_initial_list,
                               lineage_future_count,
                               lambda = 0,
                               random_initializations = 10,
                               upper_randomness = 5,
                               maxit = NA,
                               verbose = 1){
  # do some preliminary formatting
  if(!is.list(coefficient_initial_list)) coefficient_initial_list <- list(coefficient_initial_list)
  list_len <- length(coefficient_initial_list)
  
  # some cleanup
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count,
                          verbose = verbose)
  cell_features <- tmp$cell_features
  cell_lineage <- tmp$cell_lineage
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  uniq_lineages <- tmp$uniq_lineages
  coefficient_initial_list <- .append_intercept_term(coefficient_initial_list)
  p <- ncol(cell_features)

  # BFGS searches a p-dimensional space, so its iteration budget has to grow with
  # p. optim()'s default of 100 truncated a fifth of the fits at p+1 = 61 with
  # L = 100 -- and it truncated them unevenly: the small-lambda fits are the
  # ill-conditioned ones that need the most iterations, so the low end of the CV
  # curve was penalised by an optimizer artefact rather than by generalisation.
  # This is a cap, not a target: a fit meeting `reltol` stops earlier regardless.
  if(length(maxit) != 1 || (!is.na(maxit) && (!is.numeric(maxit) || maxit < 1))){
    stop("`maxit` must be a single positive number, or NA to scale it with the feature count")
  }
  if(is.na(maxit)) maxit <- max(100, 10*p)

  stopifnot(setequal(unique(cell_lineage), names(lineage_future_count)),
            is.matrix(cell_features), nrow(cell_features) == length(cell_lineage),
            all(sapply(coefficient_initial_list, length) == ncol(cell_features)),
            sum(is.na(cell_features)) == 0,
            sum(is.na(cell_lineage)) == 0,
            sum(is.na(lineage_future_count)) == 0)
  for(i in 1:list_len){
    if(length(names(coefficient_initial_list[[i]])) != 0){
      stopifnot(all(names(coefficient_initial_list[[i]]) == colnames(cell_features)))
    } else {
      names(coefficient_initial_list[[i]]) <- colnames(cell_features)
    }
  }
  
  # rearrange arguments
  optim_fn <- function(coefficient_vec,
                       cell_features,
                       cell_lineage,
                       cell_lineage_idx_list,
                       lambda,
                       lineage_future_count){
    # Inside the line search an overflowing trial point is routine: hand `optim`
    # a non-finite value and let it back off, as it did before
    # `.lineage_objective()` grew its guard. A non-finite *start* is a different
    # matter and is reported below, before `optim` is ever called.
    tryCatch(
      .lineage_objective(cell_features = cell_features,
                         cell_lineage = cell_lineage,
                         cell_lineage_idx_list = cell_lineage_idx_list,
                         coefficient_vec = coefficient_vec,
                         lambda = lambda,
                         lineage_future_count = lineage_future_count),
      error = function(e){Inf}
    )
  }
  
  optim_gr <- function(coefficient_vec,
                       cell_features,
                       cell_lineage,
                       cell_lineage_idx_list,
                       lambda,
                       lineage_future_count){
    .lineage_gradient(cell_features = cell_features,
                      cell_lineage = cell_lineage,
                      cell_lineage_idx_list = cell_lineage_idx_list,
                      coefficient_vec = coefficient_vec,
                      lambda = lambda,
                      lineage_future_count = lineage_future_count)
  }
  
  res_list <- vector("list", length = list_len+random_initializations)
  
  for(i in 1:list_len){
    if(verbose > 0) print(paste0("On provided initialization ", i))
    # let `.lineage_objective()`'s diagnostic reach the caller when the supplied
    # start already overflows -- `optim` would otherwise only report
    # "initial value in 'vmmin' is not finite"
    .lineage_objective(cell_features = cell_features,
                       cell_lineage = cell_lineage,
                       cell_lineage_idx_list = cell_lineage_idx_list,
                       coefficient_vec = coefficient_initial_list[[i]],
                       lambda = lambda,
                       lineage_future_count = lineage_future_count)
    res <- stats::optim(
      par = coefficient_initial_list[[i]],
      fn = optim_fn,
      gr = optim_gr,
      method = "BFGS",
      control = list(maxit = maxit),
      cell_features = cell_features,
      cell_lineage = cell_lineage,
      cell_lineage_idx_list = cell_lineage_idx_list,
      lambda = lambda,
      lineage_future_count = lineage_future_count
    )
    
    res_vec <- res$par
    names(res_vec) <- colnames(cell_features)
    
    res_list[[i]] <- list(coefficient_initial = coefficient_initial_list[[i]],
                          coefficient_vec = res_vec,
                          convergence = res$convergence,
                          lambda = lambda,
                          objective_val = res$value)
  }
  
  if(random_initializations > 0){
    max_feature <- stats::quantile(abs(cell_features), probs = 0.95)
    num_cells_per_lineage <- sapply(cell_lineage_idx_list, length)
    names(num_cells_per_lineage) <- uniq_lineages
    max_count_ratio <- max(lineage_future_count[uniq_lineages]/num_cells_per_lineage[uniq_lineages])
    if(max_count_ratio <= 0){
      # log(0) = -Inf would make every draw NaN, and optim would then report
      # "non-finite value supplied by optim" without naming the cause
      max_limit <- 0
      min_value <- -10
    } else {
      # The budget is on the *total* linear predictor: keep |x_i'beta| within
      # about 2*log(max growth ratio) so exp() cannot overflow. Turning that into
      # a per-coefficient allowance needs a model of how p contributions add.
      # With mixed signs they partly cancel and the sum grows like sqrt(p), so
      # the allowance scales as 1/sqrt(p) -- the same reasoning behind
      # Xavier/Glorot and He initialization. The previous 1/p came from a
      # triangle-inequality worst case that effectively never occurs, and it
      # collapsed the range to roughly [0, 0.09] at p = 60: ten restarts inside
      # one small box are one search repeated ten times, at ten times the cost,
      # on a non-convex objective.
      max_limit <- 2*log(max_count_ratio)/(sqrt(p)*max_feature)
      # the threshold, not just the sign: at max_count_ratio == 1 the limit is
      # exactly 0, and a `< 0` test would leave min_value == max_limit, so all
      # the restarts would collapse onto the same zero vector
      min_value <- ifelse(max_limit <= 1e-6, 2*(max_limit-1), 0)
    }
    for(i in 1:random_initializations){
      if(verbose > 0) print(paste0("On random initialization ", i))
      coef_vec <- pmin(stats::runif(p, min = min_value, max = max_limit), upper_randomness)
      names(coef_vec) <- colnames(cell_features)
      
      res <- stats::optim(
        par = coef_vec,
        fn = optim_fn,
        gr = optim_gr,
        method = "BFGS",
        control = list(maxit = maxit),
        cell_features = cell_features,
        cell_lineage = cell_lineage,
        cell_lineage_idx_list = cell_lineage_idx_list,
        lambda = lambda,
        lineage_future_count = lineage_future_count
      )
      
      res_vec <- res$par
      names(res_vec) <- colnames(cell_features)
      
      res_list[[i+list_len]] <- list(coefficient_initial = coef_vec,
                                     coefficient_vec = res_vec,
                                     convergence = res$convergence,
                                     lambda = lambda,
                                     objective_val = res$value)
    }
  }
  
  obj_vec <- sapply(res_list, function(lis){lis$objective_val})
  if(verbose > 1){
    print("Quantile of all the objective scores")
    print(stats::quantile(obj_vec))
  }

  # `optim`'s convergence code was stored on every fit and read nowhere, so a fit
  # that ran out of iterations was indistinguishable from one that converged --
  # and it is the returned coefficient vector, wherever BFGS happened to be when
  # the budget ran out. Report it: the selected fit is the one that propagates,
  # and the count across restarts says whether the budget is systematically tight.
  best_idx <- which.min(obj_vec)
  conv_vec <- sapply(res_list, function(lis){lis$convergence})
  if(res_list[[best_idx]]$convergence != 0){
    warning("the selected fit did not converge (optim code ",
            res_list[[best_idx]]$convergence,
            if(res_list[[best_idx]]$convergence == 1) paste0("; hit maxit = ", maxit) else "",
            "; lambda = ", signif(lambda, 4), "). Its coefficients are where the ",
            "optimizer stopped, not an optimum.")
  } else if(verbose > 0 && any(conv_vec != 0)){
    print(paste0(sum(conv_vec != 0), " of ", length(conv_vec),
                 " initializations did not converge (the selected one did)"))
  }

  structure(list(fit =  res_list[[best_idx]],
                 res_list = res_list),
            class = "lineage_imputation")
}

#' Evaluate the CYFER objective at a given coefficient vector
#'
#' Runs \code{.lineage_cleanup()} and then \code{.lineage_objective()}. This is
#' how \code{cyfer()} scores a fold: fit on the training lineages, then call
#' this on the held-out lineages with \code{lambda = 0} to get an unpenalized
#' held-out value.
#'
#' What comes back is the \emph{negative} penalized log-likelihood, averaged
#' over lineages --- so \bold{lower is better}, and it is not on the scale of any
#' literal log-likelihood. It is the same quantity stored as \code{train_loglik}
#' and \code{test_loglik} by \code{cyfer()}.
#'
#' @inheritParams cyfer
#' @param coefficient_vec A named numeric vector of coefficients on the
#'   \bold{natural-log} scale, including an \code{Intercept} entry, with names
#'   matching \code{colnames(cell_features)} after the intercept is prepended.
#'   Typically \code{fit$coefficient_vec} from a \code{lineage_imputation()}
#'   fit.
#' @param lambda Ridge penalty weight on the non-intercept coefficients. Default
#'   \code{0}, which is what makes the value comparable across the lambda path;
#'   pass a non-zero value only when the penalized objective is wanted.
#'
#' @returns A single numeric. Errors, rather than returning a non-finite value,
#'   if \code{exp(cell_features \%*\% coefficient_vec)} overflows --- see
#'   \code{.lineage_objective()}.
#'
#' @noRd
evaluate_nll <- function(cell_features,
                                   cell_lineage,
                                   coefficient_vec,
                                   lineage_future_count,
                                   lambda = 0){
  
  tmp <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)
  cell_features <- tmp$cell_features
  cell_lineage <- tmp$cell_lineage
  cell_lineage_idx_list <- tmp$cell_lineage_idx_list
  lineage_future_count <- tmp$lineage_future_count
  uniq_lineages <- tmp$uniq_lineages
  
  .lineage_objective(cell_features = cell_features,
                     cell_lineage = cell_lineage,
                     cell_lineage_idx_list = cell_lineage_idx_list,
                     coefficient_vec = coefficient_vec,
                     lambda = lambda,
                     lineage_future_count = lineage_future_count)
}

#################################

#' Align the three inputs and prepend the intercept
#'
#' The single entry point through which every estimation function normalizes its
#' inputs, so that all of them agree on lineage ordering, on the intercept
#' column, and on \code{cell_lineage} being character. Called by
#' \code{lineage_imputation()}, \code{evaluate_nll()},
#' \code{.compute_initial_parameters()}, and the test fixture
#' \code{.construct_lineage_data()}.
#'
#' Three things happen, in this order:
#' \enumerate{
#'   \item \code{cell_lineage} is coerced to character. Indexing a named vector
#'     by a factor uses the factor's integer \emph{codes} rather than its
#'     labels, which is correct by accident on full data and wrong inside a CV
#'     fold. This coercion is the fix for that bug.
#'   \item If the lineages in \code{lineage_future_count} and in
#'     \code{cell_lineage} disagree, both are cut down to the intersection ---
#'     \bold{silently unless \code{verbose > 0}}. Cells whose lineage has no
#'     future count are dropped along with their rows of \code{cell_features}.
#'   \item An \code{Intercept} column of ones is prepended to
#'     \code{cell_features}, unless a column of that name already exists. This
#'     is why callers must not supply one.
#' }
#'
#' Lineages come back in \code{sort()} order, and every returned object is
#' ordered consistently with that.
#'
#' @inheritParams cyfer
#' @param verbose A numeric. At \code{0} (the default here, note, unlike the
#'   user-facing functions) the intersection in step 2 is taken without comment;
#'   above \code{0} it raises a \code{warning()}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{cell_features}}{the matrix with \code{Intercept} as its first
#'       column, and rows for dropped cells removed.}
#'     \item{\code{cell_lineage}}{character vector, row-aligned with
#'       \code{cell_features}.}
#'     \item{\code{cell_lineage_idx_list}}{list named by lineage, each element
#'       the integer row positions of that lineage's cells. Precomputed because
#'       the objective and gradient both need it at every optimizer step.}
#'     \item{\code{lineage_future_count}}{the named vector, subset and reordered
#'       to \code{uniq_lineages}.}
#'     \item{\code{uniq_lineages}}{sorted character vector of the surviving
#'       lineage names.}
#'   }
#'
#' @noRd
.lineage_cleanup <- function(cell_features,
                             cell_lineage,
                             lineage_future_count,
                             verbose = 0){
  stopifnot(length(colnames(cell_features)) == ncol(cell_features))

  # `cell_lineage` is matched to `cell_features` by POSITION, never by name:
  # `which(cell_lineage == lineage)` below returns positions, and those integers
  # index rows of `cell_features` when the objective forms the per-lineage sum.
  # A permuted `cell_lineage` therefore assembles every clone out of the wrong
  # cells and still fits, silently. Check it here, the one place every estimation
  # path funnels through. Names are only available to check if the caller kept
  # them -- `as.character()` drops them, so the entry points coerce name-preserving.
  if(length(cell_lineage) != nrow(cell_features)){
    stop("`cell_lineage` has length ", length(cell_lineage),
         " but `cell_features` has ", nrow(cell_features), " rows")
  }
  if(!is.null(names(cell_lineage))){
    if(is.null(rownames(cell_features))){
      stop("`cell_lineage` is named but `cell_features` has no row names, so the ",
           "two cannot be checked for alignment. Supply row names, or drop the ",
           "names from `cell_lineage` to assert that it is already row-aligned.")
    }
    if(!identical(names(cell_lineage), rownames(cell_features))){
      if(setequal(names(cell_lineage), rownames(cell_features))){
        stop("`cell_lineage` and `cell_features` name the same cells in a ",
             "DIFFERENT ORDER. They are matched by position, so this would fit ",
             "the wrong cells to every lineage. Reorder one to match the other, ",
             "e.g. `cell_lineage <- cell_lineage[rownames(cell_features)]`.")
      }
      stop("`names(cell_lineage)` and `rownames(cell_features)` name different ",
           "cells: ", length(setdiff(names(cell_lineage), rownames(cell_features))),
           " only in `cell_lineage`, ",
           length(setdiff(rownames(cell_features), names(cell_lineage))),
           " only in `cell_features`.")
    }
  }

  # `cell_lineage` is documented as character or factor. Coerce once, here, so
  # that nothing downstream can index a named vector by a factor's integer codes
  # (which silently returns the wrong element, or NA, rather than erroring).
  # This drops the names, which have served their purpose above; everything
  # downstream is positional.
  cell_lineage <- as.character(cell_lineage)

  # `lineage_future_count` is the response. Nothing downstream validates it: the
  # objective just sums `mu - N*log(mu)`, so a negative or non-finite N produces a
  # finite-looking or NaN objective rather than an error, and a duplicated name
  # silently resolves to whichever entry comes first. Check it here, where every
  # estimation path passes, and where each fold's subset is checked as well.
  # These are all element-wise properties, so a subset of a valid vector is valid.
  if(!is.numeric(lineage_future_count)){
    stop("`lineage_future_count` must be numeric, not ", class(lineage_future_count)[1])
  }
  if(is.null(names(lineage_future_count))){
    stop("`lineage_future_count` must be named, with the lineage IDs as names")
  }
  if(anyNA(names(lineage_future_count)) || any(names(lineage_future_count) == "")){
    stop("`lineage_future_count` has missing or empty lineage names")
  }
  if(anyDuplicated(names(lineage_future_count)) != 0){
    dup <- unique(names(lineage_future_count)[duplicated(names(lineage_future_count))])
    stop("`lineage_future_count` has duplicated lineage names (",
         paste0(utils::head(dup, 5), collapse = ", "),
         if(length(dup) > 5) ", ..." else "",
         "). Indexing by name would silently take the first of each.")
  }
  bad <- which(!is.finite(lineage_future_count))
  if(length(bad) > 0){
    stop("`lineage_future_count` has ", length(bad), " non-finite value(s) (NA/NaN/Inf), ",
         "e.g. lineage ", names(lineage_future_count)[bad[1]], ". ",
         "An Inf here surfaces later as an overflow blamed on `cell_features`.")
  }
  bad <- which(lineage_future_count < 0)
  if(length(bad) > 0){
    stop("`lineage_future_count` has ", length(bad), " negative value(s), ",
         "e.g. lineage ", names(lineage_future_count)[bad[1]], " = ",
         lineage_future_count[bad[1]], ". These are counts of cells at the ",
         "future time point.")
  }

  # some cleanup
  if(!setequal(names(lineage_future_count), unique(cell_lineage))){
    if(verbose > 0) warning("Lineages in `lineage_future_count` are not the same as those in `cell_lineage`")
    
    uniq_lineages <- sort(intersect(unique(names(lineage_future_count)), unique(cell_lineage)))
    lineage_future_count <- lineage_future_count[names(lineage_future_count) %in% uniq_lineages]
    rm_cell_idx <- which(!cell_lineage %in% uniq_lineages)
    if(length(rm_cell_idx) > 0){
      cell_lineage <- cell_lineage[-rm_cell_idx]
      cell_features <- cell_features[-rm_cell_idx,,drop=F]
    }
  }
  
  # reorganize everything to be in the same order
  uniq_lineages <- sort(unique(names(lineage_future_count)))
  lineage_future_count <- lineage_future_count[uniq_lineages]
  cell_lineage_idx_list <- lapply(uniq_lineages, function(lineage){
    which(cell_lineage == lineage)
  })
  names(cell_lineage_idx_list) <- uniq_lineages
  
  
  # add intercept to cell_features
  if(!"Intercept" %in% colnames(cell_features)){
    cell_features <- cbind(1, cell_features)
    colnames(cell_features)[1] <- "Intercept"
  }
  
  list(cell_features = cell_features,
       cell_lineage = cell_lineage,
       cell_lineage_idx_list = cell_lineage_idx_list,
       lineage_future_count = lineage_future_count,
       uniq_lineages = uniq_lineages)
}

#' Prepend a zero intercept to each starting coefficient vector
#'
#' The counterpart of the intercept column \code{.lineage_cleanup()} adds to
#' \code{cell_features}: a caller who supplies starting coefficients for the
#' features alone would otherwise hand \code{optim()} a vector one shorter than
#' the design matrix. The starting intercept is \code{0}, i.e. one expected
#' progeny per cell before any feature contribution.
#'
#' A vector that already carries an \code{Intercept} name is left alone, so this
#' is safe to apply to a warm start taken from a previous fit.
#'
#' @param coefficient_initial_list A list of named numeric vectors. Must already
#'   be a list --- \code{lineage_imputation()} wraps a bare vector before
#'   calling.
#'
#' @returns The same list, each element having \code{Intercept} as its first
#'   entry.
#'
#' @noRd
.append_intercept_term <- function(coefficient_initial_list){
  stopifnot(is.list(coefficient_initial_list))
  
  for(i in 1:length(coefficient_initial_list)){
    if(!"Intercept" %in% names(coefficient_initial_list[[i]])){
      coefficient_initial_list[[i]] <- c(0, coefficient_initial_list[[i]])
      names(coefficient_initial_list[[i]])[1] <- "Intercept"
    }
  }
  
  coefficient_initial_list
}

#' The penalized CYFER objective
#'
#' Computes
#' \deqn{\frac{1}{L}\sum_{\ell=1}^{L}\left[\Big(\sum_{i \in \ell} e^{X_{i,\cdot}^\top \beta}\Big) - y_\ell \log\Big(\sum_{i \in \ell} e^{X_{i,\cdot}^\top \beta}\Big)\right] + \lambda\|\beta_{-0}\|_2^2}{(1/L) * sum_l [ (sum_{i in l} exp(x_i'beta)) - y_l * log(sum_{i in l} exp(x_i'beta)) ] + lambda * ||beta_{-0}||^2}
#' the negative Poisson log-likelihood of the per-lineage future counts, dropping
#' terms free of \eqn{\beta}, averaged over lineages, plus a ridge penalty. The
#' intercept is excluded from the penalty. \bold{Lower is better}; this is a
#' minimization objective, and is non-convex in \eqn{\beta}, which is why
#' \code{lineage_imputation()} uses random restarts.
#'
#' Dividing by the number of lineages is what makes values comparable between
#' the training and held-out folds, which contain different numbers of lineages.
#' It also means \code{lambda} is on a per-lineage scale.
#'
#' @param cell_features Numeric matrix with the \code{Intercept} column already
#'   present, rows = cells.
#' @param cell_lineage Character vector, row-aligned with \code{cell_features}.
#'   Accepted for signature symmetry with \code{.lineage_gradient()}; not used.
#' @param cell_lineage_idx_list List named by lineage, giving each lineage's row
#'   positions, as built by \code{.lineage_cleanup()}.
#' @param coefficient_vec Named numeric vector, aligned with
#'   \code{colnames(cell_features)}, on the natural-log scale.
#' @param lambda Ridge penalty weight.
#' @param lineage_future_count Named numeric vector of future counts. Its
#'   \emph{names determine the lineage ordering} used here.
#'
#' @returns A single numeric.
#'
#'   Errors on a non-finite result rather than returning it. That guard is
#'   deliberate and load-bearing in two places: \code{optim()} evaluates
#'   \code{fn} before \code{gr}, so without it a caller with unscaled features
#'   sees \code{"initial value in 'vmmin' is not finite"} instead of an
#'   actionable message; and a non-finite \code{test_loglik} is silently skipped
#'   by \code{which.min()} at lambda-selection time. \code{lineage_imputation()}
#'   wraps this in a \code{tryCatch} for the \emph{trial} points \code{optim()}
#'   probes during its line search, where overflow is routine and backing off is
#'   correct --- but evaluates supplied \emph{starting} points through this
#'   function directly, so the diagnostic reaches the caller. Removing either
#'   half breaks something.
#'
#' @noRd
.lineage_objective <- function(cell_features,
                               cell_lineage,
                               cell_lineage_idx_list,
                               coefficient_vec,
                               lambda,
                               lineage_future_count){
  uniq_lineages <- names(lineage_future_count)
  num_lineages <- length(uniq_lineages)
  cell_names <- rownames(cell_features)
  
  scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec))
  names(scalar1) <- cell_names
  scalar2 <- sapply(uniq_lineages, function(lineage){
    log(sum(scalar1[cell_lineage_idx_list[[lineage]]]))
  })
  
  idx_notintercept <- which(names(coefficient_vec) != "Intercept")
  scalar3 <- .l2norm(coefficient_vec[idx_notintercept])

  res <- (sum(scalar1) - sum(lineage_future_count*scalar2))/num_lineages + lambda*scalar3^2

  # `optim` evaluates fn before gr, so without this guard the gradient's
  # diagnostic below can never reach a caller whose features overflow: they see
  # "initial value in 'vmmin' is not finite" instead. A non-finite value also
  # travels silently into `test_loglik`, where `which.min()` skips it.
  if(!is.finite(res)){
    stop("`.lineage_objective()` produced a non-finite value: exp(cell_features %*% coefficient_vec) ",
         "overflowed or underflowed. Scale `cell_features`, or use a larger lambda.")
  }

  res
}

#' Analytical gradient of the penalized CYFER objective
#'
#' The exact gradient of \code{.lineage_objective()} with respect to
#' \code{coefficient_vec}, handed to \code{optim()} as \code{gr}. Writing it out
#' rather than letting BFGS difference the objective matters here: the objective
#' is non-convex and evaluated at every one of the random restarts, so a
#' finite-difference gradient would multiply the cost by \code{p + 1}.
#'
#' Each cell contributes
#' \code{exp(x_i'beta) * (1 - y_{l(i)} / sum_{j in l(i)} exp(x_j'beta))} to a
#' weight, which is then applied to its feature row. The intercept entry is the
#' plain sum of those weights; the feature entries add \code{2*lambda*beta}.
#'
#' @param cell_features Numeric matrix with the \code{Intercept} column already
#'   present, rows = cells.
#' @param cell_lineage Character vector, row-aligned with \code{cell_features}.
#'   Used --- unlike in \code{.lineage_objective()} --- to broadcast per-lineage
#'   quantities back out to cells. \bold{Must be character}: indexing
#'   \code{lineage_future_count} or \code{denom_vec} by a factor uses the
#'   factor's integer codes, not its labels. Coerced defensively here as well as
#'   in \code{.lineage_cleanup()}.
#' @param cell_lineage_idx_list List named by lineage, giving each lineage's row
#'   positions. Its \emph{names determine the lineage ordering} used here ---
#'   note this differs from \code{.lineage_objective()}, which orders by
#'   \code{names(lineage_future_count)}.
#' @param coefficient_vec Named numeric vector, aligned with
#'   \code{colnames(cell_features)}, on the natural-log scale.
#' @param lambda Ridge penalty weight.
#' @param lineage_future_count Named numeric vector of future counts.
#'
#' @returns A named numeric vector the same length as \code{coefficient_vec} and
#'   in the same order, with \code{Intercept} first.
#'
#'   Errors rather than returning \code{NA}. An \code{NA} gradient is invisible
#'   to \code{optim()}: BFGS returns its starting value and still reports
#'   \code{convergence = 0}, so a fit that never moved looks like a fit that
#'   converged immediately. The two causes are distinguished in the message ---
#'   misaligned \code{cell_lineage} against \code{lineage_future_count}, versus
#'   \code{exp()} overflow from unscaled features.
#'
#' @noRd
.lineage_gradient <- function(cell_features,
                              cell_lineage,
                              cell_lineage_idx_list,
                              coefficient_vec,
                              lambda,
                              lineage_future_count){
  stopifnot(colnames(cell_features) == names(coefficient_vec))

  uniq_lineages <- names(cell_lineage_idx_list)
  num_lineages <- length(uniq_lineages)
  cell_names <- rownames(cell_features)
  # as.character() is required: indexing a named vector by a factor uses the
  # factor's integer codes, not its labels.
  cell_lineage <- as.character(cell_lineage)
  lineage_future_count_full <- lineage_future_count[cell_lineage]
  
  # keep track of colnames(cell_features) that is not the intercept
  colname_vec <- colnames(cell_features)
  colname_vec <- colname_vec[colname_vec != "Intercept"]
  
  # construct scalar_vec, which is a vector with length of nrow(cell_features)
  scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec)) 
  names(scalar1) <- cell_names
  scalar2a <- lineage_future_count_full * scalar1 
  denom_vec <- sapply(uniq_lineages, function(lineage){
    sum(scalar1[cell_lineage_idx_list[[lineage]]]) 
  })
  names(denom_vec) <- uniq_lineages
  scalar2b <- denom_vec[cell_lineage]
  scalar_vec <- scalar1 - scalar2a/scalar2b

  # An NA gradient is not recoverable by BFGS: optim() silently returns its
  # starting value while still reporting convergence = 0. Fail loudly instead.
  # Two distinct causes reach here, so name the right one.
  if(anyNA(scalar_vec)){
    if(!all(cell_lineage %in% names(lineage_future_count))){
      stop("`.lineage_gradient()` produced NA: `cell_lineage` and `lineage_future_count` are misaligned")
    }
    stop("`.lineage_gradient()` produced a non-finite value: exp(cell_features %*% coefficient_vec) ",
         "overflowed or underflowed. Scale `cell_features`, or use a larger lambda.")
  }
  
  # gradient of the intercept
  res1 <- sum(scalar_vec)/num_lineages
  
  # gradient of the other terms
  weighted_features <- sweep(cell_features[,colname_vec,drop=F], 
                             MARGIN = 1, 
                             STATS = scalar_vec, 
                             FUN = "*")
  res2 <- Matrix::colSums(weighted_features)/num_lineages + 2*lambda*coefficient_vec[colname_vec]
  
  res <- c(res1, res2)
  names(res) <- colnames(cell_features)
  
  res
}

#' Euclidean norm
#'
#' Note this is the norm itself, not its square --- callers that want the ridge
#' penalty square the result (\code{.l2norm(...)^2}).
#'
#' @param x A numeric vector.
#'
#' @returns A single numeric.
#'
#' @noRd
.l2norm <- function(x){sqrt(sum(x^2))}