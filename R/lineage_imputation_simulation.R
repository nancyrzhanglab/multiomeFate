#' Generate a small dataset from the CYFER model, for tests
#'
#' Draws a lineage
#' centre from \code{N(1, variance_across_lineage * I)}, then \code{n_each} cells
#' around that centre from \code{N(centre, variance_within_lineage * I)}, and
#' finally one Poisson future count per lineage with mean
#' \code{sum_i exp(x_i' beta)} --- that is, exactly the generative model
#' \code{cyfer()} fits, so a fit on this data should recover
#' \code{coefficient_vec}.
#'
#' Two things to know before using it:
#' \itemize{
#'   \item The returned \code{cell_features} has already been through
#'     \code{.lineage_cleanup()}, so it \bold{carries an \code{Intercept}
#'     column}. The user-facing functions add their own, so strip it first:
#'     \code{cell_features[, setdiff(colnames(cell_features), "Intercept"),
#'     drop = FALSE]}.
#'   \item Until 2026-08-06 the count draw evaluated the \emph{outer} product
#'     rather than the inner one, so it did not generate from the CYFER model at
#'     all for \code{p >= 2}. Any characterization number recorded against this
#'     fixture before that date was computed against the wrong data and should be
#'     re-derived rather than trusted.
#' }
#'
#' The true intercept is \code{0} --- the returned \code{coefficient_vec} has a
#' zero prepended --- so no growth offset is built into the data.
#'
#' @param coefficient_vec True feature coefficients on the natural-log scale,
#'   length \code{p}, \bold{excluding} the intercept. Default \code{c(1, 1)}.
#' @param L Number of lineages. Default \code{10}.
#' @param n_each Number of cells per lineage; every lineage gets the same
#'   number. Default \code{5}.
#' @param p Number of features. Must equal \code{length(coefficient_vec)};
#'   \code{NULL} takes it from there. Default \code{2}.
#' @param variance_across_lineage Isotropic variance of the lineage centres
#'   about \code{rep(1, p)}. Larger values separate the lineages. Default
#'   \code{1}.
#' @param variance_within_lineage Isotropic variance of cells about their
#'   lineage centre --- the intra-clonal heterogeneity. Default \code{0.3}.
#' @param seed Seed set at the top of the function. Default \code{10}. Cannot
#'   be \code{NULL}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{cell_features}}{numeric matrix, \code{L * n_each} rows named
#'       \code{"c:1"}..., and \code{p + 1} columns: \code{Intercept} followed by
#'       \code{"p:1"}...}
#'     \item{\code{cell_lineage}}{character vector of lineage names
#'       \code{"lin:1"}..., row-aligned with \code{cell_features}.}
#'     \item{\code{cell_lineage_idx_list}}{list named by lineage of row
#'       positions.}
#'     \item{\code{coefficient_vec}}{the true coefficients with a leading
#'       \code{Intercept = 0}, so it aligns with \code{colnames(cell_features)}
#'       and can be passed straight to \code{.lineage_objective()}.}
#'     \item{\code{lineage_future_count}}{named numeric vector of the Poisson
#'       draws, one per lineage.}
#'   }
#'
#' @noRd
.construct_lineage_data <- function(coefficient_vec = c(1,1),
                                    L = 10,
                                    n_each = 5,
                                    p = 2,
                                    variance_across_lineage = 1,
                                    variance_within_lineage = 0.3,
                                    seed = 10){
  if(is.null(p)) p <- length(coefficient_vec)
  stopifnot(length(coefficient_vec) == p)
  set.seed(seed)
  
  # construct feature matrix
  tmp <- lapply(1:L, function(lineage){
    center <- MASS::mvrnorm(n = 1, mu = rep(1,p), Sigma = variance_across_lineage*diag(p))
    mat <- MASS::mvrnorm(n = n_each, mu = center, Sigma = variance_within_lineage*diag(p))
    if(n_each == 1){
      mat <- matrix(mat, nrow = 1, ncol = p)
    }
    rownames(mat) <- paste0("c:", (lineage-1)*n_each+1:nrow(mat))
    mat
  })
  cell_features <- do.call(rbind, tmp)
  
  # construct lineage
  uniq_lineages <- paste0("lin:", 1:L)
  cell_lineage <- rep(uniq_lineages, each = n_each)
  names(cell_lineage) <- rownames(cell_features)
  
  # construct future lineage counts
  lineage_future_count <- sapply(uniq_lineages, function(lineage){
    idx <- which(cell_lineage == lineage)
    lambda <- sum(sapply(idx, function(i){
      # inner product: the other order conforms the bare length-p vector as p x 1
      # against a 1 x p matrix and silently computes the p x p outer product
      exp(cell_features[i,,drop=F] %*% coefficient_vec)
    }))
    stats::rpois(1, lambda = lambda)
  })
  
  # name things
  colnames(cell_features) <- paste0("p:", 1:p)
  names(coefficient_vec) <- colnames(cell_features)
  names(lineage_future_count) <- uniq_lineages
  
  # do some other coding
  uniq_lineages <- sort(unique(names(lineage_future_count)))
  cell_lineage_idx_list <- lapply(uniq_lineages, function(lineage){
    which(cell_lineage == lineage)
  })
  names(cell_lineage_idx_list) <- uniq_lineages
  
  # cleanup
  res <- .lineage_cleanup(cell_features = cell_features,
                          cell_lineage = cell_lineage,
                          lineage_future_count = lineage_future_count)
  cell_features <- res$cell_features
  cell_lineage <- res$cell_lineage
  lineage_future_count <- res$lineage_future_count
  coefficient_vec <- c(0, coefficient_vec)
  names(coefficient_vec)[1] <- "Intercept"
  
  list(cell_features = cell_features,
       cell_lineage = cell_lineage,
       cell_lineage_idx_list = cell_lineage_idx_list,
       coefficient_vec = coefficient_vec,
       lineage_future_count = lineage_future_count)
}
