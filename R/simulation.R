#' Simulate lineages under pure selection ("priming")
#'
#' Generates the \bold{priming} regime: a cell's fate potential is fixed by its
#' position in the embedding, and lineages are spatially localized clumps of
#' cells. Because cells near each other share both their embedding and their
#' potential, a lineage's cells are alike and the between-lineage variation in
#' expansion is large --- the signature of selection acting on pre-existing
#' state. Contrast \code{generate_simulation_plastic()}, which builds lineages
#' that are deliberately \emph{heterogeneous} in potential.
#'
#' This is how \code{data/priming_simulation.rda} was produced; see
#' \code{?priming_simulation}.
#'
#' The construction, in the order the \code{verbose} messages report it:
#' \enumerate{
#'   \item \code{stats::kmeans()} with \code{2 * num_lineages} centres, of which
#'     one random cell from each of the first \code{num_lineages} clusters
#'     becomes a lineage seed. Over-clustering makes the seeds tighter and more
#'     spread out than \code{num_lineages} centres would.
#'   \item Each seed anchors an isotropic Gaussian whose per-dimension variance
#'     is \code{lineage_spread} times the observed variance of that embedding
#'     dimension.
#'   \item Every cell gets a posterior over lineages from those Gaussians and
#'     \code{lineage_prior}, and is \bold{sampled} from it --- so lineages are
#'     soft, overlapping clumps, and their realized sizes are random rather than
#'     equal to \code{n * lineage_prior}.
#'   \item Each cell's expected progeny count is
#'     \code{exp(coefficient_intercept + x_i' beta)}, optionally Poisson-drawn,
#'     and summed within lineage to give the future size.
#' }
#'
#' \bold{Stochastic with no \code{seed_number} argument} --- steps 1, 3, and 4
#' all draw. Set a seed before calling.
#'
#' @param embedding_mat Numeric matrix of the current time point's cells, rows =
#'   cells and columns = embedding dimensions (at least 2). Row names, if
#'   present, are carried onto the returned vectors. Supplied by the caller
#'   rather than simulated, so that simulations sit on a real embedding.
#' @param bool_add_randomness Whether to draw the realized progeny counts from
#'   \code{stats::rpois()} around the expected counts. Default \code{TRUE}.
#'   \code{FALSE} gives the noiseless best case, useful for isolating estimation
#'   error from sampling noise.
#' @param coefficient_intercept The intercept, on the natural-log scale. Sets
#'   the overall growth level; \code{0} means one expected progeny per cell
#'   before any feature contribution. Default \code{0}.
#' @param embedding_coefficient_vec True coefficients on the embedding, length
#'   \code{ncol(embedding_mat)}, natural-log scale. Default all ones.
#' @param fatefeatures_coefficient_vec True coefficients on the extra fate
#'   features, length \code{ncol(fatefeatures_mat)}. Default \code{NULL}.
#' @param fatefeatures_mat Optional numeric matrix of features that drive fate
#'   but are \bold{not} part of the embedding, same rows as
#'   \code{embedding_mat}. This is how a simulation withholds signal from the
#'   estimator: CYFER fitted on \code{embedding_mat} alone cannot recover this
#'   contribution, which is the point when studying misspecification. Default
#'   \code{NULL}.
#' @param lineage_spread Multiplier on each embedding dimension's variance when
#'   forming the lineage Gaussians. Default \code{1}. Larger values make
#'   lineages broader and more overlapping, weakening the selection signal.
#' @param lineage_prior Numeric vector of length \code{num_lineages} of prior
#'   lineage probabilities; renormalized to sum to 1. Default \code{NA}, meaning
#'   uniform. Any names it carries are overwritten (with a warning) by
#'   \code{"lineage:1"}...
#' @param num_lineages Number of lineages. Default \code{10}.
#' @param tol Tolerance for the "prior sums to 1" check. Default \code{1e-06}.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns An object of class \code{"multiomeFate_simulation_vanilla"}, a list
#'   with:
#'   \describe{
#'     \item{\code{cell_fate_potential}}{named numeric, \code{log10(realized
#'       progeny + 1)} per cell. The \code{+1} keeps a zero-progeny cell finite,
#'       so this is \emph{not} directly comparable to a
#'       \code{cell_imputed_score}.}
#'     \item{\code{cell_fate_potential_truth}}{named numeric,
#'       \code{log10(expected progeny)} per cell, with no \code{+1}. \bold{This}
#'       is the ground truth to score \code{cell_imputed_score} against.}
#'     \item{\code{coefficient_intercept}}{as supplied.}
#'     \item{\code{embedding_mat}}{as supplied.}
#'     \item{\code{fatefeatures_coefficient_vec}, \code{fatefeatures_mat}}{as
#'       supplied.}
#'     \item{\code{gaussian_list}}{list of \code{num_lineages} objects of class
#'       \code{"gaussian"}, each with \code{mean} and \code{cov}.}
#'     \item{\code{lineage_assignment}}{factor of length \code{nrow(embedding_mat)}
#'       with levels \code{"lineage:1"}..., named by cell. Pass
#'       \code{as.character()} of this as \code{cell_lineage}.}
#'     \item{\code{lineage_future_size}}{named numeric, the future count per
#'       lineage. This is \code{lineage_future_count}.}
#'     \item{\code{prob_mat}}{cells-by-lineages posterior matrix used for the
#'       assignment draw.}
#'     \item{\code{summary_mat}}{5-by-\code{num_lineages} matrix with rows
#'       \code{mean}, \code{median}, \code{sd}, \code{range} of the true log10
#'       potential within each lineage, and \code{future_size}. The \code{sd}
#'       and \code{range} rows are the intra-clonal heterogeneity, and are what
#'       distinguish this regime from the plastic one.}
#'   }
#'
#' @export
generate_simulation <- function(embedding_mat,
                                bool_add_randomness = TRUE, 
                                coefficient_intercept = 0, 
                                embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
                                fatefeatures_coefficient_vec = NULL,
                                fatefeatures_mat = NULL, 
                                lineage_spread = 1, 
                                lineage_prior = NA, 
                                num_lineages = 10, 
                                tol = 1e-06, 
                                verbose = 0) 
{
  if (all(is.na(lineage_prior))) 
    lineage_prior <- rep(1/num_lineages, length = num_lineages)
  lineage_prior <- lineage_prior/sum(lineage_prior)
  K <- num_lineages
  n <- nrow(embedding_mat)
  d <- ncol(embedding_mat)
  if(all(!is.null(fatefeatures_mat))){
    d2 <- ncol(fatefeatures_mat)
    stopifnot(nrow(fatefeatures_mat) == nrow(embedding_mat),
              length(fatefeatures_coefficient_vec) == d2)
  } else {
    d2 <- 0
  }
  
  rho <- lineage_spread
  if (length(names(lineage_prior)) > 0) {
    warning("Overwriting names in lineage_prior")
  } 
  
  names(lineage_prior) <- paste0("lineage:", 1:K)
  stopifnot(d > 1, 
            length(embedding_coefficient_vec) == d, 
            length(lineage_prior) == K, 
            all(lineage_prior >= 0), 
            abs(sum(lineage_prior) - 1) <= tol, 
            rho >= 0)
  
  if (verbose > 0) 
    print("Step 1: Selecting seed cells")
  cluster_idxs <- .kmeans_seed(embedding_mat = embedding_mat, 
                               K = K)
  if (verbose > 0) 
    print("Step 2: Computing Gaussian distributions")
  gaussian_list <- .form_gaussian_distributions(cluster_idxs = cluster_idxs, 
                                                embedding_mat = embedding_mat, 
                                                rho = rho)
  if (verbose > 0) 
    print("Step 3: Computing posterior distributions")
  prob_mat <- .compute_posteriors(embedding_mat = embedding_mat, 
                                  gaussian_list = gaussian_list, 
                                  lineage_prior = lineage_prior, 
                                  verbose = verbose - 1)
  if (verbose > 0) 
    print("Step 4: Sampling lineages")
  lineage_assignment <- sapply(1:n, function(i) {
    sample(1:K, size = 1, prob = prob_mat[i, ])
  })
  lineage_assignment <- factor(paste0("lineage:", lineage_assignment), 
                               levels = colnames(prob_mat))
  if (length(rownames(embedding_mat)) > 0) 
    names(lineage_assignment) <- rownames(embedding_mat)
  
  if (verbose > 0) 
    print("Step 5: Computing future lineage size")
  cell_contribution <- as.numeric(embedding_mat %*% embedding_coefficient_vec)
  if(d2 > 0){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution_truth <- exp(cell_contribution + coefficient_intercept)
  if (length(rownames(embedding_mat)) > 0) 
    names(cell_contribution_truth) <- rownames(embedding_mat)
  cell_contribution_random <- cell_contribution_truth
  
  if (bool_add_randomness) {
    if (verbose > 0) 
      print("Step 5b: (Optional) Adding randomness")
    cell_contribution_random <- stats::rpois(n = length(cell_contribution_random), 
                                             lambda = cell_contribution_random)
    if (length(rownames(embedding_mat)) > 0) 
      names(cell_contribution_random) <- rownames(embedding_mat)
  }
  
  lineage_future_size <- sapply(levels(lineage_assignment), 
                                function(lev) {
                                  idx <- which(lineage_assignment == lev)
                                  round(sum(cell_contribution_random[idx]))
                                })
  names(lineage_future_size) <- levels(lineage_assignment)
  
  summary_mat <- .compute_summary_lineages(cell_fate_potential_truth = log10(cell_contribution_truth),
                                           lineage_assignment = lineage_assignment,
                                           lineage_future_size = lineage_future_size)
  summary_mat
  
  if (verbose > 0) 
    print("Step 6: Outputting")
  return(
    structure(list(cell_fate_potential = log10(cell_contribution_random + 1), 
                   cell_fate_potential_truth = log10(cell_contribution_truth), 
                   coefficient_intercept = coefficient_intercept,
                   embedding_mat = embedding_mat,
                   fatefeatures_coefficient_vec = fatefeatures_coefficient_vec, 
                   fatefeatures_mat = fatefeatures_mat,
                   gaussian_list = gaussian_list, 
                   lineage_assignment = lineage_assignment, 
                   lineage_future_size = lineage_future_size, 
                   prob_mat = prob_mat,
                   summary_mat = summary_mat),
              class = "multiomeFate_simulation_vanilla")
  )
}


#############

#' Pick one seed cell per lineage
#'
#' Runs k-means with \code{2 * K} centres and takes one random cell from each of
#' clusters \code{1..K}. The over-clustering is deliberate: with exactly \code{K}
#' centres the seeds would be near the cluster centroids and roughly evenly
#' spaced, whereas doubling the centres and using half of them gives seeds that
#' are tighter and less regularly arranged, which is closer to how real clones
#' sit in an embedding. Which half is kept depends on k-means' arbitrary cluster
#' numbering, so it is effectively a random half.
#'
#' \bold{Stochastic}: both the k-means initialization and the within-cluster
#' draw. \code{stats::kmeans()} warnings (typically non-convergence) are
#' suppressed.
#'
#' @param embedding_mat Numeric matrix, rows = cells.
#' @param K Number of lineages, i.e. of seeds to return.
#'
#' @returns An integer vector of length \code{K} of row positions into
#'   \code{embedding_mat}, named \code{"lineage:1"}...
#'
#' @noRd
.kmeans_seed <- function(
    embedding_mat,
    K
){
  kmeans_res <- suppressWarnings(stats::kmeans(embedding_mat, centers = 2*K))
  cluster_idxs <- sapply(1:K, function(k){
    sample(which(kmeans_res$cluster == k), 1)
  })
  names(cluster_idxs) <- paste0("lineage:", 1:K)
  
  return(cluster_idxs)
}

#' Construct a Gaussian distribution object
#'
#' A two-slot container, so that \code{gaussian_list} in the simulation output
#' is self-describing rather than a bare pair of matrices.
#'
#' @param cov_mat Square covariance matrix.
#' @param mean_vec Mean vector, length \code{nrow(cov_mat)}.
#'
#' @returns A list of class \code{"gaussian"} with elements \code{cov} and
#'   \code{mean}. Note the element names drop the type suffixes the arguments
#'   carry.
#'
#' @noRd
.gaussian <- function(cov_mat,
                      mean_vec){
  stopifnot(nrow(cov_mat) == ncol(cov_mat),
            nrow(cov_mat) == length(mean_vec))
  
  return(structure(list(cov = cov_mat, 
                        mean = mean_vec),
                   class = "gaussian"))
}

#' Build one lineage's Gaussian around its seed cell
#'
#' The mean is the seed cell's own embedding coordinates; the covariance is
#' diagonal with entry \code{rho * sd_j^2}, where \code{sd_j} is the standard
#' deviation of embedding dimension \code{j} across \emph{all} cells. So
#' \code{rho == 1} makes each lineage as wide as the whole dataset in every
#' direction, and lineages are localized only for \code{rho} well below 1.
#'
#' Note the covariance is built from squared standard deviations, i.e. actual
#' variances --- unlike \code{.compute_previous_to_future_mapping()} in
#' \code{R/simulation_attach-future.R}, which passes unsquared standard
#' deviations.
#'
#' @param cluster_idx A single row position into \code{embedding_mat}, the seed
#'   cell.
#' @param embedding_mat Numeric matrix, rows = cells.
#' @param rho Variance multiplier; the \code{lineage_spread} argument of
#'   \code{generate_simulation()}.
#'
#' @returns An object of class \code{"gaussian"}.
#'
#' @noRd
.form_gaussian_distribution <- function(
    cluster_idx,
    embedding_mat,
    rho
){
  stopifnot(length(cluster_idx) == 1,
            cluster_idx <= nrow(embedding_mat),
            cluster_idx > 0,
            cluster_idx %% 1 == 0)
  
  mean_vec <- embedding_mat[cluster_idx,]
  sd_vec <- apply(embedding_mat, 2, stats::sd)
  
  cov_mat <- diag(rho*sd_vec^2)
  
  return(.gaussian(cov_mat = cov_mat, 
                   mean_vec = mean_vec))
}

#' Build every lineage's Gaussian
#'
#' Maps \code{.form_gaussian_distribution()} over the seeds. All lineages share
#' the same \code{rho}, so they differ only in location, not in shape.
#'
#' @param cluster_idxs Integer vector of seed row positions, one per lineage.
#' @param embedding_mat Numeric matrix, rows = cells.
#' @param rho Variance multiplier, shared across lineages.
#'
#' @returns A list of \code{"gaussian"} objects, named \code{"lineage:1"}...
#'
#' @noRd
.form_gaussian_distributions <- function(
    cluster_idxs,
    embedding_mat,
    rho
){
  K <- length(cluster_idxs)
  gaussian_list <- lapply(1:K, function(k){
    .form_gaussian_distribution(
      cluster_idx = cluster_idxs[k],
      embedding_mat = embedding_mat,
      rho = rho
    )
  })
  names(gaussian_list) <- paste0("lineage:", 1:K)
  
  return(gaussian_list)
}

#' Posterior over lineages for every cell
#'
#' Bayes' rule with \code{lineage_prior} against the per-lineage Gaussian
#' densities, computed on the log scale and normalized by
#' \code{.log_sum_exp_normalization()}. A cell's posterior concentrates on
#' whichever lineage seed it sits nearest, so \code{lineage_spread} controls how
#' sharply.
#'
#' @param embedding_mat Numeric matrix, rows = cells.
#' @param gaussian_list List of \code{"gaussian"} objects, one per lineage.
#' @param lineage_prior Numeric vector of prior probabilities, named by lineage;
#'   its names become the column names of the result.
#' @param verbose A numeric; above \code{0} prints a progress tick. Default
#'   \code{0}.
#'
#' @returns A cells-by-lineages numeric matrix whose rows are probability
#'   vectors summing to 1, with \code{rownames(embedding_mat)} and
#'   \code{names(lineage_prior)} as dimnames.
#'
#' @noRd
.compute_posteriors <- function(
    embedding_mat,
    gaussian_list,
    lineage_prior,
    verbose = 0
){
  n <- nrow(embedding_mat)
  K <- length(lineage_prior)
  
  # all the calculations are done on the log scale
  # we use the log-sum-exp trick: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  prob_mat <- matrix(NA, nrow = n, ncol = K)
  if(length(rownames(embedding_mat)) > 0) 
    rownames(prob_mat) <- rownames(embedding_mat)
  colnames(prob_mat) <- names(lineage_prior)
  
  for(i in 1:n){
    if(verbose > 0 && n > 10 && i %% floor(n/10) == 0) cat('*')
    d_vec <- sapply(1:K, function(k){
      .dmvnorm(x = embedding_mat[i,],
               mean = gaussian_list[[k]]$mean, 
               sigma = gaussian_list[[k]]$cov, 
               log = TRUE, 
               checkSymmetry = FALSE)
    })
    stopifnot(length(lineage_prior) == length(d_vec))
    log_vec <- log(lineage_prior) + d_vec
    
    prob_mat[i,] <- .log_sum_exp_normalization(log_vec)
  }
  
  return(prob_mat)
}

#' Normalize log-scale weights into a probability vector
#'
#' we use the log-sum-exp trick: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
#' vec is the vector of log(un-normalized probabilities),
#' and we wish to output the probabilities (non-negative, sums to 1)
#'
#' Subtracting the maximum before exponentiating keeps the largest term at
#' \code{exp(0) = 1}, so no entry can overflow and at least one is non-zero.
#'
#' \code{NA} and infinite entries are \emph{propagated} rather than dropped, so
#' either one trips the closing assertion instead of being silently ignored.
#'
#' @param x A numeric vector of log-scale un-normalized weights.
#' @param tol Tolerance for the closing assertion that the result is
#'   non-negative and sums to 1. Default \code{1e-6}.
#'
#' @returns A numeric probability vector the same length as \code{x}.
#'
#' @noRd
.log_sum_exp_normalization <- function(x, tol = 1e-6){
  c <- max(x)
  y <- c + log(sum(exp(x-c)))
  res <- exp(x-c)
  res <- res/sum(res)
  
  stopifnot(all(res >= -tol, abs(sum(res)-1) <= tol))
  
  return(res)
}

#' Multivariate normal density
#'
#' from the mvtnorm package: https://github.com/cran/mvtnorm/blob/master/R/mvnorm.R
#'
#' Vendored rather than depended upon, so that \code{mvtnorm} need not appear in
#' \code{Imports}. Kept byte-compatible with the upstream implementation ---
#' treat it as third-party code and re-vendor rather than edit if it ever needs
#' updating.
#'
#' Evaluates via the Cholesky factor and \code{backsolve()}, avoiding an
#' explicit inverse. A non-positive-definite \code{sigma} does not error: it
#' returns \code{Inf} at \code{x == mean} and \code{-Inf} elsewhere, matching
#' \code{stats::dnorm()}'s behaviour at zero variance.
#'
#' @param x A numeric vector (treated as one observation) or a matrix with one
#'   observation per row.
#' @param mean Mean vector, length \code{ncol(x)}. Default the zero vector.
#' @param sigma Covariance matrix. Default the identity.
#' @param log Whether to return the log density. Default \code{FALSE}; the
#'   callers here all pass \code{TRUE}.
#' @param checkSymmetry Whether to verify \code{sigma} is symmetric. Default
#'   \code{TRUE}; \code{.compute_posteriors()} passes \code{FALSE} because it
#'   builds diagonal covariances itself and the check is per-cell hot-path cost.
#'
#' @returns A numeric vector with one entry per row of \code{x}.
#'
#' @noRd
.dmvnorm <- function (x,
                      mean = rep(0, p), 
                      sigma = diag(p), 
                      log = FALSE, 
                      checkSymmetry = TRUE)
{
  
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  
  if(!missing(mean)) {
    if(!is.null(dim(mean))) dim(mean) <- NULL
    if (length(mean) != p)
      stop("x and mean have non-conforming size")
  }
  if(!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (checkSymmetry && !Matrix::isSymmetric(
      sigma, 
      tol = sqrt(.Machine$double.eps), 
      check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  
  ## <faster code contributed by Matteo Fasiolo mf364 at bath.ac.uk
  dec <- tryCatch(base::chol(sigma), error=function(e)e)
  if (inherits(dec, "error")) {
    ## warning("cannot compute chol(sigma)"); return(NaN)
    ## behave the same as dnorm(): return Inf or 0
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf # and all other f(.) == 0
  } else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp ^ 2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if(log) return(logretval) else return(exp(logretval))
}
