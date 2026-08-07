#' Attach a real future time point to a simulated current one
#'
#' The other two simulators stop at a future \emph{count} per lineage. This one
#' goes further and hands each individual future cell a parent: given a real
#' matrix of future-time-point cells, it decides which current cell each of them
#' descended from. That is what makes it possible to simulate quantities defined
#' across the bottleneck --- the adaptation index in particular, which needs the
#' embedding shift from parent to progeny and so cannot be computed from counts
#' alone.
#'
#' The construction:
#' \enumerate{
#'   \item The intercept is \bold{rescaled} so the expected progeny total equals
#'     the number of future cells actually supplied. The intercept the caller
#'     passes therefore does not survive; the adjusted one is returned. Current
#'     cells whose rounded contribution is 0 are dropped outright.
#'   \item A push-forward map from current to future embedding space is fitted:
#'     a single shared scale \code{a} and a per-dimension offset \code{b},
#'     estimated by median regression over subsamples. It is deliberately this
#'     rigid --- a flexible map would absorb the very state change the
#'     adaptation index is meant to detect.
#'   \item Each current cell is pushed forward and scored against every future
#'     cell by a Gaussian density, giving a soft parent-child mapping.
#'   \item Each future cell samples a parent from that mapping, with a parent
#'     removed once it has produced its allotted number of progeny.
#' }
#'
#' \bold{Stochastic with no \code{seed_number} argument}; steps 2 to 4 all draw.
#'
#' @param coefficient_intercept The starting intercept, natural-log scale. See
#'   step 1 --- it is rescaled, so this only sets where the search begins.
#' @param embedding_coefficient_vec True coefficients on the embedding, length
#'   \code{ncol(previous_cell_embedding_mat)}.
#' @param future_cell_embedding_mat Numeric matrix of the future time point's
#'   \bold{real} cells, rows = cells and columns = embedding dimensions (at
#'   least 2, and the same dimensions as the previous matrix). Row names are
#'   required for the assignment step.
#' @param lineage_assignment A \bold{factor} of lineage membership for the
#'   current cells (asserted), row-aligned with
#'   \code{previous_cell_embedding_mat}. Usually taken from a
#'   \code{generate_simulation()} run, so that this function attaches a future
#'   to an already-simulated present. Dropped cells are removed with
#'   \code{droplevels()}.
#' @param previous_cell_embedding_mat Numeric matrix of the current time point's
#'   cells, rows = cells. Row names required.
#' @param fatefeatures_coefficient_vec,fatefeatures_mat Optional fate-driving
#'   signal outside the embedding, as in \code{generate_simulation()}. Default
#'   \code{NULL}.
#' @param lineage_spread Multiplier on the covariance used to score future cells
#'   against a pushed-forward current cell. Default \code{1}. Larger values make
#'   parentage more diffuse.
#' @param num_pushforward_training_iter Number of random restarts when fitting
#'   the push-forward map; the restart with the lowest squared error wins.
#'   Default \code{20}.
#' @param num_subsamples Number of (current, future) pairs drawn per restart when
#'   fitting the map. Default \code{200}. Current cells are sampled with
#'   probability proportional to their expected progeny count, so the map is
#'   fitted where the descendants actually come from.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns An object of class \code{"multiomeFate_simulation_future"}, a list
#'   with:
#'   \describe{
#'     \item{\code{cell_fate_potential}}{named numeric, \code{log10(expected
#'       progeny)} per surviving current cell, under the \emph{adjusted}
#'       intercept.}
#'     \item{\code{coefficient_intercept}}{the adjusted intercept, not the one
#'       passed in.}
#'     \item{\code{future_cell_assignment}}{character vector of length
#'       \code{nrow(future_cell_embedding_mat)}, named by future cell ID, giving
#'       the current cell ID it descends from. Compose with
#'       \code{lineage_assignment} to get each future cell's lineage.}
#'     \item{\code{future_lineage_size}}{named numeric, the realized number of
#'       progeny per lineage --- counted from the assignment, so it sums exactly
#'       to the number of future cells.}
#'     \item{\code{mapping_mat}}{the soft parent-child probability matrix,
#'       current cells by future cells, \bold{scaled by 1e3 and rounded} to keep
#'       it storable. Divide by 1e3 to recover probabilities. Columns are
#'       reordered by column sum, so they are not in the input order.}
#'     \item{\code{prev_cell_num_progenitor}}{named numeric, how many future
#'       cells each current cell parented. Despite the name this counts
#'       \emph{progeny}, not progenitors.}
#'   }
#'
#' @export
generate_simulation_attachFuture <- function(
    coefficient_intercept,
    embedding_coefficient_vec,
    future_cell_embedding_mat,
    lineage_assignment,
    previous_cell_embedding_mat,
    fatefeatures_coefficient_vec = NULL,
    fatefeatures_mat = NULL, 
    lineage_spread = 1,
    num_pushforward_training_iter = 20,
    num_subsamples = 200,
    verbose = 0
){
  stopifnot(ncol(future_cell_embedding_mat) > 1,
            is.factor(lineage_assignment))
  
  cell_contribution <- coefficient_intercept + as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) 
  if(!all(is.null(fatefeatures_coefficient_vec))){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution <- exp(cell_contribution)
  
  if(verbose > 0) print("Step 1: Adjusting all the ingredients")
  num_future_cells <- nrow(future_cell_embedding_mat)
  potential_sum <- sum(cell_contribution)
  
  new_coefficient_intercept <- .adjust_coefficient_intercept(
    cell_contribution = cell_contribution,
    coefficient_intercept = coefficient_intercept,
    num_future_cells = num_future_cells
  )
  cell_contribution <- new_coefficient_intercept + as.numeric(previous_cell_embedding_mat %*% embedding_coefficient_vec) 
  names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  if(!all(is.null(fatefeatures_coefficient_vec))){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution <-  exp(cell_contribution)
  cell_contribution_rounded <- round(cell_contribution)
  
  # remove all the cells that don't contribute to the future timepoint
  non_zero_idx <- which(cell_contribution_rounded > 0)
  stopifnot(length(non_zero_idx) > 1)
  cell_contribution_rounded <- cell_contribution_rounded[non_zero_idx]
  lineage_assignment <- droplevels(lineage_assignment[non_zero_idx])
  previous_cell_embedding_mat <- previous_cell_embedding_mat[non_zero_idx,]
  if(verbose > 0) print(paste0("There are ", num_future_cells, " future cells, and the sum of potentials is ", sum(cell_contribution_rounded)))
  
  if(verbose > 0) print("Step 2: Recompute all the adjusted ingredients from generate_simulation")
  cell_fate_potential <- log10(cell_contribution)
  if(length(rownames(previous_cell_embedding_mat)) > 0) 
    names(cell_contribution) <- rownames(previous_cell_embedding_mat)
  future_lineage_size <- sapply(levels(lineage_assignment), function(lev){
    idx <- which(lineage_assignment == lev)
    round(sum(cell_contribution[idx]))
  })
  names(future_lineage_size) <- levels(lineage_assignment)
  
  if(verbose > 0) print("Step 3: Determining the push-forward function")
  pushforward_res <- .compute_pushforward(
    cell_contribution = cell_contribution_rounded,
    future_cell_embedding_mat = future_cell_embedding_mat,
    num_pushforward_training_iter = num_pushforward_training_iter,
    num_subsamples = num_subsamples,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    verbose = verbose - 1
  )
  
  if(verbose > 0) print("Step 4: Computing the mapping of each previous cell to a future cell")
  sd_vec <- apply(future_cell_embedding_mat, 2, stats::sd)
  mapping_mat <- .compute_previous_to_future_mapping(
    future_cell_embedding_mat = future_cell_embedding_mat,
    lineage_spread = lineage_spread,
    previous_cell_embedding_mat = previous_cell_embedding_mat,
    pushforward_func = pushforward_res$pushforward_func,
    sd_vec = sd_vec,
    verbose = verbose - 1
  )
  
  if(verbose > 0) print("Step 5: Assigning future cells to a lineage")
  tmp <- .assign_future_to_previous(
    mapping_mat = mapping_mat,
    previous_cell_contribution = cell_contribution_rounded,
    verbose = verbose - 1
  )
  future_cell_assignment <- tmp$future_cell_assignment
  prev_cell_num_progenitor <- tmp$prev_cell_num_progenitor
  
  stopifnot(length(prev_cell_num_progenitor) == length(lineage_assignment))
  future_lineage_size <- sapply(levels(lineage_assignment), function(lev){
    idx <- which(lineage_assignment == lev)
    sum(prev_cell_num_progenitor[idx])
  })
  names(future_lineage_size) <- levels(lineage_assignment)
  
  return(
    structure(list(cell_fate_potential = cell_fate_potential,
                   coefficient_intercept = new_coefficient_intercept,
                   future_cell_assignment = future_cell_assignment,
                   future_lineage_size = future_lineage_size,
                   mapping_mat = round(mapping_mat*1e3),
                   prev_cell_num_progenitor = prev_cell_num_progenitor),
              class = "multiomeFate_simulation_future")
  )
}

#######################################
#######################################

#' Rescale the intercept so expected progeny matches the future cell count
#'
#' The simulation must produce at least as many progeny as there are real future
#' cells to hand out, or the assignment step runs out of parents. Shifting the
#' intercept by \code{log(num_future_cells) - log(sum(cell_contribution))} makes
#' the expected total match exactly.
#'
#' Rounding then breaks that: each cell's contribution is rounded down or up
#' independently, and the rounded total can land below the target. The loop adds
#' \code{interval_add} to the intercept until the \emph{rounded} total clears the
#' target, so the shift is generally a little larger than the exact one.
#'
#' @param cell_contribution Numeric vector of expected progeny counts per
#'   current cell, computed at the caller's original intercept.
#' @param coefficient_intercept The original intercept, natural-log scale.
#' @param num_future_cells Number of real future cells to be assigned.
#' @param interval_add Step size added to the intercept each iteration. Default
#'   \code{0.1}.
#' @param max_iter Iteration cap; exceeding it errors. Default \code{100}.
#'
#' @returns A single numeric, the adjusted intercept on the natural-log scale.
#'
#' @noRd
.adjust_coefficient_intercept <- function(
    cell_contribution,
    coefficient_intercept,
    num_future_cells,
    interval_add = 0.1,
    max_iter = 100
){
  tmp <- log(num_future_cells) - log(sum(cell_contribution))
  new_cell_contribution <- cell_contribution * exp(tmp)
  
  # now make sure when rounded, it still is larger
  iter <- 0
  interval_value <- 0
  while(sum(round(new_cell_contribution)) < sum(num_future_cells)){
    iter <- iter + 1
    interval_value <- interval_value + interval_add
    new_cell_contribution <- cell_contribution * exp(tmp+interval_value)
    
    if(iter > max_iter) stop("Error with computing intercept")
  }
  
  return(coefficient_intercept + tmp + interval_value)
}

########

#' Fit the current-to-future embedding map, keeping the best restart
#'
#' Subsamples current cells \bold{with probability proportional to their
#' expected progeny count}, so the map is fitted on the cells that actually
#' contribute descendants rather than on the population at large, then fits
#' \code{.compute_pushforward_fit()} against an independent subsample of future
#' cells and keeps whichever restart minimizes total squared error.
#'
#' Pairing two independent subsamples means the fitted map captures only the
#' \emph{aggregate} shift between the two point clouds, which is the intent: any
#' per-cell correspondence at this stage would be circular, since establishing
#' correspondence is what the map is for.
#'
#' @param cell_contribution Numeric vector of expected progeny counts per
#'   current cell; used as sampling weights.
#' @param future_cell_embedding_mat Numeric matrix of future cells, rows =
#'   cells.
#' @param num_pushforward_training_iter Number of random restarts.
#' @param num_subsamples Number of cells drawn per subsample, with replacement,
#'   on both sides.
#' @param previous_cell_embedding_mat Numeric matrix of current cells, rows =
#'   cells.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns The \code{pushforward_res} element of the winning restart: a list
#'   with \code{a} (scalar scale), \code{b} (per-dimension offset vector), and
#'   \code{pushforward_func} (the closure mapping one current cell's coordinates
#'   to future space).
#'
#' @noRd
.compute_pushforward <- function(
    cell_contribution,
    future_cell_embedding_mat,
    num_pushforward_training_iter,
    num_subsamples,
    previous_cell_embedding_mat,
    verbose = 0
){
  n <- nrow(previous_cell_embedding_mat)
  m <- nrow(future_cell_embedding_mat)
  previous_idx <- sample(seq_len(n),
                         size = num_subsamples,
                         prob = cell_contribution,
                         replace = TRUE)
  previous_mat <- previous_cell_embedding_mat[previous_idx,]

  if(verbose > 0) print("Compute the pushforward functions")
  pushforward_list <- lapply(seq_len(num_pushforward_training_iter), function(kk){
    future_idx <- sample(1:m, size = num_subsamples, replace = TRUE)
    future_mat <- future_cell_embedding_mat[future_idx,]
    pushforward_res <- .compute_pushforward_fit(
      future_mat = future_mat,
      previous_mat = previous_mat
    )
    
    tmp <- sapply(1:num_subsamples, function(i){
      vec <- pushforward_res$pushforward_func(previous_mat[i,])
      .l2norm(vec - future_mat[i,])^2
    })
    fit <- sum(tmp)
    
    list(fit = fit,
         pushforward_res = pushforward_res)
  })
  fit_vec <- sapply(pushforward_list, function(x){x$fit})
  
  return(pushforward_list[[which.min(fit_vec)]]$pushforward_res)
}

#' Build the push-forward closure
#'
#' Returns \code{function(vec) a * vec + b}: a single isotropic scale plus a
#' per-dimension translation. Deliberately the most rigid affine map that can
#' still move one point cloud onto another --- no rotation, no per-dimension
#' scaling.
#'
#' @param a A single numeric scale factor, shared across dimensions.
#' @param b A numeric offset vector, one entry per embedding dimension. Its
#'   length fixes the dimension the returned function accepts.
#'
#' @returns A function taking a numeric vector of length \code{length(b)} and
#'   returning one of the same length.
#'
#' @noRd
.pushforward_func_constructor <- function(a, b){
  stopifnot(length(a) == 1)
  
  return(
    function(vec){
      stopifnot(length(vec) == length(b))
      a * vec + b
    }
  )
}

#' Estimate the shared scale and per-dimension offset
#'
#' Fits one univariate \code{stats::lm()} per embedding dimension, takes the
#' \bold{median} of the per-dimension slopes as the single shared scale \code{a},
#' then re-derives each offset as \code{median(future_j - a * previous_j)} at
#' that fixed \code{a} rather than reusing the per-dimension intercepts.
#'
#' Medians rather than means throughout, because the two subsamples are paired
#' arbitrarily (row \code{i} of one against row \code{i} of the other), so a
#' large fraction of the pairs are uninformative and a mean would follow them.
#'
#' @param future_mat Numeric matrix of future cells, rows = cells.
#' @param previous_mat Numeric matrix of current cells; must have exactly the
#'   same dimensions as \code{future_mat} (asserted), since the two are treated
#'   as row-paired.
#'
#' @returns A list with \code{a} (scalar), \code{b} (numeric vector, one entry
#'   per dimension), and \code{pushforward_func} (the closure built from them).
#'
#' @noRd
.compute_pushforward_fit <- function(
    future_mat,
    previous_mat
){
  stopifnot(all(dim(future_mat) == dim(previous_mat)))
  d <- ncol(future_mat)
  n <- nrow(future_mat)
  
  coef_mat <- sapply(1:d, function(j){
    df <- data.frame(
      x = previous_mat[,j],
      y = future_mat[,j]
    )
    lm_res <- stats::lm(y ~ ., data = df)
    stats::coef(lm_res)
  })
  coef_mat <- t(coef_mat)
  colnames(coef_mat) <- c("b", "a")
  
  a <- stats::median(coef_mat[,"a"])
  b <- sapply(1:d, function(j){
    stats::median(future_mat[,j] - a*previous_mat[,j])
  })
  
  return(
    list(a = a,
         b = b,
         pushforward_func = .pushforward_func_constructor(a = a, b = b))
  )
}

########

#' Soft parent-child probabilities between current and future cells
#'
#' Pushes each current cell forward and evaluates a Gaussian density at every
#' future cell, then normalizes \bold{down each column}. The normalization
#' direction is the point: column \code{j} becomes a probability distribution
#' over \emph{candidate parents} for future cell \code{j}, which is exactly what
#' \code{.assign_future_to_previous()} samples from.
#'
#' Columns are first reordered by decreasing column sum, so future cells with a
#' clear parent are assigned before ambiguous ones. Since the capacity-limited
#' assignment closes off parents as they fill, this gives the confident matches
#' first claim.
#'
#' The kernel is a diagonal Gaussian with covariance
#' \code{lineage_spread * diag(sd_vec^2)}, matching the convention
#' \code{.form_gaussian_distribution()} uses in \code{R/simulation.R}.
#'
#' @param future_cell_embedding_mat Numeric matrix of future cells, rows =
#'   cells; row names become the result's column names.
#' @param lineage_spread Multiplier on the covariance.
#' @param previous_cell_embedding_mat Numeric matrix of current cells, rows =
#'   cells; row names become the result's row names. Must have the same number
#'   of columns as the future matrix (asserted).
#' @param pushforward_func The closure from \code{.compute_pushforward()}.
#' @param sd_vec Per-dimension standard deviations of the future embedding.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns A current-by-future numeric matrix whose \bold{columns} sum to 1.
#'   Columns are in decreasing order of their pre-normalization mass, not in the
#'   input order.
#'
#' @noRd
.compute_previous_to_future_mapping <- function(
    future_cell_embedding_mat,
    lineage_spread,
    previous_cell_embedding_mat,
    pushforward_func,
    sd_vec,
    verbose = 0
){
  stopifnot(ncol(previous_cell_embedding_mat) == ncol(future_cell_embedding_mat))
  
  n <- nrow(previous_cell_embedding_mat)
  m <- nrow(future_cell_embedding_mat)
  mapping_mat <- matrix(NA, nrow = n, ncol = m)
  if(length(rownames(future_cell_embedding_mat)) > 0){
    rownames(mapping_mat) <- rownames(previous_cell_embedding_mat)
  }
  if(length(rownames(previous_cell_embedding_mat)) > 0){
    colnames(mapping_mat) <- rownames(future_cell_embedding_mat)
  }
  
  for(i in 1:n){
    if(verbose == 1 && n > 10 && i %% floor(n/10) == 0) cat('*')
    if(verbose > 1) print(paste0("i:", i))
    mean_vec <- pushforward_func(previous_cell_embedding_mat[i,])
    mapping_mat[i,] <- .dmvnorm_log_many_samples(
      mean = mean_vec,
      sigma = lineage_spread*diag(sd_vec^2),
      x_mat = future_cell_embedding_mat
    )
  }
  
  # rearrange the columns
  col_sum <- colSums(mapping_mat)
  mapping_mat <- mapping_mat[,order(col_sum, decreasing = TRUE)]
  
  for(j in 1:m){
    mapping_mat[,j] <- .log_sum_exp_normalization(mapping_mat[,j])
  }
  
  return(mapping_mat)
}

#' Multivariate normal log density at many points, one mean
#'
#' returns a vector of length nrow(x_mat), the log density for each row of x_mat
#'
#' A separate implementation from \code{.dmvnorm()} in \code{R/simulation.R},
#' specialized for the opposite hot path: there, many means against one point;
#' here, one mean against many points. Factoring \code{sigma} once and reusing it
#' across all rows is what makes the cells-by-cells mapping matrix tractable.
#' Uses an explicit \code{solve()} and \code{determinant(logarithm = TRUE)}
#' rather than a Cholesky factor.
#'
#' No degenerate-covariance fallback --- \code{sigma} must be full rank, and this
#' is asserted rather than handled.
#'
#' @param mean Numeric mean vector, length \code{ncol(sigma)}.
#' @param sigma Square, full-rank covariance matrix (asserted).
#' @param x_mat Numeric matrix of evaluation points, rows = points.
#'
#' @returns A numeric vector of length \code{nrow(x_mat)} of log densities.
#'   Always log scale; there is no \code{log} argument.
#'
#' @noRd
.dmvnorm_log_many_samples <- function(mean,
                                      sigma,
                                      x_mat){
  stopifnot(ncol(sigma) == nrow(sigma),
            Matrix::rankMatrix(sigma) == ncol(sigma))
  
  p <- ncol(sigma)
  n <- nrow(x_mat)
  sigma_inv <- solve(sigma)
  determinant_value <- determinant(sigma,
                                   logarithm = TRUE)$modulus
  determinant_value <- as.numeric(determinant_value)
  
  x_mat <- sweep(x_mat, 
                 MARGIN = 2,
                 STATS = mean, 
                 FUN = "-")
  lhs <- x_mat %*% sigma_inv
  rss <- sapply(1:n, function(i){
    lhs[i,] %*% x_mat[i,]
  })
  
  return(- 0.5 * determinant_value - 0.5 * p * log(2 * pi) - 0.5 * rss)
}


########

#' Give every future cell a parent, respecting each parent's quota
#'
#' Walks the columns of \code{mapping_mat} --- future cells, already ordered
#' most-confident-first --- and samples a parent for each from that column's
#' probabilities. A current cell that reaches its allotted progeny count
#' (\code{previous_cell_contribution}) is deleted from the matrix, so it cannot
#' parent any further future cells. Without that quota, the highest-potential
#' cells would take nearly every future cell and the simulated expansion would
#' not match the fate potentials that generated it.
#'
#' This is what makes the simulation's ground truth self-consistent: the realized
#' progeny counts respect the model's expected counts by construction, so
#' \code{future_lineage_size} can be recomputed from the assignment rather than
#' from the potentials.
#'
#' \bold{Stochastic}; the caller seeds. Requires
#' \code{sum(previous_cell_contribution) >= ncol(mapping_mat)} (asserted) ---
#' otherwise the parents run out mid-loop, which is why
#' \code{.adjust_coefficient_intercept()} runs first.
#'
#' @param mapping_mat Current-by-future probability matrix with both dimnames
#'   set (asserted), columns summing to 1. Columns are consumed left to right.
#' @param previous_cell_contribution Named numeric vector of each current cell's
#'   progeny quota, names matching \code{rownames(mapping_mat)}.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{future_cell_assignment}}{character vector of length
#'       \code{ncol(mapping_mat)}, named by future cell, giving its parent's
#'       cell ID.}
#'     \item{\code{prev_cell_num_progenitor}}{numeric vector of length
#'       \code{nrow(mapping_mat)}, named by current cell, counting progeny.
#'       Zero for cells that parented none. Despite the name, this counts
#'       progeny, not progenitors.}
#'   }
#'
#' @noRd
.assign_future_to_previous <- function(
    mapping_mat,
    previous_cell_contribution, 
    verbose = 0
){
  stopifnot(sum(previous_cell_contribution) >= ncol(mapping_mat),
            length(rownames(mapping_mat)) > 0,
            length(colnames(mapping_mat)) > 0)
  
  m <- ncol(mapping_mat)
  n <- nrow(mapping_mat)
  prev_cell_num_progenitor <- rep(0, length = n)
  names(prev_cell_num_progenitor) <- rownames(mapping_mat)
  future_cell_assignment <- rep(NA, length = m)
  names(future_cell_assignment) <- colnames(mapping_mat)
  
  for(j in 1:m){
    if(verbose > 0 && m > 10 && j %% floor(m/10) == 0) cat('*')
    
    prob_vec <- mapping_mat[,j]
    sample_prev_name <- sample(rownames(mapping_mat), size = 1, prob = prob_vec)
    future_cell_assignment[colnames(mapping_mat)[j]] <- sample_prev_name
    prev_cell_num_progenitor[sample_prev_name] <- prev_cell_num_progenitor[sample_prev_name] + 1
    
    if(prev_cell_num_progenitor[sample_prev_name] >= previous_cell_contribution[sample_prev_name]){
      rm_idx <- which(rownames(mapping_mat) == sample_prev_name)
      mapping_mat <- mapping_mat[-rm_idx,,drop = FALSE]
    }
  }
  
  return(
    list(future_cell_assignment = future_cell_assignment,
         prev_cell_num_progenitor = prev_cell_num_progenitor)
  )
}




