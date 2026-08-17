#' Simulate lineages that differ in heterogeneity, not in mean ("plastic")
#'
#' The counterpart to \code{generate_simulation()}. There, lineages are spatial
#' clumps and so differ in \emph{mean} fate potential --- the selection regime.
#' Here the assignment is driven by fate potential directly rather than by
#' position, and lineages are built to have the \bold{same mean potential but
#' different spreads}: lineage 1 collects cells from the extremes of the
#' potential distribution, the last lineage collects cells near the middle.
#' Between-lineage variation in expansion is then small while within-lineage
#' variation is large, which is what a plastic (adaptive) population looks like
#' and what a method that only compares clone sizes would miss.
#'
#' This is how \code{data/plastic_simulation.rda} was produced; see
#' \code{?plastic_simulation}.
#'
#' Order of construction, and note it is the reverse of the priming simulation
#' --- \bold{potentials are computed first and lineages assigned from them},
#' rather than lineages first and potentials from position:
#' \enumerate{
#'   \item Each cell's expected progeny count is
#'     \code{exp(coefficient_intercept + x_i^T beta)}, optionally Poisson-drawn.
#'   \item Each lineage gets a Gaussian over \emph{log potential}, all with the
#'     same mean and with standard deviations interpolating from
#'     \code{sd * rho} down to \code{sd / rho}; cells are scored against those.
#'   \item Cells are sampled into lineages from that posterior, sequentially,
#'     with a lineage removed from contention once it reaches
#'     \code{ceiling(n / num_lineages)} cells --- so lineages come out
#'     near-equal in size.
#' }
#'
#' \bold{Stochastic with no \code{seed_number} argument}; set a seed before
#' calling.
#'
#' @param embedding_mat Numeric matrix of the current time point's cells, rows =
#'   cells and columns = embedding dimensions (at least 2). Row names are
#'   generated as \code{"cell:1"}... if absent.
#' @param bool_add_randomness Whether to Poisson-draw the realized progeny
#'   counts. Default \code{TRUE}.
#' @param coefficient_intercept The intercept, natural-log scale. Default
#'   \code{0}.
#' @param embedding_coefficient_vec True coefficients on the embedding, length
#'   \code{ncol(embedding_mat)}, natural-log scale. Default all ones.
#' @param fatefeatures_coefficient_vec True coefficients on the extra fate
#'   features. Default \code{NULL}.
#' @param fatefeatures_mat Optional matrix of fate-driving features outside the
#'   embedding, same rows as \code{embedding_mat}; signal deliberately withheld
#'   from an estimator fitted on the embedding alone. Default \code{NULL}.
#' @param lineage_mean_spread Controls whether lineage \emph{means} are allowed
#'   to differ. Only two values are honoured: \code{1} (the default) holds every
#'   lineage's mean at the population mean, which is the plastic regime and also
#'   turns on the equal-size constraint in step 3; \code{NA} spreads the means
#'   across quantiles of the potential distribution, shrinks the working
#'   standard deviation to a quarter, and drops the equal-size constraint. Any
#'   other numeric warns and is treated as \code{1}.
#' @param lineage_sd_spread The ratio \code{rho} defining the spread ladder:
#'   lineage 1 gets standard deviation \code{sd * rho} and the last gets
#'   \code{sd / rho}, interpolated linearly in between, so values above 1 give
#'   the intended high-to-low ordering. Default \code{NA}, meaning derive it
#'   from the data as \code{max(|log potential - mean|) / sd / 2}. The value
#'   actually used comes back in the output.
#' @param num_lineages Number of lineages. Default \code{10}.
#' @param tol Unused; retained for signature compatibility with
#'   \code{generate_simulation()}. Default \code{1e-06}.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns An object of class \code{"multiomeFate_simulation_plastic"}, a list
#'   with the same elements as \code{generate_simulation()} except that
#'   \code{gaussian_list} is absent and \code{lineage_sd_spread} (the realized
#'   \code{rho}) is present:
#'   \describe{
#'     \item{\code{cell_fate_potential}}{named numeric, \code{log10(realized
#'       progeny + 1)} per cell.}
#'     \item{\code{cell_fate_potential_truth}}{named numeric,
#'       \code{log10(expected progeny)} per cell --- the ground truth for
#'       scoring \code{cell_imputed_score}.}
#'     \item{\code{coefficient_intercept}, \code{embedding_mat},
#'       \code{fatefeatures_coefficient_vec}, \code{fatefeatures_mat}}{as
#'       supplied.}
#'     \item{\code{lineage_assignment}}{factor named by cell, levels
#'       \code{"lineage:1"}..., reordered to match
#'       \code{cell_fate_potential_truth}.}
#'     \item{\code{lineage_future_size}}{named numeric, the future count per
#'       lineage.}
#'     \item{\code{lineage_sd_spread}}{the \code{rho} used, whether supplied or
#'       derived.}
#'     \item{\code{prob_mat}}{cells-by-lineages posterior matrix. Its rows are
#'       in the internal reordered cell order, not the input order.}
#'     \item{\code{summary_mat}}{5-by-\code{num_lineages} matrix, rows
#'       \code{mean}, \code{median}, \code{sd}, \code{range}, \code{future_size}.
#'       The check that the simulation did what it claims: the \code{mean} row
#'       should be near-flat across lineages while \code{sd} and \code{range}
#'       decrease.}
#'   }
#'
#' @export
generate_simulation_plastic <- function(embedding_mat,
                                        bool_add_randomness = TRUE, 
                                        coefficient_intercept = 0, 
                                        embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
                                        fatefeatures_coefficient_vec = NULL,
                                        fatefeatures_mat = NULL, 
                                        lineage_mean_spread = 1, # NA or a value 1 or larger. "1" means no spread
                                        lineage_sd_spread = NA, # NA or a numeric. Lineage 1 is lineage_sd_spread, and Lineage num_lineage is 1/lineage_sd_spread.
                                        num_lineages = 10, 
                                        tol = 1e-06, 
                                        verbose = 0) {
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
  
  gamma <- lineage_mean_spread
  rho <- lineage_sd_spread
  if (length(rownames(embedding_mat)) == 0){
    rownames(embedding_mat) <- paste0("cell:", 1:n)
  }
  
  stopifnot(d > 1, 
            length(embedding_coefficient_vec) == d)
  
  if (verbose > 0) 
    print("Step 1: Computing fate potential of all the cells")
  cell_contribution <- as.numeric(embedding_mat %*% embedding_coefficient_vec)
  if(d2 > 0){
    cell_contribution <- cell_contribution + as.numeric(fatefeatures_mat %*% fatefeatures_coefficient_vec)
  }
  cell_contribution_truth <- exp(cell_contribution + coefficient_intercept)
  names(cell_contribution_truth) <- rownames(embedding_mat) 
  cell_contribution_random <- cell_contribution_truth
  
  if (bool_add_randomness) {
    if (verbose > 0) 
      print("Step 1b: (Optional) Adding randomness")
    cell_contribution_random <- stats::rpois(n = length(cell_contribution_random), 
                                             lambda = cell_contribution_random)
    names(cell_contribution_random) <- rownames(embedding_mat) 
  }
  
  if (verbose > 0) 
    print("Step 2: Computing probability of cells in each lineage")
  tmp <- .compute_plastic_probabilities(
    cell_contribution_truth = cell_contribution_truth,
    num_lineages = num_lineages,
    gamma = gamma,
    rho = rho,
    verbose = verbose - 1
  )
  prob_mat <- tmp$prob_mat; rho <- tmp$rho
  
  if (verbose > 0) 
    print("Step 3: Assigning cells to lineages")
  enforce_equal_size <- !is.na(gamma)
  lineage_assignment <- .assign_plastic_lineages(enforce_equal_size = enforce_equal_size,
                                                 prob_mat = prob_mat)
  # reorder the lineage assignment
  lineage_assignment <- lineage_assignment[names(cell_contribution_truth)]
  
  lineage_future_size <- sapply(levels(lineage_assignment), 
                                function(lev) {
                                  idx <- which(lineage_assignment == lev)
                                  round(sum(cell_contribution_random[idx]))
                                })
  names(lineage_future_size) <- levels(lineage_assignment)
  
  summary_mat <- .compute_summary_lineages(cell_fate_potential_truth = log10(cell_contribution_truth),
                                           lineage_assignment = lineage_assignment,
                                           lineage_future_size = lineage_future_size)
  
  if (verbose > 0) 
    print("Step 4: Outputting")
  return(
    structure(list(cell_fate_potential = log10(cell_contribution_random + 1), 
                   cell_fate_potential_truth = log10(cell_contribution_truth), 
                   coefficient_intercept = coefficient_intercept,
                   embedding_mat = embedding_mat,
                   fatefeatures_coefficient_vec = fatefeatures_coefficient_vec, 
                   fatefeatures_mat = fatefeatures_mat,
                   lineage_assignment = lineage_assignment, 
                   lineage_future_size = lineage_future_size,
                   lineage_sd_spread = rho,
                   prob_mat = prob_mat,
                   summary_mat = summary_mat),
              class = "multiomeFate_simulation_plastic")
  )
}

#################

#' Score every cell against each lineage's log-potential Gaussian
#'
#' Builds the ladder of lineage distributions and evaluates each cell's log
#' potential under all of them. All the plastic regime's structure lives here:
#' every lineage shares a mean (when \code{gamma} is not \code{NA}) and they
#' differ only in standard deviation, running from \code{sd * rho} for lineage 1
#' down to \code{sd / rho} for lineage \code{K}. A wide lineage therefore favours
#' cells far from the population mean in \emph{either} direction, and a narrow
#' one favours typical cells.
#'
#' Cells are reordered by \code{.reorder_by_contribution()} before scoring, which
#' is what lets the sequential capacity-limited assignment in
#' \code{.assign_plastic_lineages()} fill the wide lineages with genuine
#' extremes rather than whatever happened to come first.
#'
#' @param cell_contribution_truth Named numeric vector of expected progeny
#'   counts per cell, strictly positive (asserted). Logged internally, so the
#'   Gaussians are over log potential.
#' @param num_lineages Number of lineages.
#' @param gamma The \code{lineage_mean_spread} argument. \code{1} holds all
#'   means equal; \code{NA} spreads them over quantiles and quarters the working
#'   \code{sd}. Any other value warns and is treated as \code{1}. \code{NA} for
#'   both \code{gamma} and \code{rho} warns, since lineages then vary in mean
#'   \emph{and} spread together, confounding the two regimes.
#' @param rho The \code{lineage_sd_spread} ratio, or \code{NA} to derive it from
#'   the data.
#' @param verbose A numeric; above \code{1} prints the realized lineage means
#'   and standard deviations, which is the quickest way to confirm the ladder
#'   came out as intended. Default \code{0}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{lineage_sd_vec}}{named numeric of the standard deviations,
#'       one per lineage.}
#'     \item{\code{prob_mat}}{cells-by-lineages matrix whose rows are
#'       probability vectors summing to 1. \bold{Rows are in the reordered cell
#'       order}, not the input order.}
#'     \item{\code{rho}}{the ratio used, whether supplied or derived.}
#'   }
#'
#' @noRd
.compute_plastic_probabilities <- function(
    cell_contribution_truth,
    num_lineages,
    gamma,
    rho,
    verbose = 0
){
  
  stopifnot(all(cell_contribution_truth > 0))
  
  # Spreading the means and the variances at once confounds the two regimes this
  # simulation exists to separate, so it is not a supported configuration.
  if(is.na(rho) && is.na(gamma)){
    stop("`lineage_mean_spread` and `lineage_sd_spread` cannot both be NA: ",
         "that gives lineages large means *and* large variances together, ",
         "which confounds priming with plasticity. Set one of them.")
  }
  if(!is.na(gamma) && abs(gamma - 1) > 1e-4){
    warning("lineage_mean_spread can only handle NA or 1. Defaulting to 1")
  }
  
  cell_contribution_truth <- log(cell_contribution_truth)
  n <- length(cell_contribution_truth)
  K <- num_lineages
  mean_val <- mean(cell_contribution_truth)
  sd_val <- stats::sd(cell_contribution_truth)

  # compute gamma if it's NA
  if(is.na(gamma)){
    lineage_mean_vec <- stats::quantile(cell_contribution_truth, 
                                        probs = seq(1, 0, length.out = num_lineages))
    sd_val <- sd_val/4
  } else {
    lineage_mean_vec <- rep(mean_val, length = num_lineages)
  }
  
  # compute rho if it's NA
  if(is.na(rho)){
    rho <- (max(abs(cell_contribution_truth - mean_val))/sd_val)/2
  }
  lineage_sd_vec <- seq(sd_val*rho, sd_val/rho, length.out = K)
  
  
  names(lineage_mean_vec) <- paste0("lineage:", 1:K) 
  names(lineage_sd_vec) <- names(lineage_mean_vec)
  
  if(verbose > 1){
    print("Lineage means: ")
    print(round(lineage_mean_vec,2))
    print("Lineage std's: ")
    print(lineage_sd_vec)
  }

  
  # reorder the cell contribution
  idx <- .reorder_by_contribution(abs(cell_contribution_truth - mean_val))
  cell_contribution_truth <- cell_contribution_truth[idx]
  
  prob_mat <- matrix(0, nrow = n, ncol = K)
  rownames(prob_mat) <- names(cell_contribution_truth)
  colnames(prob_mat) <- names(lineage_sd_vec)
  
  for(j in 1:K){
    if(verbose > 0 && K > 10 && j %% floor(K/10) == 0) cat('*')
    lineage <- colnames(prob_mat)[j]
    
    prob_mat[,lineage] <- stats::dnorm(
      x = cell_contribution_truth,
      mean = lineage_mean_vec[lineage],
      sd = lineage_sd_vec[lineage],
      log = TRUE
    )
  }
  
  for(i in 1:n){
    prob_mat[i,] <- .log_sum_exp_normalization(prob_mat[i,])
  }
  
  return(
    list(lineage_sd_vec = lineage_sd_vec,
         prob_mat = prob_mat,
         rho = rho)
  )
}

#' Interleave the extremes of a vector with its centre
#'
#' Returns positions ordered as: smallest, largest, second smallest, second
#' largest, and so on, de-duplicated at the meeting point.
#'
#' The purpose is to make the capacity-limited assignment in
#' \code{.assign_plastic_lineages()} fair. That loop walks cells in order and
#' closes a lineage once it is full, so whichever cells come last get only the
#' lineages nobody wanted. Alternating between the two tails means the extreme
#' cells --- the ones the wide lineages exist to collect --- are placed early
#' and interleaved, rather than all arriving after the wide lineages have
#' filled.
#'
#' @param vec A numeric vector; callers pass \code{abs(x - mean(x))}, so
#' "smallest" means nearest the mean.
#'
#' @returns An integer permutation of \code{seq_along(vec)}.
#'
#' @noRd
.reorder_by_contribution <- function(vec){
  n <- length(vec)
  order_dec <- order(vec, decreasing = TRUE)
  order_inc <- order(vec, decreasing = FALSE)
  
  vec <- as.numeric(rbind(order_inc, order_dec))
  vec <- vec[!duplicated(vec)]
  
  return(vec)
}

#' Sample cells into lineages, optionally capping lineage size
#'
#' Walks the rows of \code{prob_mat} in order and draws one lineage per cell
#' from that row. When \code{enforce_equal_size} is \code{TRUE}, a lineage that
#' reaches \code{ceiling(n / K)} cells has its column dropped from
#' \code{prob_mat}, removing it from every subsequent draw. Sizes then come out
#' near-equal, which keeps lineage size from itself carrying the signal --- in
#' the plastic regime the lineages are supposed to differ in composition, not in
#' how many cells they have.
#'
#' Because the cap is enforced by deletion rather than by renormalizing, the
#' order of rows matters: cells drawn late choose among whatever remains. See
#' \code{.reorder_by_contribution()} for why the caller interleaves the extremes
#' first.
#'
#' A row whose probabilities are all at or below \code{1e-6} falls back to a
#' uniform draw over the remaining lineages, so a cell in the far tail of every
#' lineage's Gaussian is still placed rather than erroring.
#'
#' \bold{Stochastic}; the caller seeds.
#'
#' @param enforce_equal_size Whether to cap lineage sizes. \code{TRUE} when
#'   \code{lineage_mean_spread} is not \code{NA}.
#' @param prob_mat Cells-by-lineages matrix of probabilities, row names being
#'   cell IDs and column names lineage names.
#'
#' @returns A factor of length \code{nrow(prob_mat)}, named by cell, with all
#'   \code{K} lineages as levels even if a lineage received no cells.
#'
#' @noRd
.assign_plastic_lineages <- function(enforce_equal_size,
                                     prob_mat){
  n <- nrow(prob_mat)
  K <- ncol(prob_mat)
  
  maximum_lineage_size <- ceiling(n/K)
  current_size <- rep(0, K)
  names(current_size) <- colnames(prob_mat)
  
  lineage_assignment <- rep(NA, n)
  for(i in 1:n) {
    stopifnot(is.matrix(prob_mat), ncol(prob_mat) >= 1)
    
    prob_vec <- prob_mat[i,]
    if(all(prob_vec <= 1e-6)) 
      prob_vec <- rep(1/ncol(prob_mat), length = ncol(prob_mat))
    
    lineage <- colnames(prob_mat)[sample(1:ncol(prob_mat), size = 1, prob = prob_vec)]
    current_size[lineage] <- current_size[lineage] + 1
    lineage_assignment[i] <- lineage
    
    if(enforce_equal_size & current_size[lineage] >= maximum_lineage_size){
      prob_mat <- prob_mat[,-which(colnames(prob_mat) == lineage), drop = FALSE]
    }
  }
  
  lineage_assignment <- factor(lineage_assignment, 
                               levels = names(current_size))
  names(lineage_assignment) <- rownames(prob_mat)
  
  return(lineage_assignment)
}

#' Per-lineage summary of the true fate potentials
#'
#' The table that makes a simulation's regime legible at a glance, and the one
#' place the priming and plastic simulations can be compared directly (both call
#' this). Read the \code{mean} row against the \code{sd} and \code{range} rows:
#' priming spreads the means and keeps the spreads small, plastic holds the
#' means flat and varies the spreads.
#'
#' @param cell_fate_potential_truth Named numeric vector of true log10
#'   potentials per cell. Must be name-aligned with \code{lineage_assignment};
#'   asserted.
#' @param lineage_assignment Factor of lineage membership, named by cell. Its
#'   \code{levels()} set the column order.
#' @param lineage_future_size Named numeric of future counts per lineage;
#'   reordered to the factor levels internally.
#'
#' @returns A 5-by-\code{num_lineages} numeric matrix with row names
#'   \code{"mean"}, \code{"median"}, \code{"sd"}, \code{"range"},
#'   \code{"future_size"} and column names the lineage levels. A lineage holding
#'   one cell gives \code{NA} for \code{sd} and \code{0} for \code{range}.
#'
#' @noRd
.compute_summary_lineages <- function(cell_fate_potential_truth,
                                      lineage_assignment,
                                      lineage_future_size){
  stopifnot(all(names(cell_fate_potential_truth) == names(lineage_assignment)))
  stopifnot(all(sort(unique(lineage_assignment)) == names(lineage_future_size)))
  
  mean_val <- sapply(levels(lineage_assignment), function(lineage){
    mean(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  med_val <- sapply(levels(lineage_assignment), function(lineage){
    stats::median(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  sd_val <- sapply(levels(lineage_assignment), function(lineage){
    stats::sd(cell_fate_potential_truth[lineage_assignment == lineage])
  })
  range_val <- sapply(levels(lineage_assignment), function(lineage){
    diff(range(cell_fate_potential_truth[lineage_assignment == lineage]))
  })
  
  lineage_future_size <- lineage_future_size[levels(lineage_assignment)]
  
  mat <- rbind(mean_val, med_val, sd_val, range_val, lineage_future_size)
  rownames(mat) <- c("mean", "median", "sd", "range", "future_size")
  colnames(mat) <- levels(lineage_assignment)
  
  return(mat)
}