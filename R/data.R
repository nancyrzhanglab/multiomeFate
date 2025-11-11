#' Priming simulation dataset (15 lineages)
#'
#' @title Priming simulation dataset
#'
#' @description
#' A small, fully synthetic dataset used to demonstrate and test lineage-fate
#' imputation in a **priming** setting, where cells at an early time point
#' possess heterogeneous growth potential (fate propensity) and downstream
#' lineage sizes scale with \eqn{\exp(\alpha + X\beta)}.
#'
#' The dataset was reduced to **15 lineages** for a fast example.
#'
#' @usage
#' data(priming_simulation)
#'
#' @format
#' A named list with four elements:
#' \describe{
#'   \item{cell_features}{Numeric matrix of size \eqn{n \times p}.
#'     Rows are cells (rownames formatted as \code{"cell:1"}, \code{"cell:2"}, ...).
#'     Columns are preprocessed expression-derived features (no intercept column;
#'     modeling functions add it internally). No missing values.}
#'
#'   \item{cell_lineage}{Character vector of length \eqn{n} giving each cell's
#'     lineage label used for training/evaluation.}
#'
#'   \item{lineage_future_count}{Named numeric vector (length = 15). Names are
#'     lineage IDs matching \code{unique(cell_lineage)}. Each value is an integer
#'     equal to the rounded sum of per-cell expected progenies
#'     (i.e., \eqn{\sum_{i \in \ell} \exp(\alpha + x_i^\top \beta)}).}
#'
#'   \item{tab_mat}{Integer matrix with 15 rows and 2 columns \code{c("now","future")};
#'     rownames are lineage IDs. Column \code{"now"} is the number of current
#'     cells observed in that lineage (in this subset); column \code{"future"}
#'     equals \code{lineage_future_count}.}
#' }
#'
#' @details
#' \strong{How this dataset was generated (summary).}
#'
#' \enumerate{
#' \item Loaded example embeddings via \code{multiomeFate:::data_loader("fasttopics")};
#'   used the \code{"fasttopic.COCL2"} cell embeddings as an early-timepoint RNA
#'   feature space and preprocessed them with an internal helper
#'   \code{.preprocess_rna(…, "day10_COCL2")}.
#'
#' \item Estimated a priming direction and intercept using internal helpers:
#'   \code{.search_for_priming_parameters()} returned \code{coefficient_vec}
#'   and \code{coefficient_intercept}. The intercept was decreased in small steps
#'   until the expected total number of descendants across all cells matched a
#'   target size (chosen from a later timepoint, e.g. \code{"week5_COCL2"}), using
#'   \code{.compute_mean_total_cells()}.
#'
#' \item Assigned cells to lineages with \code{.assign_lineages_priming()},
#'   simulating exponential clonal expansion over multiple rounds. The per-cell
#'   fate potential was \eqn{\log_{10} \exp(\alpha + x^\top \beta)}.
#'   Cells were grouped into lineages by quantiles of this log10 potential;
#'   lineage future sizes were the rounded sums of \eqn{\exp(\alpha + x^\top \beta)}
#'   within each lineage.
#'
#' \item For packaging, lineages were ranked by future size; 15 lineages were
#'   retained by taking equally spaced ranks along this ordering. Cells not
#'   belonging to the retained lineages were dropped. Cell rownames were reset to
#'   \code{"cell:1..n"} for clarity. The four objects in \code{priming_simulation}
#'   were then assembled as below and saved with \code{usethis::use_data()}:
#'
#'   \preformatted{
#'   priming_simulation <- list(
#'     cell_features = embedding_mat_subset,
#'     cell_lineage = as.character(lineage_assignment_subset),
#'     lineage_future_count = lineage_future_size_subset,
#'     tab_mat = cbind(now = table(lineage_assignment_subset)[names(lineage_future_size_subset)],
#'                     future = lineage_future_size_subset)
#'   )
#'   }
#' }
#'
#' \strong{Intended use.}
#' \itemize{
#'   \item Quick-start examples for \code{\link{lineage_cv}},
#'     \code{\link{lineage_cv_finalize}}, and \code{\link{lineage_imputation_sequence}}.
#'   \item Unit tests where a small but structured dataset is helpful.
#' }
#'
#' \strong{Notes.}
#' \itemize{
#'   \item An intercept column is \emph{not} present in \code{cell_features}; modeling
#'         functions add it automatically where required.
#'   \item The data are \emph{fully synthetic} and derived from internal simulation helpers;
#'         no real subject-level outcomes are included.
#' }
#'
#' @seealso
#' \code{\link{lineage_cv}}, \code{\link{lineage_cv_finalize}},
#' \code{\link{lineage_imputation_sequence}}
#'
#' @examples
#' data(priming_simulation)
#' str(priming_simulation)
#'
#' \donttest{
#' # Minimal end-to-end example:
#' with(priming_simulation, {
#'   set.seed(10)
#'   cv <- lineage_cv(
#'     cell_features = cell_features,
#'     cell_lineage  = cell_lineage,
#'     future_timepoint = "future",
#'     lineage_future_count = lineage_future_count,
#'     lambda_initial = NA,
#'     lambda_sequence_length = 20,
#'     tab_mat = tab_mat,
#'     num_folds = 5,
#'     verbose = 1
#'   )
#'   fit <- lineage_cv_finalize(
#'     cell_features = cell_features,
#'     cell_lineage  = cell_lineage,
#'     fit_res = cv,
#'     lineage_future_count = lineage_future_count
#'   )
#'
#'   fit$lambda
#'   head(fit$cell_imputed_score)
#'   head(fit$lineage_imputed_count)
#' })
#' }
#'
#' @keywords datasets
"priming_simulation"

#' Plastic simulation dataset (15 lineages)
#'
#' @title Plastic simulation dataset
#'
#' @description
#' A small, fully synthetic dataset used to demonstrate and test lineage-fate
#' imputation in a **plastic** setting, where each lineage has roughly the
#' same median growth potential, but certain lineages have a few cells with
#' extremely high growth potential (fate propensity). Downstream
#' lineage sizes scale with \eqn{\exp(\alpha + X\beta)}.
#'
#' The dataset was reduced to **15 lineages** for a fast example.
#'
#' @usage
#' data(plastic_simulation)
#'
#' @format
#' A named list with four elements:
#' \describe{
#'   \item{cell_features}{Numeric matrix of size \eqn{n \times p}.
#'     Rows are cells (rownames formatted as \code{"cell:1"}, \code{"cell:2"}, ...).
#'     Columns are preprocessed expression-derived features (no intercept column;
#'     modeling functions add it internally). No missing values.}
#'
#'   \item{cell_lineage}{Character vector of length \eqn{n} giving each cell's
#'     lineage label used for training/evaluation.}
#'
#'   \item{lineage_future_count}{Named numeric vector (length = 15). Names are
#'     lineage IDs matching \code{unique(cell_lineage)}. Each value is an integer
#'     equal to the rounded sum of per-cell expected progenies
#'     (i.e., \eqn{\sum_{i \in \ell} \exp(\alpha + x_i^\top \beta)}).}
#'
#'   \item{tab_mat}{Integer matrix with 15 rows and 2 columns \code{c("now","future")};
#'     rownames are lineage IDs. Column \code{"now"} is the number of current
#'     cells observed in that lineage (in this subset); column \code{"future"}
#'     equals \code{lineage_future_count}.}
#' }
#'
#' @details
#' \strong{How this dataset was generated (summary).}
#'
#' \enumerate{
#' \item Loaded example embeddings via \code{multiomeFate:::data_loader("fasttopics")};
#'   used the \code{"fasttopic.COCL2"} cell embeddings as an early-timepoint RNA
#'   feature space and preprocessed them with an internal helper
#'   \code{.preprocess_rna(…, "day10_COCL2")}.
#'
#' \item Estimated a plastic direction and intercept using internal helpers:
#'   \code{.generate_simulation_plastic()} returned \code{coefficient_vec}
#'   and \code{coefficient_intercept}. This function also assigns cells to lineages.
#'   The per-cell fate potential was \eqn{\log_{10} \exp(\alpha + x^\top \beta)}.
#'   Cells were probabilistically assigned to lineages via the internal functions
#'   \code{.compute_plastic_probabilities()} and \code{.assign_plastic_lineages()}.
#'
#' \item For packaging, lineages were ranked by future size; 15 lineages were
#'   retained by taking equally spaced ranks along this ordering. Cells not
#'   belonging to the retained lineages were dropped. Cell rownames were reset to
#'   \code{"cell:1..n"} for clarity. The four objects in \code{plastic_simulation}
#'   were then assembled as below and saved with \code{usethis::use_data()}:
#'
#'   \preformatted{
#'   plastic_simulation <- list(
#'     cell_features = embedding_mat_subset,
#'     cell_lineage = as.character(lineage_assignment_subset),
#'     lineage_future_count = lineage_future_size_subset,
#'     tab_mat = cbind(now = table(lineage_assignment_subset)[names(lineage_future_size_subset)],
#'                     future = lineage_future_size_subset)
#'   )
#'   }
#' }
#'
#' \strong{Intended use.}
#' \itemize{
#'   \item Quick-start examples for \code{\link{lineage_cv}},
#'     \code{\link{lineage_cv_finalize}}, and \code{\link{lineage_imputation_sequence}}.
#'   \item Unit tests where a small but structured dataset is helpful.
#' }
#'
#' \strong{Notes.}
#' \itemize{
#'   \item An intercept column is \emph{not} present in \code{cell_features}; modeling
#'         functions add it automatically where required.
#'   \item The data are \emph{fully synthetic} and derived from internal simulation helpers;
#'         no real subject-level outcomes are included.
#' }
#'
#' @seealso
#' \code{\link{lineage_cv}}, \code{\link{lineage_cv_finalize}},
#' \code{\link{lineage_imputation_sequence}}
#'
#' @examples
#' data(plastic_simulation)
#' str(plastic_simulation)
#'
#' \donttest{
#' # Minimal end-to-end example:
#' with(plastic_simulation, {
#'   set.seed(10)
#'   cv <- lineage_cv(
#'     cell_features = cell_features,
#'     cell_lineage  = cell_lineage,
#'     future_timepoint = "future",
#'     lineage_future_count = lineage_future_count,
#'     lambda_initial = NA,
#'     lambda_sequence_length = 20,
#'     tab_mat = tab_mat,
#'     num_folds = 5,
#'     verbose = 1
#'   )
#'   fit <- lineage_cv_finalize(
#'     cell_features = cell_features,
#'     cell_lineage  = cell_lineage,
#'     fit_res = cv,
#'     lineage_future_count = lineage_future_count
#'   )
#'
#'   fit$lambda
#'   head(fit$cell_imputed_score)
#'   head(fit$lineage_imputed_count)
#' })
#' }
#'
#' @keywords datasets
"plastic_simulation"
