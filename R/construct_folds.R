#' Partition lineages into cross-validation folds
#'
#' Splits the \emph{lineages} --- not the cells --- into \code{num_folds}
#' groups, and reports both the lineage membership and the corresponding cell
#' row indices for each fold. Cross-validation in CYFER has to hold out whole
#' lineages, because a lineage's future count is a single observation shared by
#' all of its cells; splitting cells across folds would leak the held-out
#' response into training.
#'
#' Lineages are sorted by descending \code{lineage_future_count}, shuffled
#' within contiguous blocks, then dealt round-robin. Sorting first is what keeps
#' the folds comparable: dealing at random would let one fold collect several of
#' the largest clones and dominate the held-out objective. Because the deal is
#' round-robin, fold sizes are \code{floor()} or \code{ceiling()} of
#' \code{num_lineages / num_folds}, never a short final fold.
#'
#' \bold{This function calls \code{sample()} and takes no \code{seed_number}
#' argument} --- the caller seeds the stream. \code{cyfer()} does so immediately
#' before calling this.
#'
#' @param cell_lineage A character or factor vector, element \code{i} naming the
#'   lineage of cell \code{i}. Coerced to character internally. Its
#'   \emph{positions} are what \code{cv_cell_list} reports, so it must be
#'   row-aligned with the \code{cell_features} the caller will subset.
#' @param lineage_future_count A named numeric vector: names are lineage IDs,
#'   values are cell counts at the future time point. Only the ordering induced
#'   by the values is used here. Names must be unique and non-missing. Lineages
#'   named here but absent from \code{cell_lineage} still occupy a fold slot ---
#'   \code{cyfer()} drops them beforehand for exactly this reason.
#' @param num_folds Number of folds. Must be at least 2 (one fold would leave no
#'   training data) and at most the number of lineages (an empty fold would make
#'   the training set the whole dataset). Default is \code{10}.
#'
#' @returns A list with two elements, both named \code{"fold:1"} ...
#'   \code{"fold:num_folds"} and both of length \code{num_folds}:
#'   \describe{
#'     \item{\code{cv_cell_list}}{a list of integer vectors, each holding the
#'       \emph{positions in \code{cell_lineage}} of the cells belonging to that
#'       fold's lineages. These index rows of \code{cell_features}. An element
#'       is \code{NULL} rather than \code{integer(0)} if none of the fold's
#'       lineages has any cell --- which \code{cyfer()} prevents by dropping
#'       cell-less lineages before calling.}
#'     \item{\code{fold_lineage_list}}{a list of character vectors giving the
#'       lineage names assigned to each fold.}
#'   }
#'   The folds partition the lineages exactly: every name in
#'   \code{names(lineage_future_count)} appears in exactly one fold, and no fold
#'   is empty. This is asserted before returning.
#'
#' @noRd
construct_folds <- function(cell_lineage,
                            lineage_future_count,
                            num_folds = 10){
  cell_lineage <- as.character(cell_lineage)

  lineage_names <- names(lineage_future_count)
  num_lineages <- length(lineage_names)

  if(num_lineages == 0 || is.null(lineage_names) || any(is.na(lineage_names))){
    stop("`lineage_future_count` must be a non-empty named vector")
  }
  if(anyDuplicated(lineage_names) != 0){
    stop("`lineage_future_count` has duplicated lineage names")
  }
  if(length(num_folds) != 1 || is.na(num_folds) || num_folds < 2){
    stop("`num_folds` must be a single number >= 2")
  }
  # Every fold must be non-empty: `cyfer()` trains on the complement of a fold,
  # and an empty fold would leave the training set equal to the whole dataset.
  if(num_folds > num_lineages){
    stop("`num_folds` (", num_folds, ") exceeds the number of lineages (",
         num_lineages, ")")
  }

  lineages_ordered <- lineage_names[order(lineage_future_count, decreasing = TRUE)]
  # Misnomer kept for continuity: this is the *maximum* fold size, and below it
  # doubles as the shuffle block width. Because the deal is round-robin, actual
  # fold sizes are floor() or ceiling() of num_lineages/num_folds -- e.g. 7
  # lineages into 3 folds gives sizes (3, 2, 2), not (3, 3, 1).
  num_per_fold <- ceiling(num_lineages/num_folds)

  # Shuffle within each contiguous block of `num_per_fold` lineages, so the
  # round-robin deal below does not always assign the same rank to the same fold.
  # There are `num_blocks` such blocks -- not `num_per_fold` of them.
  #
  # Note: the deal cycles every `num_folds`, so blocks of width `num_per_fold`
  # let lineages move between deal cycles. Using `num_folds`-wide blocks instead
  # would balance the per-fold future-count totals roughly 2-3x more tightly.
  # Left as-is deliberately: that is a change to fold *composition*, not a bug,
  # and it would shift every existing result again.
  num_blocks <- ceiling(num_lineages/num_per_fold)
  for(i in seq_len(num_blocks)){
    idx_vec <- ((i-1)*num_per_fold+1):min(i*num_per_fold, num_lineages)
    # `sample(x)` on a length-1 vector permutes seq_len(x) instead of returning x,
    # so permute the positions rather than the values.
    lineages_ordered[idx_vec[sample(length(idx_vec))]] <- lineages_ordered[idx_vec]
  }

  # Deal the sorted lineages round-robin into folds. Indices past the end are
  # dropped; clamping them instead would assign the last lineage to every fold.
  fold_lineage_list <- lapply(seq_len(num_folds), function(i){
    idx_vec <- i + (0:num_per_fold)*num_folds
    lineages_ordered[idx_vec[idx_vec <= num_lineages]]
  })
  names(fold_lineage_list) <- paste0("fold:", seq_len(num_folds))

  cv_cell_list <- lapply(fold_lineage_list, function(lineages){
    unlist(lapply(lineages, function(lineage){
      which(cell_lineage == lineage)
    }))
  })
  names(cv_cell_list) <- names(fold_lineage_list)

  # The folds must partition the lineages: each appears in exactly one fold.
  # Compare against the lineage names as they arrived, not against the shuffled
  # ordering, so that a shuffle which corrupted the set is also caught.
  assigned_lineages <- unlist(fold_lineage_list, use.names = FALSE)
  stopifnot(!anyNA(assigned_lineages),
            anyDuplicated(assigned_lineages) == 0,
            setequal(assigned_lineages, lineage_names),
            lengths(fold_lineage_list) > 0)

  list(cv_cell_list = cv_cell_list,
       fold_lineage_list = fold_lineage_list)
}
