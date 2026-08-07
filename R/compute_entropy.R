#' Assemble the per-cell fate-composition table used by the simplex plots
#'
#' Turns a matrix of per-cell fate potentials --- one column per candidate
#' \emph{fate} (in the paper's LARRY analysis: Monocyte, Neutrophil,
#' Undifferentiated) --- into the data frame that \code{\link{plot_simplex}}
#' consumes. Each row is normalized to a composition summing to 1, giving the
#' cell's position in the simplex, and is then annotated with its current cell
#' type, its total predicted progeny count, and two \emph{lineage}-level
#' quantities computed from what that cell's lineage actually became at the
#' later time point: the dominant observed fate and the Shannon entropy of the
#' observed fate distribution.
#'
#' The columns of \code{cell_imputation_mat} are separate CYFER fits --- one
#' \code{cyfer_finalize()} run per fate --- placed side by side. Nothing here
#' checks that, so the caller is responsible for the columns being on a common
#' scale.
#'
#' Note the asymmetry the name does not convey: \code{entropy} and
#' \code{dominant_fate} are properties of the \emph{lineage}, computed from
#' observed cell types at \code{later_timepoint}, so every cell of a lineage
#' carries the same value. The composition columns are per-cell predictions.
#'
#' @param cell_imputation_mat A numeric matrix, rows = cells (row names must be
#'   cell IDs present in the Seurat object) and columns = candidate fates
#'   (column names required, and used as the composition column names of the
#'   result). Typically the \code{cell_imputed_score} vectors from one
#'   \code{cyfer_finalize()} fit per fate, column-bound.
#' @param later_timepoint The value of \code{variable_timepoint} identifying the
#'   future time point whose observed cell types define the dominant fate and
#'   entropy.
#' @param seurat_object A Seurat object whose \code{meta.data} supplies the cell
#'   type, lineage, and time point annotations. Must contain every row name of
#'   \code{cell_imputation_mat}, and also the later-time-point cells, which are
#'   generally \emph{not} rows of \code{cell_imputation_mat}.
#' @param variable_celltype Name of the \code{meta.data} column holding the cell
#'   type annotation.
#' @param variable_lineage Name of the \code{meta.data} column holding the
#'   lineage assignment (usually \code{"assigned_lineage"}).
#' @param variable_timepoint Name of the \code{meta.data} column holding the
#'   time point.
#' @param bool_10_power Whether to apply \code{10^} to
#'   \code{cell_imputation_mat} before doing anything else. Default \code{TRUE},
#'   which is correct for \code{cell_imputed_score}, since that is on the log10
#'   scale --- see the "Scales" section of \code{\link{cyfer_finalize}}. Set to
#'   \code{FALSE} only if the caller has already exponentiated.
#' @param bool_jitter Whether to add \code{Uniform(0, min_jitter)} noise to each
#'   composition before renormalizing. Default \code{TRUE}. This exists to pull
#'   points off the edges and vertices of the simplex, where a cell with a
#'   single non-zero fate would otherwise overplot, and it means the function is
#'   \bold{stochastic}: two calls give different coordinates unless the caller
#'   sets a seed beforehand. There is no \code{seed_number} argument.
#' @param entropy_bump A constant added to every entropy value at the end.
#'   Default \code{0.01}. Purely cosmetic: entropy is mapped to point size, and
#'   a lineage with a single observed fate has entropy exactly 0 and would
#'   otherwise be drawn invisibly small.
#' @param min_imputation Cells whose row of \code{cell_imputation_mat} sums to
#'   at most this are dropped, since their composition would be the ratio of two
#'   near-zero numbers. Default \code{0.01}. Applied \emph{after}
#'   \code{bool_10_power}, so it is a threshold on predicted progeny count, not
#'   on the log10 score.
#' @param min_jitter Upper bound of the jitter draw. Default \code{0.1}. Large
#'   relative to a composition that sums to 1, so the jitter is a visible
#'   perturbation rather than a nudge.
#'
#' @returns A data frame with one row per surviving cell, row names being cell
#'   IDs. Columns: one numeric column per fate (named as in
#'   \code{colnames(cell_imputation_mat)}) holding the normalized, jittered
#'   composition; \code{celltype} (factor, the cell's own annotation);
#'   \code{cellsize} (numeric, the cell's total predicted progeny count
#'   \emph{before} normalization); \code{lineage} (the cell's lineage);
#'   \code{dominant_fate} (factor, the most common observed cell type among that
#'   lineage's cells at \code{later_timepoint}, with the literal string
#'   \code{"NA"} for lineages having no such cells); and \code{entropy}
#'   (numeric, Shannon entropy in bits of that same distribution, plus
#'   \code{entropy_bump}; \code{NA} for lineages with no later-time-point
#'   cells).
#'
#' @export
compute_entropy <- function(cell_imputation_mat,
                            later_timepoint,
                            seurat_object,
                            variable_celltype,
                            variable_lineage,
                            variable_timepoint,
                            bool_10_power = TRUE,
                            bool_jitter = TRUE,
                            entropy_bump = 0.01,
                            min_imputation = 0.01,
                            min_jitter = 0.1){
  stopifnot(
    length(rownames(cell_imputation_mat)) > 0,
    length(colnames(cell_imputation_mat)) == ncol(cell_imputation_mat)
  )
  
  k <- ncol(cell_imputation_mat)
  if(bool_10_power) {
    cell_imputation_mat <- 10^cell_imputation_mat
  }
  
  cellsize <- rowSums(cell_imputation_mat)
  n <- nrow(cell_imputation_mat)
  
  for(i in 1:n){
    tmp <- cell_imputation_mat[i,]
    if(sum(tmp) <= min_imputation){
      cell_imputation_mat[i,] <- NA
    } else {
      cell_imputation_mat[i,] <- tmp/sum(tmp)
    }
  }
  
  idx <- unique(unlist(apply(cell_imputation_mat, 2, function(x){
    which(is.na(x))
  })))
  if(length(idx) > 0) {
    cell_imputation_mat <- cell_imputation_mat[-idx,,drop = FALSE]
    cellsize <- cellsize[-idx]
  }
  
  if(bool_jitter) {
    n <- nrow(cell_imputation_mat)
    for(i in 1:n){
      if(any(is.na(cell_imputation_mat[i,]))) next()
      cell_imputation_mat[i,] <- cell_imputation_mat[i,] + stats::runif(k, min = 0, max = min_jitter)
      cell_imputation_mat[i,] <- cell_imputation_mat[i,]/sum(cell_imputation_mat[i,])
    }
  }
  
  cell_imputation_mat <- cell_imputation_mat[which(!is.na(cell_imputation_mat[,1])),,drop = FALSE]

  df <- as.data.frame(cell_imputation_mat)
  metadata <- seurat_object@meta.data
  df$celltype <- metadata[rownames(cell_imputation_mat), variable_celltype]
  df$celltype <- factor(df$celltype)
  df$cellsize <- cellsize[rownames(cell_imputation_mat)]
  
  # compute the dominant fate
  df$lineage <- metadata[rownames(df),variable_lineage]
  df$dominant_fate <- rep(NA, nrow(df))
  df$entropy <- rep(NA, nrow(df))
  
  for(lineage in unique(df$lineage)){
    df_idx <- which(df$lineage == lineage)
    seurat_idx <- intersect(which(metadata[,variable_timepoint] == later_timepoint),
                            which(metadata[,variable_lineage] == lineage))
    tab_vec <- table(metadata[seurat_idx, variable_celltype])
    if(length(tab_vec) == 0) next()
    df$dominant_fate[df_idx] <- names(tab_vec)[which.max(tab_vec)]
    df$entropy[df_idx] <- .shannon_entropy(tab_vec)
  }
  df$dominant_fate[which(is.na(df$dominant_fate))] <- "NA"
  df$dominant_fate <- factor(df$dominant_fate)
  
  df$entropy <- df$entropy + entropy_bump
  
  df
}

#' Shannon entropy of a count vector, in bits
#'
#' Normalizes \code{x} to a probability vector and returns
#' \code{-sum(p * log2(p))}, with the \code{0 * log2(0) = 0} convention applied
#' explicitly (\code{log2(0)} is \code{-Inf}, so the product would otherwise be
#' \code{NaN}). Zero entries therefore contribute nothing, which makes the
#' result invariant to padding a table with unobserved categories.
#'
#' @param x A numeric vector of non-negative counts, typically a
#'   \code{table()}. A length-1 input short-circuits to 0.
#'
#' @returns A single numeric between 0 and \code{log2(length(x))}: 0 when all
#'   the mass sits on one category, the maximum when it is spread uniformly.
#'
#' @noRd
.shannon_entropy <- function(x){
  if(length(x) == 1) return(0)
  x <- x/sum(x)
  y <- log2(x)
  y[x == 0] <- 0
  -sum(x*y)
}