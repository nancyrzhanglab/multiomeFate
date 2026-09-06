#' UMAP coloured by per-cell fate potential
#'
#' Puts \code{cell_imputed_score} onto an existing embedding, so that fate
#' potential can be read against transcriptional state. Cells absent from
#' \code{cell_imputed_score} --- typically the future time point's cells, which
#' were never scored --- are drawn in \code{na_color}.
#'
#' Two thresholds are applied to the colour scale, both silently, and both
#' matter when comparing panels:
#' \itemize{
#'   \item Scores are \bold{winsorized at their 99th percentile}, so a handful
#'     of very high cells cannot flatten the rest of the scale. The top 1\% of
#'     cells are therefore drawn at an indistinguishable colour.
#'   \item Cells below the \bold{5th percentile} are handed to
#'     \code{na_cutoff}, so they take \code{na_color} rather than the low end of
#'     the gradient --- the low tail is greyed out, not coloured.
#' }
#' Both cutoffs are recomputed per call from the supplied scores, so two panels
#' built from different fits do not share a scale.
#'
#' The scores are written into the object as the metadata columns
#' \code{cell_imputed_score} (unwinsorized) and \code{tmp} (winsorized, the one
#' actually plotted). The object is modified only locally.
#'
#' @param seurat_object A Seurat object carrying the reduction named by
#'   \code{reduction}.
#' @param cell_imputed_score A named numeric vector of per-cell fate potentials,
#'   names being cell IDs --- the \code{cell_imputed_score} element of
#'   \code{\link{cyfer_finalize}}, on the log10 scale. Need not cover every cell.
#' @param colors_use List of colours passed to
#'   \code{scCustomize::FeaturePlot_scCustom()} as the gradient, low to high.
#'   Default \code{list("blue", "lightgray", "red")}.
#' @param na_color Colour for unscored cells and for cells below the 5th
#'   percentile. Default \code{"bisque"}.
#' @param reduction Name of the reduction to plot. Default \code{"umap"}.
#' @param title Plot title. Default \code{""}.
#' @param order Whether to draw high-value cells on top of low-value ones.
#'   Default \code{TRUE}, which makes high-potential cells visible in dense
#'   regions but overstates how many of them there are.
#'
#' @returns A \code{ggplot} object.
#'
#' @examples
#' set.seed(10)
#' count_mat <- matrix(stats::rpois(50 * 120, lambda = 3), nrow = 50,
#'                     dimnames = list(paste0("gene", 1:50), paste0("cell", 1:120)))
#' seurat_object <- suppressWarnings(Seurat::CreateSeuratObject(counts = count_mat))
#' seurat_object$assigned_lineage <- rep(paste0("L", 1:12), length.out = 120)
#' seurat_object$time_celltype <- rep(c("day0", "day7"), each = 60)
#' # scores for the day0 cells only, as cyfer_finalize() would return them
#' score_vec <- stats::rnorm(60)
#' names(score_vec) <- colnames(seurat_object)[1:60]
#' # any 2-D embedding stored under the name passed as `reduction`
#' umap_mat <- matrix(stats::rnorm(120 * 2), nrow = 120, ncol = 2,
#'                    dimnames = list(colnames(seurat_object), c("UMAP_1", "UMAP_2")))
#' seurat_object[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat,
#'                                                          key = "UMAP_",
#'                                                          assay = "RNA")
#' plot_cellGrowthUmap(seurat_object = seurat_object,
#'                     cell_imputed_score = score_vec,
#'                     reduction = "umap")
#' @export
plot_cellGrowthUmap <- function(seurat_object,
                                cell_imputed_score,
                                colors_use = list("blue", "lightgray", "red"),
                                na_color = "bisque",
                                reduction = "umap",
                                title = "",
                                order = TRUE){
  
  cell_imputed_score_full <- rep(NA, ncol(seurat_object))
  names(cell_imputed_score_full) <- colnames(seurat_object)
  cell_imputed_score_full[names(cell_imputed_score)] <- cell_imputed_score
  seurat_object$cell_imputed_score <- cell_imputed_score_full
  
  max_val <- stats::quantile(cell_imputed_score_full, 
                             probs = 0.99, 
                             na.rm = TRUE)
  cell_imputed_score_thres <- pmin(cell_imputed_score_full, max_val)
  
  seurat_object$tmp <- cell_imputed_score_thres
  
  na_cutoff <- stats::quantile(cell_imputed_score_thres, 
                               probs = 0.05, 
                               na.rm = TRUE)
  plot1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                             colors_use = colors_use,
                                             na_cutoff = na_cutoff,
                                             na_color = na_color,
                                             reduction = reduction, 
                                             features = "tmp",
                                             order = order)
  plot1 <- plot1 + ggplot2::ggtitle(title)
  plot1
}