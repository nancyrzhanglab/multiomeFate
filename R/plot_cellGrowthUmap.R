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
  
  na_cutoff <- quantile(cell_imputed_score_thres, 
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