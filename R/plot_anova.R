plot_anova <- function(seurat_object,
                       cell_imputed_score,
                       assigned_lineage_variable,
                       time_celltype_variable,
                       day_later,
                       bool_add_future_size = TRUE,
                       bool_anova = TRUE,
                       bool_mark_mean = TRUE,
                       bool_mark_max = FALSE,
                       col_all_lineages = "#E69F00",
                       col_lineages = "#999999",
                       min_lineage_size = 2,
                       num_lineages_top = 10,
                       num_lineages_bottom = 10,
                       ylab = "",
                       ylim = NA){
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  
  time_celltype <- seurat_object@meta.data[,time_celltype_variable]
  names(time_celltype) <- Seurat::Cells(seurat_object)
  stopifnot(day_later %in% time_celltype)
  
  # determine which lineages qualify to be in the plot
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_mat <- table(assigned_lineage, time_celltype)
  lineage_future_size <- tab_mat[, day_later]
  names(lineage_future_size) <- rownames(tab_mat)
  
  .plot_anova_helper(seurat_object = seurat_object,
                     cell_imputed_score = cell_imputed_score,
                     assigned_lineage_variable = assigned_lineage_variable,
                     lineage_future_size = lineage_future_size,
                     bool_anova = bool_anova,
                     bool_mark_mean = bool_mark_mean,
                     bool_mark_max = bool_mark_max,
                     col_all_lineages = col_all_lineages,
                     col_lineages = col_lineages,
                     min_lineage_size = min_lineage_size,
                     num_lineages_top = num_lineages_top,
                     num_lineages_bottom = num_lineages_bottom,
                     ylab = ylab,
                     ylim = ylim)
}

#########################

.plot_anova_helper <- function(seurat_object,
                               cell_imputed_score,
                               assigned_lineage_variable,
                               lineage_future_size,
                               bool_add_future_size = TRUE,
                               bool_anova = TRUE,
                               bool_mark_mean = TRUE,
                               bool_mark_max = TRUE,
                               col = "#E69F00",
                               min_lineage_size = 2,
                               num_lineages_top = 10,
                               num_lineages_bottom = 10,
                               ylab = "",
                               ylim = NA){
  stopifnot(length(names(cell_imputed_score)) == length(cell_imputed_score))
  
  if(any(is.na(cell_imputed_score))){
    cell_imputed_score <- cell_imputed_score[!is.na(cell_imputed_score)]
  }
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  assigned_lineage <- assigned_lineage[names(cell_imputed_score)]
  
  # filter out lineages that too small
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_vec <- table(assigned_lineage)
  tab_vec <- tab_vec[tab_vec >= min_lineage_size] # current size needs to be big enough
  passed_lineages <- names(tab_vec)
  passed_cells_names <- names(assigned_lineage)[which(assigned_lineage %in% passed_lineages)]
  cell_imputed_score <- cell_imputed_score[passed_cells_names]
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  lineage_future_size <- lineage_future_size[which(names(lineage_future_size) %in% passed_lineages)]
  
  # determine which lineages qualify to be in the plot
  lineage_names_ordered <- names(lineage_future_size)[order(lineage_future_size, decreasing = TRUE)]
  lineage_names_top <- lineage_names_ordered[1:num_lineages_top]
  lineage_names_bottom <- lineage_names_ordered[(length(lineage_names_ordered)-num_lineages_bottom+1):length(lineage_names_ordered)]
  lineage_names <- unique(c(lineage_names_top, lineage_names_bottom))
  idx <- which(lineage_vec %in% lineage_names)
  
  if(bool_add_future_size){
    lineage_names_new <- sapply(lineage_names, function(lineage_name){
      paste0(lineage_name, " (", lineage_future_size[lineage_name], ")")
    })
    lineage_vec <- plyr::mapvalues(lineage_vec, from = lineage_names, to = lineage_names_new)
    lineage_names <- lineage_names_new
  }
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score[idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  if(bool_anova) anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score)
  df <- rbind(df, df2)
  
  # compute percentage
  if(bool_anova){
    lineage_effect <- .anova_percentage(
      df = df_tmp,
      lineage_variable = "lineage",
      value_variable = "imputed_count"
    )
  }
 
  names(col_vec) <- c(lineage_names, "All")
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage, y = imputed_count))
  plot1 <- plot1 + ggplot2::geom_violin(trim = TRUE, 
                                        scale = "width", 
                                        ggplot2::aes(fill=lineage))
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = "lightgray") 
  plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                        position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                             guide = ggplot2::guide_axis(angle = 45))
  plot1 <- plot1 + ggplot2::ylab(ylab)
  
  if(!all(is.na(ylim))){
    plot1 <- plot1 + ggplot2::ylim(ylim[1], ylim[2])
  }
  
  if(bool_mark_mean) 
    plot1 <- plot1 + ggplot2::stat_summary(fun = median, geom = "crossbar", 
                                           width = 0.75, color = col)
  
  if(bool_mark_max) 
    plot1 <- plot1 + ggplot2::stat_summary(fun = max, geom = "crossbar", 
                                           width = 0.75, color = col)
 
  if(bool_anova){
    plot1 <- plot1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), 
                                             ", Lineage effect = ", lineage_effect, "%"))
  }
 
  plot1
}

##################

.anova_percentage <- function(df,
                              lineage_variable,
                              value_variable){
  stopifnot(is.factor(df[,lineage_variable]))
  imputed_count <- df[,value_variable]
  lineage <- df[,lineage_variable]
  
  total_std <- sum((imputed_count - mean(imputed_count))^2)
  
  within_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    sum((imputed_count[idx] - mean(imputed_count[idx]))^2)
  }))
  
  across_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    mean_val <- mean(imputed_count[idx])
    length(idx) * (mean_val - mean(imputed_count))^2 
  }))
  
  lineage_effect <- round(across_lineage_std/total_std*100,1)
  lineage_effect
}