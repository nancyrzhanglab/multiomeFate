#' @importFrom rlang .data
plot_lineageScatterplot <- function(lineage_future_count,
                                    lineage_imputed_count,
                                    num_lineage = 10,
                                    threshold_x = 1.5,
                                    threshold_y = 1.5,
                                    title = ""){
  
  stopifnot(length(names(lineage_imputed_count)) > 0,
            length(names(lineage_future_count)) > 0,
            all(names(lineage_imputed_count) == names(lineage_future_count)))
  
  lineage_imputed_count2 <- log10(lineage_imputed_count+1)
  lineage_future_count2 <- log10(lineage_future_count+1)
  
  labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
  if(!is.na(threshold_y)){
    labeling_vec[intersect(which(lineage_imputed_count2 >= threshold_y),
                           order(lineage_imputed_count2, decreasing = TRUE)[1:num_lineage])] <- TRUE
  }
  if(!is.na(threshold_x)){
    labeling_vec[intersect(which(lineage_future_count2 >= threshold_x),
                           order(lineage_future_count2, decreasing = TRUE)[1:num_lineage])] <- TRUE
  }
  
  n <- length(lineage_imputed_count2)
  df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                   lineage_future_count = log10(lineage_future_count + 1 + 
                                                  stats::runif(n, min = 0, max = 0.5)),
                   name = names(lineage_imputed_count2),
                   labeling = labeling_vec)
  # put all the labeling == TRUE on bottom
  if(any(df$labeling)) df <- df[c(which(!df$labeling), which(df$labeling)),]
  
  df$labeling <- factor(df$labeling, levels = c(FALSE, TRUE))
  lims <- range(c(df$lineage_imputed_count, df$lineage_future_count), na.rm = TRUE)
  
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lineage_future_count, 
                                            y = .data$lineage_imputed_count))
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = .data$labeling))
  plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
  plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, .data$labeling == TRUE),
                                            ggplot2::aes(label = .data$name, 
                                                         color = .data$labeling),
                                            box.padding = ggplot2::unit(0.5, 'lines'),
                                            point.padding = ggplot2::unit(1.6, 'lines'),
                                            max.overlaps = 50)
  plot1 <- plot1 + ggplot2::ggtitle(paste0(
    title,
    "\nCorr:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
  )
  plot1 <- plot1 + ggplot2::xlab("Observed lineage count (Log10, jittered)") + ggplot2::ylab("Predicted lineage count (Log10)")
  plot1 <- plot1 + Seurat::NoLegend() + ggplot2::coord_fixed(ratio = 1, xlim = lims, ylim = lims)
  
  plot1
}