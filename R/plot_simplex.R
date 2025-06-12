#' @importFrom rlang .data
plot_simplex <- function(df,
                         aes_formula, # ggplot2::aes(x = Monocyte, y = Neutrophil, z = Undifferentiated, color = celltype, size = size)
                         col_palette = NULL,
                         xlab = "x",
                         ylab = "y",
                         zlab = "z",
                         title = "title"){
  plot1 <- ggtern::ggtern(data = df,
                          mapping = aes_formula)
  plot1 <- plot1 + ggplot2::geom_point()
  
  if(!all(is.null(col_palette))){
    plot1 <- plot1 + ggplot2::scale_color_manual(values = col_palette)
  }
  
  plot1 <- plot1 + ggtern::theme_showarrows() 
  plot1 <- plot1 + ggplot2::scale_size_area() 
  plot1 <- plot1 + ggplot2::labs(x = xlab, 
                  y = ylab,
                  z = zlab,
                  title = title)
  
  plot1
}