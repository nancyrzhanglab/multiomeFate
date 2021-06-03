plot_igraph <- function(nn_g, color_by = NA, bool_default_par = T,
                        bool_continuous = F){
  if(bool_default_par){
    graphics::par(mfrow = c(1,1), mar = rep(0.5, 4))
  }
  
  if(is.na(color_by)){
    graphics::plot(nn_g, vertex.label = NA)
  } else {
    attr_vec <- igraph::vertex_attr(nn_g, name = color_by) 
    
    if(!bool_continuous){
      nlevels <- length(unique(attr_vec))
      color_palette <- scales::hue_pal()(nlevels)
      color_vec <- color_palette[as.numeric(as.factor(attr_vec))]
    } else {
      attr_vec <- as.numeric(attr_vec)
      color_palette <- grDevices::colorRampPalette(c("red", "blue"))(10)
      quantile_vec <- stats::quantile(attr_vec, probs = seq(0, 1, length.out = 10))
      color_vec <- color_palette[sapply(attr_vec, function(x){which.min(abs(x - quantile_vec))})]
    }
    
    graphics::plot(nn_g, vertex.label = NA, vertex.color = color_vec)
  }
 
  invisible()
}