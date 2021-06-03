assign_graph_attributes <- function(nn_g, df_cell){
  attr_vec <- colnames(df_cell)
  attr_vec <- attr_vec[!attr_vec %in% "name"]
  
  for(attr in attr_vec){
    idx <- which(colnames(df_cell) == attr)
    nn_g <- igraph::set_vertex_attr(nn_g, name = attr, value = df_cell[,idx])
  }
 
  nn_g
}