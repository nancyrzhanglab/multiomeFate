nearest_neighbor <- function(mat, nn_options){
  if(nn_options[["method"]] == "annoy"){
    res <- .nearest_neighbor_annoy(mat, nn_options)
  } else {
    stop("Nearest neighbor method not found")
  }
  
  res
}

################3

.nearest_neighbor_annoy <- function(mat, nn_options){
  p <- ncol(mat); n <- nrow(mat)
  if(nn_options$metric == "euclidean"){
    nn_obj <- new(RcppAnnoy::AnnoyEuclidean, p)
  } else if(nn_options$metric == "cosine"){
    nn_obj <- new(RcppAnnoy::AnnoyAngular, p)
  } else if(nn_options$metric == "manhattan"){
    nn_obj <- new(RcppAnnoy::AnnoyManhattan, p)
  } else if(nn_options$metric == "hamming"){
    nn_obj <- new(RcppAnnoy::AnnoyHamming, p)
  } 
  
  for(i in 1:n){
    nn_obj$addItem(i-1, mat[i,])
  }
  
  nn_obj$build(nn_options$ntrees)
  
  nn_obj
}

###################

.query_nn <- function(nn_obj, nn_options){
  n <- nn_obj$getNItems()
  
  my_sapply <- ifelse(
    test = !nn_options$parallel && future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  
  nn_mat <- my_sapply(1:n, function(i){nn_obj$getNNsByItem(i-1, nn_options$nn+1)})
  nn_mat <- t(nn_mat[-1,]+1) # since RcppAnnoy starts indexing with 0, and the entry itself is always included
  
  nn_mat
}