#' Plot UMAP embedding (for \code{mf_simul})
#' 
#' Here, red denotes younger cells and blue denotes older cells.
#'
#' @param obj object of class \code{mf_simul}, from the output of \code{generate_data}
#' @param mode_x boolean, where \code{TRUE} means the data from Modality 1 is involved in the UMAP
#' @param mode_y boolean, where \code{TRUE} means the data from Modality 2 is involved in the UMAP
#' @param noiseless boolean, where \code{TRUE} means information from 
#' \code{dat$true_x} and/or \code{dat$true_y} is plotted, and \code{FALSE}
#' means information from \code{dat$obs_x} and/or \code{dat$obs_y} is plotted
#' @param k positive numeric for the number of principal scores extracted
#' @param ghost_neighbor \code{NA} or a \code{hash} object from \code{chromatin_potential}
#' that denotes which cells are matched to which other cells
#' @param col_vec \code{NA} or pre-specified vector of colors
#' @param reorder boolean. If \code{TRUE} and \code{all(is.na(col_vec))}, color the cells by their pseudotime.
#' If \code{FALSE}, color the cells by the order they appear in simulation
#' @param num_col positive numeric where if \code{all(is.na(col_vec))}, 
#' denotes the number of distinct colors used
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_umap.mf_simul <- function(obj, mode_x = T, mode_y = T, noiseless = F, k = 10,
                               ghost_neighbor = NA,
                               col_vec = NA, reorder = T, num_col = 10, ...){
  stopifnot(class(obj) == "mf_simul", mode_x | mode_y)
  
  # extract color
  if(all(is.na(col_vec))){
    if(reorder){
      vec <- obj$df_info$time
    } else {
      vec <- obj$df_info$counter
    }
    col_palette <- grDevices::colorRampPalette(c("red", "blue"))(num_col) 
    vec_val <- seq(min(vec), max(vec), length.out = 10)
    col_vec <- sapply(vec, function(x){
      col_palette[which.min(abs(x - vec_val))]
    })
  } else {
    stopifnot(length(col_vec) == nrow(obj$obs_x))
  }
  
  if(noiseless){
    mat_x <- obj$true_x; mat_y <- obj$true_y
  } else {
    mat_x <- obj$obs_x; mat_y <- obj$obs_y
  }
  
  if(class(ghost_neighbor) == "hash"){
    tmp <- .ghost_matrix(mat_x, mat_y, ghost_neighbor)
    mat_x <- tmp$mat_x; mat_y <- tmp$mat_y
  }
  
  #embedding <- .umap_embedding(mat_x, mat_y, mode_x, mode_y, k)
  embedding= umap_coord

  graphics::plot(embedding[,1], embedding[,2], pch = 16, col = col_vec, ...)
  invisible()
}

#' Plot UMAP embedding (for \code{chromatin_potential})
#'
#' @param obj object of class \code{chromatin_potential}, from the output of \code{chromatin_potential}
#' @param mode_x boolean, where \code{TRUE} means the data from Modality 1 is involved in the UMAP
#' @param mode_y boolean, where \code{TRUE} means the data from Modality 2 is involved in the UMAP
#' @param k positive numeric for the number of principal scores extracted
#' @param col_vec \code{NA} or pre-specified vector of colors for the cells
#' @param multiple_to string (\code{"ghost"} where we make "ghost" samples
#' prior to using the UMAP, or \code{"umap_avg"} where we simply average
#' among the resulting UMAP coordinates)
#' @param percent_arrows numeric
#' @param arrow_length numeric
#' @param col_arrows_by string, either \code{"order_rec"} to color arrows
#' by the recruitment order, or \code{"direction"} to color arrows 
#' based on whether they are point forward in time or backward in time
#' @param num_col_arrows if \code{NA}, color all the arrows black. If numeric,
#' then color from red to blue, where Here, red denotes matches made in the later-iterations
#' of the algorithm and blue denotes matches made in the earlier-iterations.
#' This is only used if \code{col_arrows_by="order_rec"}
#' @param vec_time \code{NA} or numeric vector of pseudotimes, with length \code{nrow(obj$mat_x)}.
#' This is only used if \code{col_arrows_by="direction"} (in which case, this input cannot be \code{NA})
#' @param ... additional graphical parameters
#'
#' @return nothing. A plot is made
#' @export
plot_umap.chromatin_potential <- function(obj, mode_x = T, mode_y = T, k = 10, 
                                          multiple_to = "ghost", col_vec = NA,
                                          percent_arrows = 0.3,
                                          arrow_length = 0.05,
                                          col_arrows_by = "order_rec",
                                          num_col_arrows = NA, 
                                          vec_time = NA, ...){
  stopifnot(class(obj) == "chromatin_potential", mode_x | mode_y, percent_arrows >= 0,
            percent_arrows <= 1, multiple_to %in% c("ghost", "umap_avg"),
            col_arrows_by %in% c("order_rec", "direction"))
  if(col_arrows_by == "direction") stopifnot(length(vec_time) %in% nrow(obj$mat_x))
  
  if(all(is.na(col_vec))){
    col_vec <- grDevices::rgb(0.5, 0.5, 0.5, 0.5)
  } else {
    stopifnot(length(col_vec) == nrow(obj$mat_x))
  }
  
  mat_x <- obj$mat_x; mat_y <- obj$mat_y; n <- nrow(mat_x)
  if(multiple_to == "ghost"){
    tmp <- .ghost_matrix(mat_x, mat_y, obj$ht_neighbor)
    mat_x <- tmp$mat_x; mat_y <- tmp$mat_y
  }
  
  embedding <- .umap_embedding(mat_x, mat_y, mode_x, mode_y, k)
  graphics::plot(embedding[1:n,1], embedding[1:n,2], pch = 16, col = col_vec, ...)
  
  if(!is.na(num_col_arrows)){
    col_palette <- grDevices::colorRampPalette(c("blue", "red"))(num_col_arrows)
    max_order <- max(obj$df_res$order_rec)
  }
  
  # add arrows
  if(percent_arrows > 0){
    n <- nrow(obj$mat_x)
    if(percent_arrows != 1){
      idx <- sample(1:n, size = round(n * percent_arrows))
    } else {
      idx <- 1:n
    }
    
    for(i in idx){
      vec_from <- embedding[i,]
      if(multiple_to == "ghost"){
        vec_to <- embedding[n+i,]
      } else {
        neigh <- obj$ht_neighbor[[as.character(i)]]
        if(length(neigh) == 1){
          vec_to <- embedding[neigh[1],]
        } else {
          vec_to <- colMeans(embedding[neigh,])
        }
      }
      
      if(col_arrows_by == "order_rec"){
        if(is.na(num_col_arrows)){
          col <- "black"
        } else {
          tmp <- ceiling(obj$df_res$order_rec[i]/max_order*num_col_arrows)
          tmp <- max(tmp, 1) # in order to not output 0
          col <- col_palette[tmp]
        }
      } else {
        time_from <- vec_time[i]
        time_to <- mean(vec_time[obj$ht_neighbor[[as.character(i)]]])
        
        col <- ifelse(time_from < time_to, "black", "green")
      }
     
      suppressWarnings(graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2], col = col,
                       length = arrow_length))
    }
  }
  
  invisible()
}

############################

.umap_embedding <- function(mat_x, mat_y, mode_x, mode_y, k){
  mat <- numeric(0)
  if(mode_x){
    if(k > ncol(mat_x)){
      tmp <- scale(mat_x)
    } else {
      tmp <- stats::prcomp(mat_x)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  if(mode_y){
    if(k > ncol(mat_y)){
      tmp <- scale(mat_y)
    } else {
      tmp <- stats::prcomp(mat_y)$x[,1:k]
    }
    
    mat <- cbind(mat, tmp)
  }
  
  # extract embedding
  embedding <- suppressWarnings(Seurat::RunUMAP(mat, verbose = F, metric="euclidean" )@cell.embeddings)
}

.ghost_matrix <- function(mat_x, mat_y, ht_neighbor){
  n <- nrow(mat_x)
  
  mat_x <- rbind(mat_x, matrix(NA, n, ncol(mat_x)))
  mat_y <- rbind(mat_y, matrix(NA, n, ncol(mat_y)))
  
  for(i in 1:n){
    neigh <- ht_neighbor[[as.character(i)]]
    mat_x[n+i,] <- colMeans(mat_x[neigh,,drop = F])
    mat_y[n+i,] <- colMeans(mat_y[neigh,,drop = F])
  }
  
  list(mat_x = mat_x, mat_y = mat_y)
}