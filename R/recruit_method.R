.recruit_next <- function(mat_x, vec_cand, mat_y1, idx1, res_g, rec_options){
  stopifnot(all(idx1 <= nrow(mat_x)), length(idx1) == nrow(mat_y1), length(idx1) == length(unique(idx1)),
            all(idx1 %% 1 == 0), all(idx1 > 0), all(idx1 <= nrow(mat_x)))
  stopifnot(all(vec_cand %% 1 == 0), all(vec_cand > 0), all(vec_cand <= nrow(mat_x)),
            length(vec_cand) == length(unique(vec_cand)))
  stopifnot(!any(vec_cand %in% idx1))
  
  if(rec_options[["method"]] == "nn"){
    res <- .recruit_next_nn(mat_x, vec_cand, mat_y1, idx1, res_g, rec_options)
  } else {
    stop("Recruit method not found")
  }

  res
}

###################

.recruit_next_nn <- function(mat_x, vec_cand, mat_y1, idx1, res_g, 
                                    rec_options){
  num_rec <- min(rec_options$num_rec, length(vec_cand))
  nn <- min(c(rec_options$nn, ceiling(nrow(mat_y1)/2)))
  
  # apply mat_g to mat_x
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g)
  
  # see which prediction is closest to mat_y1
  res <- RANN::nn2(mat_y1, query = pred_y, k = nn)
  idx <- order(apply(res$nn.dist, 1, mean), decreasing = F)[1:num_rec]
  
  vec_from <- vec_cand[idx]
  list_to <- lapply(idx, function(i){idx1[res$nn.idx[i,]]})
  list(vec_from = vec_from, list_to = list_to)
}

#########################

# [[note to self: I'm not sure about this function name, also, Poisson hard-coded right now]]
.predict_yfromx <- function(mat_x, res_g){
  stopifnot(c("vec_g", "mat_g") %in% names(res_g))
  
  p2 <- ncol(res_g$mat_g)
  nat_param <- mat_x %*% res_g$mat_g
  
  #[[note to self: There's gotta be a cleaner way to do this]]
  for(j in 1:p2){
    nat_param[,j] <- nat_param[,j] + res_g$vec_g[j]
  }
  
  exp(nat_param)
}