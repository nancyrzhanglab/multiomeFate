generate_data <- function(df_x, df_y, list_xnoise, list_ynoise, 
                           df_cell, blueprint_resolution = 500, verbose = T){
  
  # checks
  .check_data_inputs(df_x, df_y, list_xnoise, list_ynoise, df_cell)
  
  # setup
  n <- nrow(df_cell)
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  num_branches <- .extract_num_branches(df_y)
  stopifnot(num_branches == length(unique(df_cell$branch)))
  vec_time <- seq(min(df_cell$time), max(df_cell$time), length.out = blueprint_resolution)
  stopifnot(all(diff(vec_time) > 0))
  
  # create x's blueprint that is (vec_time) x (num. of variables) for each branch
  blueprint_x <- .create_blueprint_x(df_x, list_xnoise, num_branches, vec_time)
  # set the amount of biological noise in each cell
  noise_x <- .generate_noise_x(df_x, df_cell)
  # based on x's blueprint, start generating the mean expression for each cell
  mean_x <- .generate_mean_x(df_cell, blueprint_x, noise_x)
  
  # based on the amount of biological noise and blueprint for x, generate the mean for y
  mean_y <- .generate_mean_y(df_x, df_y, df_cell, blueprint_x, noise_x, list_ynoise)
  
  # generate the observed x and y
  obs_x <- .generate_obs(df_x, mean_x)
  obs_y <- .generate_obs(df_y, mean_y)
  
  list(obs_x = obs_x, obs_y = obs_y, mean_x = mean_x, mean_y = mean_y,
       blueprint_x = blueprint_x)
}

########################

.extract_num_branches <- function(df_x){
  max(as.numeric(unlist(lapply(df_x$branch, function(x){strsplit(x, split = ",")[[1]]}))))
}

# create a list (equal to the number of branches) where each element is 
# (vec_time) x (number of peaks), denoting the mean expression across time
# for each branch
.create_blueprint_x <- function(df_x, list_xnoise, num_branches, vec_time){
  p <- nrow(df_x)
  
  list_blueprint <- lapply(1:num_branches, function(branch){
    mat_blueprint <- sapply(1:p, function(j){
      tmp <- rep(df_x$exp_baseline[j], length(vec_time))
    })
    rownames(mat_blueprint) <- vec_time
    colnames(mat_blueprint) <- df_x$name
    
    mat_blueprint
  })
  names(list_blueprint) <- paste0("blueprint_", 1:num_branches)
  
  
  for(j in 1:p){
    if(df_x$branch[j] == "0"){
      idx <- .extract_unrelated_idx(vec_time, list_xnoise[[df_x$name[j]]])
      for(branch in 1:length(list_blueprint)){
        list_blueprint[[branch]][idx,j] <- df_x$exp_max[j]
      }
      
    } else {
      branches <- as.numeric(strsplit(df_x$branch[j], split = ",")[[1]])
      tmp <- .compute_expression_values(vec_time, df_x[j,])
      
      for(branch in branches){
        list_blueprint[[branch]][tmp$idx,j] <- tmp$val
      }
    }
  }
  
  list_blueprint
}

# create a template for the biological noise in each cell, where this is a 
# matrix that is (number of cells) x (number of peaks)
.generate_noise_x <- function(df_x, df_cell){
  n <- nrow(df_cell); p <- nrow(df_x)
  
  mat <- sapply(1:p, function(j){
    stats::rnorm(n, mean = 0, sd = df_x$sd_biological)
  })
  
  colnames(mat) <- df_x$name
  rownames(mat) <- df_cell$name
  
  mat
}

# take the blueprint_x (dictating for this branch, at this time, what is the supposed-expression)
# and add the noise for that cell (given in noise_x). the output is a 
# a matrix that is (number of cells) x (number of peaks)
.generate_mean_x <- function(df_cell, blueprint_x, noise_x){
  n <- nrow(df_cell); p <- ncol(noise_x)
  
  # [[note to self: hard coded to be the same across all elements in blueprint_x]]
  vec_time <- as.numeric(rownames(blueprint_x[[1]]))
  
  mat <- t(sapply(1:n, function(i){
    branch <- df_cell$branch[i]
    idx <- which.min(abs(vec_time - df_cell$time[i]))
    pmax(blueprint_x[[branch]][idx,] + noise_x[i,], 0)
  }))
  
  colnames(mat) <- colnames(noise_x)
  rownames(mat) <- df_cell$name
  
  mat
}

# create a matrix that is (number of cells) x (number of y variables).
# to do this, for a particular cell, look at it's time and branch. 
# Based on df_y, determine if a gene is informative for this time/branch.
# Then, look up in df_x which peaks are associated with said gene and what the 
# time-offset is. then, for the peaks' offset time, look at blueprint_x for
# what that peak's supposed-expression is at said time, and add on the cell's
# biological noise. This gives this cell's mean-x vector "earlier in time." 
# Then, based on the peak's coefficients in df_x, compute the mean-y vector
.generate_mean_y <- function(df_x, df_y, df_cell, blueprint_x, noise_x, 
                             list_ynoise){
  n <- nrow(df_cell); p1 <- nrow(df_x); p2 <- nrow(df_y)
  vec_time <- df_cell$time
  
  mat <- sapply(1:p2, function(j){
    vec <- rep(df_y$exp_baseline[j], n)
    
    if(df_y$branch[j] == "0"){
      idx <- .extract_unrelated_idx(vec_time, list_ynoise[[df_y$name[j]]])
      vec[idx] <- df_y$exp_max[j]
      
    } else {
      branches <- as.numeric(strsplit(df_y$branch[j], split = ",")[[1]])
      x_idx <- which(df_x$gene == df_y$name[j])
      
      # handle cells in the gene's non-informative branches
      non_idx <- which(!df_cell$branch %in% branches)
      if(length(non_idx) > 0){
        mat <- sapply(1:length(x_idx), function(j){
          df_x$exp_baseline[x_idx[j]] + noise_x[non_idx, x_idx[j]]
        })
        mat <- pmax(mat, 0)
        vec[non_idx] <- as.numeric(mat %*% df_x$coef[x_idx])
      }
   
      # handle cells in the gene's informative branch
      for(k in 1:length(branches)){
        branch <- branches[k]
        idx <- which(df_cell$branch %in% branch)
        val <- .compute_y_expression(time_cell = df_cell$time[idx],
                                     mat_noise = noise_x[idx, x_idx],
                                     df_param = df_x[x_idx,],
                                     blueprint_matx = blueprint_x[[branch]][,x_idx,drop = F])
        
        vec[idx] <- val
      }
    }
    
    vec
  })
  
  colnames(mat) <- df_y$name
  rownames(mat) <- df_cell$name
  
  mat
}

.generate_obs <- function(df, mat){
  n <- nrow(mat); p <- nrow(df)
  
  mat_obs <- sapply(1:p, function(j){
    mat[,j] + stats::rnorm(n, mean = 0, sd = df$sd_technical[j])
  })
  
  colnames(mat_obs) <- colnames(mat)
  rownames(mat_obs) <- rownames(mat)
  
  mat_obs <- pmax(mat_obs, 0)
  
  mat_obs
}

#########################

.check_data_inputs <- function(df_x, df_y, list_xnoise, list_ynoise, df_cell){
  stopifnot(all(c("name", "gene", "time_start", "time_end", "time_windup", "time_lag",
                  "exp_baseline", "exp_max", "branch", "sd_biological", 
                  "sd_technical", "coef") %in% colnames(df_x)))
  stopifnot(all(c("name", "branch", "sd_technical", "exp_max", "exp_baseline") %in% colnames(df_y)))
  stopifnot(all(c("name", "branch", "time") %in% colnames(df_cell)))
  stopifnot(all(df_x$gene[!is.na(df_x$gene)] %in% df_y$name))
  stopifnot(length(which(df_x$branch == "0")) == length(list_xnoise),
            length(which(df_y$branch == "0")) == length(list_ynoise))
  stopifnot(all(df_x$time_start[df_x$branch != "0"] <= df_x$time_end[df_x$branch != "0"]), 
            all(df_x$exp_baseline[df_x$branch != "0"] <= df_x$exp_max[df_x$branch != "0"]))
  for(i in 1:length(list_xnoise)){
    tmp <- list_xnoise[[i]]; len <- ncol(tmp)
    stopifnot(all(tmp["time_start",] <= tmp["time_end",]))
    if(len > 1) stopifnot(all(tmp["time_start",-1] >= tmp["time_end",-len]))
  }
  for(i in 1:length(list_ynoise)){
    tmp <- list_ynoise[[i]]; len <- ncol(tmp)
    stopifnot(all(tmp["time_start",] <= tmp["time_end",]))
    if(len > 1) stopifnot(all(tmp["time_start",-1] >= tmp["time_end",-len]))
  }
  
  invisible()
}

.extract_unrelated_idx <- function(vec_time, mat_noiseparam){
  unlist(lapply(1:ncol(mat_noiseparam), function(j){
    intersect(which(vec_time >= mat_noiseparam["time_start",j]),
              which(vec_time <= mat_noiseparam["time_end",j]))
  }))
}

.compute_expression_values <- function(vec_time, vec_param){
  idx_high <- intersect(which(vec_time >= vec_param$time_start), 
                        which(vec_time <= vec_param$time_end))
  val_high <- rep(vec_param$exp_max, length(idx_high))
  
  # [[note to self: assumes a constant change in vec_time]]
  idx_spacing <- ceiling(vec_param$time_windup/min(abs(diff(vec_time))))
  val_change <- (vec_param$exp_max - vec_param$exp_baseline)/idx_spacing
    
  idx_up <- intersect(which(vec_time >= vec_param$time_start - vec_param$time_windup),
                      which(vec_time < vec_param$time_start))
  if(length(idx_up) > 0){
    val_up <- vec_param$exp_max-c(length(idx_up):1)*val_change
  } else {
    val_up <- numeric(0)
  }
 
  idx_down <- intersect(which(vec_time >= vec_param$time_end),
                        which(vec_time < vec_param$time_end + vec_param$time_windup))
  if(length(idx_down) > 0){
    val_down <- vec_param$exp_max-c(1:length(idx_down))*val_change
  } else {
    val_down <- numeric(0)
  }
  
  list(idx = c(idx_up, idx_high, idx_down), val = c(val_up, val_high, val_down))
}

# [[note to self: assumes the rownames in blueprint_matx are a fixed spacing]]
.compute_y_expression <- function(time_cell, mat_noise, df_param, blueprint_matx){
  stopifnot(nrow(df_param) == ncol(blueprint_matx), nrow(df_param) == ncol(mat_noise),
            length(time_cell) == nrow(mat_noise))
  
  time_blueprint <- as.numeric(rownames(blueprint_matx))
  time_diff <- min(diff(time_blueprint))
  time_lag <- df_param$time_lag
  idx_diff <- round(time_lag/time_diff)
  
  idx_cell <- sapply(time_cell, function(x){
    which.min(abs(x - time_blueprint))
  })
  
  mat <- sapply(1:ncol(blueprint_matx), function(j){
    blueprint_matx[pmax(idx_cell-idx_diff[j],1), j]
  })
  
  mat <- mat + mat_noise
  mat <- pmax(mat, 0)
  
  as.numeric(mat %*% df_param$coef)
}
