generate_data2 <- function(df_x, df_y, list_xnoise, list_ynoise, 
                           df_cell, blueprint_resolution = 500, verbose = T){
  
  # checks
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
  
  # setup
  n <- nrow(df_cell)
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  num_branches <- .extract_num_branches(df_y)
  stopifnot(num_branches == length(unique(df_cell$branch)))
  vec_time <- seq(min(df_cell$time), max(df_cell$time), length.out = blueprint_resolution)
  stopifnot(all(diff(vec_time) > 0))
  
  # create x's blueprint that is (vec_time) x (num. of variables) for each branch
  blueprint_x <- .create_blueprint_x(df_x, list_xnoise, num_branches, vec_time, verbose = T)
  # set the amount of biological noise in each cell
  noise_x <- .generate_noise_x(df_x, df_cell)
  # based on x's blueprint, start generating the mean expression for each cell
  mean_x <- .generate_mean_x(df_cell, blueprint_x, noise_x)
  
  # based on the amount of biological noise and blueprint for x, generate the mean for y
  mean_y <- .generate_mean_y(df_x, df_y, df_cell, blueprint_x, noise_x)
  
  # generate the observed x and y
  obs_x <- .generate_obs(df_x, mean_x)
  obs_y <- .generate_obs(df_y, mean_y)
  
  list(mat_x = res_x$mat_x, mat_y = res_y$mat_y,
       mat_meanx = res_x$mat_meanx, mat_meany = res_y$mat_meany)
}

########################

.extract_num_branches <- function(df_x){
  max(as.numeric(unlist(lapply(df_x$branch, function(x){strsplit(x, split = ",")[[1]]}))))
}

# create a list (equal to the number of branches) where each element is 
# (vec_time) x (number of peaks), denoting the mean expression across time
# for each branch
.create_blueprint_x <- function(df_x, list_xnoise, num_branches, vec_time,
                                verbose = T){
  p <- nrow(df_x)
  
  list_blueprint <- lapply(1:num_branches, function(branch){
    if(verbose) print(paste0("branch ", branch))
    mat_blueprint <- t(sapply(1:length(vec_time), function(idx){
      if(verbose && length(vec_time) >= 10 && idx %% floor(length(vec_time)/10) == 0) cat('*')
      time <- vec_time[idx]
      
      # [[note to self: there's definitely a better way to code this...]]
      sapply(1:p, function(j){
        if(df_x$branch[j] == "0"){
          .extract_exp_unrelated(time, df_x[j,], list_xnoise[[df_x$name[j]]])
          
          # informative gene
        } else if(branch %in% as.numeric(strsplit(df_x$branch[j], split = ",")[[1]])){
          .extract_exp_informative(time, df_x[j,])
          
          # gene not informative for this gene, so generate from its baseline value
        } else {
          df_x$exp_baseline[j]
        }
      })
    }))
    
    rownames(mat_blueprint) <- vec_time
    colnames(mat_blueprint) <- df_x$name
    
    mat_blueprint
  })
  
  names(list_blueprint) <- paste0("blueprint", 1:k)
  
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
  
  mat <- t(sapply(1:n, function(i){
    branch <- df_cell$branch[i]
    .generate_cell_x(time = df_cell$time[i], noise_vec = noise_x[i,], 
                     blueprint_x = blueprint_x[[branch]])
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
.generate_mean_y <- function(df_x, df_y, df_cell, blueprint_x, noise_x){
  n <- nrow(df_cell)
  
  mat <- t(sapply(1:n, function(i){
    branch <- df_cell$branch[i]
    .generate_cell_y(df_x, df_y, time = df_cell$time[i], noise_vec = noise_x[i,], 
                     blueprint_x = blueprint_x[[branch]])
  }))
  
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
  
  mat_obs
}

#########################

# [[note to self: no windup time used here]]
.extract_exp_unrelated <- function(time, vec_param, mat_interval){
  idx <- intersect(which(time >= mat_interval["time_start",]), 
                   which(time <= mat_interval["time_end",]))
  
  if(length(idx) == 1){
    vec_param["exp_max"]
  } else {
    vec_param["exp_baseline"]
  }
}

.extract_exp_informative <- function(time, vec_param){
  # determine if time is in the high region
  if(time >= vec_param["time_start"] & time <= vec_param["time_end"]){
    return(vec_param["exp_max"])
    
  # determine if time is in the windup period
  } else if(time >= vec_param["time_start"]-vec_param["time_windup"] & 
            time <= vec_param["time_end"]+vec_param["time_windup"]){
    
    if(time >= vec_param["time_start"]-vec_param["time_windup"] & time <= vec_param["time_start"]){
      frac <- abs(time - vec_param["time_start"])/vec_param["time_windup"]
      val <- (1-frac)*vec_param["exp_max"] + frac*vec_param["exp_baseline"]
    } else {
      frac <- abs(time - vec_param["time_end"])/vec_param["time_windup"]
      val <- (1-frac)*vec_param["exp_max"] + frac*vec_param["exp_baseline"]
    }
    return(val)
    
  } else {
    return(vec_param["exp_baseline"])
  }
}

.generate_cell_x <- function(time, noise_vec, blueprint_matx){
  vec_time <- as.numeric(rownames(blueprint_matx))
  
  idx_lower <- which.max(sapply(vec_time, function(x){ifelse(x <= time, x, NA)}))
  idx_upper <- which.min(sapply(vec_time, function(x){ifelse(x >= time, x, NA)}))
  
  if(length(idx_lower) == 0 || length(idx_upper) == 0 || idx_lower == idx_upper){
    idx <- ifelse(length(idx_lower) == 0, idx_upper, idx_lower)
    pmax(blueprint_matx[idx,] + noise_vec, 0)
    
  } else {
    weight <- (x-vec_time[idx_lower])/(vec_time[idx_upper] - vec_time[idx_lower])
    vec <- weight*blueprint_matx[idx_upper,] + (1-weight)*blueprint_matx[idx_lower,]
    
    pmax(vec + noise_vec, 0)
  }
}

.generate_cell_y <- function(df_x, df_y, time, noise_vec, blueprint_x){
  stopifnot(ncol(noise_vec) == ncol(blueprint_x))
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  
  sapply(1:p2, function(j){
    if(df_y$branch[j] == "0"){
      .extract_exp_unrelated(time, vec_param = df_y[j,], 
                             mat_interval = list_ynoise[[df_y$name[j]]])
      
    } else {
      # [[note to self: there could be a check on the branches here?]]
      idx_x <- which(df_x$gene == df_y$name[j])
      time_lag <- df_x$time_lag[idx_x]
      vec_time <- pmax(time - time_lag, 0)
      
      vec_x <- sapply(1:length(idx_x), function(i){
        as.numeric(.generate_cell_x(time = vec_time[i], 
                                    noise_vec = noise_vec[idx_x[i]],
                                    blueprint_matx = blueprint_x[,idx_x[i],drop=F]))
      })
      
      as.numeric(df_x$coef[idx_x] %*% vec_x)
    }
  })
}
