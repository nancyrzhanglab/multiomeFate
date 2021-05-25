# generate how many genes there are per wave that belong to each branch combination
generate_combn_wave_mat <- function(g, idx_root, num_waves = 20, num_per_wave = 20){
  stopifnot(class(g) == "igraph")
  
  idx_leaves <- which(igraph::degree(g) == 1)
  stopifnot(length(unique(idx_leaves)) == length(idx_leaves),
            all(idx_leaves > 0), all(idx_leaves %% 1 == 0), 
            max(idx_leaves) == length(idx_leaves),
            !idx_root %in% idx_leaves)
  
  res <- .enumerate_combn(g, idx_root)
  stopifnot(all(res$vec_lag < num_waves), all(res$vec_lag >= 0))
  ncombn <- length(res$combn_list)
  
  # [[note to self: i need to test to make sure the order in res$combn_list
  # and res$comparison_list is sensible]]
  mat <- matrix(0, nrow = ncombn, ncol = num_waves)
  rownames(mat) <- unlist(res$combn_list)
  colnames(mat) <- paste0("w", 1:num_waves)
  mat[res$combn_list[[1]],] <- num_per_wave
  
  for(i in 1:length(res$comparison_list)){
    # [[note to self: most likely buggy for more than 10 clusters]]
    common_row <- which(rownames(mat) == paste0(sort(res$comparison_list[[i]]), 
                                                collapse = ","))
    distinct_row_vec <- sapply(res$comparison_list[[i]], function(x){
      tmp <- which(rownames(mat) == x)
      tmp
    })
    names(distinct_row_vec) <- NULL
    
    col_idx <- (res$vec_lag[i]+1):ncol(mat)
    for(j in 1:length(col_idx)){
      val <- mat[common_row, col_idx[j]]
      percentage <- j/length(col_idx)
      mat[common_row, col_idx[j]] <- round((1-percentage)*val)
      mat[distinct_row_vec, col_idx[j]] <- round(percentage*val)
    }
  }
  
  mat
}

generate_data_input <- function(combn_wave_mat, num_x_per_y = 5,
                                genome_length = 10000,
                                time_max = 100, time_on = 2*time_max/ncol(combn_wave_mat), 
                                time_windup = time_on/2, 
                                max_lag = time_on/2, min_lag = 0,
                                x_exp_baseline = 0, x_exp_max = 1,
                                x_sd_biological = 0.1, x_sd_technical = 0.01,
                                x_coef = 1, 
                                y_exp_baseline = 0, y_sd_technical = 0.01,
                                num_unrelated_x = 10, num_unrelated_y = 10,
                                x_unrelated_baseline = 0, x_unrelated_max = 1,
                                x_unrelated_intervals = 5,
                                y_unrelated_baseline = 0, y_unrelated_max = 1,
                                y_unrelated_intervals = 5){
  
  # construct the gene information
  p2 <- sum(combn_wave_mat)
  branch_vec <- unlist(lapply(1:ncol(combn_wave_mat), function(i){
    rep(rownames(combn_wave_mat), times = combn_wave_mat[,i])
  }))
  wave_length <- time_max/ncol(combn_wave_mat)
  time_start_scaffold <- unlist(lapply(1:ncol(combn_wave_mat), function(j){
    unlist(lapply(1:nrow(combn_wave_mat), function(i){
      seq((j-1)*wave_length, j*wave_length, length.out = combn_wave_mat[i,j])
    }))
  }))
  df_y <- data.frame(name = paste0("vary_",1:p2), location = round(seq(1, genome_length, length.out = p2)),
                     branch = branch_vec, time_start_scaffold = round(time_start_scaffold), 
                     time_end_scaffold = round(pmin(time_start_scaffold + time_on, time_max)),
                     exp_baseline = y_exp_baseline, exp_max = NA,
                     sd_technical = y_sd_technical)

  # construct the peak information
  spacing <- min(diff(df_y$location))
  p1 <- num_x_per_y*p2
  df_x <- data.frame(name = paste0("varx_", 1:p1), location = NA, gene = NA,
                     branch = NA,
                     time_start = NA, time_end = NA, time_lag = NA,
                     time_windup = time_windup, exp_baseline = x_exp_baseline,
                     exp_max = x_exp_max, sd_technical = x_sd_technical,
                     sd_biological = x_sd_biological, coef = x_coef)
  for(i in 1:p2){
    idx <- (num_x_per_y*(i-1)+1):(num_x_per_y*i)
    df_x$location[idx] <- round(df_y$location[i]+seq(-spacing/2, spacing/2, length.out = num_x_per_y))
    df_x$location[idx] <- pmin(pmax(df_x$location[idx], 0), genome_length)
    
    df_x$gene[idx] <- rep(df_y$name[i], num_x_per_y)
    df_x$branch[idx] <- rep(df_y$branch[i], num_x_per_y)
    df_x$time_lag[idx] <- round(seq(max_lag, min_lag, length.out = num_x_per_y))
    
    df_x$time_start[idx] <- round(pmax(df_y$time_start_scaffold[i] - df_x$time_lag[idx], 0))
    df_x$time_end[idx] <- round(pmax(df_y$time_end_scaffold[i] - df_x$time_lag[idx], time_on))
  }
  
  # add in the unrelated genes
  if(num_unrelated_y > 0){
    tmp_y <- data.frame(name = paste0("vary_",(p2+1):(p2+num_unrelated_y)), 
                        location = round(seq(1, genome_length, length.out = num_unrelated_y)),
                        branch = "0", time_start_scaffold = NA, 
                        time_end_scaffold = NA,
                        exp_baseline = y_unrelated_baseline, exp_max = y_unrelated_max,
                        sd_technical = y_sd_technical)
    df_y <- rbind(df_y, tmp_y)
    list_ynoise <- lapply(1:num_unrelated_y, function(i){
      vec <- sort(sample(1:time_max, size = 2*y_unrelated_intervals, replace = F))
      mat <- matrix(vec, nrow = 2, ncol = length(vec)/2)
      rownames(mat) <- c("time_start", "time_end")
      
      mat
    })
    names(list_ynoise) <- paste0("vary_", (p2+1):(p2+num_unrelated_y))
  } else {
    list_ynoise <- list()
  }
  
  # add in the unrelated peaks
  if(num_unrelated_x > 0){
    tmp_x <- data.frame(name = paste0("varx_", (p1+1):(p1+num_unrelated_x)), 
                        location = round(seq(1, genome_length, length.out = num_unrelated_x)),
                        gene = NA, branch = "0",
                        time_start = NA, time_end = NA, time_lag = NA,
                        time_windup = NA, exp_baseline = x_unrelated_baseline,
                        exp_max = x_unrelated_max, sd_technical = x_sd_technical,
                        sd_biological = x_sd_biological, coef = NA)
    df_x <- rbind(df_x, tmp_x)
    list_xnoise <- lapply(1:num_unrelated_x, function(i){
      vec <- sort(sample(1:time_max, size = 2*x_unrelated_intervals, replace = F))
      mat <- matrix(vec, nrow = 2, ncol = length(vec)/2)
      rownames(mat) <- c("time_start", "time_end")
      
      mat
    })
    names(list_xnoise) <- paste0("varx_", (p1+1):(p1+num_unrelated_x))
  } else {
    list_xnoise <- list()
  }
  
  list(df_x = df_x, df_y = df_y, list_xnoise = list_xnoise, 
       list_ynoise = list_ynoise)
}

generate_df_cell <- function(n, time_max, num_branch){
  df <- expand.grid(round(seq(0, time_max, length.out = ceiling(n/num_branch))),
                    1:num_branch)
  if(nrow(df) >= n) df <- df[1:n,]
  colnames(df) <- c("time", "branch")
  df$name <- paste0("cell_", 1:nrow(df))
  df <- df[,c("name", "time", "branch")]
  df
}

generate_df_cell_random <- function(n_list, func_list, time_max){
  list_time <- lapply(1:length(n_list), function(branch){
    sapply(1:n_list[[branch]], function(i){
      func_list[[i]]()*time_max
    })
  })
  
  df <- do.call(rbind, lapply(1:length(list_time), function(i){
    data.frame(time = list_time[[i]], branch = rep(i, length(list_time[[i]])))
  }))
  df$name <- paste0("cell_", 1:nrow(df))
  df <- df[,c("name", "time", "branch")]
  df
}

###############

.enumerate_combn <- function(g, idx_root){
  # set up queue
  igraph::V(g)$name <- paste0("v", 1:igraph::vcount(g))
  q <- dequer::queue()
  dequer::pushback(q, list(g = g, root = paste0("v", idx_root)))
  
  # set up combn_list
  combn_list <- list()
  comparison_list <- list(); vec_lag <- numeric(0)
  
  while(length(q) > 0){
    # pop and determine all the leaves
    res <- dequer::pop(q)
    idx_leaves <- names(which(igraph::degree(res$g) == 1))
    if(length(idx_leaves) == 0){
      idx_leaves <- igraph::V(res$g)$name[1]
    }
    combn_list[[length(combn_list) + 1]] <- paste0(idx_leaves, collapse = ",")
    
    # find all the next roots
    if(length(idx_leaves) != 1){
      next_roots <- names(igraph::neighbors(res$g, v = res$root))
      tmp <- igraph::vertex_attr(res$g, name = "lag", index = res$root)
      vec_lag <- c(vec_lag, ifelse(is.na(tmp), 0, tmp))
      g_next <- igraph::delete_vertices(res$g, v = res$root)
      comp <- igraph::components(g_next)
      vec <- numeric(0)
      
      for(i in 1:comp$no){
        idx <- names(comp$membership[which(comp$membership == i)])
        g_sub <-igraph::induced_subgraph(res$g, idx)
        
        # [[note to self: ugly coding]]
        idx_leaves <- names(which(igraph::degree(g_sub) == 1))
        if(length(idx_leaves) == 0){
          idx_leaves <- names(which(igraph::degree(g_sub) == 0))
        }
        vec <- c(vec, paste0(idx_leaves, collapse = ","))
        
        tmp_root <- next_roots[which(next_roots %in% igraph::V(g_sub)$name)]
        stopifnot(length(tmp_root) == 1)
        dequer::pushback(q, list(g = g_sub, root = tmp_root))
      }
      
      comparison_list[[length(comparison_list)+1]] <- vec
    }
  }
  
  # [[note to self: test things still work if there's more than 10 clusters]]
  combn_list <- lapply(combn_list, function(x){gsub("v", "", x)})
  for(i in 1:length(comparison_list)){
    comparison_list[[i]] <- sapply(comparison_list[[i]], function(x){
      gsub("v", "", x)
    })
    names(comparison_list[[i]]) <- NULL
  }
  
  list(combn_list = combn_list, comparison_list = comparison_list,
       vec_lag = vec_lag)
}
