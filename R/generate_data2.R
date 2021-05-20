# starting to port over nancy's simulator

# one hash-table: for each gene, which peaks are nearby

# one gene dataframe: one row for each gene.
# the attributes: TSS for gene, which branch (string, ex: "1" or "2,3", or "0"), which wave,
# 5 values: interval-start, interval-end, windup, low-exp, high-exp
# recall: the waves simply help us translate "which group of genes is this" to "what times are this gene maxed out"

# one peak dataframe: one row for each peak.
# the attributes: peak start and end, which branch (string, ex: "1" or "2,3" or "0"), which wave,
# 5 values: interval-start, interval-end, windup, low-exp, high-exp
## this is made after gene dataframe, so the construction of this df requires chromatin lead/lag

# one noise gene list: length of number of noise genes, named by the gene
# the elements: matrix of start&ends, vector of 4 (windup, low-exp, high-exp, prob of expressed)

# one noise peak list: length of number of noise peak, named by the peak
# the elements: matrix of start&ends, vector of 4 (windup, low-exp, high-exp, prob of expressed)

# one cell dataframe:
# one row for each cell: what its time, its branch

generate_data2 <- function(df_x, df_y, list_xnoise, list_ynoise, 
                           df_cell, blueprint_resolution = 500, options = list(),
                           verbose = T){
  n <- nrow(df_cell)
  p1 <- nrow(df_x); p2 <- nrow(df_y)
  num_branches <- .extract_num_branches(df_x)
  time_vec <- seq(min(df_cell$time), max(df_cell$time), length.out = blueprint_resolution)
  
  # create x's blueprint that is (time_vec) x (num. of variables) for each branch
  blueprint_x <- .create_blueprint_x(df_x, list_xnoise, num_branches, time_vec)
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

# create a list (equal to the number of branches) where each element is 
# (time_vec) x (number of peaks), denoting the mean expression across time
# for each branch
.create_blueprint_x <- function(df_x, list_xnoise, num_branches, time_vec){
  
}

# create a template for the biological noise in each cell, where this is a 
# matrix that is (number of cells) x (number of peaks)
.generate_noise_x <- function(df_x, df_cell){
  
}

# take the blueprint_x (dictating for this branch, at this time, what is the supposed-expression)
# and add the noise for that cell (given in noise_x). the output is a 
# a matrix that is (number of cells) x (number of peaks)
.generate_mean_x <- function(df_cell, blueprint_x, noise_x){
  
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
  
}


#########################

.generate_cell <- function(time, branch, df_var, list_noise, options){
  p <- nrow(df_var)
 
  # compute a matrix with the mean and sd across all the genes
  mat_exp <- sapply(1:p, function(j){
    # gene unrelated to trajectories
    if(df_var[j,"branch"] == 0){
      .extract_exp_unrelated(time, list_noise[[rownames(df_var)[j]]])
      
    # informative gene
    } else if(branch %in% as.numeric(strsplit(df_var[j,"branch"], split = ",")[[1]])){
      .extract_exp_informative(time, df_var[j,])
      
    # gene not informative for this gene, so generate from its baseline value
    } else {
      c(df_var[j, "baseline_exp"], df_var[j, "sd_exp"])
    }
  })
  
}

.extract_exp_unrelated <- function(time, list_param){
  
}

.extract_exp_informative <- function(time, vec_param){
  # determine if time is in the high region
  if(time >= vec_param["start_time"] & time <= vec_param["end_time"]){
    return(c(vec_param["max_exp"], vec_param["sd"]))
    
  # determine if time is in the windup period
  } else if(time >= vec_param["start_time"]-vec_param["windup"] & 
            time <= vec_param["end_time"]+vec_param["windup"]){
    
    if(time >= vec_param["start_time"]-vec_param["windup"] & time <= vec_param["start_time"]){
      frac <- abs(time - vec_param["start_time"])/vec_param["windup"]
      val <- (1-frac)*vec_param["max_exp"] + frac*vec_param["baseline_exp"]
    } else {
      frac <- abs(time - vec_param["end_time"])/vec_param["windup"]
      val <- (1-frac)*vec_param["max_exp"] + frac*vec_param["baseline_exp"]
    }
    return(c(val, vec_param["sd_exp"]))
    
  } else {
    return(c(vec_param["baseline_exp"], vec_param["sd_exp"]))
  }
}

# for an informative gene's information and the time, compute the mean expression
.extract_exp <- function(){
  
}
