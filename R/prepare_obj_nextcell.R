# some inputs: 
# - df_x (data frame of inputs for modal 1: name, location)
# - df_y (name, location, baseline expression [the intercept])
# - matrix of regression coef for g, 
# - list of mat_1 [can be proportions between [0,1]], 
# - list of mat_2 [can be proportions between [0,1]],
# - igraph for branching structure
#
# some outputs: essentially all the ingredients for the simulation:
# - df_x, df_y
# - the matrix of coefficients formed by g
# - a huge hash table that stores each unique element of "g(mat_1)"
# and contains 1) probability of which traj to use, and 2) that traj's logistic function coefficients,
# and 3) the psuedotime
# - a huge matrix of all the elements in "g(mat_1)", corresponding to the hash table 
# [note: in the future, replace this with an exposed C++ obj from RANN: https://github.com/jefferislab/RANN/blob/master/R/nn.R]
# WARNING: We'll code as if there's no branching for now. But in the future, it'll prob require putting information in the nodes of 
#  \code{branching_graph}, and we'll need fancy functions to grab the correct rows, etc.
prepare_obj_nextcell <- function(df_x, df_y, mat_g, list_x1, list_x2, 
                                 branching_graph = NA, resolution = 0.1, max_y = 1e5){
  
  stopifnot(length(list_x1) == length(list_x2), all(sapply(length(list_x1), function(i){all(dim(list_x1[[i]]) == dim(list_x2[[i]]))})))
  if(is.na(branching_graph)) stopifnot(length(list_x1) == 1) else {
    stopifnot(class(branching_graph) == "igraph", igraph::vcount(branching_graph) == length(list_x1), igraph::components(branching_graph)$no == 1)
  }
  
  # [in the future: put a function here to check branching_graph wrt list_x1 and list_x2]
  
  # compute pseudotime of each rows in list_x1
  
  # initialize hash table
  ## count how many unique rows there are
  ## record the address of the vector [which element of list_x1, and which row?]
  ## populate the nearest neighbor matrix
  
  # start estimating the logistic regression coefficients
  ## grab the relevant rows 
  ## [for now: hard-code the fact it's a linear trajectory]
  ## perform logistic regression
  ## store the values
  
  # prepare outputs
}

generate_ygivenx <- function(obj_next, x){
  # generate y from x
  
  # add the intercepts given by dat_y
}

generate_xgiveny <- function(obj_next, y){
  # find the nearest neighbor 
  
  # grab the information from the hash table
  
  # determine which traj to use
  
  # use logistic regression
  
  # prepare output (include the pseudotime from the hashtable)
}

####################################

.possion_ygivenx <- function(x, mat_g, max_val = 1e5){
  as.numeric(pmax(exp(mat_g %*% x), max_val))
}

.glmnet_logistic <- function(covariate, response_prob){
  
}