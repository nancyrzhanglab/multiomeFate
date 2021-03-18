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
prepare_obj_nextcell <- function(){

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