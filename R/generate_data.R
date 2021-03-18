# some inputs: obj_nextcell, the initial point, max_cells
# to add in future: non-uniform density of cells along time
#
# some outputs: data frame of cells [uncorrupted atac and rna] as well as its branch and psuedotime,
# true datapoint mother
# and another data of simply cells with the noisy data
generate_data <- function(obj_next, initial_x, max_n, number_runs = 1){
  # initialize the noiseless matrix
  
  # while loop
  ## generate next y, based on x
  ## based on the next y, generate that corresponding x
  ## expand the matrix if necessary
  
  # add technical noise to matrices
  
  # prepare output
}

#############################