# output: mat_g, dataframe of when each row got recruited, # of times it was a candidate, order of recruitment, and 
# hash table of who its nearest neighbors are
chromatin_potential <- function(mat_x, mat_y, df_x, df_y, vec_start, list_end,
                                form_method = "literal", est_method = "glmnet",
                                cand_method = "nn", rec_method = "singleton", 
                                options = list(),
                                verbose = T){
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), ncol(mat_y) == nrow(df_y))
  n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  stopifnot(all(vec_start > 0), all(vec_start %% 1 == 0), all(vec_start <= n))
  for(i in 1:length(list_end)){
    stopifnot(all(list_end[[i]] > 0), all(list_end[[i]] %% 1 == 0), all(list_end[[i]] <= n))
  }
  tmp <- c(vec_start, unlist(list_end))
  stopifnot(length(unique(tmp)) == length(tmp))
  
  # check all the options
  tmp <- .chrom_options(form_method, est_method, cand_method, rec_method, options)
  form_options <- tmp$form_options; est_options <- tmp$est_options
  cand_options <- tmp$cand_options; rec_options <- tmp$rec_options
  
  # initialize
  tmp <- .init_est_matrices(mat_x, mat_y, vec_start, list_end)
  mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
  df_res <- .init_chrom_df(n, vec_start, list_end, cell_name)
  ht_neighbor <- .init_chrom_ht(list_end)
  counter <- 1
  if(est_options$enforce_cis){
    est_options <- .gene_peak_map(df_x, df_y, est_options)
  }
  
  # while:
  while(length(ht) < n){
    ## estimate res_g
    res_g <- .estimate_g(mat_x1, mat_y2, df_y, est_options)
    
    ## construct candidate set
    vec_cand <- .candidate_set(mat_x, mat_y, df_res, cand_options)
    df_res <- .update_chrom_df_cand(df_res, vec_cand)
    
    ## recruit an element from the candidate set
    rec <- .recruit_next(mat_x[vec_cand,,drop = F], mat_y1, res_g, rec_options, est_options)
    
    ## update
    tmp <- .update_estimation_matrices(mat_x1, mat_y1, mat_y2, rec, form_options)
    mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
    ht_neighbor <- .update_chrom_ht(rec$vec_from, rec$list_to)
    df_res <- .update_chrom_df_rec(df_res, rec$vec_from)
    
    counter <- counter+1
  }

  # output
  list(res_g = res_g, df_res = df_res, ht_neighbor = ht_neighbor)
}

#########################

# columns: steady-state (neg for initial, pos for end, or NA)
# num times was a candidate, and order of recruitment
.init_chrom_df <- function(n, vec_start, list_end, cell_name){
  df_res <- data.frame(idx = 1:n, init_state = rep(NA, n), num_cand = rep(0, n),
                       order_rec = rep(NA, n))
  if(length(cell_name) == n) rownames(df_res) <- cell_name
  
  df_res$init_state[vec_start] <- -1
  for(i in 1:length(list_end)){
    df_res$init_state[list_end[[i]]] <- i
    df_res$order_rec[list_end[[i]]] <- 0
  }
  
  df_res
}

.init_chrom_ht <- function(list_end){
  ht_neighbor <- hash::hash()
  vec <- unlist(list_end)
  for(i in vec){
    ht_neighbor[[as.character(i)]] <- c(neighbor = i)
  }
  
  ht_neighbor
}



########################3

# ## fate prediction model
# library(glmnet)
# library(matrixStats)
# 
# # function to find the nearest neighbors in y (row) for each cell (row) in x 
# get_nn_2sets=function(x,y, k=10){
#   # x and y should have the same number of column
#   if(ncol(x) != ncol(y)){stop("x and y should have the same dimenstion.")}
#   nn_2mat=list()
#   nn_mat=matrix(nrow=nrow(x), ncol=k)
#   dist_mat=matrix(nrow=nrow(x), ncol=k)
#   for(rr in 1:nrow(x)){
#     # compute the distance. here I used Euclidian distance
#     ##! we could try other distance later
#     dd=sqrt(rowSums((matrix(rep(x[rr,], nrow(y)), nrow=nrow(y), byrow=T)-y)^2)) 
#     dist_mat[rr,]=sort(dd)[1:k]
#     nn_mat[rr,]=order(dd)[1:k]
#   }
#   rownames(nn_mat)=rownames(dist_mat)=rownames(x)
#   
#   nn_2mat$dist_mat=dist_mat
#   nn_2mat$nn_mat=nn_mat
#   return(nn_2mat)
# }
# 
# # randomly generate cell by peak matrix, cell by gene matrix and the cell identity for testing
# atac=matrix(1:10000,nrow=100)
# rna=matrix(1:5000, nrow=100)
# rownames(atac)=paste0("c", 1:nrow(atac)) # cell by peak (100*100)
# rownames(rna)=paste0("c", 1:nrow(rna)) # cell by gene (100*50)
# celltype=sample(1:5,100,replace = T)
# 
# # a vecotr of cell name for mapping
# cellname=1:nrow(rna)
# names(cellname)=rownames(rna)
# 
# # set initial values
# sd_cell=c(1,5)
# end_cell=c(5)
# lambdas <- 10^seq(2, -3, by = -.1)
# cells_rest=1:nrow(rna)
# cellmodel=which(celltype %in% sd_cell)  #which cells included in the model training
# cellind=which(celltype %in% end_cell) # which cells in the end states
# nneighbor=10 ##! could be different each time based on some threshold 
# nsmooth=10
# rna_sub=rna[cellmodel,]
# 
# # set output matrix
# dir_mat=matrix(0,nrow=length(cellname),ncol=2) # col:if the cell is identified; col2 nearest cells 
# 
# # main algorithm
# for(tt in 1:100){ #numer of iterations
#   rna_est=matrix(nrow=(length(cellname)-length(cellmodel)), ncol=ncol(rna)) # set matrix for the estimated rna values for all other cells (not in the training set)
#   for(ii in 1:ncol(rna)){ # for each gene
#     
#     # generate the sub matrices for rna and atac for training
#     ##! we need to add constraint to the peaks for each gene later
#     atac_train=atac[cellmodel,] 
#     rna_train=rna_sub[,ii]
#     
#     # cv for selecting lambda
#     reg=cv.glmnet(atac_train, rna_train, alpha=1, lambda = lambdas, standardize=T, nfold=5)
#     lambda_best <- reg$lambda.min 
#     
#     # build the model
#     ##! we will need to select a model that works for the data
#     lasso_model <- glmnet(atac_train, rna_train, alpha = 1, lambda = lambda_best, standardize = TRUE)
#     
#     # use the trained model to predict the rna values for other cells not in the training set
#     atac_test=atac[-cellmodel,, drop=F]
#     predictions_train <- predict(lasso_model, s = lambda_best, newx = atac_test)
#     rna_est[,ii]=predictions_train
#     
#   }
#   rownames(rna_est)=rownames(rna)[-cellmodel]
#   
#   # a vector indicating which cells are not in the training set
#   cells_rest=cells_rest[which(!cells_rest %in% cellmodel)]
#   
#   # Now I select k neighbors each interation for all the time points. If the number of rest cells<k, stop the iteration
#   if(length(cells_rest)<nneighbor){
#     break
#   }
#   
#   # for each cell not in the "cellmodel", find the nsmooth nearest neighbors in the "cellind"
#   rna_nn=get_nn_2sets(rna_est,rna[cellind,], k=nsmooth)
#   
#   # select the nneighbor cells not in the "cellmodel"
#   cell_dist=rowMedians(rna_nn$dist_mat) # for each cell not in the "cellmodel", get the median distance for the nsmooth nearest neighbors in the "cellind"
#   ind=order(cell_dist)[1:nneighbor]  # ind in the rna_nn list
#   celladd=cellname[rownames(rna_est)[order(cell_dist)[1:nneighbor]]] # map the nneighbor cells to the original cell list
#   
#   # save the info
#   dir_mat[celladd,1]=dir_mat[celladd,1]+1 
#   dir_mat[celladd,2]=cellind[rna_nn$nn_mat[ind,1]] # direct the cells to the nearest neighbor
#   
#   # smooth the matched rna signals (for the nsmooth neighbors in the "cellind") for the nneighbor cells selected each iteration
#   rna_smooth=matrix(ncol=ncol(rna_sub), nrow=length(celladd))
#   for(ss in 1:length(ind)){
#     rna_smooth[ss,]=colMedians(rna[cellind,][rna_nn$nn_mat[ind[ss],],])
#   }
#   
#   # combine the original rna matrix with the smooth rna matrix
#   rna_sub=rbind(rna_sub, rna_smooth)
#   
#   # update the cells in "cellind" and "cellmodel"
#   cellind=c(cellind, celladd)
#   cellmodel=c(cellmodel,celladd)
#   
# }
# 
# ##! simulations
# ##! multiple fates
# 
# 
