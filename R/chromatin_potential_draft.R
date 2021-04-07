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

### example codes for plotting segments
# umap_sub$umap1 is the starting x coordinate and umap_delta_sub$umap1 is the delta_x
# umap_sub$umap2 is the starting y coordinate and umap_delta_sub$umapy is the delta_y
# library(ggplot2)
# pp=ggplot(df, aes(x=UMAP1, y=UMAP2)) +
#   geom_point(aes( color=celltype), size=0.5, alpha=0.3)+
#   geom_segment(aes(xend = umap_sub$umap1+umap_delta_sub$umap1, yend =umap_sub$umap2+umap_delta_sub$umap2),
#                arrow = arrow(length = unit(0.1, "cm")))
# pp

## example codes for plotting arrows and streams--metR
# ggplot+
# metR::geom_vector(data = uv.se, aes(x = lon, y = lat, dx = u, dy = v), 
#                     arrow.angle = 30, arrow.type = "open", arrow.length = .5, 
#                     pivot = 0,preserve.dir = TRUE, direction = "ccw")+

# ggplot+
# metR::geom_streamline(data = uv.se, 
# aes(x = lon, y = lat, dx = u, dy = v),
# L = 1.75, res = .9, n = 40, jitter = 4)
