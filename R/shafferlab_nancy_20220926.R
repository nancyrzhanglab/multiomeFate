# ------------------------------------------------------------- #
# Builds off shafferlab_nancy_20220609.
# Assigns cells to lineages based on model 3.
# Identifies barcodes that co-occur within the same cells, merge.
# Then, assess ATAC distance between cells of the same lineage.
# converge on a barcode assignment for downstream analysis.
# ------------------------------------------------------------- #

setwd("C:/Users/Nancy R. Zhang/Dropbox (Penn)/Projects/MultiomeFate/code")

library(Seurat)
library(Signac)
library(ggplot2)

load("../ShafferLab_Data/Writeup4e_timeAll_peakmerging_complete.RData")  # takes a few minutes.

# ------------------------------------------------------------- #
# Get the X matrix
# ------------------------------------------------------------- #


### Get the lineage barcode counts matrix.
X=all_data@assays$Lineage  # rows are lineage, cols are cells
rsX=rowSums(X)   # how many reads have this lineage barcode across all cells
csX=colSums(X)   # how many barcodes are in this cell.
X=X[rsX>3, ]  # KEEP ALL CELLS! But keep lineages with at least 3 reads.

X=as.matrix(X)
X=t(X)       # rows are now cells, cols are lineages.
ncells=nrow(X)
nlineages=ncol(X)
libsize=all_data$nCount_RNA 



# ------------------------------------------------------------- #
# Co-occurring barcodes?  
# ------------------------------------------------------------- #

barcodeName=colnames(X)
Xsqrt = sqrt(X)

Xcor=cor(Xsqrt)
diag(Xcor)=0

highcor=which(Xcor>0.6, arr.ind=TRUE)
plot(highcor[,1], highcor[,2], main=paste("Correlation > 0.6, total ", length(unique(highcor[,1])), " lineages involved", sep=""))


ind=ind+1
plot(X[,highcor[ind,1]], X[,highcor[ind,2]], pch=17, cex=1.5, 
     xlab=paste("barcode ", barcodeName[highcor[ind,1]]), ylab=paste("barcode ", barcodeName[highcor[ind,1]]), 
     main=paste("Correlation=", format(Xcor[highcor[ind,1], highcor[ind,2]], digits=2)))
grid(); abline(h=0); abline(v=0); abline(0,1,col="red")





# ------------------------------------------------------------- #
# Compute Bhat matrix.
# ------------------------------------------------------------- #



### Estimate beta0, beta1, exploratory analysis of beta0, beta1.
cell_max_barcode = apply(X, 1, which.max)
cell_max_barcode_count = apply(X, 1, max)

beta0_in_nonzero=rep(NA, nlineages)
beta0=rep(NA, nlineages)
beta1=rep(NA, nlineages)
n_won = rep(NA, nlineages)
n_lost_nonzero = rep(NA, nlineages)
n_nonzero = rep(NA, nlineages)
for(b in 1:nlineages){
  if(b %% 100==0) cat(b," out of ", nlineages," done.\n", sep="")
  won_cells = which(X[,b] == cell_max_barcode_count & X[,b]>0)
  nonzero_cells = which(X[,b]>0)
  zero_cells=which(X[,b]==0)
  lost_cells = which(X[,b]>0 & X[,b] != cell_max_barcode_count)
  n_won[b] = length(won_cells)
  n_nonzero[b] = length(nonzero_cells)
  n_lost_nonzero[b] = length(lost_cells)
  if(length(lost_cells)>0) beta0_in_nonzero[b] = mean(X[lost_cells,b]/libsize[lost_cells])
  if(length(won_cells)>0) beta1[b] = mean(X[won_cells,b]/libsize[won_cells])
  beta0[b] = mean(X[c(lost_cells, zero_cells), b]/libsize[c(lost_cells, zero_cells)])
}

beta1_mean=mean(beta1[!is.na(beta1)])
beta0_in_nonzero_mean=mean(beta0_in_nonzero[!is.na(beta0_in_nonzero)])

### Various diagnostic plots of the  beta0, beta1 distributions.

par(mfrow=c(2,2))
plot(log(n_won+1), beta1,main="Does E[X_cb|X_cb>0, B_cb=1] depend on n_won?")  # does the beta1 estimate depend on n_won?  (smaller n_won gives higher variance in estimate.)
grid()
abline(h=beta1_mean, col="red")

plot(log(n_lost_nonzero+1), beta0_in_nonzero, main="Does E[X_cb|X_cb>0, B_cb=0] depend on n_lost_nonzero?")  # does the background rate among nonzeros depend on n_lost_nonzeros?
grid()
abline(h=beta0_in_nonzero_mean, col="red")

plot(n_won+1, beta0, main="Does background rate depend on n_won?"); grid()  # does the background rate depend on n_won?
plot(log(n_won+1), log(beta0),main="Does background rate depend on n_won (LOG)?"); grid()  # does the background rate depend on n_won?

par(mfrow=c(1,1))
plot(log10(beta0_in_nonzero), log10(beta0)); grid(); abline(0,1,col="red")
sum(is.na(beta0))  # are there any barcodes with no estimate of beta0?


### Compute gamma

par(mfrow=c(1,1))  # how to truncate beta0 distribution.
hist(log(beta0))
abline(v=quantile(log(beta0),0.02), col="red")
abline(v=quantile(log(beta0),0.98), col="red")

beta0_threshed = pmax(pmin(beta0, quantile(beta0, 0.98)), quantile(beta0, 0.02))

gamma=beta1_mean/beta0_threshed
lgamma = log(gamma)
hist(lgamma)
Bhat=matrix(nrow=ncells, ncol=nlineages, data=0)

for(c in 1:ncells){
  if(c%%1000==0) cat(c,"... ")
  if(max(X[c,])>10){
    # high counts, avoid overflow.
    bstarc= which.max(X[c,])
    ldeltabc = lgamma-lgamma[bstarc]
    Xdiff= X[c,] - X[c,bstarc]
    temp = exp(X[c,]*ldeltabc+Xdiff*lgamma[bstarc])
    Bhat[c,] = temp/sum(temp)
    
  } else {
    lgammaX = X[c,]*lgamma
    denom = sum(exp(lgammaX))
    Bhat[c,] = exp(lgammaX)/denom
  }
}
cat("\n")

### get the maximums of Bhat and X.
maxBhat=apply(Bhat,1,max)
argmaxX = apply(X, 1, which.max)
maxX = apply(X,1,max)
argmaxBhat = apply(Bhat, 1, which.max)


### Plots to explore Bhat

hist(maxbhat); grid()


plotXvsBhat<-function(cellnum){
  Xjittered=X[cellnum,]+rnorm(nlineages, 0, 0.03)
  plot(X[cellnum,]+rnorm(nlineages, 0, 0.03), Bhat[cellnum,], cex=2, pch=2, ylim=c(0,1), cex.lab=1.5, cex.main=2,main=paste("Cell", cellnum),xlab="Count", ylab="Posterior barcode probability")
  grid(); abline(h=1, col="red"); abline(h=0, col="red")
  ord = order(X[cellnum,], decreasing=TRUE)
  str=format(lgamma[ord[1:numprint]],digits=2)
  text(Xjittered[ord[1:numprint]], pmin(Bhat[cellnum, ord[1:numprint]]+0.05,0.95), labels=str, col="blue", cex=2)
  
}

par(mfrow=c(1,1))
plot(argmaxX, argmaxBhat)
#disagrees=which(argmaxX!=argmaxBhat)
disagrees=which(argmaxX!=argmaxBhat & maxX>5& maxBhat>0.9)
#disagrees=which(argmaxX!=argmaxBhat & maxX>5 & maxBhat<0.8)
ind=0

ind=ind+1
plotXvsBhat(disagrees[ind])

# ------------------------------------------------------------- #
# Compute embedded distances, look at lineages in embedding.
# ------------------------------------------------------------- #



### Do cells of the same lineage lie closer together in ATAC space?
dim(X)
# [1] 57661  3389
dim(all_data@reductions$lsi@cell.embeddings)
# [1] 57661    50
dim(all_data@reductions$saverpca@cell.embeddings)
# [1] 57661    50
## make sure nrow(X) and nrow(all_data) are the same!
ncells=nrow(all_data@reductions$lsi@cell.embeddings)

idx=sample(ncells, 2000)
#pairs(all_data@reductions$lsi@cell.embeddings[idx,2:10], col=as.numeric(as.factor(all_data$dataset[idx])))
#pairs(all_data@reductions$saverpca@cell.embeddings[idx,1:9], col=as.numeric(as.factor(all_data$dataset[idx])))

compute_embedded_distances<-function(reduc, which_cell, n_components){
  if(class(reduc)!="DimReduc") error("reduc needs to be a SeuratObject::DimReduc")
  ncells = nrow(reduc@cell.embeddings)
  if(which_cell>ncells) error("which_cell needs to be an integer between 1 and nrow(reduc@cell.embeddings)")
  which_U = reduc@cell.embeddings[which_cell, ]
  which_U_mat = matrix(nrow=ncells, ncol=length(which_U), data=which_U, byrow=TRUE)
  dsquare = (which_U_mat-reduc@cell.embeddings)^2  
  return(dsquare[,1:n_components]%*% (reduc@stdev[1:n_components]^2)/sum(reduc@stdev[1:n_components]^2))
}

# test the above function.
# cell-cell distances, assuming reductions have been computed.
# To test the code:
which_cell=1000

which_cell=which_cell+1000  # skip faster so that we get to all the cells.
cell_d = compute_embedded_distances(all_data@reductions$saverpca, which_cell, n_components=40)
all_data$cell_d = pmin(cell_d, quantile(cell_d,probs=0.75))
g1=FeaturePlot(all_data, reduction="saverumap", features= "cell_d")
g1+geom_point(aes(x=all_data@reductions$saverumap@cell.embeddings[which_cell,1],
                  y=all_data@reductions$saverumap@cell.embeddings[which_cell,2]), col="red", size=4, shape=18)

### Now, actually try it on the lineages we computed.

plot_reduc_name="atac.umap"

## Set parameters.
## Identify lineages that have at least min_cells_in_lineage cells in this condition.
bhat_threshold1= 0.99
min_cells_in_lineage=5
condition="day0"

# Names of output files.
pdffilename= paste("LineageAssign_", 100*bhat_threshold1,"_min",min_cells_in_lineage,"cells_",condition,".pdf", sep="")
summaryfilename= paste("LineageAssign_", 100*bhat_threshold1,"_min",min_cells_in_lineage,"cells_",condition,"_summary.png", sep="")

assignedLineage = rep(NA, ncells)  # this holds the assigned lineage to cells.
assignedLineage[which(maxBhat>bhat_threshold1)] = argmaxBhat[which(maxBhat>bhat_threshold1)]
not_assigned = (maxBhat<=bhat_threshold1)



# Get all the cells from the given condition that has > min_cells_in_lineage cells
sel_condition = all_data$dataset==condition
lineageTab = table(assignedLineage[sel_condition])
sel_lineages= as.numeric(names(lineageTab)[which(lineageTab>min_cells_in_lineage)]) # which lineage has more than ?? cells passing threshold.

ratio_lineage_overall = rep(NA, length(sel_lineages))

pdf(file=pdffilename, height=6,width=6)
cat("Total of ", length(sel_lineages), "lineages with greater than", min_cells_in_lineage,"cells.\n")
for(lineage_idx in 1:length(sel_lineages)){
  cat("Computing ", lineage_idx,"-th lineage...\n")
  sel_lineage=sel_lineages[lineage_idx]
  cells_in_lineage = which(sel_condition & assignedLineage==sel_lineage)
  cell_ds = matrix(nrow=ncells, ncol=length(cells_in_lineage), data=0)
  for(i in 1:length(cells_in_lineage)){
    cell_ds[,i]=compute_embedded_distances(all_data@reductions$lsi, cells_in_lineage[i], n_components=40)
  }
  # get lower triangular part of the matrix.
  cell_ds_in_lineage=cell_ds[cells_in_lineage,][lower.tri(cell_ds[cells_in_lineage,])]
  mean_in_lineage = mean(cell_ds_in_lineage)
  mean_overall = mean(cell_ds)
  ratio_lineage_overall[lineage_idx] = mean_in_lineage/mean_overall
  
  par(mfrow=c(1,1))
  sel_cells = which(sel_condition)
  hist(pmin(cell_ds[sel_cells,], quantile(cell_ds[sel_cells,], 0.95)), 50, 
       main="Distances to cells in lineage, same lineage in red", cex.main=1)
  points(cell_ds_in_lineage,rep(0,length(cell_ds_in_lineage)), col="red", pch=17, cex=2)
  maxy = par("usr")[4]
  abline(v=mean_in_lineage, col="red", lty=1)
  text(mean_in_lineage,maxy*0.9, pos=4,"Lineage mean", cex=1.5, col="red") 
  abline(v=mean_overall, col="blue", lty=3)
  text(mean_overall,maxy*0.9, pos=4,"Overall mean", cex=1.5, col="blue") 
  
}

dev.off()

mean_ratio = mean(ratio_lineage_overall)
prop_assigned=sum(!not_assigned)/length(not_assigned)

png(filename=summaryfilename, height=500, width=800)
par(mfrow=c(2,1))
hist(pmin(rsX[not_assigned], 100),  xlab="Total barcode count",
     main=paste("Total barcode count of ",sum(not_assigned)," un-assigned cells", sep=""), breaks=seq(-0.5,100.5,1))
hist(pmin(rsX[!not_assigned], 100),  xlab="Total barcode count",  
     main=paste("Total barcode count of ",sum(!not_assigned)," assigned cells,    assign %=", format(prop_assigned, digits=2), ", ATAC distance ratio=", format(mean_ratio, digits=2),sep=""), breaks=seq(-0.5,100.5,1))
dev.off()



## Plot a specific cell in lineage.

lineage_idx=0

lineage_idx=lineage_idx+1
sel_lineage=sel_lineages[lineage_idx]
cells_in_lineage = which(sel_condition & assignedLineage==sel_lineage)


which_cell = cells_in_lineage[1]
cell_d = compute_embedded_distances(all_data@reductions$lsi, which_cell, n_components=40)

par(mfrow=c(1,1))
sel_cells = which(sel_condition)
hist(pmin(cell_d[sel_cells], quantile(cell_d[sel_cells], 0.95)), 50, main="Distances to selected cell, same lineage in red")
points(cell_d[cells_in_lineage],rep(0, length(cells_in_lineage)), col="red", pch=17, cex=2)

all_data$cell_d = pmin(cell_d, quantile(cell_d,probs=0.6))  # threshold distances at 60%.

g1=FeaturePlot(all_data, reduction=plot_reduc_name, features= "cell_d")
g2=g1+
  geom_point(data=data.frame(x=all_data@reductions[[plot_reduc_name]]@cell.embeddings[setdiff(cells_in_lineage,which_cell),1],
                             y=all_data@reductions[[plot_reduc_name]]@cell.embeddings[setdiff(cells_in_lineage,which_cell),2]), aes(x=x,y=y), 
             col="chartreuse4", size=5, shape=17)+
  geom_point(aes(x=all_data@reductions[[plot_reduc_name]]@cell.embeddings[which_cell,1],
                 y=all_data@reductions[[plot_reduc_name]]@cell.embeddings[which_cell,2]), col="red", size=5, shape=18)+
  ggtitle("Distances to selected cell (red), same lineage in green")

g3=DimPlot(all_data, group.by="dataset", reduction=plot_reduc_name)
g3|g2


