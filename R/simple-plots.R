rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("~/nzhanglab/project/Multiome_fate/out/kevin/Writeup6b/Writeup6b_all-data.RData")

all_data

head(all_data@meta.data[,c("dataset", "assigned_lineage", "assigned_posterior")])

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[order(rowSums(tab_mat), decreasing = T),]
head(tab_mat)

quantile(tab_mat[,"day0"])
quantile(tab_mat[which(tab_mat[,"day0"] != 0),"day0"])
nrow(tab_mat)
length(which(tab_mat[,"day0"] != 0))
length(which(tab_mat[,"day0"] >= 4))

cor_mat <- stats::cor(tab_mat[,c("day0", "day10_CIS", 
                                 "day10_COCL2", "day10_DABTRAM")])
round(cor_mat,2)

day0_idx <- intersect(which(all_data$dataset == "day0"),
                      which(!is.na(all_data$assigned_lineage)))
day0_idx <- intersect(day0_idx,
                      which(all_data$assigned_posterior >= 0.5))
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[day0_idx] <- TRUE
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep == TRUE)
Seurat::DefaultAssay(all_data2) <- "Saver"
all_data2$geneActivity <- NULL
all_data2$pca <- NULL
all_data2$umap <- NULL
all_data2$lsi <- NULL
all_data2$atac.umap <- NULL
all_data2$activity.umap <- NULL
all_data2$saverpca <- NULL
all_data2$saverumap <- NULL
all_data2$fasttopic_CIS <- NULL
all_data2$fasttopic_COCL2 <- NULL
all_data2$fasttopic_DABTRAM <- NULL
all_data2$common_tcca <- NULL
all_data2$distinct1_tcca <- NULL
all_data2$distinct2_tcca <- NULL

all_data2
dim(all_data2[["Saver"]]@scale.data)
all_data2[["Saver"]]@scale.data[1:5,1:5]
dim(all_data2[["ATAC"]]@data)
all_data2[["ATAC"]]@data[1:5,1:5]
