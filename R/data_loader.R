#' Data loader designed for the paper 
#' 
#' Even if \code{remove_unassigned_cells} is \code{FALSE}, all cells have a posterior of a lineage above 0.5
#' The reason some cells do not have a lineage is because during the calculations,
#' the posterior of a lineage above 0.5 wasn't large compared to the second-highest-posterior lineage.
#'
#' @param which_files A vector of strings denoting which things you want to load in. 
#' The choices are among: 
#' "atac", "chromvar", "lineage", "rna", "saver", 
#' "fasttopics", "peakvi", "rna_dimred", "wnn".
#' @param folder_path Folder for all the preprocessed objects. Do not change.
#' @param file_atac File name for the ATAC counts. Do not change.
#' @param file_chromvar Vector of file name for the chromVar analysis on each treatment-timepoint (7 of them). Do not change.
#' @param file_empty File name for the base file. Do not change.
#' @param file_fasttopics Vector of file name for the fastTopic analysis on each treatment arms and their UMAPs (3 of them). Do not change.
#' @param file_fatepotential File name for the fate potential fits. Do not change.
#' @param file_lineage File name for the Lineage gRNA counts. Do not change.
#' @param file_peakvi Vector of file name for the peakVI analysis on each treatment arms and their UMAPs (4 of them). Do not change.
#' @param file_rna File name for the RNA counts. Do not change.
#' @param file_rna_dimred File name for the PCA and UMAP on the RNA counts. Do not change.
#' @param file_saver File name for the SAVER analysis on the RNA counts. Do not change.
#' @param file_wnn File name for the WNN analysis. Do not change.
#' @param remove_unassigned_cells boolean, on whether to include the cells that do not have an assigned lineage
#' @param verbose positive integer
#'
#' @return a Seurat object
#' @export
data_loader <- function(
    which_files = c("rna", "atac", "lineage"),
    folder_path = "~/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/",
    file_atac = "Writeup10a_data_atac.RData",
    file_chromvar = c(day0 = "Writeup10a_data_chromVar_day0.RData",
                      day10_CIS = "Writeup10a_data_chromVar_day10_CIS.RData",
                      day10_COCL2 = "Writeup10a_data_chromVar_day10_COCL2.RData",
                      day10_DABTRAM = "Writeup10a_data_chromVar_day10_DABTRAM.RData",
                      week5_CIS = "Writeup10a_data_chromVar_week5_CIS.RData",
                      week5_COCL2 = "Writeup10a_data_chromVar_week5_COCL2.RData",
                      week5_DABTRAM = "Writeup10a_data_chromVar_week5_DABTRAM.RData"),
    file_empty = "Writeup10a_data_empty.RData",
    file_fasttopics = c(CIS = "Writeup10a_data_fasttopic_CIS.RData",
                        COCL2 = "Writeup10a_data_fasttopic_COCL2.RData",
                        DABTRAM = "Writeup10a_data_fasttopic_DABTRAM.RData"),
    file_fatepotential = "Writeup10a_data_fatepotential.RData",
    file_lineage = "Writeup10a_data_lineage.RData",
    file_peakvi = c(All = "Writeup10a_data_peakVI_All.RData",
                    CIS = "Writeup10a_data_peakVI_CIS.RData",
                    COCL2 = "Writeup10a_data_peakVI_COCL2.RData",
                    DABTRAM = "Writeup10a_data_peakVI_DABTRAM.RData"),
    file_rna = "Writeup10a_data_rna.RData",
    file_rna_dimred = "Writeup10a_data_rna_dimred.RData",
    file_saver = "Writeup10a_data_saver.RData",
    file_wnn = "Writeup10a_data_wnn.RData",
    remove_unassigned_cells = TRUE,
    verbose = 1
){
  possible_assay_vec <- c("atac", "chromvar", "lineage", "rna", "saver")
  possible_dimred_vec <- c("fasttopics", "peakvi", "rna_dimred", "wnn")
  possible_misc <- c("fatepotential")
  stopifnot(all(which_files %in% c(possible_assay_vec, possible_dimred_vec, possible_misc)))
  
  # to appease R's check
  all_data_atac <- NULL
  all_data_rna <- NULL
  all_data_saver <- NULL
  all_data_saver_pca <- NULL
  all_data_saver_umap <- NULL
  all_data_lineage <- NULL
  all_data_pca <- NULL
  all_data_umap <- NULL
  all_data_wnn <- NULL
  all_data_fatepotential <- NULL
  keep <- NULL
  
  # start by loading the empty 
  load(paste0(folder_path, file_empty))
  
  # put in all the assays
  possible_assay_vec <- intersect(possible_assay_vec, which_files)
  for(assay in possible_assay_vec){
    if(assay == "atac"){
      if(print("Loading ATAC"))
        load(paste0(folder_path, file_atac))
      all_data[["ATAC"]] <- all_data_atac
    }
    
    if(assay == "rna"){
      print("Loading RNA")
      load(paste0(folder_path, file_rna))
      all_data[["RNA"]] <- all_data_rna
    }
    
    if(assay == "saver"){
      print("Loading Saver")
      load(paste0(folder_path, file_saver))
      all_data[["Saver"]] <- all_data_saver
      all_data[["Saver.pca"]] <- all_data_saver_pca
      all_data[["Saver.umap"]] <- all_data_saver_umap
    }
    
    if(assay == "lineage"){
      print("Loading Lineage")
      load(paste0(folder_path, file_lineage))
      all_data[["Lineage"]] <- all_data_lineage
    }
    
    if(assay == "chromvar"){
      print("Loading Chromvar")
      stopifnot(length(names(file_chromvar)) > 0)
      
      for(kk in 1:length(file_chromvar)){
        dataset <- names(file_chromvar)[kk]
        print(paste0("Loading Chromvar: ", dataset))
        
        load(paste0(folder_path, file_chromvar[kk]))
        all_data[[paste0("chromVar.", dataset)]] <- eval(parse(text = paste0("all_data_chromVar_", dataset)))
      }
    }
  }
  
  # put in all the dimreds
  possible_dimred_vec <- intersect(possible_dimred_vec, which_files)
  for(dimred in possible_dimred_vec){
    if(dimred == "rna_dimred"){
      print("Loading RNA dimred")
      load(paste0(folder_path, file_rna_dimred))
      all_data[["pca"]] <- all_data_pca
      all_data[["umap"]] <- all_data_umap
    }
    
    if(dimred == "wnn"){
      print("Loading WNN")
      load(paste0(folder_path, file_wnn))
      all_data[["wnn.umap"]] <- all_data_wnn
    }
    
    if(dimred == "peakvi"){
      print("Loading PeakVI")
      stopifnot(length(names(file_peakvi)) > 0)
      
      for(kk in 1:length(file_peakvi)){
        treatment <- names(file_peakvi)[kk]
        print(paste0("Loading PeakVI: ", treatment))
        
        load(paste0(folder_path, file_peakvi[kk]))
        all_data[[paste0("peakVI.", treatment)]] <- eval(parse(text = paste0("all_data_peakVI_", treatment)))
        all_data[[paste0("pVI.", treatment, ".umap")]] <- eval(parse(text = paste0("all_data_pVI_", treatment, "_umap")))
      }
    }
    
    if(dimred == "fasttopics"){
      print("Loading fastTopics")
      stopifnot(length(names(file_fasttopics)) > 0)
      
      for(kk in 1:length(file_fasttopics)){
        treatment <- names(file_fasttopics)[kk]
        print(paste0("Loading fastTopics: ", treatment))
        
        load(paste0(folder_path, file_fasttopics[kk]))
        all_data[[paste0("fasttopic.", treatment)]] <- eval(parse(text = paste0("all_data_fasttopic_", treatment)))
        all_data[[paste0("ft.", treatment, ".umap")]] <- eval(parse(text = paste0("all_data_ft_", treatment, "_umap")))
      }
    }
  }
  
  # put in all the miscs
  possible_misc <- intersect(possible_misc, which_files)
  for(misc in possible_misc){
    if(misc == "fatepotential"){
      print("Loading fate potentials")

      load(paste0(folder_path, file_fatepotential))
      all_data@misc <- all_data_fatepotential
    }
  }
  
  
  # put in the colors
  dataset_colors <- c(day0 = "gray",
                      day10_CIS = "#FBD08C",
                      day10_COCL2 = "#6DC49C",
                      day10_DABTRAM = "#9D85BE",
                      week5_CIS = "#C96D29",
                      week5_COCL2 = "#0F8241",
                      week5_DABTRAM = "#623594")
  all_data@misc$dataset_colors <- dataset_colors
  all_data@misc$fatepotential_colors <- c("red", "bisque", "blue")
  all_data@misc$fatepotential_na_colors <- "gray90"
  all_data@misc$modality_colors <- c(RNA = "#FF7F7F",
                                     chromVar = "#4682B4")
  
  # remove the empty
  assay_vec <- Seurat::Assays(all_data)
  assay_vec <- setdiff(assay_vec, "Empty")
  if(length(assay_vec) > 0){
    Seurat::DefaultAssay(all_data) <- assay_vec[1]
    all_data[["Empty"]] <- NULL
  }
  
  # remove cells with no lineage
  if(remove_unassigned_cells) {
    print("Removing cells with no assigned lineage")
    all_data$keep <- !is.na(all_data$assigned_lineage)
    if(any(!all_data$keep)){
      print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
      all_data <- subset(all_data, keep == TRUE)
    }
  }
  
  return(all_data)
}