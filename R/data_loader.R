data_loader <- function(
    which_files = c("rna", "atac", "lineage"),
    folder_path = "~/project/Multiome_fate/out/kevin/Writeup10a/",
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
    file_lineage = "Writeup10a_data_lineage.RData",
    file_metadata = "Writeup10a_data_metadata.RData",
    file_peakvi = c(All = "Writeup10a_data_peakVI_All.RData",
                    CIS = "Writeup10a_data_peakVI_CIS.RData",
                    COCL2 = "Writeup10a_data_peakVI_COCL2.RData",
                    DABTRAM = "Writeup10a_data_peakVI_DABTRAM.RData"),
    file_rna = "Writeup10a_data_rna.RData",
    file_rna_dimred = "Writeup10a_data_rna_dimred.RData",
    file_saver = "Writeup10a_data_saver.RData",
    file_wnn = "Writeup10a_data_wnn.RData",
    verbose = 1
){
  possible_assay_vec <- c("atac", "chromvar", "lineage", "rna", "saver")
  possible_dimred_vec <- c("fasttopics", "peakvi", "rna_dimred", "wnn")
  stopifnot(all(which_files %in% c(possible_assay_vec, possible_dimred_vec)))
  
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
      
      for(kk in 1:length(files_chromvar)){
        dataset <- names(files_chromvar)[kk]
        print(paste0("Loading Chromvar: ", dataset))
        
        load(paste0(folder_path, files_chromvar[kk]))
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
  
  # put in the colors
  dataset_colors <- c(day0 = "gray",
                      day10_CIS = "#FBD08C",
                      day10_COCL2 = "#6DC49C",
                      day10_DABTRAM = "#9D85BE",
                      week5_CIS = "#C96D29",
                      week5_COCL2 = "#0F8241",
                      week5_DABTRAM = "#623594")
  all_data@misc$dataset_colors <- dataset_colors
  
  # remove the empty
  assay_vec <- Seurat::Assays(all_data)
  assay_vec <- setdiff(assay_vec, "Empty")
  if(length(assay_vec) > 0){
    Seurat::DefaultAssay(all_data) <- assay_vec[1]
    all_data[["Empty"]] <- NULL
  }
  
  return(all_data)
}