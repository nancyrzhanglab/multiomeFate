extract_peaks <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-")
){
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  all_atac_peak <- object[["ATAC"]]@ranges
  overlap_res <- GenomicRanges::findOverlaps(
    query = all_atac_peak,
    subject = region
  )
  
  if(length(overlap_res) > 0){
    region_gene_peaks <- all_atac_peak[overlap_res@from]
    for(i in 1:length(region_gene_peaks)){
      region_gene_peaks[i] <- intersect(x = region_gene_peaks[i],
                                        y = region)
    }
  }
  
  ranges_obj <- region_gene_peaks@ranges
  tmp <- cbind(ranges_obj@start, ranges_obj@start + ranges_obj@width - 1)
  colnames(tmp) <- c("start", "end")
  tmp
}

extract_cutmatrix <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    cells = NULL, #NULL or names of cells
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-"),
    window = 100
){
  
  cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
  
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  # form cutmatrix
  cutmat <- Signac:::CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- (IRanges::start(x = region)):(IRanges::end(x = region))
  
  cutmat
}
