extract_peaks <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-"),
    verbose = 0
){
  if(verbose > 0) print("Looking up gene coordinates")
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  if(verbose > 0) print("Finding gene coordinates")
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  if(verbose > 0) print("Finding gene overlaps")
  all_atac_peak <- object[["ATAC"]]@ranges
  overlap_res <- GenomicRanges::findOverlaps(
    query = all_atac_peak,
    subject = region
  )
  
  if(length(overlap_res) > 0){
    if(verbose > 0) print("Intersecting regions")
    region_gene_peaks <- all_atac_peak[overlap_res@from]
    for(i in 1:length(region_gene_peaks)){
      region_gene_peaks[i] <- GenomicRanges::intersect(x = region_gene_peaks[i],
                                                       y = region)
    }
  }
  
  if(verbose > 0) print("Formatting output")
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
  
  region <- Signac::FindRegion(
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

###########################

# # taken directly from https://github.com/stuart-lab/signac/blob/master/R/utilities.R
# # here only because FindRegion is not exported
# FindRegion <- function(
    #     object,
#     region,
#     sep = c("-", "-"),
#     assay = NULL,
#     extend.upstream = 0,
#     extend.downstream = 0
# ) {
#   if (!is(object = region, class2 = "GRanges")) {
#     # first try to convert to coordinates, if not lookup gene
#     region <- tryCatch(
#       expr = suppressWarnings(
#         expr = Signac::StringToGRanges(regions = region, sep = sep)
#       ),
#       error = function(x) {
#         region <- Signac::LookupGeneCoords(
#           object = object,
#           assay = assay,
#           gene = region
#         )
#         return(region)
#       }
#     )
#     if (is.null(x = region)) {
#       stop("Gene not found")
#     }
#   }
#   region <- suppressWarnings(expr = Signac::Extend(
#     x = region,
#     upstream = extend.upstream,
#     downstream = extend.downstream
#   )
#   )
#   return(region)
# }
