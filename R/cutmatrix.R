.compute_cutmatrix <- function(
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
