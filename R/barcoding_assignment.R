#' Posterior probability that each cell carries each lineage barcode
#'
#' Step 3 of the CloneClean barcode pipeline (\code{barcode_clustering()} ->
#' \code{barcode_combine()} -> \code{barcoding_posterior()} ->
#' \code{barcoding_assignment()}). Given the raw barcode-by-cell count matrix,
#' estimates for each barcode a "signal" rate and a "background" rate, and
#' converts their ratio into a posterior over barcodes for every cell.
#'
#' The model behind it: a cell truly carrying barcode \code{b} produces reads of
#' \code{b} at library-size-normalized rate \code{beta1[b]}; a cell not carrying
#' it still produces reads at the lower ambient rate \code{beta0[b]}, from
#' free-floating barcode contamination. Neither label is observed, so the
#' estimate is bootstrapped from the argmax: a cell whose \emph{largest} count is
#' in barcode \code{b} is provisionally treated as carrying \code{b}, giving
#' \code{beta1[b]}, and every other cell contributes to \code{beta0[b]}. The
#' per-barcode enrichment is then \code{gamma[b] = mean(beta1) / beta0[b]}, and
#' a cell's posterior is the multinomial-style normalization of
#' \code{gamma^count} across barcodes.
#'
#' Note that the numerator of \code{gamma} is the \emph{global mean} of
#' \code{beta1}, not the barcode's own \code{beta1[b]}: the signal rate is
#' assumed shared across barcodes, and all the barcode-specific variation is
#' carried by the background rate. This is deliberate --- \code{beta1[b]} is
#' estimated from however many cells happened to maximize at \code{b}, which for
#' a rare barcode can be one cell or none (\code{NA}).
#'
#' @param lin_mat A barcode-by-cell count matrix, either a \code{dgCMatrix} or a
#'   base matrix: rows are barcodes (row names required --- they become the
#'   lineage names), columns are cells.
#' @param bool_force_rebase Whether to always use the max-shifted
#'   (overflow-safe) form of the normalization. Default \code{FALSE}, in which
#'   case the shift is applied only when the largest log-scale term exceeds 10.
#'   The two branches are mathematically identical; this only changes numerical
#'   conditioning.
#' @param tol Threshold below which a background rate is treated as zero when
#'   computing the winsorizing quantiles. Default \code{1e-8}.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{beta0}}{numeric vector named by barcode, the background rate
#'       per barcode.}
#'     \item{\code{beta1}}{numeric vector named by barcode, the signal rate per
#'       barcode. \code{NA} for any barcode that no cell maximized at.}
#'     \item{\code{beta1_mean}}{single numeric, \code{mean(beta1, na.rm = TRUE)}.}
#'     \item{\code{posterior_mat}}{numeric barcode-by-cell matrix with the same
#'       dimnames as \code{lin_mat}; each \emph{column} is a probability vector
#'       over barcodes summing to 1. This is what
#'       \code{barcoding_assignment()} consumes.}
#'     \item{\code{gamma}}{numeric vector of length \code{nrow(lin_mat)}, named
#'       by barcode: the enrichment \code{beta1_mean / beta0}, after
#'       winsorizing \code{beta0} to its 2nd and 98th percentiles (among
#'       barcodes with \code{beta0 > tol}) to keep the ratio from exploding on a
#'       barcode with near-zero background.}
#'     \item{\code{lineage_num_winner}}{integer vector of length
#'       \code{nrow(lin_mat)}, how many cells maximized at each barcode. A
#'       diagnostic: zeros mark the barcodes whose \code{beta1} is \code{NA}.}
#'   }
#'
#' @export
# cells as columns, lineage as rows
barcoding_posterior <- function(lin_mat, # barcode-by-cell matrix
                                bool_force_rebase = FALSE,
                                tol = 1e-8,
                                verbose = 0){
  stopifnot(inherits(lin_mat, c("matrix", "dgCMatrix")),
            length(rownames(lin_mat)) > 0)
  library_size <- Matrix::colSums(lin_mat)
  
  n <- ncol(lin_mat) 
  nlineages <- nrow(lin_mat)
  
  # to initialize the estimator, find the maximum count for each cell
  if(verbose > 0) print("Starting barcoding assignment")
  cell_max_barcode_count <- apply(lin_mat, 2, max)
  lineage_maximizing <- lapply(seq_len(nlineages), function(b){
    which(lin_mat[b,] == cell_max_barcode_count & lin_mat[b,] > 0)
  })
  lin_num_winner <- sapply(lineage_maximizing, length)
  names(lin_num_winner) <- rownames(lin_mat)

  beta0 <- rep(NA, nlineages)
  names(beta0) <- rownames(lin_mat)
  beta1 <- rep(NA, nlineages)
  names(beta1) <- rownames(lin_mat)

  for(b in seq_len(nlineages)){
    if(verbose == 1 && b %% 100==0) print(paste0(b," out of ", nlineages," done"))
    if(verbose == 2) print(paste0(b," out of ", nlineages," done"))
    
    # for barcode b, first find all the cells that have its maximum count in barcode b
    won_cells <- lineage_maximizing[[b]]
    other_cells <- setdiff(seq_len(n), lineage_maximizing[[b]])

    # beta1 estimate among maximizing cells
    if(length(won_cells) > 0) beta1[b] <- mean(lin_mat[b,won_cells] / (library_size[won_cells]+1))
    # beta0 estimate among all non-maximizing cells
    beta0[b] <- mean(lin_mat[b,other_cells] / (library_size[other_cells]+1))
  }

  # global averages
  beta1_mean <- mean(beta1, na.rm = TRUE)

  # Winsorize beta0 so a barcode with near-zero background does not send
  # gamma to Inf.
  lower_val <- stats::quantile(beta0[beta0 > tol], 0.02)
  upper_val <- stats::quantile(beta0[beta0 > tol], 0.98)
  beta0_thresh <- pmax(pmin(beta0, upper_val), lower_val)
  
  gamma <- beta1_mean/beta0_thresh # vector of length nlineages
  names(gamma) <- rownames(lin_mat)
  
  posterior_mat <- .multinomial_posterior(bool_force_rebase = bool_force_rebase,
                                          gamma = gamma,
                                          lin_mat = lin_mat)
  
  list(beta0 = beta0,
       beta1 = beta1,
       beta1_mean = beta1_mean,
       posterior_mat = posterior_mat,
       gamma = gamma,
       lineage_num_winner = lin_num_winner)
}

#' Group barcodes that behave like one lineage
#'
#' Step 1 of the CloneClean barcode pipeline. Two barcodes whose count profiles
#' across cells are highly correlated are almost certainly the same clone read
#' two ways --- a sequencing error in the barcode, or two barcodes integrated
#' into the same founder cell. This finds those groups so that
#' \code{barcode_combine()} can merge them; leaving them apart would split one
#' clone into several small ones and understate its expansion.
#'
#' Barcodes carried by fewer than \code{cell_lower_limit} cells are excluded
#' from consideration entirely, since a correlation computed over a handful of
#' cells is not informative. Excluded barcodes are simply never merged; they are
#' not dropped from the data.
#'
#' Clustering is single-linkage by construction and done incrementally: each
#' above-threshold pair either starts a new cluster, joins an existing one, or
#' --- when its two barcodes already sit in \emph{different} clusters --- forces
#' those two clusters to merge. That last case is what \code{warn_merging}
#' reports, and it is worth watching, because transitive merging can chain
#' barcodes together that were never directly correlated. Inspect
#' \code{minimum_correlation} to see how far each cluster was stretched.
#'
#' @param lin_mat A \code{dgCMatrix} barcode-by-cell count matrix: rows are
#'   barcodes (row names required), columns are cells.
#' @param cell_lower_limit Minimum number of cells with a non-zero count for a
#'   barcode to be eligible for merging. Default \code{100}.
#' @param cor_threshold Pearson correlation at or above which a pair of barcodes
#'   is considered the same lineage. Default \code{0.55}.
#' @param warn_merging Whether to \code{warning()} each time two existing
#'   clusters are joined through a shared barcode. Default \code{TRUE}. On a
#'   real dataset this can fire many times; set \code{FALSE} to silence it once
#'   the behaviour has been inspected.
#' @param verbose A numeric; larger values print more. Levels above 2 report
#'   cluster counts and size quantiles as they accumulate. Default \code{0}.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{arr_idx}}{2-column matrix of the above-threshold barcode
#'       index pairs, as returned by \code{which(..., arr.ind = TRUE)} on the
#'       upper triangle. Indices are into the \emph{filtered} correlation
#'       matrix, not into \code{lin_mat}.}
#'     \item{\code{lineage_clusters}}{named list (\code{"c1"}, \code{"c2"}, ...)
#'       of character vectors, each holding the barcode names in one cluster.
#'       This is what \code{barcode_combine()} takes as its
#'       \code{lineage_clusters}. Barcodes not in any cluster do not appear.}
#'     \item{\code{uniq_lineage}}{data frame with columns \code{Lineage}
#'       (barcode name) and \code{Cluster} (cluster label), one row per barcode
#'       involved in any merge.}
#'     \item{\code{minimum_correlation}}{numeric vector, one entry per cluster:
#'       the smallest pairwise correlation within that cluster. For a cluster
#'       built transitively this can fall well below \code{cor_threshold}, which
#'       is the main diagnostic for an over-aggressive merge.}
#'   }
#'   When no pair clears \code{cor_threshold}, all four elements are
#'   \code{NULL}; \code{barcode_combine()} handles that case by returning
#'   \code{lin_mat} unchanged.
#'
#' @export
barcode_clustering <- function(lin_mat,
                               cell_lower_limit = 100,
                               cor_threshold = 0.55,
                               warn_merging = TRUE,
                               verbose = 0){
  stopifnot(inherits(lin_mat, "dgCMatrix"))
  stopifnot(length(rownames(lin_mat)) > 0)
  
  lin_mat_t <- Matrix::t(lin_mat)
  nlineages <- nrow(lin_mat)
  num_cells <- sapply(1:nlineages, function(j){
    length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
  })
  
  if(verbose > 0) print("Compute correlation matrix")
  # discard any lineages in question (for the purposes of merging) that are too small
  lin_mat_t <- lin_mat_t[,which(num_cells >= cell_lower_limit)]
  # compute a correlation matrix (# rows/columsn = number of lineages)
  cor_mat <- .custom_correlation(lin_mat_t)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  # determine all the lineages to merge. arr_idx is a 2-column matrix
  arr_idx <- which(cor_mat >= cor_threshold, arr.ind = TRUE)
  if(verbose > 2){
    print(paste0("There are ", nrow(arr_idx), " number of highly correlated lineages to resolve."))
  }
  
  if(nrow(arr_idx) == 0){
    if(verbose > 0) print("Finishing since there is nothing to do")
    
    return(list(arr_idx = NULL,
                lineage_clusters = NULL,
                uniq_lineage = NULL,
                minimum_correlation = NULL))
  }
  
  if(verbose > 0) print("Determining how to merge lineages")
  # tabulate uniq_lineage, which is going to keep track of which lineage is 
  #  assigned to which "cluster"
  lineage_list <- vector("list", length = 0)
  uniq_lineage <- sort(unique(as.numeric(arr_idx)))
  uniq_lineage <- data.frame(Lineage = uniq_lineage, 
                             Cluster = rep(NA, length(uniq_lineage)))
  lineage_name <- rownames(cor_mat)
  
  # determine the clusters of lineages to merge
  for(i in 1:nrow(arr_idx)){
    if(verbose == 1 && nrow(arr_idx) > 10 && i %% floor(nrow(arr_idx)/10) == 0) cat('*')
    if(verbose >= 2) print(paste0("Working on lineage correlation pair ", i , " out of ", nrow(arr_idx)))
    
    # find the rows in the uniq_lineage table on which we're currently working on
    # this represents a pair of lineages
    lineage_idx <- which(uniq_lineage[,"Lineage"] %in% arr_idx[i,])
    
    if(length(lineage_list) == 0) {
      # if we have not yet merged any lineages, then this is straight-forward
      # create a new cluster
      lineage_list[[1]] <- sort(arr_idx[1,])
      names(lineage_list)[[1]] <- "c1"
      uniq_lineage[lineage_idx,"Cluster"] <- "c1"
      
    } else {
      # otherwise...
      
      if(all(is.na(uniq_lineage[lineage_idx,2]))) {
        # if this is a completely new cluster (i.e., all unassigned), also pretty straight-forward
        # create a new cluster
        
        new_list <- list(sort(arr_idx[i,])); names(new_list) <- paste0("c", length(lineage_list)+1)
        lineage_list <- c(lineage_list, new_list)
        uniq_lineage[lineage_idx,"Cluster"] <- paste0("c", length(lineage_list))
        if(verbose > 2 && length(lineage_list) %% 100 == 0) {
          print(paste0("There are currently ", length(lineage_list), " cluster of lineages"))
          
          if(verbose > 3) {
            print("The size of the lineages are: ")
            print(stats::quantile(sapply(lineage_list, length)))
          }
        }
        
      } else {
        # the difficulty is this step. The new lineage in question has a high
        #  correlation with an existing cluster of lineages
        
        # first find all the lineages in this existing cluster
        val <- uniq_lineage[lineage_idx,"Cluster"]
        val <- val[!is.na(val)]
        val <- unique(val)
        if(length(val) > 1) {
          if(warn_merging) {warning("Merging happening")}

          # A pair names two barcodes and each sits in at most one cluster, so
          # a bridging pair can join at most two clusters.
          if(length(val) != 2){
            stop("A correlated pair spans ", length(val), " clusters; it can ",
                 "span at most 2. `lineage_clusters` is inconsistent.")
          }
          val <- sort(val)
          val_pick <- val[1]
          val_dump <- val[2]
          lineage_list[[val_pick]] <- sort(unique(c(lineage_list[[val_pick]], lineage_list[[val_dump]], arr_idx[i,])))
          lineage_list[[val_dump]] <- NA
          
          ## we also need to update all the values in uniq_lineage
          uniq_lineage[lineage_idx,"Cluster"] <- val_pick
          changing_idx <- which(uniq_lineage[,"Cluster"] == val_dump)
          uniq_lineage[changing_idx,"Cluster"] <- val_pick
          
        } else {
          # just add this lineage to the existing cluster
          lineage_list[[val]] <- sort(unique(c(lineage_list[[val]], arr_idx[i,])))
          uniq_lineage[lineage_idx,"Cluster"] <- val
        }
      }
    }
  }
  
  # cleanup
  notna_idx <- which(sapply(lineage_list, function(x){!all(is.na(x))}))
  lineage_list <- lineage_list[notna_idx]
  
  # rename all the indices with their lineage name
  lineage_list_name <- lapply(lineage_list, function(vec){
    lineage_name[vec]
  })
  uniq_lineage[,"Lineage"] <- lineage_name[uniq_lineage[,"Lineage"]]
  
  # compute the minimum correlations
  min_cor_vec <- sapply(lineage_list_name, function(vec){
    min(cor_mat[vec,vec], na.rm = TRUE)
  })
  min_cor_vec <- min_cor_vec[!is.na(min_cor_vec)]
  
  list(arr_idx = arr_idx,
       lineage_clusters = lineage_list_name,
       uniq_lineage = uniq_lineage,
       minimum_correlation = min_cor_vec)
}

#' Collapse each barcode cluster into a single row
#'
#' Step 2 of the CloneClean barcode pipeline. Sums the count rows of every
#' barcode in a cluster found by \code{barcode_clustering()}, so the merged
#' clone is represented by one row carrying the pooled evidence. The surviving
#' row takes the name of the \bold{first barcode of the cluster in sorted
#' order}, which is arbitrary but stable; the other names disappear from the
#' matrix.
#'
#' Row order is not preserved: untouched barcodes come first, followed by one
#' row per cluster. Anything downstream must therefore index by row name rather
#' than position.
#'
#' @param lin_mat A barcode-by-cell count matrix, either a \code{dgCMatrix} or a
#'   base matrix. Row names required.
#' @param lineage_clusters A list of character vectors of barcode names, as
#'   returned in the \code{lineage_clusters} element of
#'   \code{barcode_clustering()}. Entries that are \code{NA} (the tombstones
#'   \code{barcode_clustering()} leaves when two clusters merge) are skipped.
#'   \code{NULL}, an empty list, or a list of nothing but tombstones all mean
#'   there is nothing to do.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns A barcode-by-cell matrix with the same columns as \code{lin_mat} and
#'   fewer rows --- one per unclustered barcode plus one per cluster. Returned
#'   unchanged when \code{lineage_clusters} is \code{NULL}.
#'
#' @export
barcode_combine <- function(lin_mat,
                            lineage_clusters,
                            verbose = 0){
  stopifnot(inherits(lin_mat, c("matrix", "dgCMatrix")),
            length(rownames(lin_mat)) > 0,
            is.null(lineage_clusters) || is.list(lineage_clusters))

  # Drop the NA tombstones before deciding there is nothing to do, so a list of
  # nothing but tombstones is treated the same as NULL.
  if(length(lineage_clusters) > 0){
    keep_idx <- which(sapply(lineage_clusters, function(x){all(!is.na(x))}))
    lineage_clusters <- lineage_clusters[keep_idx]
  }

  if(length(lineage_clusters) == 0){
    if(verbose > 0) print("No lineage clusters, so returning original matrix")
    return(lin_mat)
  }

  nlineages <- nrow(lin_mat)

  if(verbose > 0) print("Starting combination")
  lineage_included_names <- sort(unique(unlist(lineage_clusters)))
  lineage_included_idx <- which(rownames(lin_mat) %in% lineage_included_names)
  lineage_excluded_idx <- setdiff(seq_len(nlineages), lineage_included_idx)

  if(verbose > 0) print("Extracting unaffected lineages")
  # `drop = FALSE` throughout: a single surviving row would otherwise collapse
  # to a vector and rbind() would name it after the variable.
  lin_untouched <- lin_mat[lineage_excluded_idx,,drop = FALSE]

  if(verbose > 0) print("Extracting affected lineages")
  len <- length(lineage_clusters)
  lin_list <- lapply(lineage_clusters, function(vec){
    lin_idx <- which(rownames(lin_mat) %in% vec)
    lin_mat[lin_idx,,drop = FALSE]
  })

  if(verbose > 0) print("Adding lineages to be combined")
  for(i in seq_len(len)){
    if(verbose > 1 && len > 10 && i %% floor(len/10) == 0) cat('*')
    vec <- rownames(lin_list[[i]])
    colname_vec <- colnames(lin_list[[i]])
    lin_list[[i]] <- matrix(Matrix::colSums(lin_list[[i]]), ncol = ncol(lin_list[[i]]), nrow = 1)
    rownames(lin_list[[i]]) <- vec[1]
    colnames(lin_list[[i]]) <- colname_vec
  }

  if(verbose > 0) print("Formatting final matrix")
  rbind(lin_untouched, do.call(rbind, lin_list))
}

#' Assign each cell to a lineage, or to none
#'
#' Step 4 --- the last --- of the CloneClean barcode pipeline. Takes the
#' posterior matrix from \code{barcoding_posterior()} and calls a lineage for
#' each cell, but only when the call is unambiguous.
#'
#' The criterion is a \emph{margin}, not a level: the gap between the largest
#' and second-largest posterior must be at least \code{difference_val}. A cell
#' whose top posterior is 0.95 is still left unassigned if the runner-up is
#' 0.90. This is what the note on \code{data_loader()} refers to --- an
#' unassigned cell is generally not a cell with weak evidence, it is a cell with
#' two competing barcodes, typically a doublet or a cell carrying two
#' integrations.
#'
#' @param posterior_mat A numeric barcode-by-cell matrix whose columns are
#'   probability vectors over barcodes, i.e. the \code{posterior_mat} element of
#'   \code{barcoding_posterior()}. Row names are the lineage names that will be
#'   returned; column names are the cell IDs.
#' @param difference_val Minimum gap between the top two posteriors for a cell
#'   to be assigned. Default \code{0.2}. Raising it trades cell count for
#'   assignment purity.
#' @param verbose A numeric; larger values print more. Default \code{0}.
#'
#' @returns A character vector of length \code{ncol(posterior_mat)}, named by
#'   cell ID: the assigned lineage name, or \code{NA} where the margin was too
#'   small. Suitable for writing straight into
#'   \code{seurat_object$assigned_lineage}, which is the form the estimation
#'   functions expect as \code{cell_lineage}.
#'
#' @export
barcoding_assignment <- function(posterior_mat,
                                 difference_val = 0.2,
                                 verbose = 0){
  n <- ncol(posterior_mat)
  
  lineage_names <- rownames(posterior_mat)
  lineage_idx <- apply(posterior_mat, 2, which.max)
  if(verbose > 0) print("Computing difference between maximizing and second-maximizing")
  difference_vec <- apply(posterior_mat, 2, function(x){
    abs(diff(sort(x, decreasing = T)[1:2]))
  })
  
  assignment_vec <- sapply(1:n, function(i){
    if(verbose > 0 && n > 10 && i %% floor(n/10) == 0) cat('*')
    if(difference_vec[i] >= difference_val){
      return(lineage_names[lineage_idx[i]]) 
    } else {
      return(NA)
    }
  })
  
  names(assignment_vec) <- colnames(posterior_mat)
  assignment_vec
}

###################################

#' Posterior over barcodes for every cell
#'
#' Loops \code{.multinomial_posterior_vector()} over the columns (cells) of
#' \code{lin_mat}. Split out from \code{barcoding_posterior()} so the per-cell
#' arithmetic can be tested on its own.
#'
#' @param bool_force_rebase Passed through: force the overflow-safe branch.
#' @param gamma Numeric vector of per-barcode enrichments, length
#'   \code{nrow(lin_mat)} and in the same order as its rows. Logged once here
#'   rather than per cell.
#' @param lin_mat Dense barcode-by-cell count matrix.
#' @param verbose A numeric; larger values print a progress tick. Default
#'   \code{0}.
#'
#' @returns A numeric matrix with \code{dimnames(lin_mat)}, each column a
#'   probability vector over barcodes summing to 1.
#'
#' @noRd
.multinomial_posterior <- function(bool_force_rebase,
                                   gamma,
                                   lin_mat,
                                   verbose = 0){
  n <- ncol(lin_mat)
  nlineages <- nrow(lin_mat)
  
  lgamma <- log(gamma)
  Bhat <- matrix(0, ncol = n, nrow = nlineages)
  colnames(Bhat) <- colnames(lin_mat)
  rownames(Bhat) <- rownames(lin_mat)
  
  # we do calculation on the log-scale and then exponentiate
  for(i in 1:n){
    if(verbose > 0 && i %% floor(n/10) == 0) cat('*')
    Bhat[,i] <- .multinomial_posterior_vector(
      bool_force_rebase = bool_force_rebase,
      lgamma = lgamma,
      lin_count = lin_mat[,i]
    )
  }
  
  Bhat
}

#' Posterior over barcodes for one cell
#'
#' Normalizes \code{gamma[b]^count[b]} across barcodes \code{b}, computed on the
#' log scale as \code{count * log(gamma)}. A barcode with a high enrichment and
#' a high count therefore dominates multiplicatively, which is what makes the
#' posterior concentrate as sequencing depth grows.
#'
#' The two branches are algebraically identical; the first subtracts the maximum
#' before exponentiating so that a cell with large counts does not overflow to
#' \code{Inf / Inf = NaN}. The plain branch is kept for the common
#' low-count case.
#'
#' @param bool_force_rebase If \code{TRUE}, always take the max-shifted branch.
#'   Otherwise it is taken only when \code{max(count * log(gamma)) > 10}.
#' @param lgamma Numeric vector, \code{log(gamma)}, one entry per barcode.
#' @param lin_count Numeric vector of that cell's counts, same length and order
#'   as \code{lgamma}.
#'
#' @returns A numeric probability vector over barcodes, summing to 1.
#'
#' @noRd
.multinomial_posterior_vector <- function(bool_force_rebase,
                                          lgamma,
                                          lin_count){
  lgammaX <- lin_count*lgamma
  
  if(bool_force_rebase || max(lgammaX) > 10){
    # high counts, avoid overflow
    max_val <- max(lgammaX)
    tmp <- exp(lgammaX-max_val)
    vec <- tmp/sum(tmp)
    
  } else {
    denom <- sum(exp(lgammaX))
    vec <- exp(lgammaX)/denom
  }
  
  vec
}

#' Column-wise Pearson correlation of a sparse matrix
#'
#' \code{stats::cor()} densifies its input, which is not viable for a
#' cells-by-barcodes matrix. This computes the covariance from
#' \code{crossprod(x) - n * tcrossprod(colMeans(x))}, which keeps the sparse
#' product sparse and only densifies the (much smaller) barcode-by-barcode
#' result.
#'
#' from https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
#'
#' @param x A sparse matrix, \bold{cells as rows and barcodes as columns} ---
#'   note this is the transpose of the \code{lin_mat} orientation used
#'   elsewhere in this file, which is why \code{barcode_clustering()} passes
#'   \code{Matrix::t(lin_mat)}.
#'
#' @returns A dense numeric correlation matrix, \code{ncol(x)} by
#'   \code{ncol(x)}, with \code{colnames(x)} as its dimnames. A column with zero
#'   variance yields \code{NaN} in its row and column.
#'
#' @noRd
.custom_correlation <- function(x){
  n <- nrow(x)
  cMeans <- Matrix::colMeans(x)
  covmat <- (as.matrix(Matrix::crossprod(x)) - n*Matrix::tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(Matrix::diag(covmat)) 
  cormat <- covmat/Matrix::tcrossprod(sdvec)
  cormat
}