#' Extract the non-zero entries of one column of a sparse matrix
#'
#' Reads the compressed-column slots (\code{@p}, \code{@i}, \code{@x}) directly
#' rather than subsetting, which avoids materializing the column as a dense
#' vector. Only the entries stored in the sparse representation are returned, so
#' the length of the result is the number of structural non-zeros in that
#' column, not \code{nrow(mat)}.
#'
#' @param mat A \code{dgCMatrix} or \code{lgCMatrix}.
#' @param col_idx A single column index, between 1 and \code{ncol(mat)}.
#' @param bool_value If \code{TRUE}, return the stored \emph{values}; if
#'   \code{FALSE}, return the 1-based \emph{row indices} of those values.
#'
#' @returns A numeric vector, empty (\code{numeric(0)}) when the column stores
#'   no non-zeros. The two modes return vectors of the same length, aligned
#'   element-wise, so calling twice pairs each row index with its value.
#'
#' @noRd
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))

  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}
