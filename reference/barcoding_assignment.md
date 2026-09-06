# Assign each cell to a lineage, or to none

Step 4 — the last — of the CloneClean barcode pipeline. Takes the
posterior matrix from
[`barcoding_posterior()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcoding_posterior.md)
and calls a lineage for each cell, but only when the call is
unambiguous.

## Usage

``` r
barcoding_assignment(posterior_mat, difference_val = 0.2, verbose = 0)
```

## Arguments

- posterior_mat:

  A numeric barcode-by-cell matrix whose columns are probability vectors
  over barcodes, i.e. the `posterior_mat` element of
  [`barcoding_posterior()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcoding_posterior.md).
  Row names are the lineage names that will be returned; column names
  are the cell IDs.

- difference_val:

  Minimum gap between the top two posteriors for a cell to be assigned.
  Default `0.2`. Raising it trades cell count for assignment purity.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

A character vector of length `ncol(posterior_mat)`, named by cell ID:
the assigned lineage name, or `NA` where the margin was too small.
Suitable for writing straight into `seurat_object$assigned_lineage`,
which is the form the estimation functions expect as `cell_lineage`.

## Details

The criterion is a *margin*, not a level: the gap between the largest
and second-largest posterior must be at least `difference_val`. A cell
whose top posterior is 0.95 is still left unassigned if the runner-up is
0.90. This is what the note on `data_loader()` refers to — an unassigned
cell is generally not a cell with weak evidence, it is a cell with two
competing barcodes, typically a doublet or a cell carrying two
integrations.

## Examples

``` r
# Ten barcodes read across 300 cells: each cell's true barcode at a high
# rate, every other barcode at a low ambient rate.
set.seed(10)
truth_idx <- rep(1:10, length.out = 300)
lin_mat <- matrix(stats::rpois(10 * 300, lambda = 2), nrow = 10, ncol = 300)
for(i in seq_len(300)){
  lin_mat[truth_idx[i], i] <- stats::rpois(1, lambda = 200)
}
rownames(lin_mat) <- paste0("bc", 1:10)
colnames(lin_mat) <- paste0("cell", seq_len(300))
res <- barcoding_posterior(lin_mat = lin_mat)
assignment <- barcoding_assignment(posterior_mat = res$posterior_mat)
head(assignment)
#> cell1 cell2 cell3 cell4 cell5 cell6 
#> "bc1" "bc2" "bc3" "bc4" "bc5" "bc6" 
# NA marks a cell whose top two posteriors were within `difference_val`
sum(is.na(assignment))
#> [1] 0
```
