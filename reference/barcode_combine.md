# Collapse each barcode cluster into a single row

Step 2 of the CloneClean barcode pipeline. Sums the count rows of every
barcode in a cluster found by
[`barcode_clustering()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_clustering.md),
so the merged clone is represented by one row carrying the pooled
evidence. The surviving row takes the name of the **first barcode of the
cluster in sorted order**, which is arbitrary but stable; the other
names disappear from the matrix.

## Usage

``` r
barcode_combine(lin_mat, lineage_clusters, verbose = 0)
```

## Arguments

- lin_mat:

  A barcode-by-cell count matrix, either a `dgCMatrix` or a base matrix.
  Row names required.

- lineage_clusters:

  A list of character vectors of barcode names, as returned in the
  `lineage_clusters` element of
  [`barcode_clustering()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_clustering.md).
  Entries that are `NA` (the tombstones
  [`barcode_clustering()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_clustering.md)
  leaves when two clusters merge) are skipped. `NULL`, an empty list, or
  a list of nothing but tombstones all mean there is nothing to do.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

A barcode-by-cell matrix with the same columns as `lin_mat` and fewer
rows — one per unclustered barcode plus one per cluster. Returned
unchanged when `lineage_clusters` is `NULL`.

## Details

Row order is not preserved: untouched barcodes come first, followed by
one row per cluster. Anything downstream must therefore index by row
name rather than position.

## Examples

``` r
set.seed(10)
lin_mat <- matrix(stats::rpois(6 * 50, lambda = 5), nrow = 6, ncol = 50)
rownames(lin_mat) <- paste0("bc", 1:6)
colnames(lin_mat) <- paste0("cell", 1:50)
cluster_list <- list(c("bc1", "bc2"), c("bc4", "bc5", "bc6"))
combined <- barcode_combine(lin_mat = lin_mat, lineage_clusters = cluster_list)
rownames(combined)                          # one row per cluster survives
#> [1] "bc3" "bc1" "bc4"
all(colSums(combined) == colSums(lin_mat))  # counts are conserved
#> [1] TRUE
```
