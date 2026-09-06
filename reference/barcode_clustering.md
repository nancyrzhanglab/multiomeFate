# Group barcodes that behave like one lineage

Step 1 of the CloneClean barcode pipeline. Two barcodes whose count
profiles across cells are highly correlated are almost certainly the
same clone read two ways — a sequencing error in the barcode, or two
barcodes integrated into the same founder cell. This finds those groups
so that
[`barcode_combine()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_combine.md)
can merge them; leaving them apart would split one clone into several
small ones and understate its expansion.

## Usage

``` r
barcode_clustering(
  lin_mat,
  cell_lower_limit = 100,
  cor_threshold = 0.55,
  warn_merging = TRUE,
  verbose = 0
)
```

## Arguments

- lin_mat:

  A `dgCMatrix` barcode-by-cell count matrix: rows are barcodes (row
  names required), columns are cells.

- cell_lower_limit:

  Minimum number of cells with a non-zero count for a barcode to be
  eligible for merging. Default `100`.

- cor_threshold:

  Pearson correlation at or above which a pair of barcodes is considered
  the same lineage. Default `0.55`.

- warn_merging:

  Whether to [`warning()`](https://rdrr.io/r/base/warning.html) each
  time two existing clusters are joined through a shared barcode.
  Default `TRUE`. On a real dataset this can fire many times; set
  `FALSE` to silence it once the behaviour has been inspected.

- verbose:

  A numeric; larger values print more. Levels above 2 report cluster
  counts and size quantiles as they accumulate. Default `0`.

## Value

A list with:

- `arr_idx`:

  2-column matrix of the above-threshold barcode index pairs, as
  returned by `which(..., arr.ind = TRUE)` on the upper triangle.
  Indices are into the *filtered* correlation matrix, not into
  `lin_mat`.

- `lineage_clusters`:

  named list (`"c1"`, `"c2"`, ...) of character vectors, each holding
  the barcode names in one cluster. This is what
  [`barcode_combine()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_combine.md)
  takes as its `lineage_clusters`. Barcodes not in any cluster do not
  appear.

- `uniq_lineage`:

  data frame with columns `Lineage` (barcode name) and `Cluster`
  (cluster label), one row per barcode involved in any merge.

- `minimum_correlation`:

  numeric vector, one entry per cluster: the smallest pairwise
  correlation within that cluster. For a cluster built transitively this
  can fall well below `cor_threshold`, which is the main diagnostic for
  an over-aggressive merge.

When no pair clears `cor_threshold`, all four elements are `NULL`;
[`barcode_combine()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_combine.md)
handles that case by returning `lin_mat` unchanged.

## Details

Barcodes carried by fewer than `cell_lower_limit` cells are excluded
from consideration entirely, since a correlation computed over a handful
of cells is not informative. Excluded barcodes are simply never merged;
they are not dropped from the data.

Clustering is single-linkage by construction and done incrementally:
each above-threshold pair either starts a new cluster, joins an existing
one, or — when its two barcodes already sit in *different* clusters —
forces those two clusters to merge. That last case is what
`warn_merging` reports, and it is worth watching, because transitive
merging can chain barcodes together that were never directly correlated.
Inspect `minimum_correlation` to see how far each cluster was stretched.

## Examples

``` r
# Six barcodes over two disjoint blocks of cells; barcodes 1-3 share one
# per-cell rate profile and 4-6 another, so they cluster into two groups.
set.seed(10)
rate_a <- stats::runif(100, min = 5, max = 50)
rate_b <- stats::runif(100, min = 5, max = 50)
zero_vec <- rep(0, 100)
lin_mat <- rbind(c(stats::rpois(100, rate_a), zero_vec),
                 c(stats::rpois(100, rate_a), zero_vec),
                 c(stats::rpois(100, rate_a), zero_vec),
                 c(zero_vec, stats::rpois(100, rate_b)),
                 c(zero_vec, stats::rpois(100, rate_b)),
                 c(zero_vec, stats::rpois(100, rate_b)))
rownames(lin_mat) <- paste0("bc", 1:6)
colnames(lin_mat) <- paste0("cell", 1:200)
lin_mat <- Matrix::Matrix(lin_mat, sparse = TRUE)
res <- barcode_clustering(lin_mat = lin_mat, cell_lower_limit = 50)
res$lineage_clusters
#> $c1
#> [1] "bc1" "bc2" "bc3"
#> 
#> $c2
#> [1] "bc4" "bc5" "bc6"
#> 
```
