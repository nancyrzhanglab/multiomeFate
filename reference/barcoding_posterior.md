# Posterior probability that each cell carries each lineage barcode

Step 3 of the CloneClean barcode pipeline
([`barcode_clustering()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_clustering.md)
-\>
[`barcode_combine()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcode_combine.md)
-\> `barcoding_posterior()` -\>
[`barcoding_assignment()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcoding_assignment.md)).
Given the raw barcode-by-cell count matrix, estimates for each barcode a
"signal" rate and a "background" rate, and converts their ratio into a
posterior over barcodes for every cell.

## Usage

``` r
barcoding_posterior(
  lin_mat,
  bool_force_rebase = FALSE,
  tol = 1e-08,
  verbose = 0
)
```

## Arguments

- lin_mat:

  A barcode-by-cell count matrix, either a `dgCMatrix` or a base matrix:
  rows are barcodes (row names required — they become the lineage
  names), columns are cells.

- bool_force_rebase:

  Whether to always use the max-shifted (overflow-safe) form of the
  normalization. Default `FALSE`, in which case the shift is applied
  only when the largest log-scale term exceeds 10. The two branches are
  mathematically identical; this only changes numerical conditioning.

- tol:

  Threshold below which a background rate is treated as zero when
  computing the winsorizing quantiles. Default `1e-8`.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

A list with:

- `beta0`:

  numeric vector named by barcode, the background rate per barcode.

- `beta1`:

  numeric vector named by barcode, the signal rate per barcode. `NA` for
  any barcode that no cell maximized at.

- `beta1_mean`:

  single numeric, `mean(beta1, na.rm = TRUE)`.

- `posterior_mat`:

  numeric barcode-by-cell matrix with the same dimnames as `lin_mat`;
  each *column* is a probability vector over barcodes summing to 1. This
  is what
  [`barcoding_assignment()`](https://nancyrzhanglab.github.io/multiomeFate/reference/barcoding_assignment.md)
  consumes.

- `gamma`:

  numeric vector of length `nrow(lin_mat)`, named by barcode: the
  enrichment `beta1_mean / beta0`, after winsorizing `beta0` to its 2nd
  and 98th percentiles (among barcodes with `beta0 > tol`) to keep the
  ratio from exploding on a barcode with near-zero background.

- `lineage_num_winner`:

  integer vector of length `nrow(lin_mat)`, how many cells maximized at
  each barcode. A diagnostic: zeros mark the barcodes whose `beta1` is
  `NA`.

## Details

The model behind it: a cell truly carrying barcode `b` produces reads of
`b` at library-size-normalized rate `beta1[b]`; a cell not carrying it
still produces reads at the lower ambient rate `beta0[b]`, from
free-floating barcode contamination. Neither label is observed, so the
estimate is bootstrapped from the argmax: a cell whose *largest* count
is in barcode `b` is provisionally treated as carrying `b`, giving
`beta1[b]`, and every other cell contributes to `beta0[b]`. The
per-barcode enrichment is then `gamma[b] = mean(beta1) / beta0[b]`, and
a cell's posterior is the multinomial-style normalization of
`gamma^count` across barcodes.

Note that the numerator of `gamma` is the *global mean* of `beta1`, not
the barcode's own `beta1[b]`: the signal rate is assumed shared across
barcodes, and all the barcode-specific variation is carried by the
background rate. This is deliberate — `beta1[b]` is estimated from
however many cells happened to maximize at `b`, which for a rare barcode
can be one cell or none (`NA`).

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
dim(res$posterior_mat)                     # barcodes x cells, columns sum to 1
#> [1]  10 300
called <- rownames(res$posterior_mat)[apply(res$posterior_mat, 2, which.max)]
mean(called == paste0("bc", truth_idx))
#> [1] 1
```
