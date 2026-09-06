# Fit lineage imputation at a single lambda

Fit lineage imputation at a single lambda

## Usage

``` r
lineage_imputation(
  cell_features,
  cell_lineage,
  coefficient_initial_list,
  lineage_future_count,
  lambda = 0,
  random_initializations = 10,
  upper_randomness = 5,
  maxit = NA,
  verbose = 1
)
```

## Arguments

- cell_features:

  A numeric matrix where each row represents a cell, and each column
  represents a feature (for instance, the fastTopics scores). Let `n`
  denote the number of cells (rows). Row names (cell IDs) and column
  names (feature names) are required. Do **not** supply an intercept
  column, and more generally no constant column: the intercept is added
  internally (by `.lineage_cleanup()` when fitting, and by
  [`cyfer_finalize()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md)
  when scoring). A constant column you supply yourself is therefore
  duplicated, making the design collinear — and a column already named
  `Intercept` produces two columns of that name.
  [`cyfer_finalize()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md)
  errors on any constant column for this reason;
  [`cyfer()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md)
  tolerates one, so the error may not appear until the finalize step. It
  is also conventional to [`scale()`](https://rdrr.io/r/base/scale.html)
  the features before fitting, both to keep
  [`exp()`](https://rdrr.io/r/base/Log.html) away from overflow and to
  make the ridge penalty comparable across features.

- cell_lineage:

  A character or factor vector of length `n` where element `i` of
  `cell_lineage` denotes which lineage cell `i` belongs to. Factors are
  coerced to character internally, so unused factor levels are harmless.

- coefficient_initial_list:

  A numeric vector or list of numeric vectors of starting coefficients
  (names should match feature names; `Intercept` added if missing).

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future time point for each lineage.

- lambda:

  Ridge penalty weight on non-intercept coefficients.

- random_initializations:

  Number of additional random starts.

- upper_randomness:

  Upper cap for random initial coefficients.

- maxit:

  Iteration cap passed to
  [`optim()`](https://rdrr.io/r/stats/optim.html).

- verbose:

  A numeric, where numbers larger than 1 successively request more
  information to be printed out as the algorithm proceeds.

## Value

An object of class `"lineage_imputation"` with `fit` and `res_list`.
`fit` is the element of `res_list` with the smallest `objective_val`;
`res_list` holds one entry per initialization (supplied first, then
random).

`fit$coefficient_vec` is on the **natural-log** scale, so
`exp(cbind(Intercept = 1, cell_features) %*% coefficient_vec)` is the
expected number of progeny per cell. Note that
[`cyfer_finalize()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md)
returns its `cell_imputed_score` on the **log10** scale instead — see
the "Scales" section of
[`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md).

## Examples

``` r
data(priming_simulation)
# start from the null model: zero on every feature (the intercept is added)
coefficient_initial <- rep(0, ncol(priming_simulation$cell_features))
names(coefficient_initial) <- colnames(priming_simulation$cell_features)
set.seed(10)
fit <- lineage_imputation(cell_features = priming_simulation$cell_features,
                          cell_lineage = priming_simulation$cell_lineage,
                          coefficient_initial_list = coefficient_initial,
                          lineage_future_count = priming_simulation$lineage_future_count,
                          lambda = 1,
                          random_initializations = 2,
                          verbose = 0)
fit$fit$objective_val
#> [1] -430.46
head(fit$fit$coefficient_vec)
#>        Intercept fastTopicCOCL2_1 fastTopicCOCL2_2 fastTopicCOCL2_3 
#>      -0.46390210       0.07407481       0.03575851       0.06745309 
#> fastTopicCOCL2_4 fastTopicCOCL2_5 
#>      -0.03704310       0.03258583 
# one entry per initialization: the supplied start, then the random ones
length(fit$res_list)
#> [1] 3
```
