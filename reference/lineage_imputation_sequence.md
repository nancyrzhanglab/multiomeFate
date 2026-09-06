# Fit lineage imputation along a lambda path

This function calls
[`lineage_imputation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation.md)
for a sequence of lambdas.

## Usage

``` r
lineage_imputation_sequence(
  cell_features,
  cell_lineage,
  lineage_future_count,
  lambda_initial = NA,
  lambda_max = 101,
  lambda_min = 0.01,
  lambda_sequence_length = 50,
  multipler = 10000,
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

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future time point for each lineage.

- lambda_initial:

  The initial value of lambda to perform cross-validation on.

- lambda_min, lambda_max:

  Floor and cap applied to the internal data-driven `lambda_initial`
  heuristic. Only used when `lambda_initial` is `NA`.

- lambda_sequence_length:

  The number of lambdas to perform cross-validation on. The search
  starts with `lambda_initial` and then decays exponentially to 0.

- multipler:

  Scaling factor for the internal `lambda_initial` heuristic.

- verbose:

  A numeric, where numbers larger than 1 successively request more
  information to be printed out as the algorithm proceeds.

## Value

A list with `fit_list` (solution estimated by
[`lineage_imputation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation.md)
per lambda) and `lambda_sequence`. `lambda_sequence` starts at
`lambda_initial` and decays exponentially to `0`, so it is strictly
*decreasing*; `fit_list[[i]]` is the fit at `lambda_sequence[i]`,
warm-started from `fit_list[[i-1]]`. Each `coefficient_vec` is on the
**natural-log** scale — see the "Scales" section of
[`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md).

## Examples

``` r
# \donttest{
data(priming_simulation)
set.seed(10)
path <- lineage_imputation_sequence(
  cell_features = priming_simulation$cell_features,
  cell_lineage = priming_simulation$cell_lineage,
  lineage_future_count = priming_simulation$lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 4,
  verbose = 0)
path$lambda_sequence                       # decreasing, ends at 0
#> [1] 3.0000000 1.5198421 0.5874011 0.0000000
sapply(path$fit_list, function(fit){fit$objective_val})
#> [1] -430.2529 -430.4045 -430.5054 -430.5767
# }
```
