# CYFER: Clonal Fate Estimation by Exponential Regression (cross-validation)

Runs K-fold CV over a decreasing sequence of `lambda` values produced by
[`lineage_imputation_sequence()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation_sequence.md),
selecting `lambda` by held-out objective. This implements the CYFER
method described in Chen, Lin et al.

## Usage

``` r
cyfer(
  cell_features,
  cell_lineage,
  lineage_future_count,
  lambda_initial,
  lambda_sequence_length,
  num_folds = 10,
  savefile_tmp = NULL,
  seed_number = 10,
  verbose = 0
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
  errors on any constant column for this reason; `cyfer()` tolerates
  one, so the error may not appear until the finalize step. It is also
  conventional to [`scale()`](https://rdrr.io/r/base/scale.html) the
  features before fitting, both to keep
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

- lambda_sequence_length:

  The number of lambdas to perform cross-validation on. The search
  starts with `lambda_initial` and then decays exponentially to 0.

- num_folds:

  Number of folds to do cross-validation on. Default is `10`. Must be at
  least 2 and at most the number of distinct lineages.

- savefile_tmp:

  Filepath to save files to. Default is `NULL` (no temporary save
  files).

- seed_number:

  Seed value for reproducibility reasons. Default is `10`. Governs both
  the fold assignment and the optimizer's random restarts.

- verbose:

  A numeric, where numbers larger than 1 successively request more
  information to be printed out as the algorithm proceeds.

## Value

A list with the number of elements corresponding to `num_folds`. Each
element of this list contains: `test_loglik` (the negative
log-likelihood on the held-out lineages), `train_loglik` (the negative
log-likelihood on the trained lineages), and `train_fit` (the actual
fit, after using the
[`lineage_imputation_sequence()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation_sequence.md)).

`test_loglik` and `train_loglik` are the *unpenalized* objective
evaluated at each lambda along `train_fit$lambda_sequence` (lower is
better), not literal log-likelihoods. Pass the result to
[`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md)
to select lambda and score the cells; note that the per-cell scores it
returns are on the log10 scale while the coefficients here are on the
natural-log scale.

Lineages named in `lineage_future_count` that have no cells in
`cell_lineage` are dropped before the folds are built, so they neither
occupy a fold nor count towards `num_folds`.

## Examples

``` r
# \donttest{
data(priming_simulation)
cv <- cyfer(cell_features = priming_simulation$cell_features,
            cell_lineage = priming_simulation$cell_lineage,
            lineage_future_count = priming_simulation$lineage_future_count,
            lambda_initial = 3,
            lambda_sequence_length = 5,
            num_folds = 5,
            verbose = 0)
class(cv)
#> [1] "cyfer"
names(cv)
#> [1] "fold:1" "fold:2" "fold:3" "fold:4" "fold:5"
# held-out objective along the lambda path, one row per fold
sapply(cv, function(fold){fold$test_loglik})
#>         fold:1    fold:2    fold:3    fold:4    fold:5
#> [1,] -405.4184 -440.5515 -419.8813 -482.7639 -403.8150
#> [2,] -405.4332 -440.5615 -419.8933 -482.8811 -403.8197
#> [3,] -405.4417 -440.5670 -419.9029 -482.9528 -403.8235
#> [4,] -405.4464 -440.5706 -419.9120 -483.0174 -403.8278
#> [5,] -405.4477 -440.5703 -419.9228 -483.0921 -403.8339
# }
```
