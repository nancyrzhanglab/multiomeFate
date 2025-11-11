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
  lambda_min = 101,
  lambda_sequence_length = 50,
  multipler = 10000,
  verbose = 1
)
```

## Arguments

- cell_features:

  A numeric matrix where each row represents a cell, and each column
  represents a feature (for instance, the fastTopics scores). Let `n`
  denote the number of cells (rows). Please ensure there are row names
  for `cell_features` (denoting the cell IDs), and column names
  (denoting the feature names).

- cell_lineage:

  A character or factor vector of length `n` where element `i` of
  `cell_lineage` denotes which lineage cell `i` belongs to.

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future timepoint for each lineage.

- lambda_initial:

  The initial value of lambda to perform cross-validation on.

- lambda_max, lambda_min:

  Bounds used when computing an internal `lambda_initial`.

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
per lambda) and `lambda_sequence`.
