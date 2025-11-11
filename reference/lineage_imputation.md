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

- coefficient_initial_list:

  A numeric vector or list of numeric vectors of starting coefficients
  (names should match feature names; `Intercept` added if missing).

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future timepoint for each lineage.

- lambda:

  Ridge penalty weight on non-intercept coefficients.

- random_initializations:

  Number of additional random starts.

- upper_randomness:

  Upper cap for random initial coefficients.

- verbose:

  A numeric, where numbers larger than 1 successively request more
  information to be printed out as the algorithm proceeds.

## Value

An object of class `"lineage_imputation"` with `fit` and `res_list`.
