# Lineage imputation with cross validation

Runs K-fold CV over a decreasing sequence of `lambda` values produced by
[`lineage_imputation_sequence()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation_sequence.md),
selecting `lambda` by held-out objective.

## Usage

``` r
lineage_cv(
  cell_features,
  cell_lineage,
  future_timepoint,
  lineage_future_count,
  lambda_initial,
  lambda_sequence_length,
  tab_mat,
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
  denote the number of cells (rows). Please ensure there are row names
  for `cell_features` (denoting the cell IDs), and column names
  (denoting the feature names).

- cell_lineage:

  A character or factor vector of length `n` where element `i` of
  `cell_lineage` denotes which lineage cell `i` belongs to.

- future_timepoint:

  A character that is one of the column names of `tab_mat` to denote
  which column denotes which timepoint is the "future" time point .

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future timepoint for each lineage.

- lambda_initial:

  The initial value of lambda to perform cross-validation on.

- lambda_sequence_length:

  The number of lambdas to perform cross-validation on. The search
  starts with `lambda_initial` and then decays exponentially to 0.

- tab_mat:

  A matrix where the rows are named and are of each lineage (which
  appeared in `cell_lineage`). There are two columns, one for the number
  of cells in the current timepoint (corresponding to `cell_lineage`).
  The other is the number of cells in the future timepoint
  (corresponding to `lineage_future_count`).

- num_folds:

  Number of folds to do cross-validation on. Default is `10`.

- savefile_tmp:

  Filepath to save files to. Default is `NULL` (no temporary save
  files).

- seed_number:

  Seed value for reproducibility reasons. Default is `10`.

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
