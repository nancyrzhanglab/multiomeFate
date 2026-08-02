# CYFER: Cell Fate via Exponential Regression (cross-validation)

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
  names (feature names) are required.

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
