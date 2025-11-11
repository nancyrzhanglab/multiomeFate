# Finalize lineage imputation after cross validation

This function is used after running
[`multiomeFate::lineage_cv()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_cv.md).
Chooses `lambda` by minimizing the median held-out objective across
folds, then refits once on all cells at the chosen `lambda`.

## Usage

``` r
lineage_cv_finalize(cell_features, cell_lineage, fit_res, lineage_future_count)
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

- fit_res:

  This is the output of
  [`lineage_cv()`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_cv.md).

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future timepoint for each lineage.

## Value

A list with the following elements: `cell_imputed_score` (a vector of
length `nrow(cell_features)`) that denotes the predicted progenies
spawning from each particular cell, `coefficient_vec` (the coefficient
vector of lenght `ncol(cell_features)+1`) that denotes the coefficients
in the GLM, `lambda` (the chosen parameter after cross-validation), and
`lineage_imputed_count` (the vector of length `lineage_future_count`
that denotes the number of predicted cells at the future timepoint in
each lineage).
