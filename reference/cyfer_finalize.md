# Finalize CYFER after cross-validation

This function is used after running
[`multiomeFate::cyfer()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md).
Chooses `lambda` by minimizing the median held-out objective across
folds, then refits once on all cells at the chosen `lambda`.

## Usage

``` r
cyfer_finalize(
  cell_features,
  cell_lineage,
  fit_res,
  lineage_future_count,
  seed_number = 10
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
  `cyfer_finalize()` when scoring). A constant column you supply
  yourself is therefore duplicated, making the design collinear — and a
  column already named `Intercept` produces two columns of that name.
  `cyfer_finalize()` errors on any constant column for this reason;
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

- fit_res:

  This is the output of
  [`cyfer()`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md).

- lineage_future_count:

  A named numeric vector (where the names are the lineage names that
  appeared in `cell_lineage`) that denotes the number of cells at the
  future time point for each lineage.

- seed_number:

  Seed value for reproducibility reasons. Default is `10`. Governs the
  optimizer's random restarts in the final refit. Set to `NULL` to leave
  the random number stream untouched.

## Value

A list with the following elements: `cell_imputed_score` (a named vector
of length `nrow(cell_features)`) that denotes the predicted progenies
spawning from each particular cell, `coefficient_vec` (the coefficient
vector of length `ncol(cell_features)+1`) that denotes the coefficients
in the GLM, `lambda` (the chosen parameter after cross-validation), and
`lineage_imputed_count` (the vector of length `lineage_future_count`
that denotes the number of predicted cells at the future time point in
each lineage).

Every cell in `cell_features` is scored and named by its row name, even
cells whose lineage is absent from `lineage_future_count` or is `NA` (an
unassigned cell). Such cells are dropped from the refit (they carry no
future count to fit against) but the fitted coefficients still apply to
them, so the returned score covers every row that was passed in. A
message reports how many `NA`-lineage cells were excluded.

## Scales — exp() versus 10^()

Three scales are in play here and they are easy to confuse. The GLM's
linear predictor is on the **natural-log** scale,

\$\$Z_i = \beta_0 + X\_{i,\cdot}^\top \beta,\$\$

so `coefficient_vec` is on that scale, and cell `i`'s expected number of
progeny at the future time point is `exp(Z_i)` — with
[`exp()`](https://rdrr.io/r/base/Log.html), not `10^`.

The two vectors this function returns then sit on *different* scales,
deliberately:

- `cell_imputed_score`:

  **log10** of the expected progeny count, that is `log10(exp(Z_i))`.
  This is the per-cell fate potential the paper reports, and what the
  plotting helpers and the downstream selection/adaptation indices
  expect.

- `lineage_imputed_count`:

  a count on the **natural** scale, `sum(exp(Z_i))` over the cells of
  the lineage, directly comparable to `lineage_future_count`.

Because the cell scores have already been converted to base 10, the
identity linking the two uses `10^` and not
[`exp()`](https://rdrr.io/r/base/Log.html):

    lineage_imputed_count[l] == sum(10^cell_imputed_score[cells in l])

Applying [`exp()`](https://rdrr.io/r/base/Log.html) to
`cell_imputed_score` is the mistake this note exists to prevent: it
computes `exp(log10(exp(Z_i)))`, which is a count on no scale at all and
is silently plausible-looking. To go back: `10^cell_imputed_score`
recovers the expected progeny count, and `log(10)*cell_imputed_score`
recovers the linear predictor `Z_i`.

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
fit <- cyfer_finalize(cell_features = priming_simulation$cell_features,
                      cell_lineage = priming_simulation$cell_lineage,
                      fit_res = cv,
                      lineage_future_count = priming_simulation$lineage_future_count)
fit$lambda
#> [1] 0
head(fit$cell_imputed_score)                 # log10 scale
#>     cell:1     cell:2     cell:3     cell:4     cell:5     cell:6 
#> -0.2578001 -0.1488415 -0.1168919 -0.2955593  0.2184627 -0.1935971 
head(fit$lineage_imputed_count)              # natural-scale counts
#>  lineage:1 lineage:10 lineage:11 lineage:12 lineage:13 lineage:14 
#>   341.9535   171.9057   143.5195   135.6339   157.0583   126.2237 
# the identity linking the two scales
lineage <- names(fit$lineage_imputed_count)[1]
idx <- which(priming_simulation$cell_lineage == lineage)
all.equal(sum(10^fit$cell_imputed_score[idx]),
          unname(fit$lineage_imputed_count[lineage]))
#> [1] TRUE
# }
```
