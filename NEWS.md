# multiomeFate 1.0.2.000

Cross-validation bug fixes. **Any `cyfer()` result produced by an earlier version
should be re-run.** Reported by Raymond Ng (Shaffer lab) via Emilia Chen, July 2026.

## Bug fixes

* `cyfer()` now returns correct results when `cell_lineage` is a **factor**.
  `.lineage_gradient()` indexed the named vectors `lineage_future_count` and
  `denom_vec` with `cell_lineage` directly; when that is a factor, R indexes by
  the factor's integer codes rather than its labels. On a full-data fit the level
  order coincides with the sorted count vector, so the result was correct by
  accident. Inside a cross-validation fold it was not: subsetting a factor
  retains its unused levels, while `lineage_future_count` is trimmed to the fold,
  so the lookup returned `NA` and the whole gradient became `NaN`.
  `optim(method = "BFGS")` then returned its starting value while still
  reporting `convergence = 0`. The visible symptom was that the selected lambda
  was always the largest value in the path, no matter how large `lambda_initial`
  was set — held-out log-likelihood was bit-identical at every lambda, and
  `which.min()` on an exact tie returns the first (largest) element.
  `cell_lineage` is now coerced to character at every entry point, and
  `.lineage_gradient()` raises an error rather than returning a `NaN` gradient.

* `construct_folds()` now assigns every lineage to exactly one fold. It
  previously had three independent defects: the block-shuffling loop iterated
  `num_per_fold` times instead of once per block (silently `NA`-padding past the
  end of the vector); out-of-range fold indices were clamped with `pmin()`
  instead of dropped (placing the last lineage in *every* fold); and
  `sample(idx_vec)` on a single-element block permuted `seq_len(idx_vec)` rather
  than returning the index. Depending on the lineage/fold counts this dropped a
  lineage from held-out evaluation entirely or duplicated one across all folds.
  A post-condition now asserts that the folds partition the lineages.

* `cyfer()` now sets the seed *before* calling `construct_folds()`, so
  `seed_number` makes the fold assignment reproducible. Previously the folds
  depended on whatever RNG state the caller happened to leave behind.

* `.lineage_cleanup()` now reconciles *partial* mismatches between
  `cell_lineage` and `lineage_future_count`. The guard read `all(a != b)`, which
  is `TRUE` only when every element differs, so the reconciliation branch almost
  never fired. It now uses `setequal()`, as does the corresponding `stopifnot()`
  in `lineage_imputation()`.

* The `lambda_min` / `lambda_max` clamp on the auto-computed `lambda_initial` is
  now functional. Both defaulted to `101` and were applied in the wrong order, so
  `min(max(lambda_initial, lambda_max), lambda_min)` collapsed to the constant
  `101` for any input. `lambda_min` now defaults to `0.01` and the clamp floors
  at `lambda_min` and caps at `lambda_max`. In practice the heuristic still
  saturates at the `lambda_max = 101` cap on typical data, so this changes
  results only when the heuristic lands below 101 (previously it was raised to
  101 regardless). Affects only callers who leave `lambda_initial = NA`.

* `construct_folds()` now validates its inputs instead of failing obscurely
  downstream: `num_folds` must be between 2 and the number of lineages, and
  `lineage_future_count` must be a non-empty named vector without duplicate
  names. Previously `num_folds > num_lineages` was masked by the `pmin()` clamp
  above; with that clamp correctly removed it would have produced empty folds and
  an "invalid argument to unary operator" error from `cyfer()`.

## Testing

* New `tests/testthat/test_construct_folds.R`, including a partition invariant
  checked across a grid of lineage counts and fold counts.
* New regression tests covering factor-valued `cell_lineage` end to end, the
  held-out log-likelihood varying across the lambda path, the selected lambda not
  pinning to the top of the path, and fold reproducibility under `seed_number`.
