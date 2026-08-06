# multiomeFate 1.0.2.001

Fixes for the defects pinned by the test suite added in
`additional_context/test-plan_2026-08-06_kevin.md`.

## Bug fixes

* `.construct_lineage_data()` (the internal test-data generator) now draws counts
  from the CYFER model. `exp(coefficient_vec %*% cell_features[i,,drop=F])`
  conformed the bare length-`p` vector as `p x 1` against a `1 x p` matrix and so
  evaluated the `p x p` **outer** product, summing `exp()` over `p^2` entries.
  This was correct only at `p == 1`. Test-fixture defect only — no estimation
  code and no published result is affected — but no fixture-based test written
  before this fix may be read as evidence about estimation accuracy.

* `.lineage_objective()` now errors on a non-finite value, as
  `.lineage_gradient()` already did. Because `optim()` evaluates `fn` before
  `gr`, a caller with unscaled features previously got `"initial value in 'vmmin'
  is not finite"` rather than the gradient's actionable diagnostic; and at CV
  evaluation time a non-finite `test_loglik` was silently skipped by
  `which.min()`. A non-finite *trial* point inside the BFGS line search is still
  routine and `optim()` backs off from it as before.

* `lineage_imputation()`'s random-restart range no longer collapses. When
  `max_count_ratio == 1` (no net growth) `max_limit` and `min_value` were both
  `0`, so all ten restarts were the same zero vector; `min_value` is now
  `2*(max_limit-1)` whenever `max_limit <= 1e-6`. When `max_count_ratio == 0`
  (all future counts zero) `log(0) = -Inf` made every draw `NaN` and `optim()`
  died with `"non-finite value supplied by optim"`; that case is now hard-coded
  to `[-10, 0]`. The growth regime (`max_limit > 1e-6`) keeps `min_value = 0` and
  the shrinking regime keeps a wholly negative range, both unchanged.

* `.compute_initial_parameters()` smooths the growth ratio by one count, so an
  all-zero `lineage_future_count` no longer yields `lambda_initial = NaN` and
  `coefficient_initial["Intercept"] = -Inf` silently. This shifts
  `lambda_initial` and the intercept start on *every* input, not only the
  degenerate one; in practice the heuristic still saturates at the
  `lambda_max = 101` cap on typical data.

* `.compute_initial_parameters()`'s `multipler` now defaults to `1e4`, matching
  what its only caller has always passed.

* `cyfer()` drops lineages with no cells at the first time point *before*
  building folds. Such a lineage was dealt into a fold where it contributed
  nothing, and a fold made up entirely of them left `cv_cell_list[[fold]]` as
  `NULL`, so `cell_features[-NULL,,drop=F]` trained on zero rows.

* `cyfer_finalize()` supports a length-1 lambda path. `sapply()` returned a
  vector rather than a matrix, so `apply(test_vec, 1, median)` failed with
  `"dim(X) must have a positive length"`.

* `cyfer_finalize()` errors when `cell_features` carries a constant column. It
  prepends its own intercept unconditionally, so such a column produced a
  collinear design — and one already named `Intercept` produced two columns of
  that name, caught only by a bare `stopifnot()`. It also now requires row and
  column names, as `cyfer()` does.

## New features

* `cyfer_finalize()` gains a `seed_number` argument (default `10`, `NULL` to
  leave the stream untouched). The final refit uses ten random restarts, which
  were previously unseeded even though `cyfer()`'s folds were reproducible.

## Documentation

* `?cyfer_finalize` gains a **"Scales — exp() versus 10^()"** section spelling
  out the three scales in play: the linear predictor `Z_i` and `coefficient_vec`
  are on the natural-log scale, so expected progeny per cell is `exp(Z_i)`;
  `cell_imputed_score` is `log10(exp(Z_i))`; and `lineage_imputed_count` is a
  natural-scale count. The identity linking the last two is
  `lineage_imputed_count[l] == sum(10^cell_imputed_score[cells in l])` — with
  `10^`, not `exp()`. Calling `exp()` on `cell_imputed_score` is the specific
  mistake the section exists to prevent, since the result looks plausible and is
  a count on no scale at all.
* The shared `cell_features` documentation (defined on `cyfer()` and inherited by
  `lineage_imputation()`, `lineage_imputation_sequence()`, and
  `cyfer_finalize()`) now states that no intercept column and no constant column
  of any kind may be supplied, explains why (the intercept is added internally,
  so a supplied constant column is duplicated and the design becomes collinear),
  and records that `cyfer_finalize()` errors on one while `cyfer()` tolerates it.
* `?cyfer_finalize` documents that every cell is scored and named, including
  cells whose lineage is absent from `lineage_future_count` and therefore
  excluded from the refit.
* `?lineage_imputation` and `?lineage_imputation_sequence` state that their
  `coefficient_vec` is on the natural-log scale and cross-reference the Scales
  section; the latter also pins that `lambda_sequence` is strictly decreasing
  from `lambda_initial` to `0` and that each fit warm-starts from the previous.
* `?cyfer` states that `test_loglik` / `train_loglik` are the unpenalized
  objective (lower is better) rather than literal log-likelihoods, and that
  cell-less lineages are dropped before the folds are built.

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
