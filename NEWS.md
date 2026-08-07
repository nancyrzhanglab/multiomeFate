# multiomeFate 1.0.2.002

The API changes and defect fixes agreed in
`additional_context/test-plan-full_2026-08-06_kevin.md`, which extended the test
suite from the estimation core to the barcoding pipeline, the simulators, the
plotting functions, `compute_entropy()` and the numerical utilities. Writing
those tests is what surfaced most of the defects below.

## API changes

* **The public API grows from 6 functions to 17.** Newly exported: the four
  plotting functions (`plot_anova()`, `plot_trainTest()`,
  `plot_cellGrowthUmap()`, `plot_lineageScatterplot()`), all three simulators
  (`generate_simulation()`, `generate_simulation_plastic()`,
  `generate_simulation_attachFuture()`), `compute_entropy()`, and the four
  barcoding functions (`barcoding_posterior()`, `barcode_clustering()`,
  `barcode_combine()`, `barcoding_assignment()`). All were already user-facing in
  practice, reached with `:::` from the analysis code, and now have help pages.

* **`data_loader()` is no longer exported.** It reads `.RData` from a hardcoded
  lab path and is a helper for the package authors rather than public API. Still
  available as `multiomeFate:::data_loader()`.

* **`evaluate_loglikelihood()` is renamed `evaluate_nll()`**, because it returns
  the *negative* penalized log-likelihood -- lower is better -- and the old name
  said the opposite. **There is deliberately no deprecated alias**: calls to the
  old name error rather than silently return a value whose sign the caller may
  have misread. It remains unexported. Analysis code calling
  `multiomeFate:::evaluate_loglikelihood()` must be updated.

* **Four unused helpers were removed from `R/util.R`**: `.mult_vec_mat()`,
  `.mult_mat_vec()`, `.log_sum_exp()` and `.exp_ratio()`. None was called from
  anywhere in `R/`. `.log_sum_exp()` in particular invited confusion with
  `.log_sum_exp_normalization()` in `R/simulation.R`, which is a different
  function (vector out, not scalar, and it propagates `NA` rather than dropping
  it); that one is unaffected.

* `stats`, `utils`, `graphics` and `methods` are now declared in `Imports:`.
  They are used throughout `R/` and were previously undeclared, which `R CMD
  check --as-cran` flags.

## Bug fixes -- plotting

* **`plot_trainTest()` and `plot_lineageScatterplot()` errored on every call.**
  Three call sites passed the `rlang` `.data` pronoun to `base::subset()`, which
  evaluates outside a data mask; current `rlang` aborts with `"Can't subset
  .data outside of a data mask context"` rather than resolving it. Both
  functions now use plain `df$column`. (`.data` is still used inside `aes()`,
  where it is required.)

* **`plot_anova()` errored whenever fewer than nine lineages qualified.** The
  bottom slice `(length(x) - num_lineages_bottom + 1):length(x)` went negative
  at the default `num_lineages_bottom = 10`, and R refuses to mix negative and
  positive subscripts. Both ends of the selection are now clamped, so asking for
  more lineages than exist falls back to all of them. The top end previously ran
  off the other side and put a phantom `"NA (NA)"` category on the x axis.

* **`plot_anova()` never forwarded `bool_add_future_size`** to
  `.plot_anova_helper()`, so `FALSE` still drew the `"lineage (n)"` labels.

* `.plot_trainTest_helper()` now passes `linewidth` rather than the `size`
  aesthetic, which `ggplot2` deprecated for lines in 3.4.0.

* `plot_simplex()` now drops rows whose three values sum to zero, with a warning
  naming how many went, instead of mapping them to `NaN`.

* `plot_trainTest()` now requires `quantile_vec[2] == 0.5` and errors otherwise.
  The centre curve is what the marked lambda is chosen by, and
  `cyfer_finalize()` always refits at the median, so any other centre drew a
  dashed line at a lambda the refit never used. The two outer entries are still
  free.

## Bug fixes -- barcoding

* **`barcode_combine(lin_mat, NULL)` errored** at
  `stopifnot(is.list(lineage_clusters))`, because `is.list(NULL)` is `FALSE`.
  The documented "nothing to do" contract was therefore unreachable -- and
  `barcode_clustering()` returns exactly `NULL` when no pair clears
  `cor_threshold`, so piping the two together crashed on any dataset with no
  correlated barcodes. `NULL`, an empty list, and a list of nothing but `NA`
  tombstones now all return `lin_mat` unchanged.

* **`barcode_combine()` silently renamed a barcode** when exactly one was left
  unclustered: the single row collapsed to a vector and `rbind()` named it after
  the variable, so the output row was called `"lin_untouched"`. A cluster naming
  a single barcode errored for the same reason. Both fixed with `drop = FALSE`.

* **`barcoding_posterior()` returned `beta0` and `beta1` unnamed** -- the code
  assigned `names(lin_mat)`, which is `NULL` for a matrix. They are now named by
  barcode, as `gamma` already was, along with `lineage_num_winner`.

* **`barcoding_posterior()` now accepts a `dgCMatrix`**, so the whole pipeline
  works in sparse form and callers no longer need to bridge
  `barcode_clustering()` and `barcoding_posterior()` with `as.matrix()`. Row
  names are now required rather than assumed.

* `barcode_clustering()`'s cluster-merge assertion is now an explicit `stop()`
  naming the invariant, rather than a bare `stopifnot()`.

## Bug fixes -- simulation

* **`generate_simulation_attachFuture()` fitted its push-forward map from a
  single restart**, never `num_pushforward_training_iter` of them: `lapply()`
  iterated over the *scalar* rather than `seq_len()` of it, so the `which.min()`
  selecting the best restart was a no-op.

* **`.compute_previous_to_future_mapping()` built its covariance from unsquared
  standard deviations**, `lineage_spread * diag(sd_vec)`, where
  `.form_gaussian_distribution()` uses variances. The parent-child kernel was
  therefore too tight, and `lineage_spread` did not mean the same thing in the
  two files. Now `lineage_spread * diag(sd_vec^2)`.

  Together with the restart fix above this changes every parent-child
  assignment: `mapping_mat`, `future_cell_assignment`,
  `prev_cell_num_progenitor` and `future_lineage_size` all move. The
  `cell_fate_potential` and `coefficient_intercept` elements are computed before
  the push-forward step and are unaffected.

* `generate_simulation_plastic()` now **errors** when `lineage_mean_spread` and
  `lineage_sd_spread` are both `NA`, rather than warning. That configuration
  gives lineages large means *and* large variances together, which confounds the
  two regimes the simulation exists to separate.

* `.compute_pushforward()` no longer declares an unused
  `previous_cell_potential` argument.

## Bug fixes -- other

* `compute_entropy()` now subsets with `drop = FALSE` near the end. With exactly
  one surviving cell the matrix collapsed to a vector, `rownames()` became
  `NULL`, and the metadata lookup that follows silently produced garbage.

* `.nonzero_col()` was defined twice, byte-identically, in `R/util.R` and
  `R/barcoding_assignment.R`. The duplicate is removed.

## Testing

* The suite grows from 637 to 1149 passing tests, covering the barcoding
  pipeline, the three simulators, the plotting functions, `compute_entropy()`
  and `data_loader()`'s pure helper -- all of which previously had little or no
  coverage. Plotting had none at all, which is why the two defects that broke
  `plot_trainTest()` and `plot_lineageScatterplot()` outright had gone
  unnoticed.

* New: a test that the two simulators produce the two regimes they are named
  for (priming spreads lineage *means*, plastic spreads within-lineage
  *variances*); an end-to-end test that `cyfer()` recovers the simulated fate
  potential; and a test pinning the exported-name list, so an accidental
  `@export` cannot slip in unnoticed.

* The five `_claude`-suffixed test files are folded into their unsuffixed
  counterparts.

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
