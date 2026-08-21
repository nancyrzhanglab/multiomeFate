# multiomeFate (development version)

Tentative, unreleased. Changes on branch `emilia_review` (emiliac), from the
code review in `Reports/CYFER_code_review_2026-08-14.pdf`; each entry is recorded
in full in `Reports/CYFER_fix_log_2026-08-17.pdf`.

## Bug fixes

* **`cyfer()` now resolves `lambda_initial` once, on the full data, when it is
  passed as `NA`**. Previously the `NA` fallback was taken separately
  inside each fold, by `.compute_initial_parameters()` acting on that fold's
  *training* lineages — and every input to that heuristic (`future_total`,
  `current_total`, `term2`, `num_lineages`) is a sum or count over the lineages it
  is handed. Each fold therefore built a **different `lambda_sequence`**, while
  `cyfer_finalize()` stacks the folds' `test_loglik` vectors by *position*,
  medians down each row, and reads the winning lambda off `fit_res[[1]]`'s grid
  alone. Row `k` was thus an average of held-out error measured at `k` different
  penalties, and the reported lambda was fold 1's value for that position. Nothing
  errored, because all the paths have the same *length*.

* **`cyfer_finalize()` now asserts that every fold was fit on the same
  `lambda_sequence`**, instead of assuming it and reading the path off fold 1.
  This cannot fire on a fresh run — the fix above guarantees the shared grid —
  but `fit_res` objects saved *before* it, including `savefile_tmp`
  checkpoints, do carry a different grid per fold and are still on disk.
  A `fit_res` that is not `cyfer()` output at all is also rejected rather than
  reaching `fit_res[[1]]`.

* `cyfer()` now validates `lambda_initial`: it must be a single number,
  or `NA`. `c(1, 2)`, `"3"` and `NULL` previously reached `is.na()` and behaved
  unpredictably from there.

* **A broken `\link{data_loader}` cross-reference is fixed** in
  `barcoding_assignment`'s help page. `data_loader()` was unexported in
  1.0.2.001, which left the link with no target, so `R CMD check` reported
  `Missing link or links in Rd file`.

* **`man/lineage_imputation.Rd` regenerated** to document `maxit`, which was
  added to the roxygen but never propagated.

* **`.lineage_cleanup()` now checks that `cell_lineage` is row-aligned with
  `cell_features`.

  The **length** of `cell_lineage` against `nrow(cell_features)` is
  always verified. And **when `cell_lineage` carries names**, they must equal 
  `rownames(cell_features)` exactly.

  **An unnamed `cell_lineage` passes**, exactly as before. Passing no names is
  the caller asserting the vector is already row-aligned. 

  `cyfer()` and `cyfer_finalize()` coerce with
  `setNames(as.character(cell_lineage), names(cell_lineage))` rather than plain
  `as.character()`, which drops names. 

* **`cyfer()` now refuses a fit whose unpenalized endpoint is not identified.**
  CYFER's effective sample size is the number of *lineages*.

  `cyfer()` now checks the number of training lineages against `ncol(cell_features)
  + 1` after the folds are built, comparing against the *largest* fold since that
  leaves the smallest training set. The identifiability floor is `L(k-1)/k >= p+1`, 
  which is much weaker than one clone per feature — at `p = 60`, `k = 10` it needs 
  only `L >= 68`.

* **The random restarts in `lineage_imputation()` are drawn from a range that
  scales as `1/sqrt(p)` rather than `1/p`.** The range is a per-coefficient
  allowance carved out of a budget on the *total* linear predictor — keep
  `|x_i' beta|` within about `2*log(max_count_ratio)` so `exp()` cannot overflow
  — and converting the budget into a per-coefficient figure needs a model of how
  `p` contributions add. The old `1/p` came from a triangle-inequality worst
  case, `|x' beta| <= p * max|x| * max|beta|`, which requires every term to have
  the same sign and the same magnitude. With mixed signs they partly cancel and
  the sum grows like `sqrt(p)`, which is the reasoning behind Xavier/Glorot and
  He initialization.

  The practical effect was that the range collapsed as features were added. On
  scaled features with a maximum growth ratio of order 250:

  ```
  p+1    old (1/p)   new (1/sqrt p)   ratio
    3       1.8405           3.1878    1.7x
   11       0.5020           1.6648    3.3x
   31       0.1781           0.9917    5.6x
   61       0.0905           0.7070    7.8x
  ```

  At `p+1 = 61` all ten restarts drew from a range about 0.09 wide, so they were
  one search of a non-convex objective repeated ten times at ten times the cost.
  The package's own fixtures use 2 features, where the range is wide and the
  mechanism works as intended, which is why no test caught this.

  **This changes fitted numbers**: the restarts start elsewhere, so `optim()`
  follows different paths. On a synthetic `L = 100`, `p = 60`, 800-cell fit the
  wider range did disperse the search — the spread of objective values across
  restarts rose from 31.1 to 48.9 — but the best objective found was identical to
  six decimal places at both `lambda = 1` and `lambda = 0`. The change buys
  coverage of the parameter space; it is not yet demonstrated to find better
  optima on real data.

  The range is still floored at zero rather than centred on it, so restarts still
  begin in the non-negative orthant. That is a separate decision and is
  deliberately left unchanged. `upper_randomness` still does not bind at
  realistic feature counts (0.71 against a cap of 5 at `p+1 = 61`), though it can
  now bind at small `p` with a large growth ratio.

* **`.lineage_cleanup()` now validates `lineage_future_count`.** The response was
  unchecked for no negative count, no `Inf`, no `NA` / `NaN`, no duplicated names,
  and no unnamed or character.

* **`lineage_imputation()` gains a `maxit` argument, defaulting to
  `min(500, max(100, 10*p))`, and `optim()`'s convergence code is now read.** BFGS
  searches a `p`-dimensional space, so its iteration budget has to grow with `p`.
  `optim()`'s own default of 100 was left in place, and the `convergence` code it
  returns was stored on every fit and read nowhere in `R/` — so a fit that ran out
  of iterations was indistinguishable from one that converged, and its
  coefficients are wherever BFGS happened to be when the budget ran out.

  Measured on synthetic data, the fraction of stored fits hitting the cap:

  ```
  L=100, p=60,  800 cells   20.0%      (L/p ~ 1.6)
  L=400, p=60, 3200 cells    2.0%      (L/p ~ 6.6)
  L=100, p=10,  800 cells    0.0%      (L/p ~ 9.1)
  ```

  The convergence code is now surfaced. If the **selected** fit did not converge,
  `lineage_imputation()` warns, naming the code, the `maxit` it hit and the
  `lambda`, and saying that its coefficients are where the optimizer stopped
  rather than an optimum. If only some restarts failed, that is reported at
  `verbose > 0`. Warning on every non-converged restart was rejected deliberately:
  a `cyfer()` run makes roughly 500 `lineage_imputation()` calls, and only the
  selected fit propagates into the result.

## Continuous integration

* **`R CMD check` now runs on every push and pull request**
  (`.github/workflows/R-CMD-check.yaml`). Previously `.github/workflows/` held
  only `pkgdown.yaml`, so nothing ran the test suite or `R CMD check`. The
  malformed `.Rd` that made 1.0.2.002 uninstallable would have been caught on
  the commit that introduced it, since `tools::parse_Rd()` runs inside the check.

  The matrix is ubuntu-latest `release` and `oldrel-1`. macOS and Windows are
  left out because `Seurat` and `scCustomize` resolve quickly only from Linux
  binaries; a comment in the workflow marks where to add them.

  Running the check before adding the workflow turned up three problems that
  would have failed the first CI run — the undocumented `maxit`, the broken
  `data_loader` link, and the test below. All three are fixed, and the check is
  clean apart from one pre-existing NOTE (`graphics` and `methods` are declared
  in `Imports` but never imported from).

* **`test_data_loader.R` no longer assumes the lab data is absent.** It called
  the real `data_loader()` and asserted the error was `cannot open|No such file`,
  which is true on CI but false on the lab machines, where the calls instead
  load multi-gigabyte `.RData` files one after another — enough to be killed by
  a 24 GB cgroup limit. It now skips when the data directory is present.

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
