# Simulate lineages that differ in heterogeneity, not in mean ("plastic")

The counterpart to
[`generate_simulation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation.md).
There, lineages are spatial clumps and so differ in *mean* fate
potential — the selection regime. Here the assignment is driven by fate
potential directly rather than by position, and lineages are built to
have the **same mean potential but different spreads**: lineage 1
collects cells from the extremes of the potential distribution, the last
lineage collects cells near the middle. Between-lineage variation in
expansion is then small while within-lineage variation is large, which
is what a plastic (adaptive) population looks like and what a method
that only compares clone sizes would miss.

## Usage

``` r
generate_simulation_plastic(
  embedding_mat,
  bool_add_randomness = TRUE,
  coefficient_intercept = 0,
  embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
  fatefeatures_coefficient_vec = NULL,
  fatefeatures_mat = NULL,
  lineage_mean_spread = 1,
  lineage_sd_spread = NA,
  num_lineages = 10,
  tol = 1e-06,
  verbose = 0
)
```

## Arguments

- embedding_mat:

  Numeric matrix of the current time point's cells, rows = cells and
  columns = embedding dimensions (at least 2). Row names are generated
  as `"cell:1"`... if absent.

- bool_add_randomness:

  Whether to Poisson-draw the realized progeny counts. Default `TRUE`.

- coefficient_intercept:

  The intercept, natural-log scale. Default `0`.

- embedding_coefficient_vec:

  True coefficients on the embedding, length `ncol(embedding_mat)`,
  natural-log scale. Default all ones.

- fatefeatures_coefficient_vec:

  True coefficients on the extra fate features. Default `NULL`.

- fatefeatures_mat:

  Optional matrix of fate-driving features outside the embedding, same
  rows as `embedding_mat`; signal deliberately withheld from an
  estimator fitted on the embedding alone. Default `NULL`.

- lineage_mean_spread:

  Controls whether lineage *means* are allowed to differ. Only two
  values are honoured: `1` (the default) holds every lineage's mean at
  the population mean, which is the plastic regime and also turns on the
  equal-size constraint in step 3; `NA` spreads the means across
  quantiles of the potential distribution, shrinks the working standard
  deviation to a quarter, and drops the equal-size constraint. Any other
  numeric warns and is treated as `1`.

- lineage_sd_spread:

  The ratio `rho` defining the spread ladder: lineage 1 gets standard
  deviation `sd * rho` and the last gets `sd / rho`, interpolated
  linearly in between, so values above 1 give the intended high-to-low
  ordering. Default `NA`, meaning derive it from the data as
  `max(|log potential - mean|) / sd / 2`. The value actually used comes
  back in the output.

- num_lineages:

  Number of lineages. Default `10`.

- tol:

  Unused; retained for signature compatibility with
  [`generate_simulation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation.md).
  Default `1e-06`.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

An object of class `"multiomeFate_simulation_plastic"`, a list with the
same elements as
[`generate_simulation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation.md)
except that `gaussian_list` is absent and `lineage_sd_spread` (the
realized `rho`) is present:

- `cell_fate_potential`:

  named numeric, `log10(realized progeny + 1)` per cell.

- `cell_fate_potential_truth`:

  named numeric, `log10(expected progeny)` per cell — the ground truth
  for scoring `cell_imputed_score`.

- `coefficient_intercept`, `embedding_mat`,
  `fatefeatures_coefficient_vec`, `fatefeatures_mat`:

  as supplied.

- `lineage_assignment`:

  factor named by cell, levels `"lineage:1"`..., reordered to match
  `cell_fate_potential_truth`.

- `lineage_future_size`:

  named numeric, the future count per lineage.

- `lineage_sd_spread`:

  the `rho` used, whether supplied or derived.

- `prob_mat`:

  cells-by-lineages posterior matrix. Its rows are in the internal
  reordered cell order, not the input order.

- `summary_mat`:

  5-by-`num_lineages` matrix, rows `mean`, `median`, `sd`, `range`,
  `future_size`. The check that the simulation did what it claims: the
  `mean` row should be near-flat across lineages while `sd` and `range`
  decrease.

## Details

This is how `data/plastic_simulation.rda` was produced; see
[`?plastic_simulation`](https://nancyrzhanglab.github.io/multiomeFate/reference/plastic_simulation.md).

Order of construction, and note it is the reverse of the priming
simulation — **potentials are computed first and lineages assigned from
them**, rather than lineages first and potentials from position:

1.  Each cell's expected progeny count is
    `exp(coefficient_intercept + x_i^T beta)`, optionally Poisson-drawn.

2.  Each lineage gets a Gaussian over *log potential*, all with the same
    mean and with standard deviations interpolating from `sd * rho` down
    to `sd / rho`; cells are scored against those.

3.  Cells are sampled into lineages from that posterior, sequentially,
    with a lineage removed from contention once it reaches
    `ceiling(n / num_lineages)` cells — so lineages come out near-equal
    in size.

**Stochastic with no `seed_number` argument**; set a seed before
calling.

## Examples

``` r
set.seed(10)
embedding_mat <- matrix(stats::rnorm(300 * 5), nrow = 300, ncol = 5)
rownames(embedding_mat) <- paste0("cell:", seq_len(300))
sim <- generate_simulation_plastic(embedding_mat = embedding_mat,
                                   coefficient_intercept = -1,
                                   embedding_coefficient_vec = c(1, 0.5, 0, 0, 0),
                                   num_lineages = 5)
table(sim$lineage_assignment)
#> 
#> lineage:1 lineage:2 lineage:3 lineage:4 lineage:5 
#>        60        60        60        60        60 
sim$lineage_future_size
#> lineage:1 lineage:2 lineage:3 lineage:4 lineage:5 
#>        49        31        45        33        27 
# per-lineage spread of the log potential, which is what "plastic" means here
sim$summary_mat
#>              lineage:1  lineage:2  lineage:3  lineage:4  lineage:5
#> mean        -0.4472834 -0.5211972 -0.4083935 -0.4329289 -0.4295606
#> median      -0.4980247 -0.4833333 -0.4147407 -0.4850237 -0.4034181
#> sd           0.5557908  0.5471189  0.4772796  0.3891652  0.3573510
#> range        2.8342158  2.1998376  2.2885444  1.7628528  1.7348575
#> future_size 49.0000000 31.0000000 45.0000000 33.0000000 27.0000000
```
