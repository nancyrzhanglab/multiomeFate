# Simulate lineages under pure selection ("priming")

Generates the **priming** regime: a cell's fate potential is fixed by
its position in the embedding, and lineages are spatially localized
clumps of cells. Because cells near each other share both their
embedding and their potential, a lineage's cells are alike and the
between-lineage variation in expansion is large — the signature of
selection acting on pre-existing state. Contrast
[`generate_simulation_plastic()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation_plastic.md),
which builds lineages that are deliberately *heterogeneous* in
potential.

## Usage

``` r
generate_simulation(
  embedding_mat,
  bool_add_randomness = TRUE,
  coefficient_intercept = 0,
  embedding_coefficient_vec = rep(1, ncol(embedding_mat)),
  fatefeatures_coefficient_vec = NULL,
  fatefeatures_mat = NULL,
  lineage_spread = 1,
  lineage_prior = NA,
  num_lineages = 10,
  tol = 1e-06,
  verbose = 0
)
```

## Arguments

- embedding_mat:

  Numeric matrix of the current time point's cells, rows = cells and
  columns = embedding dimensions (at least 2). Row names, if present,
  are carried onto the returned vectors. Supplied by the caller rather
  than simulated, so that simulations sit on a real embedding.

- bool_add_randomness:

  Whether to draw the realized progeny counts from
  [`stats::rpois()`](https://rdrr.io/r/stats/Poisson.html) around the
  expected counts. Default `TRUE`. `FALSE` gives the noiseless best
  case, useful for isolating estimation error from sampling noise.

- coefficient_intercept:

  The intercept, on the natural-log scale. Sets the overall growth
  level; `0` means one expected progeny per cell before any feature
  contribution. Default `0`.

- embedding_coefficient_vec:

  True coefficients on the embedding, length `ncol(embedding_mat)`,
  natural-log scale. Default all ones.

- fatefeatures_coefficient_vec:

  True coefficients on the extra fate features, length
  `ncol(fatefeatures_mat)`. Default `NULL`.

- fatefeatures_mat:

  Optional numeric matrix of features that drive fate but are **not**
  part of the embedding, same rows as `embedding_mat`. This is how a
  simulation withholds signal from the estimator: CYFER fitted on
  `embedding_mat` alone cannot recover this contribution, which is the
  point when studying misspecification. Default `NULL`.

- lineage_spread:

  Multiplier on each embedding dimension's variance when forming the
  lineage Gaussians. Default `1`. Larger values make lineages broader
  and more overlapping, weakening the selection signal.

- lineage_prior:

  Numeric vector of length `num_lineages` of prior lineage
  probabilities; renormalized to sum to 1. Default `NA`, meaning
  uniform. Any names it carries are overwritten (with a warning) by
  `"lineage:1"`...

- num_lineages:

  Number of lineages. Default `10`.

- tol:

  Tolerance for the "prior sums to 1" check. Default `1e-06`.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

An object of class `"multiomeFate_simulation_vanilla"`, a list with:

- `cell_fate_potential`:

  named numeric, `log10(realized progeny + 1)` per cell. The `+1` keeps
  a zero-progeny cell finite, so this is *not* directly comparable to a
  `cell_imputed_score`.

- `cell_fate_potential_truth`:

  named numeric, `log10(expected progeny)` per cell, with no `+1`.
  **This** is the ground truth to score `cell_imputed_score` against.

- `coefficient_intercept`:

  as supplied.

- `embedding_mat`:

  as supplied.

- `fatefeatures_coefficient_vec`, `fatefeatures_mat`:

  as supplied.

- `gaussian_list`:

  list of `num_lineages` objects of class `"gaussian"`, each with `mean`
  and `cov`.

- `lineage_assignment`:

  factor of length `nrow(embedding_mat)` with levels `"lineage:1"`...,
  named by cell. Pass
  [`as.character()`](https://rdrr.io/r/base/character.html) of this as
  `cell_lineage`.

- `lineage_future_size`:

  named numeric, the future count per lineage. This is
  `lineage_future_count`.

- `prob_mat`:

  cells-by-lineages posterior matrix used for the assignment draw.

- `summary_mat`:

  5-by-`num_lineages` matrix with rows `mean`, `median`, `sd`, `range`
  of the true log10 potential within each lineage, and `future_size`.
  The `sd` and `range` rows are the intra-clonal heterogeneity, and are
  what distinguish this regime from the plastic one.

## Details

This is how `data/priming_simulation.rda` was produced; see
[`?priming_simulation`](https://nancyrzhanglab.github.io/multiomeFate/reference/priming_simulation.md).

The construction, in the order the `verbose` messages report it:

1.  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) with
    `2 * num_lineages` centres, of which one random cell from each of
    the first `num_lineages` clusters becomes a lineage seed.
    Over-clustering makes the seeds tighter and more spread out than
    `num_lineages` centres would.

2.  Each seed anchors an isotropic Gaussian whose per-dimension variance
    is `lineage_spread` times the observed variance of that embedding
    dimension.

3.  Every cell gets a posterior over lineages from those Gaussians and
    `lineage_prior`, and is **sampled** from it — so lineages are soft,
    overlapping clumps, and their realized sizes are random rather than
    equal to `n * lineage_prior`.

4.  Each cell's expected progeny count is
    `exp(coefficient_intercept + x_i^T beta)`, optionally Poisson-drawn,
    and summed within lineage to give the future size.

**Stochastic with no `seed_number` argument** — steps 1, 3, and 4 all
draw. Set a seed before calling.

## Examples

``` r
set.seed(10)
embedding_mat <- matrix(stats::rnorm(300 * 5), nrow = 300, ncol = 5)
rownames(embedding_mat) <- paste0("cell:", seq_len(300))
sim <- generate_simulation(embedding_mat = embedding_mat,
                           coefficient_intercept = -1,
                           embedding_coefficient_vec = c(1, 0.5, 0, 0, 0),
                           lineage_prior = rep(0.2, 5),
                           num_lineages = 5)
table(sim$lineage_assignment)
#> 
#> lineage:1 lineage:2 lineage:3 lineage:4 lineage:5 
#>        55        55        82        32        76 
sim$lineage_future_size
#> lineage:1 lineage:2 lineage:3 lineage:4 lineage:5 
#>         9        42        37        57        53 
head(sim$cell_fate_potential)
#>  cell:1  cell:2  cell:3  cell:4  cell:5  cell:6 
#> 0.00000 0.00000 0.00000 0.30103 0.30103 0.00000 
```
