# Attach a real future time point to a simulated current one

The other two simulators stop at a future *count* per lineage. This one
goes further and hands each individual future cell a parent: given a
real matrix of future-time-point cells, it decides which current cell
each of them descended from. That is what makes it possible to simulate
quantities defined across the bottleneck — the adaptation index in
particular, which needs the embedding shift from parent to progeny and
so cannot be computed from counts alone.

## Usage

``` r
generate_simulation_attachFuture(
  coefficient_intercept,
  embedding_coefficient_vec,
  future_cell_embedding_mat,
  lineage_assignment,
  previous_cell_embedding_mat,
  fatefeatures_coefficient_vec = NULL,
  fatefeatures_mat = NULL,
  lineage_spread = 1,
  num_pushforward_training_iter = 20,
  num_subsamples = 200,
  verbose = 0
)
```

## Arguments

- coefficient_intercept:

  The starting intercept, natural-log scale. See step 1 — it is
  rescaled, so this only sets where the search begins.

- embedding_coefficient_vec:

  True coefficients on the embedding, length
  `ncol(previous_cell_embedding_mat)`.

- future_cell_embedding_mat:

  Numeric matrix of the future time point's **real** cells, rows = cells
  and columns = embedding dimensions (at least 2, and the same
  dimensions as the previous matrix). Row names are required for the
  assignment step.

- lineage_assignment:

  A **factor** of lineage membership for the current cells (asserted),
  row-aligned with `previous_cell_embedding_mat`. Usually taken from a
  [`generate_simulation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation.md)
  run, so that this function attaches a future to an already-simulated
  present. Dropped cells are removed with
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html).

- previous_cell_embedding_mat:

  Numeric matrix of the current time point's cells, rows = cells. Row
  names required.

- fatefeatures_coefficient_vec, fatefeatures_mat:

  Optional fate-driving signal outside the embedding, as in
  [`generate_simulation()`](https://nancyrzhanglab.github.io/multiomeFate/reference/generate_simulation.md).
  Default `NULL`.

- lineage_spread:

  Multiplier on the covariance used to score future cells against a
  pushed-forward current cell. Default `1`. Larger values make parentage
  more diffuse.

- num_pushforward_training_iter:

  Number of random restarts when fitting the push-forward map; the
  restart with the lowest squared error wins. Default `20`.

- num_subsamples:

  Number of (current, future) pairs drawn per restart when fitting the
  map. Default `200`. Current cells are sampled with probability
  proportional to their expected progeny count, so the map is fitted
  where the descendants actually come from.

- verbose:

  A numeric; larger values print more. Default `0`.

## Value

An object of class `"multiomeFate_simulation_future"`, a list with:

- `cell_fate_potential`:

  named numeric, `log10(expected progeny)` per surviving current cell,
  under the *adjusted* intercept.

- `coefficient_intercept`:

  the adjusted intercept, not the one passed in.

- `future_cell_assignment`:

  character vector of length `nrow(future_cell_embedding_mat)`, named by
  future cell ID, giving the current cell ID it descends from. Compose
  with `lineage_assignment` to get each future cell's lineage.

- `future_lineage_size`:

  named numeric, the realized number of progeny per lineage — counted
  from the assignment, so it sums exactly to the number of future cells.

- `mapping_mat`:

  the soft parent-child probability matrix, current cells by future
  cells, **scaled by 1e3 and rounded** to keep it storable. Divide by
  1e3 to recover probabilities. Columns are reordered by column sum, so
  they are not in the input order.

- `prev_cell_num_progenitor`:

  named numeric, how many future cells each current cell parented.
  Despite the name this counts *progeny*, not progenitors.

## Details

The construction:

1.  The intercept is **rescaled** so the expected progeny total equals
    the number of future cells actually supplied. The intercept the
    caller passes therefore does not survive; the adjusted one is
    returned. Current cells whose rounded contribution is 0 are dropped
    outright.

2.  A push-forward map from current to future embedding space is fitted:
    a single shared scale `a` and a per-dimension offset `b`, estimated
    by median regression over subsamples. It is deliberately this rigid
    — a flexible map would absorb the very state change the adaptation
    index is meant to detect.

3.  Each current cell is pushed forward and scored against every future
    cell by a Gaussian density, giving a soft parent-child mapping.

4.  Each future cell samples a parent from that mapping, with a parent
    removed once it has produced its allotted number of progeny.

**Stochastic with no `seed_number` argument**; steps 2 to 4 all draw.

## Examples

``` r
set.seed(10)
d <- 5
previous_cell_embedding_mat <- matrix(stats::rnorm(100 * d), nrow = 100, ncol = d)
rownames(previous_cell_embedding_mat) <- paste0("prev:", seq_len(100))
future_cell_embedding_mat <- matrix(stats::rnorm(200 * d), nrow = 200, ncol = d)
rownames(future_cell_embedding_mat) <- paste0("fut:", seq_len(200))
lineage_assignment <- factor(sample(paste0("lineage:", 1:5), size = 100,
                                    replace = TRUE))
names(lineage_assignment) <- rownames(previous_cell_embedding_mat)
res <- generate_simulation_attachFuture(
  coefficient_intercept = 0,
  embedding_coefficient_vec = rep(1, d),
  future_cell_embedding_mat = future_cell_embedding_mat,
  lineage_assignment = lineage_assignment,
  previous_cell_embedding_mat = previous_cell_embedding_mat,
  verbose = 0)
names(res)
#> [1] "cell_fate_potential"      "coefficient_intercept"   
#> [3] "future_cell_assignment"   "future_lineage_size"     
#> [5] "mapping_mat"              "prev_cell_num_progenitor"
# each future cell is attached to one previous cell, hence one lineage
table(res$future_lineage_assignment)
#> < table of extent 0 >
```
