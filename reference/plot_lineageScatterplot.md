# Predicted versus observed lineage size

The calibration check on a CYFER fit: each point is one lineage, its
observed future count on the x axis against the count CYFER predicts on
the y axis, both on a log10 scale with a fixed 1:1 aspect ratio so that
departures from the diagonal are readable. The largest lineages on
either axis are labelled. The title carries the Pearson correlation
between the two log10 vectors.

## Usage

``` r
plot_lineageScatterplot(
  lineage_future_count,
  lineage_imputed_count,
  num_lineage = 10,
  threshold_x = 1.5,
  threshold_y = 1.5,
  title = ""
)
```

## Arguments

- lineage_future_count:

  Named numeric vector of observed future counts per lineage.

- lineage_imputed_count:

  Named numeric vector of predicted counts per lineage — the
  `lineage_imputed_count` element of
  [`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md),
  which is on the natural count scale, not log10. The two vectors must
  be **in the same order with identical names**; this is asserted by
  position, not matched by name.

- num_lineage:

  Maximum number of lineages to label on each axis, taken in decreasing
  order of that axis. Default `10`. Up to `2 * num_lineage` labels can
  therefore appear.

- threshold_x, threshold_y:

  Minimum log10 count for a lineage to be eligible for labelling on the
  observed and predicted axis respectively. Default `1.5` each, i.e.
  about 32 cells. `NA` disables labelling from that axis.

- title:

  Plot title; the correlation is appended on a second line. Default
  `""`.

## Value

A `ggplot` object.

## Details

**The x coordinate is jittered and the plot is stochastic.** Observed
counts are small integers, so many lineages would otherwise stack on the
same vertical line; `Uniform(0, 0.5)` is added *before* the log10, which
spreads the small counts more than the large ones. There is no
`seed_number` argument, so set a seed beforehand for a reproducible
figure. The reported correlation is computed on the *unjittered* values,
so it does not move between calls even though the points do.

## Examples

``` r
lineage_future_count <- c(L1 = 5, L2 = 40, L3 = 100, L4 = 3, L5 = 60)
lineage_imputed_count <- c(L1 = 4, L2 = 52, L3 = 88, L4 = 6, L5 = 45)
set.seed(10)
plot_lineageScatterplot(lineage_future_count = lineage_future_count,
                        lineage_imputed_count = lineage_imputed_count,
                        num_lineage = 3)
```
