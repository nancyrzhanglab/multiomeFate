# Plastic simulation dataset

A small, fully synthetic dataset used to demonstrate and test
lineage-fate imputation in a \*\*plastic\*\* setting, where each lineage
has roughly the same median growth potential, but certain lineages have
a few cells with extremely high growth potential (fate propensity).
Downstream lineage sizes scale with \\\exp(\alpha + X\beta)\\.

All \*\*50 lineages\*\* of the simulation are included (7,940 cells).

## Usage

``` r
data(plastic_simulation)
```

## Format

A named list with four elements:

- cell_features:

  Numeric matrix of size \\n \times p\\. Rows are cells (rownames
  formatted as `"cell:1"`, `"cell:2"`, ...). Columns are preprocessed
  expression-derived features (no intercept column; modeling functions
  add it internally). No missing values.

- cell_lineage:

  Character vector of length \\n\\ giving each cell's lineage label used
  for training/evaluation.

- lineage_future_count:

  Named numeric vector (length = 50). Names are lineage IDs matching
  `unique(cell_lineage)`. Each value is an integer equal to the rounded
  sum of per-cell expected progenies (i.e., \\\sum\_{i \in \ell}
  \exp(\alpha + x_i^\top \beta)\\).

- tab_mat:

  Integer matrix with 50 rows and 2 columns `c("now","future")`;
  rownames are lineage IDs. Column `"now"` is the number of current
  cells observed in that lineage (in this subset); column `"future"`
  equals `lineage_future_count`.

## Details

Plastic simulation dataset (50 lineages)

**How this dataset was generated (summary).**

1.  Loaded example embeddings via
    `multiomeFate:::data_loader("fasttopics")`; used the
    `"fasttopic.COCL2"` cell embeddings as an early-timepoint RNA
    feature space and preprocessed them with an internal helper
    `.preprocess_rna(..., "day10_COCL2")`.

2.  Estimated a plastic direction and intercept using internal helpers:
    `.generate_simulation_plastic()` returned `coefficient_vec` and
    `coefficient_intercept`. This function also assigns cells to
    lineages. The per-cell fate potential was \\\log\_{10} \exp(\alpha +
    x^\top \beta)\\. Cells were probabilistically assigned to lineages
    via the internal functions `.compute_plastic_probabilities()` and
    `.assign_plastic_lineages()`.

3.  For packaging, lineages were ordered by decreasing future size and
    all 50 retained (an earlier release shipped 15, too few to identify
    a 30-feature fit; see
    [`cyfer`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md)).
    Cell rownames were reset to `"cell:1..n"` for clarity. The four
    objects in `plastic_simulation` were then assembled as below and
    saved with `usethis::use_data()`:


          plastic_simulation <- list(
            cell_features = embedding_mat_subset,
            cell_lineage = as.character(lineage_assignment_subset),
            lineage_future_count = lineage_future_size_subset,
            tab_mat = cbind(now = table(lineage_assignment_subset)[names(lineage_future_size_subset)],
                            future = lineage_future_size_subset)
          )
          

**Intended use.**

- Quick-start examples for
  [`cyfer`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md),
  [`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md),
  and
  [`lineage_imputation_sequence`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation_sequence.md).

- Unit tests where a small but structured dataset is helpful.

**Notes.**

- An intercept column is *not* present in `cell_features`; modeling
  functions add it automatically where required.

- The data are *fully synthetic* and derived from internal simulation
  helpers; no real subject-level outcomes are included.

## See also

[`cyfer`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer.md),
[`cyfer_finalize`](https://nancyrzhanglab.github.io/multiomeFate/reference/cyfer_finalize.md),
[`lineage_imputation_sequence`](https://nancyrzhanglab.github.io/multiomeFate/reference/lineage_imputation_sequence.md)

## Examples

``` r
data(plastic_simulation)
str(plastic_simulation)
#> List of 4
#>  $ cell_features       : num [1:7940, 1:30] 0.066 0.798 -0.134 -0.818 1.016 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7940] "cell:1" "cell:2" "cell:3" "cell:4" ...
#>   .. ..$ : chr [1:30] "fastTopicCOCL2_1" "fastTopicCOCL2_2" "fastTopicCOCL2_3" "fastTopicCOCL2_4" ...
#>  $ cell_lineage        : chr [1:7940] "lineage:26" "lineage:44" "lineage:49" "lineage:26" ...
#>  $ lineage_future_count: Named num [1:50] 532 256 232 191 187 173 160 156 143 139 ...
#>   ..- attr(*, "names")= chr [1:50] "lineage:16" "lineage:9" "lineage:10" "lineage:14" ...
#>  $ tab_mat             : num [1:50, 1:2] 159 159 159 159 159 150 159 159 159 159 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:50] "lineage:16" "lineage:9" "lineage:10" "lineage:14" ...
#>   .. ..$ : chr [1:2] "now" "future"

# \donttest{
# Minimal end-to-end example:
with(plastic_simulation, {
  set.seed(10)
  cv <- cyfer(
    cell_features = cell_features,
    cell_lineage  = cell_lineage,
    lineage_future_count = lineage_future_count,
    lambda_initial = NA,
    lambda_sequence_length = 10,
    num_folds = 5,
    verbose = 0
  )
  fit <- cyfer_finalize(
    cell_features = cell_features,
    cell_lineage  = cell_lineage,
    fit_res = cv,
    lineage_future_count = lineage_future_count
  )

  fit$lambda
  head(fit$cell_imputed_score)
  head(fit$lineage_imputed_count)
})
#>  lineage:1 lineage:10 lineage:11 lineage:12 lineage:13 lineage:14 
#>  135.92109  229.71790   93.84451  134.29372  128.30764  196.10378 
# }
```
