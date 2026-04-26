# multiomeFate — CLAUDE.md

## What this package is

`multiomeFate` implements **CYFER** (Cell Fate via Exponential Regression), the statistical framework from Chen, Lin et al. "Resolution of Selection Versus Adaptation in Cellular Evolution." CYFER estimates per-cell fate potentials from lineage-traced single-cell multiomic data by fitting an exponential Poisson GLM that explicitly models clonal expansion.

The package name is `multiomeFate`; the method name used in the paper and documentation is **CYFER**.

---

## Key API

| Function | Role |
|---|---|
| `cyfer()` | K-fold cross-validation over a lambda path to select regularization |
| `cyfer_finalize()` | Refit at the chosen lambda, return per-cell imputed scores |
| `lineage_imputation_sequence()` | Fit the GLM along a full lambda path |
| `lineage_imputation()` | Fit the GLM at a single lambda |
| `data_loader()` | Load the paper's preprocessed Seurat objects |
| `plot_simplex()` | Ternary/simplex plot (pure ggplot2, no ggtern) |

### Input conventions for `cyfer` / `cyfer_finalize`

- `cell_features`: numeric matrix, **rows = cells, columns = features**. Must have both row names (cell IDs) and column names (feature names). Do **not** include an intercept column — functions add it internally via `.lineage_cleanup()`.
- `cell_lineage`: character vector mapping each cell to its lineage.
- `lineage_future_count`: named numeric vector — names are lineage IDs, values are cell counts at the future time point.
- No `tab_mat` or `future_timepoint` arguments. These were removed; fold stratification now uses `lineage_future_count` values directly.

### Class names

- `cyfer()` returns an object of class `"cyfer"` (a list of per-fold results).
- `lineage_imputation()` returns class `"lineage_imputation"`.
- `plot_trainTest()` expects a `"cyfer"`-class object.

---

## Design decisions and history

### `tab_mat` removed (2025 revision)
`cyfer()` previously required a `tab_mat` matrix (rows = lineages, columns = timepoints) and a `future_timepoint` string to identify which column held future counts. This was redundant because `lineage_future_count` already encodes that information. Both arguments were dropped; `construct_folds()` now orders lineages by `lineage_future_count` directly. The stored datasets (`priming_simulation`, `plastic_simulation`) still contain a `tab_mat` element for historical reference but it is not used by any current function.

### `lineage_cv` → `cyfer` rename (2025 revision)
`lineage_cv()` and `lineage_cv_finalize()` were renamed to `cyfer()` and `cyfer_finalize()` to match the method name used in the paper. The S3 class changed from `"lineage_cv"` to `"cyfer"`.

### `data_loader` DimReduc assay fix
Seurat v5 errors when a `DimReduc` object references an assay (via `@assay.used`) that is not present in the object. The fastTopics, peakVI, and other DimReduc files were computed with `assay.used = "RNA"`, so loading them without loading the RNA assay triggered `"Cannot find assay 'RNA'"`. The internal helper `.fix_dimreduc_assay()` updates `@assay.used` to the first available non-Empty assay before adding any DimReduc to the Seurat object.

### `plot_simplex` without ggtern
`ggtern` overrides core ggplot2 methods when loaded, breaking other plots in the same session. `plot_simplex()` is implemented using a manual ternary-to-Cartesian coordinate transform and plain `ggplot2` — no ggtern dependency. The API uses `x_col`, `y_col`, `z_col` character arguments (column names in `df`) instead of ggtern's `aes()` formula.

---

## Documentation conventions

- Use **"time point"** (two words) in all narrative documentation and prose.
- Parameter *names* like `future_timepoint` stay as-is (they are code identifiers, not prose).
- Do not add an intercept column to `cell_features` in examples or docs; functions add it internally.

---

## Testing

Tests live in `tests/testthat/`. Run with `devtools::test()`.

- The internal helper `.construct_lineage_data()` (in `R/lineage_imputation_simulation.R`) returns `cell_features` **with** an Intercept column already prepended (it calls `.lineage_cleanup()` internally). Strip it before passing to user-facing functions: `cell_features[, setdiff(colnames(cell_features), "Intercept"), drop=FALSE]`.
- `numDeriv` is used in gradient tests to verify the analytical gradient against numerical differentiation.

---

## Package structure

```
R/
  lineage_cv.R              # cyfer() — main CV function
  lineage_cv_finalize.R     # cyfer_finalize()
  construct_folds.R         # fold construction (internal)
  lineage_imputation.R      # single-lambda fit + internals
  lineage_imputation_sequence.R  # lambda path fit
  lineage_imputation_simulation.R  # .construct_lineage_data() test helper
  data_loader.R             # Seurat object loader for paper data
  plot_simplex.R            # ternary plot (no ggtern)
  plot_trainTest.R          # CV diagnostic plots
  plot_anova.R              # ANOVA-style plots
  plot_cellGrowthUmap.R     # UMAP coloured by fate potential
  plot_lineageScatterplot.R # lineage scatter plots
  simulation*.R             # priming / plastic / attachFuture simulations
  barcoding_assignment.R    # CloneClean barcode assignment
  compute_entropy.R         # Shannon entropy of fate predictions
  util.R                    # sparse matrix / numerical stability helpers
  data.R                    # roxygen docs for bundled datasets
```
