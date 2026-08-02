# CLAUDE.md — multiomeFate (CYFER R package)

## Workflow Instructions
1. **Always enter plan mode** before starting any non-trivial task.
2. **Use superpower skills** where relevant: `/code-review` for code changes.
3. **After every prompt**, run `/project-state`: refresh the current-state sections of the individual contributor's `CLAUDE_[name].md` in place, and append a dated entry to their `HISTORY_[name].md`. Write only non-obvious things; skip anything already in the code or git history.

## Project Context (High Level)
**Paper/Project**: `multiomeFate` is the R package implementing **CYFER** (Clonal Fate Estimation by Exponential Regression), the statistical framework from Chen, Lin, Schaff et al., *"Resolution of Selection Versus Adaptation in Cellular Evolution."* CYFER estimates a per-cell **fate potential** — the expected number of progeny a cell contributes to a designated future population — from lineage-barcoded single-cell multiomic data, by fitting a Poisson GLM with an exponential (log-linear) mean that explicitly models clonal expansion and intra-clonal heterogeneity.

**Authors** (manuscript): Xinyi E. Chen, Kevin Z. Lin, Dylan Schaff (co-first); Robert Vander Velde, Christopher Cote, Sijia Huang; Andy J. Minn, Sydney M. Shaffer\* and Nancy R. Zhang\* (co-corresponding).
**Package authors**: Kevin Z. Lin (`aut`, `cre`, kzlin@uw.edu), Xinyi Chen (`aut`).

**Submission status**: being **resubmitted to Nature Methods**. The manuscript was previously reviewed at Nature Genetics; that reviewer report and response letter sit in the paper repo's `additional_context/`, and the revision simulations in the analysis repo (`kevin/Writeup_Simulations/sim1`–`sim7`) were written to answer those critiques. Note that `paper.tex` still uses the Science-family template and `paper_nbt.tex` is a second manuscript source — neither has been reformatted for Nature Methods.

**Goal**: Distinguish **selection** (pre-existing high-potential cells preferentially expand) from **adaptation** (surviving cells change state) in cellular evolution under stress, and link each force to molecular features. The package provides the estimation machinery; the analysis repository applies it to the paper's data; the paper repository holds the manuscript.

**The name split**: the package is `multiomeFate`; the method is **CYFER**. Use CYFER in prose and documentation.

## Repository Layout

This is one of **three sibling repositories** under
`/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/`:

| Repo | Role | GitHub |
|---|---|---|
| `multiomeFate/` (**this repo**) | The R package: CYFER estimation, plotting, simulation helpers, vignette, pkgdown site | `nancyrzhanglab/multiomeFate` |
| `multiomeFate_analysis/` | All analyses reproducing the paper (branch `kevin`); `kevin/Writeup*/` directories, `kevin/Writeup_Simulations/` for the revision simulations, `csv/` exports | `nancyrzhanglab/multiomeFate_analysis` |
| `multiome_fate_paper/` | Manuscript LaTeX (`paper.tex`, `paper_nbt.tex`), `references.bib`, `fig/`, and `additional_context/` holding `Nature-genetics-review.pdf` + `Nature-genetics-response.docx` | (local / Dropbox) |

Each sibling repo has its own `CLAUDE.md`. Read the analysis repo's `CLAUDE.md` when a question concerns how CYFER is *applied*; read this one for how it is *implemented*.

**This repository's internals:**

```
R/
  lineage_cv.R                     # cyfer() — K-fold CV over the lambda path
  lineage_cv_finalize.R            # cyfer_finalize() — refit at chosen lambda
  construct_folds.R                # fold construction (internal)
  lineage_imputation.R             # single-lambda fit + optimizer internals
  lineage_imputation_sequence.R    # lambda-path fit
  lineage_imputation_simulation.R  # .construct_lineage_data() test helper
  data_loader.R                    # Seurat object loader for the paper's data
  plot_simplex.R                   # ternary plot (no ggtern)
  plot_trainTest.R                 # CV diagnostic plots
  plot_anova.R                     # ANOVA-style plots
  plot_cellGrowthUmap.R            # UMAP coloured by fate potential
  plot_lineageScatterplot.R        # lineage scatter plots
  simulation*.R                    # priming / plastic / attachFuture simulations
  barcoding_assignment.R           # CloneClean barcode assignment
  compute_entropy.R                # Shannon entropy of fate predictions
  util.R                           # sparse-matrix / numerical-stability helpers
  data.R                           # roxygen docs for bundled datasets
data/                              # priming_simulation.rda, plastic_simulation.rda
man/                               # roxygen-generated .Rd
tests/testthat/                    # devtools::test()
vignettes/simulation.rmd           # the CYFER walkthrough vignette
docs/                              # pkgdown output (gitignored; built by .github/workflows/pkgdown.yaml)
```

## Who Is Using This Session?
**Detect the current user** by running: `echo $USER`. This table maps each login to that person's **first-name** context file; it is the source of truth that `/brainstorm` and `/project-state` use to resolve `$USER` to the right filename (so the login `kevinlin` maps to `CLAUDE_kevin.md`, never `CLAUDE_kevinlin.md`).

| Username (login) | Current-state file (first name) | History archive |
|---|---|---|
| `kevinlin` | `CLAUDE_kevin.md` | `HISTORY_kevin.md` |

**File ownership.** Each row above names one person's files, and **only that person's session writes them.** Once `$USER` resolves to a first name, that is the only suffix you may create, edit, append to, rename, or delete — every other collaborator's `CLAUDE_[name].md`, `HISTORY_[name].md`, and `brainstorming_[name].md` is read-only. Read them for context when useful; never modify them, not even to fix a typo or add a cross-reference. This repository is shared, so an edit lands in the owner's working copy immediately and can overwrite state they wrote from their own machine. If a collaborator's file looks wrong or stale, say so instead of editing it. This master `CLAUDE.md` is the exception: it is shared and any collaborator may update it. The single override is the user, in the current turn, *directing you to write* that exact file ("add this to `CLAUDE_sarah.md`") — confirm once, then write it. Merely mentioning a collaborator, referring to their file, or being away does not qualify, and attribution does not launder the edit — a change stamped with your name and confined to two lines is still a write to a file you do not own. Record the decision here in the master `CLAUDE.md` and tell the user to contact the owner directly.

After detecting the user, **read that person's `CLAUDE_[name].md` immediately** before doing any other work — it holds the current project state and restores context in ~30 seconds. Do **not** read `HISTORY_[name].md` at startup; it is the append-only session log, consulted only on demand when deep history is needed. If no match is found, ask who the user is and create a new `CLAUDE_[name].md` using `/project-state`.

## Post-Prompt Update Instructions
After completing each user prompt, run `/project-state`. It will:
- **Refresh in place** the current-state sections of `CLAUDE_[name].md` (Project Status, Key Methodological Details, Open Questions / Next Steps).
- **Append a dated entry at the bottom** of `HISTORY_[name].md` recording new decisions, resolved/open questions, non-obvious code or LaTeX rationale, and empirical findings.

Do NOT record: things already in the LaTeX/code, git history, or reproducible from code.

---

# Package Reference

## The CYFER model

For cell $i$ sampled at time $t_1$ with feature vector $X_{i,\cdot} \in \mathbb{R}^d$, fate potential is linear in the features:

$$Z_i = \beta_0 + X_{i,\cdot}^\top \beta.$$

The observed clone size at $t_2$ is the sum of independent Poisson draws over the clone's $t_1$ cells:

$$y_\ell = \sum_{i \in \text{clone } \ell} N_i, \qquad N_i \sim \text{Poisson}\!\left(e^{\beta_0 + X_{i,\cdot}^\top \beta}\right).$$

Fitting minimizes the penalized negative log-likelihood

$$\sum_{\ell=1}^{L}\left[\Big(\sum_{i \in \ell} e^{\beta_0 + X_{i,\cdot}^\top \beta}\Big) - y_\ell \log\Big(\sum_{i \in \ell} e^{\beta_0 + X_{i,\cdot}^\top \beta}\Big)\right] + \lambda\|\beta\|_2^2$$

via `stats::optim` (BFGS) with an analytical gradient and $M$ random restarts (default 10), because the objective is non-convex.

**Crucially, features are needed only at $t_1$** — the later time point contributes only per-clone counts. That is what makes the method work across a selection bottleneck.

**Typical features**: column-wise concatenation of each cell's fastTopics embedding (RNA) and PeakVI embedding (ATAC), not raw genes/peaks — the low-dimensional representation is far more stable. To interpret, correlate a gene/peak against the estimated fate potential across cells; $\hat\beta$ itself is not directly interpretable.

**Downstream indices** (computed in the analysis repo, not here):
- `fate_potential` = log10 expected future progeny per cell (`cell_imputed_score`).
- `selection_index` = Gini coefficient of estimated fate potentials.
- `adaptation_index` = fate-potential-weighted average embedding shift $t_1 \to t_2$, normalized by a "no adaptation" distance.

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

### Fold construction and lambda selection

`construct_folds()` sorts clones in descending order of `lineage_future_count` and deals them round-robin into folds, so no fold dominates the others in cell number. `cyfer_finalize()` picks lambda by minimizing the **median across folds** of `test_loglik` at each lambda (lower is better; this is the unpenalized objective value, not a literal log-likelihood), then refits once on all clones at that lambda.

### Preprocessing conventions used by callers

These live in the analysis scripts, not in the package, but every caller repeats them:

- Filter to clones with `lineage_future_count > 0`.
- Filter to clones with `>= 2` cells for stable CV.
- `num_folds = min(3, length(unique(cell_lineage)) - 1)` for small simulations; the paper used `k = 20` folds on real data.
- `scale()` the feature matrix before fitting.

## Design decisions and history

### `tab_mat` removed (2025 revision)
`cyfer()` previously required a `tab_mat` matrix (rows = lineages, columns = timepoints) and a `future_timepoint` string to identify which column held future counts. This was redundant because `lineage_future_count` already encodes that information. Both arguments were dropped; `construct_folds()` now orders lineages by `lineage_future_count` directly. The stored datasets (`priming_simulation`, `plastic_simulation`) still contain a `tab_mat` element for historical reference but it is not used by any current function.

### `lineage_cv` → `cyfer` rename (2025 revision)
`lineage_cv()` and `lineage_cv_finalize()` were renamed to `cyfer()` and `cyfer_finalize()` to match the method name used in the paper. The S3 class changed from `"lineage_cv"` to `"cyfer"`. The **file** names (`R/lineage_cv.R`, `R/lineage_cv_finalize.R`) were deliberately not renamed, so grepping for `cyfer` will not find them by filename.

### `data_loader` DimReduc assay fix
Seurat v5 errors when a `DimReduc` object references an assay (via `@assay.used`) that is not present in the object. The fastTopics, peakVI, and other DimReduc files were computed with `assay.used = "RNA"`, so loading them without loading the RNA assay triggered `"Cannot find assay 'RNA'"`. The internal helper `.fix_dimreduc_assay()` updates `@assay.used` to the first available non-Empty assay before adding any DimReduc to the Seurat object.

### `plot_simplex` without ggtern
`ggtern` overrides core ggplot2 methods when loaded, breaking other plots in the same session. `plot_simplex()` is implemented using a manual ternary-to-Cartesian coordinate transform and plain `ggplot2` — no ggtern dependency. The API uses `x_col`, `y_col`, `z_col` character arguments (column names in `df`) instead of ggtern's `aes()` formula.

## Documentation conventions

- Use **"time point"** (two words) in all narrative documentation and prose.
- Parameter *names* like `future_timepoint` stay as-is (they are code identifiers, not prose).
- Do not add an intercept column to `cell_features` in examples or docs; functions add it internally.
- New R files drafted by Claude get a `_claude` suffix so the human can review before integrating.

## Testing

Tests live in `tests/testthat/`. Run with `devtools::test()`.

- The internal helper `.construct_lineage_data()` (in `R/lineage_imputation_simulation.R`) returns `cell_features` **with** an Intercept column already prepended (it calls `.lineage_cleanup()` internally). Strip it before passing to user-facing functions: `cell_features[, setdiff(colnames(cell_features), "Intercept"), drop=FALSE]`.
- `numDeriv` is used in gradient tests to verify the analytical gradient against numerical differentiation.

## Build and release

- Roxygen 7.3.3; regenerate `man/` and `NAMESPACE` with `devtools::document()`.
- pkgdown site at https://nancyrzhanglab.github.io/multiomeFate/, built by `.github/workflows/pkgdown.yaml`. `docs/` is gitignored.
- `.Rbuildignore` excludes `CLAUDE*.md`, `HISTORY*.md`, `brainstorming*.md`, `additional_context/`, `.claude/`, `.githooks/`, `.github/`, `docs/`, and `_pkgdown.yml` so `R CMD check` stays clean.
- Developed and tested primarily on R 4.3.2 (macOS).

## Git hooks
`.githooks/pre-commit` blocks any staged file ≥ 50 MB. Git does not run tracked hooks automatically — each collaborator runs `git config core.hooksPath .githooks` once per clone. See `.githooks/README.md`.

## Additional context
`additional_context/summary.md` indexes reference material for this project. The Nature Genetics reviewer report and response letter live in the **paper** repo's `additional_context/`, not here.
