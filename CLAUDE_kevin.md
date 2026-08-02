# CLAUDE_kevin.md — Kevin's Context

> **Current state only.** Every section below is updated *in place* each session — overwrite, don't append. The append-only dated log lives in `HISTORY_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session writes this file and `HISTORY_kevin.md`. Other collaborators may read it for context but must not edit it.

## About Kevin
- Role in project: co-first author; statistical methodology and software — author and maintainer (`cre`) of the `multiomeFate` R package
- Background: statistics / biostatistics, Department of Biostatistics, University of Washington
- Email: kzlin@uw.edu

## Environment (local paths)
All three repos are siblings under
`/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/`:

- Package root (this repo): `multiomeFate/` — branch `master`, remote `nancyrzhanglab/multiomeFate`
- Analysis root: `multiomeFate_analysis/` — branch `kevin`, remote `nancyrzhanglab/multiomeFate_analysis`
- LaTeX root: `multiome_fate_paper/` — `paper.tex` (Science-family template) and `paper_nbt.tex`
- Simulation outputs: `/Users/kevinlin/.../archive/Nancy/multiomeFate/out/Writeup_Simulations/` (RDS), exported to `multiomeFate_analysis/csv/kevin/Writeup_Simulations/`
- Paper PDF and reviewer notes: `.../archive/Nancy/multiomeFate/papers/`

## Project Status (as of 2026-08-01)

The paper is **Chen, Lin, Schaff et al., "Resolution of Selection Versus Adaptation in Cellular Evolution"**, now being **resubmitted to Nature Methods**. It was previously reviewed at Nature Genetics; that reviewer report (`Nature-genetics-review.pdf`) and response letter (`Nature-genetics-response.docx`) are in the paper repo's `additional_context/`, and `sim1`–`sim7` in the analysis repo were written against those critiques.

**Package (this repo):** stable at version 1.0.1.000. Six exported functions (`cyfer`, `cyfer_finalize`, `lineage_imputation`, `lineage_imputation_sequence`, `data_loader`, `plot_simplex`); the 2025 revision renamed `lineage_cv*` → `cyfer*` and dropped the `tab_mat` / `future_timepoint` arguments. pkgdown site builds from `.github/workflows/pkgdown.yaml`; one vignette (`vignettes/simulation.rmd`) walks through the API. Eight testthat files. Working tree currently holds uncommitted scaffolding, README fixes, and `.DS_Store` deletions.

**Analysis repo:** the active work is `kevin/Writeup_Simulations/` — seven revision simulation scripts (`sim1`–`sim7`) answering specific reviewer critiques (growth-model misspecification, rare resistance, heritability grid, power analysis, barcode dropout, adaptation-index bias under cell death, sensitivity/specificity). Each writes an RDS with its own ad-hoc schema; `make_csvs.R` flattens them to CSV.

**Project-file scaffolding:** as of this session the package repo carries the standard three-file pattern (master `CLAUDE.md`, `CLAUDE_kevin.md`, `HISTORY_kevin.md`) plus `.gitignore`, `.githooks/`, and `additional_context/`.

## Key Methodological Details

- **CYFER = Clonal Fate Estimation by Exponential Regression** (per `paper.tex`). The package's older `CLAUDE.md` expanded it as "Cell Fate via Exponential Regression" — that was wrong and is now corrected.
- **Current title** is "Resolution of Selection Versus Adaptation in Cellular Evolution." The earlier title, "Temporal and Clonal Resolution of Cellular Evolution Under Stress," survives only as a comment in `paper.tex`.
- **Features are only needed at $t_1$.** The future time point contributes nothing but per-clone counts. This is the property that lets CYFER work through a selection bottleneck and is the main differentiator from trajectory/optimal-transport methods.
- **Feature matrix in practice** = fastTopics (RNA) embedding concatenated with PeakVI (ATAC) embedding, not raw genes/peaks. $\hat\beta$ is not interpretable at the dimension level; interpretation is done by correlating a gene or peak against estimated fate potential across cells.
- **Objective is non-convex** → 10 random restarts by default, with a data-driven initialization range derived from `max_ratio` (largest $t_2/t_1$ clone-size ratio) and the 95% quantile of `|X|`.
- **Folds** are dealt round-robin after sorting clones by descending future count, so no fold dominates in cell number. Paper used $k=20$ on real data; simulations use `min(3, n_clones - 1)`.
- **Lambda choice** in `cyfer_finalize()` minimizes the *median* across folds of the held-out unpenalized objective — median, not mean, for robustness to a fold blowing up.
- **`cell_features` must not carry an Intercept column** — `.lineage_cleanup()` adds it. But the test helper `.construct_lineage_data()` returns one *with* Intercept already prepended, so tests must strip it.
- **`plot_simplex()` avoids `ggtern` deliberately** — loading ggtern overrides ggplot2 methods session-wide and breaks every other plot.
- **File names lag function names**: `cyfer()` lives in `R/lineage_cv.R`, `cyfer_finalize()` in `R/lineage_cv_finalize.R`.

## Open Questions / Next Steps

1. **Open: neither manuscript source is formatted for Nature Methods.** `paper.tex` uses the Science-family template (with its live abstract commented "for submission to Nature", 200 words); `paper_nbt.tex` is a second source. Decide which becomes the Nature Methods submission, then reformat — Nature Methods wants a ~150-word editorial summary and abstract, Methods-style article structure, and its own reference style.
2. **Open: `sim1`–`sim7` were written against the Nature Genetics reviewer report.** Check which of those analyses still belong in the Nature Methods submission, and whether the response letter needs to be reframed as new-submission cover material rather than a point-by-point reply.
3. **Open: the three commits across the repos are uncommitted.** Package repo and analysis repo both carry staged `.DS_Store` deletions; the paper repo has a new untracked `.gitignore`.
4. **Open: no `brainstorming_kevin.md` yet.** Run `/brainstorm` if research directions need capturing.
5. **Open: the revision simulations have no uniform RDS schema.** If another review round arrives, consider unifying before adding an eighth script.
