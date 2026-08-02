# HISTORY_kevin.md — Kevin's Session Log

> **Append-only, ascending chronological order** (oldest at top, newest at the bottom). Add each session's dated entry to the END of this file. Never read at session startup — consulted only on demand for deep history. Current project state lives in `CLAUDE_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session appends to this file. Other collaborators may read it but must not edit or rewrite entries.

---

### [2026-08-01] (Session 1 — project-setup / project-state scaffolding for the package repo)
- Read all three sibling repos for context: this package, `multiomeFate_analysis` (branch `kevin`), and `multiome_fate_paper` (`paper.tex`, Methods §"Overview of CYFER" / "Estimation procedure for CYFER").
- Applied `/project-setup`: rewrote the master `CLAUDE.md` on the standard template (workflow instructions, "Who Is Using This Session?" table, file-ownership rules, post-prompt update instructions), merged the `.gitignore` template into the existing R-only ignore list, added `.githooks/` with the ≥ 50 MB pre-commit guard, and enabled it via `git config core.hooksPath .githooks`.
- Applied `/project-state`: created `CLAUDE_kevin.md` and this file from the templates.
- **Resolved: CYFER acronym was wrong in the old package `CLAUDE.md`.** It read "Cell Fate via Exponential Regression"; `paper.tex` line 176 defines it as **Clonal Fate Estimation by Exponential Regression**. Corrected in the master `CLAUDE.md`.
- **Resolved: paper title.** The live title is "Resolution of Selection Versus Adaptation in Cellular Evolution"; the older "Temporal and Clonal Resolution of Cellular Evolution Under Stress" is commented out in `paper.tex` and is only still asserted by this package's `README.md`.
- Preserved the package-specific knowledge that was already in the old `CLAUDE.md` (input conventions, class names, the `tab_mat` removal, the `lineage_cv`→`cyfer` rename, the `data_loader` DimReduc/`assay.used` fix, the deliberate no-ggtern `plot_simplex` implementation, the `.construct_lineage_data()` Intercept gotcha) rather than replacing it — the template only covers workflow scaffolding, and that material is the actual institutional memory of the package.
- Added a "The CYFER model" section to the master `CLAUDE.md` transcribing the model, objective, and estimation procedure from the paper's Methods, so package sessions do not have to re-derive them from `paper.tex`.
- Extended `.Rbuildignore` to cover `CLAUDE*.md`, `HISTORY*.md`, `brainstorming*.md`, `additional_context/`, and `.githooks/` — otherwise `R CMD check` flags the new scaffolding as non-standard top-level files.
- Created `additional_context/summary.md` as a stub index pointing at the paper repo's `additional_context/` (which holds `Nature-genetics-review.pdf` and `Nature-genetics-response.docx`); no reference PDFs are duplicated into the package repo.
- **Open: target journal.** Reviewer materials are named "Nature-genetics-*", but `paper.tex` uses the Science-family template and its live abstract is commented "for submission to Nature" (200-word limit), with a separate `paper_nbt.tex` also present. Did not attempt to resolve.
- **Open: stale `README.md`.** Names the old title and still lists `ggtern` as a dependency, which `DESCRIPTION` does not and `plot_simplex()` deliberately avoids. Left unedited — outside the scope of this session's request.
- **Open: six tracked `.DS_Store` files** (root, `R/`, `data/`, `man/`, `tests/`, `tests/testthat/`). Now matched by `.gitignore`, but gitignore does not untrack; they need `git rm --cached`. Not run — it touches the index, which already had staged changes.
- **Open: analysis repo `CLAUDE.md` absolute paths are stale** — they omit the `archive/` segment the tree acquired when it moved under `Collaboration-and-People/archive/Nancy/`. That file is a shared master, so it is editable, but it belongs to a different repo and was left alone this session.

### [2026-08-01] (Session 2 — README fix, .DS_Store purge, Nature Methods retarget)
- **Decision: the paper is being resubmitted to Nature Methods** (previously reviewed at Nature Genetics). This resolves Session 1's open question about the target journal. Recorded in this package's master `CLAUDE.md`, in the analysis repo's `CLAUDE.md`, and in `CLAUDE_kevin.md`.
- Fixed `README.md`: paper title updated to the live one, CYFER spelled out on first use, and the dependency list corrected against `DESCRIPTION` (dropped `ggtern`, added `MASS`/`Matrix`/`rlang`). Also repaired a dangling pointer — the old text sent readers to "the last section of this README" for package sources, but no such section exists; the file ends with a `session_info()` dump. All Imports are on CRAN, so that is now what the sentence says. The `session_info()` block itself was left verbatim: it is a record of a tested environment (including the `ggtern` that was loaded at the time), not a claim about current dependencies.
- **Purged `.DS_Store` everywhere.** 6 tracked in the package repo and 20 in the analysis repo were `git rm --cached`'d, and every on-disk copy across all three repos was deleted. The analysis repo had *already* ignored `.DS_Store` — the files predated the ignore rule, which is exactly the case gitignore cannot fix on its own.
- Added a `.gitignore` to the paper repo, which had none at all (OS cruft + LaTeX build artifacts, keeping `.tex`/`.bib`/`fig/` tracked). Without it `.DS_Store` would immediately reappear as untracked noise there.
- Fixed the analysis repo's three stale absolute paths by inserting the `archive/` segment; verified `papers/` and `out/` both exist at the corrected locations before editing.
- Also corrected the CYFER acronym in the analysis repo's `CLAUDE.md` ("Cell Fate via Exponential Regression" → "Clonal Fate Estimation by Exponential Regression"), the same error found in the package's file last session — it had been copied between the two.
- **Open: neither `paper.tex` (Science-family template) nor `paper_nbt.tex` is formatted for Nature Methods.** No LaTeX was touched this session; reformatting was not requested.
- **Open: `sim1`–`sim7` were written to answer Nature Genetics reviewers.** Worth deciding which still belong in a fresh Nature Methods submission, and whether the response letter should become cover-letter material instead of a point-by-point reply.
- Nothing committed in any of the three repos.
