context("Testing the plotting functions")

## Plot tests do not compare rendered images. They assert the data the plot is
## built from, and that construction does not error -- which is enough to catch
## a ggplot2 deprecation, a broken coordinate transform, or an argument that is
## silently ignored.

## Fixtures -------------------------------------------------------------------

.plot_cyfer_fit <- function(num_folds = 3,
                            seed_number = 10){
  set.seed(seed_number)
  res <- .construct_lineage_data()
  cell_features <- res$cell_features[
    , setdiff(colnames(res$cell_features), "Intercept"), drop = FALSE]

  list(cell_features = cell_features,
       cell_lineage = res$cell_lineage,
       cv_fit_list = cyfer(cell_features = cell_features,
                           cell_lineage = res$cell_lineage,
                           lineage_future_count = res$lineage_future_count,
                           lambda_initial = NA,
                           lambda_sequence_length = 10,
                           num_folds = num_folds,
                           seed_number = 1,
                           verbose = 0),
       lineage_future_count = res$lineage_future_count)
}

## A minimal Seurat object -- counts, a lineage column and a time column are all
## `plot_anova()` reads. Deliberately not a realistic object; the point is to
## exercise the selection and labelling arithmetic without a fixture file.
## Note the default of 12 lineages rather than a handful: `plot_anova()` errors
## outright when fewer lineages survive `min_lineage_size` than
## `num_lineages_bottom` (default 10) -- see defect D10 below.
.plot_seurat_object <- function(num_cells = 120,
                                num_lineages = 12,
                                seed_number = 10){
  set.seed(seed_number)
  count_mat <- matrix(stats::rpois(50 * num_cells, lambda = 3),
                      nrow = 50,
                      dimnames = list(paste0("gene", seq_len(50)),
                                      paste0("cell", seq_len(num_cells))))
  seurat_object <- suppressWarnings(
    Seurat::CreateSeuratObject(counts = count_mat))
  seurat_object$assigned_lineage <- rep(paste0("L", seq_len(num_lineages)),
                                        length.out = num_cells)
  seurat_object$time_celltype <- rep(c("day0", "day7"),
                                     each = num_cells / 2)

  seurat_object
}

##############################################################################
## PL1 / PL2 -- smoke tests and the class contract
##############################################################################

## PL1: there is not even a smoke test today, so a ggplot2 deprecation can break
## every figure in the analysis repo and stay silent until someone re-runs a
## script by hand.
test_that("every plotting function returns a ggplot object (PL1)", {
  simplex_df <- data.frame(a = c(1, 0, 0, 1), b = c(0, 1, 0, 1),
                           c = c(0, 0, 1, 1))
  expect_true(inherits(plot_simplex(df = simplex_df, x_col = "a",
                                    y_col = "b", z_col = "c"), "ggplot"))

  seurat_object <- .plot_seurat_object()
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  expect_true(inherits(plot_anova(
    seurat_object = seurat_object,
    cell_imputed_score = score_vec,
    assigned_lineage_variable = "assigned_lineage",
    time_celltype_variable = "time_celltype",
    day_later = "day7"), "ggplot"))
})

test_that("plot_trainTest and plot_lineageScatterplot build a plot (PL1)", {
  fit_list <- .plot_cyfer_fit()
  expect_true(inherits(plot_trainTest(fit_list$cv_fit_list), "ggplot"))

  count_vec <- stats::setNames(c(5, 40, 100, 3, 60), paste0("L", seq_len(5)))
  imputed_vec <- stats::setNames(c(4, 52, 88, 6, 45), paste0("L", seq_len(5)))
  set.seed(10)
  expect_true(inherits(plot_lineageScatterplot(
    lineage_future_count = count_vec,
    lineage_imputed_count = imputed_vec), "ggplot"))
})

## PL2: pins the class contract the lineage_cv -> cyfer rename established.
test_that("plot_trainTest requires a cyfer object (PL2)", {
  fit_list <- .plot_cyfer_fit()

  expect_error(plot_trainTest(unclass(fit_list$cv_fit_list)))
  expect_error(plot_trainTest(list(a = 1)))

  wrong_class_list <- fit_list$cv_fit_list
  class(wrong_class_list) <- "lineage_cv"
  expect_error(plot_trainTest(wrong_class_list))
})

##############################################################################
## PL3 / PL4 / PL5 -- the lambda the plot marks
##############################################################################

## PL3
test_that(".prepare_trainTest_data returns the documented shape (PL3)", {
  fit_list <- .plot_cyfer_fit()
  lambda_sequence <- fit_list$cv_fit_list[[1]]$train_fit$lambda_sequence

  for(what in c("train", "test")){
    label <- paste0("what = ", what)
    res <- .prepare_trainTest_data(fit_list$cv_fit_list,
                                   quantile_vec = c(0.1, 0.5, 0.9),
                                   what = what)

    expect_true(nrow(res$df) == 3 * length(lambda_sequence), info = label)
    expect_true(all(sort(unique(res$df$quantile_str)) ==
                      c("lower", "median", "upper")), info = label)
    expect_true(all(table(res$df$quantile_str) == length(lambda_sequence)),
                info = label)
    # The `lambda` column is already shifted by +1 for the log axis.
    expect_true(max(abs(sort(unique(res$df$lambda)) -
                          sort(lambda_sequence + 1))) <= 1e-10, info = label)
    expect_true(res$lambda %in% lambda_sequence, info = label)
  }

  expect_error(.prepare_trainTest_data(fit_list$cv_fit_list,
                                       quantile_vec = c(0.1, 0.5),
                                       what = "test"))
  expect_error(.prepare_trainTest_data(fit_list$cv_fit_list,
                                       quantile_vec = c(0.5, 0.1, 0.9),
                                       what = "test"))
  expect_error(.prepare_trainTest_data(fit_list$cv_fit_list,
                                       quantile_vec = c(0.1, 0.5, 0.9),
                                       what = "validation"))
})

## PL4: two independent code paths compute the selected lambda -- the plot's
## dashed line and `cyfer_finalize()`'s refit point -- and nothing checks they
## agree. At the default `quantile_vec` they must.
test_that("the marked lambda equals the one cyfer_finalize picks (PL4)", {
  for(num_folds in c(3, 4)){
    label <- paste0("num_folds = ", num_folds)
    fit_list <- .plot_cyfer_fit(num_folds = num_folds)

    plot_lambda <- .prepare_trainTest_data(fit_list$cv_fit_list,
                                           quantile_vec = c(0.1, 0.5, 0.9),
                                           what = "test")$lambda
    final_res <- cyfer_finalize(
      cell_features = fit_list$cell_features,
      cell_lineage = fit_list$cell_lineage,
      fit_res = fit_list$cv_fit_list,
      lineage_future_count = fit_list$lineage_future_count)

    expect_true(abs(plot_lambda - final_res$lambda) <= 1e-10, info = label)
  }
})

## PL5
test_that("plot_trainTest always marks the median lambda (PL5)", {
  fit_list <- .plot_cyfer_fit()
  expect_error(plot_trainTest(fit_list$cv_fit_list,
                              quantile_vec = c(0.1, 0.25, 0.9)))
})

## The rejection happens in `.prepare_trainTest_data()`, so it fires whichever
## curve is being summarized, and the outer two quantiles stay free.
test_that("only the middle quantile is pinned to 0.5 (PL5)", {
  fit_list <- .plot_cyfer_fit()

  for(what in c("train", "test")){
    expect_error(.prepare_trainTest_data(fit_list$cv_fit_list,
                                         quantile_vec = c(0.1, 0.25, 0.9),
                                         what = what),
                 "must be 0.5", info = what)
  }

  # Moving the band edges is fine and does not move the marked lambda.
  narrow_lambda <- .prepare_trainTest_data(fit_list$cv_fit_list,
                                           quantile_vec = c(0.25, 0.5, 0.75),
                                           what = "test")$lambda
  wide_lambda <- .prepare_trainTest_data(fit_list$cv_fit_list,
                                         quantile_vec = c(0.05, 0.5, 0.95),
                                         what = "test")$lambda
  expect_true(abs(narrow_lambda - wide_lambda) <= 1e-10)
  expect_true(inherits(plot_trainTest(fit_list$cv_fit_list,
                                      quantile_vec = c(0.25, 0.5, 0.75)),
                       "ggplot"))
})

##############################################################################
## PL6 / PL7 -- .anova_percentage()
##############################################################################

## PL6
test_that(".anova_percentage matches a hand computation (PL6)", {
  df <- data.frame(lineage = factor(c("a", "a", "b", "b")),
                   value = c(1, 3, 6, 10))
  # grand mean 5; group means 2 and 8.
  # across = 2*(2-5)^2 + 2*(8-5)^2 = 36
  # total  = 16 + 4 + 1 + 25 = 46, so the answer is 78.3.
  expect_true(abs(.anova_percentage(df = df,
                                    lineage_variable = "lineage",
                                    value_variable = "value") - 78.3) <= 1e-10)
  expect_true(abs(.anova_percentage(df = df,
                                    lineage_variable = "lineage",
                                    value_variable = "value") -
                    round(36 / 46 * 100, 1)) <= 1e-10)

  # Identical group means: no across-lineage variance at all.
  df_flat <- data.frame(lineage = factor(c("a", "a", "b", "b")),
                        value = c(1, 3, 1, 3))
  expect_true(.anova_percentage(df = df_flat,
                                lineage_variable = "lineage",
                                value_variable = "value") == 0)

  # Zero within-lineage variance: everything is across-lineage.
  df_tight <- data.frame(lineage = factor(c("a", "a", "b", "b")),
                         value = c(2, 2, 8, 8))
  expect_true(.anova_percentage(df = df_tight,
                                lineage_variable = "lineage",
                                value_variable = "value") == 100)

  expect_error(.anova_percentage(
    df = data.frame(lineage = c("a", "a", "b", "b"), value = c(1, 3, 6, 10)),
    lineage_variable = "lineage",
    value_variable = "value"))
})

## PL7: `.anova_percentage()` iterates `levels()`, not observed values, so an
## unused level contributes `mean(numeric(0))` = NaN. `.plot_anova_helper()`
## calls `droplevels()` first, which is the only thing standing between a stray
## level and a NaN in a figure title.
test_that("unused factor levels give NaN, and droplevels prevents it (PL7)", {
  df_unused <- data.frame(
    lineage = factor(c("a", "a", "b", "b"), levels = c("a", "b", "ghost")),
    value = c(1, 3, 6, 10))

  expect_true(is.nan(.anova_percentage(df = df_unused,
                                       lineage_variable = "lineage",
                                       value_variable = "value")))

  df_dropped <- df_unused
  df_dropped$lineage <- droplevels(df_dropped$lineage)
  expect_true(!is.nan(.anova_percentage(df = df_dropped,
                                        lineage_variable = "lineage",
                                        value_variable = "value")))

  # And the helper does apply droplevels(), so a stray level cannot reach it.
  seurat_object <- .plot_seurat_object()
  seurat_object$assigned_lineage <- factor(
    seurat_object$assigned_lineage,
    levels = c(paste0("L", seq_len(12)), "ghost"))
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = score_vec,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = "day7")
  expect_true(!grepl("NaN", plot1$labels$title))
})

##############################################################################
## PL8 / PL9 / PL10 -- plot_anova() selection and labelling
##############################################################################

## PL8
test_that("bool_add_future_size = FALSE drops the ' (n)' suffix (PL8)", {
  seurat_object <- .plot_seurat_object()
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = score_vec,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = "day7",
                      bool_add_future_size = FALSE)

  expect_true(!any(grepl("\\(", unique(as.character(plot1$data$lineage)))))
})

## The helper honours the argument correctly; it is only the forwarding that is
## missing. Pinning it here means the D2 fix is a one-line change with a test
## already in place on both sides.
test_that(".plot_anova_helper does honour bool_add_future_size (D2)", {
  seurat_object <- .plot_seurat_object()
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  size_vec <- stats::setNames(rep(10, 12), paste0("L", seq_len(12)))

  plot_with <- .plot_anova_helper(
    seurat_object = seurat_object,
    cell_imputed_score = score_vec,
    assigned_lineage_variable = "assigned_lineage",
    lineage_future_size = size_vec,
    bool_add_future_size = TRUE)
  plot_without <- .plot_anova_helper(
    seurat_object = seurat_object,
    cell_imputed_score = score_vec,
    assigned_lineage_variable = "assigned_lineage",
    lineage_future_size = size_vec,
    bool_add_future_size = FALSE)

  expect_true(any(grepl(" \\(10\\)",
                        unique(as.character(plot_with$data$lineage)))))
  expect_true(!any(grepl(" \\(10\\)",
                         unique(as.character(plot_without$data$lineage)))))
})

## PL9: `unique(c(top, bottom))` is what keeps an overlapping lineage from being
## drawn twice. Off-by-one country -- with 4 eligible lineages and 3 requested
## from each end, two lineages are in both sets.
test_that("plot_anova draws an overlapping lineage once (PL9)", {
  seurat_object <- .plot_seurat_object(num_lineages = 4)
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))

  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = score_vec,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = "day7",
                      num_lineages_top = 3,
                      num_lineages_bottom = 3)

  lineage_vec <- unique(as.character(plot1$data$lineage))
  lineage_vec <- setdiff(lineage_vec, "All")
  expect_true(length(lineage_vec) == 4)
  expect_true(anyDuplicated(lineage_vec) == 0)
})

## Asking for more lineages than exist, from either end, must fall back to all
## of them. The bottom slice used to run off the front of the vector and error
## ("only 0's may be mixed with negative subscripts") whenever fewer lineages
## qualified than `num_lineages_bottom`, which is every small simulation at the
## default of 10; the top slice ran off the back and produced NA names.
test_that("plot_anova clamps both ends of the lineage selection (D10)", {
  for(num_lineages in c(2, 3, 4, 6, 8, 9, 12)){
    label <- paste0("num_lineages = ", num_lineages)
    seurat_object <- .plot_seurat_object(num_lineages = num_lineages)
    score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                                 colnames(seurat_object))

    plot1 <- plot_anova(seurat_object = seurat_object,
                        cell_imputed_score = score_vec,
                        assigned_lineage_variable = "assigned_lineage",
                        time_celltype_variable = "time_celltype",
                        day_later = "day7")

    expect_true(inherits(plot1, "ggplot"), info = label)
    lineage_vec <- setdiff(unique(as.character(plot1$data$lineage)), "All")
    expect_true(!any(is.na(lineage_vec)), info = label)
    expect_true(length(lineage_vec) == num_lineages, info = label)

    # The selected names also set the x axis limits, so an unclamped top slice
    # shows up as a phantom "NA (NA)" category with no data behind it.
    limit_vec <- ggplot2::layer_scales(plot1)$x$limits
    expect_true(!any(is.na(limit_vec)), info = label)
    expect_true(!any(grepl("NA", limit_vec)), info = label)
    expect_true(setequal(limit_vec, c(lineage_vec, "All")), info = label)
  }

  # Asking for far more than exist from one end only is the same story.
  seurat_object <- .plot_seurat_object(num_lineages = 3)
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  for(count_pair in list(c(100, 1), c(1, 100), c(100, 100))){
    label <- paste0("top = ", count_pair[1], ", bottom = ", count_pair[2])
    plot1 <- plot_anova(seurat_object = seurat_object,
                        cell_imputed_score = score_vec,
                        assigned_lineage_variable = "assigned_lineage",
                        time_celltype_variable = "time_celltype",
                        day_later = "day7",
                        num_lineages_top = count_pair[1],
                        num_lineages_bottom = count_pair[2])
    lineage_vec <- setdiff(unique(as.character(plot1$data$lineage)), "All")
    expect_true(!any(is.na(lineage_vec)), info = label)
    expect_true(length(lineage_vec) == 3, info = label)
  }
})

## PL10
test_that("plot_anova drops NA scores and filters on current size (PL10)", {
  seurat_object <- .plot_seurat_object(num_lineages = 6)
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  # L1's cells are all NA, so it must disappear entirely.
  na_idx <- which(seurat_object$assigned_lineage == "L1")
  score_vec[na_idx] <- NA

  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = score_vec,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = "day7",
                      num_lineages_top = 3,
                      num_lineages_bottom = 2)

  expect_true(!any(is.na(plot1$data$imputed_count)))
  expect_true(!any(grepl("^L1", as.character(plot1$data$lineage))))

  ## `min_lineage_size` counts cells at the *current* time point -- scored cells
  ## -- not the future size. The two are easy to get backwards on a later edit,
  ## so give a lineage a large future size and a single current cell and assert
  ## it is dropped anyway.
  seurat_object <- .plot_seurat_object(num_lineages = 6)
  lineage_vec <- seurat_object$assigned_lineage
  lineage_vec[lineage_vec == "L6"] <- "L5"
  lineage_vec[1] <- "L6"
  seurat_object$assigned_lineage <- lineage_vec
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))

  plot_strict <- .plot_anova_helper(
    seurat_object = seurat_object,
    cell_imputed_score = score_vec,
    assigned_lineage_variable = "assigned_lineage",
    lineage_future_size = stats::setNames(c(rep(1, 5), 9999),
                                          paste0("L", seq_len(6))),
    min_lineage_size = 2,
    num_lineages_top = 3,
    num_lineages_bottom = 2)

  # L6 has the largest future size by far, and would top the ordering -- but it
  # holds one current cell, so it never reaches the plot.
  expect_true(!any(grepl("^L6", as.character(plot_strict$data$lineage))))
})

##############################################################################
## PL11 / PL12 -- plot_lineageScatterplot()
##############################################################################

## PL11: the two inputs are asserted equal *by position*, not matched by name.
## Two same-length vectors with different names must error rather than silently
## pair the wrong lineages.
test_that("plot_lineageScatterplot aligns by position, not name (PL11)", {
  count_vec <- stats::setNames(c(5, 40, 100), paste0("L", seq_len(3)))
  imputed_vec <- stats::setNames(c(4, 52, 88), paste0("M", seq_len(3)))

  expect_error(plot_lineageScatterplot(lineage_future_count = count_vec,
                                       lineage_imputed_count = imputed_vec))

  # Same names in a different order is also a mismatch, by the same assertion.
  shuffled_vec <- imputed_vec
  names(shuffled_vec) <- paste0("L", c(2, 1, 3))
  expect_error(plot_lineageScatterplot(lineage_future_count = count_vec,
                                       lineage_imputed_count = shuffled_vec))

  # Unnamed input is rejected too.
  expect_error(plot_lineageScatterplot(
    lineage_future_count = as.numeric(count_vec),
    lineage_imputed_count = imputed_vec))
})

## PL12: the x axis is jittered but the reported correlation is computed on the
## unjittered values, so the title is stable while the points move. A reader who
## assumed the correlation described the plotted points would be wrong.
test_that("the reported correlation ignores the jitter (PL12)", {
  count_vec <- stats::setNames(c(5, 40, 100, 3, 60, 12, 250, 8),
                               paste0("L", seq_len(8)))
  imputed_vec <- stats::setNames(c(4, 52, 88, 6, 45, 20, 190, 11),
                                 paste0("L", seq_len(8)))

  set.seed(1)
  plot_a <- plot_lineageScatterplot(lineage_future_count = count_vec,
                                    lineage_imputed_count = imputed_vec)
  set.seed(2)
  plot_b <- plot_lineageScatterplot(lineage_future_count = count_vec,
                                    lineage_imputed_count = imputed_vec)

  # Different jitter draws move the points ...
  expect_true(max(abs(plot_a$data$lineage_future_count -
                        plot_b$data$lineage_future_count)) > 1e-6)
  # ... but the title is identical.
  expect_true(identical(plot_a$labels$title, plot_b$labels$title))

  expected_val <- round(stats::cor(log10(imputed_vec + 1),
                                   log10(count_vec + 1)), 2)
  expect_true(grepl(paste0("Corr:", expected_val), plot_a$labels$title,
                    fixed = TRUE))
})

##############################################################################
## PL14 / PL15 -- plot_simplex()
##############################################################################

## PL14: the coordinate transform is the entire reason this function exists --
## it is what lets the package avoid `ggtern`, which overrides ggplot2 methods
## session-wide. It deserves an actual numerical test rather than a smoke test.
test_that("plot_simplex maps the ternary corners correctly (PL14)", {
  simplex_df <- data.frame(a = c(1, 0, 0, 1, 2),
                           b = c(0, 1, 0, 1, 2),
                           c = c(0, 0, 1, 1, 2))
  plot1 <- plot_simplex(df = simplex_df, x_col = "a", y_col = "b", z_col = "c")

  s3_val <- sqrt(3) / 2
  expected_x <- c(0, 1, 0.5, 0.5, 0.5)
  expected_y <- c(0, 0, s3_val, s3_val / 3, s3_val / 3)

  expect_true(max(abs(plot1$data$.cx - expected_x)) <= 1e-10)
  expect_true(max(abs(plot1$data$.cy - expected_y)) <= 1e-10)

  # Rows are normalized to sum to 1, so scaling a row does not move its point:
  # rows 4 and 5 are (1,1,1) and (2,2,2) and land on the same centroid.
  expect_true(abs(plot1$data$.cx[4] - plot1$data$.cx[5]) <= 1e-10)
  expect_true(abs(plot1$data$.cy[4] - plot1$data$.cy[5]) <= 1e-10)

  # Every point must sit inside the triangle.
  expect_true(all(plot1$data$.cy >= -1e-10))
  expect_true(all(plot1$data$.cy <= s3_val + 1e-10))
})

## PL15
test_that("plot_simplex drops rows summing to zero (PL15)", {
  simplex_df <- data.frame(a = c(1, 0, 0), b = c(0, 1, 0), c = c(0, 0, 0))
  expect_warning(plot1 <- plot_simplex(df = simplex_df, x_col = "a",
                                       y_col = "b", z_col = "c"))
  expect_true(nrow(plot1$data) == 2)
  expect_true(all(!is.nan(plot1$data$.cx)))
})

##############################################################################
## Defects found while writing these tests -- see the plan file, section 0
##############################################################################

## D10: `lineage_names_ordered[(length(x) - num_lineages_bottom + 1):length(x)]`
## goes negative as soon as fewer lineages survive `min_lineage_size` than
## `num_lineages_bottom` (default 10), and R refuses to mix negative and
## positive subscripts. So `plot_anova()` errors outright on any dataset with
## fewer than ten qualifying lineages -- which is most simulations.
test_that("plot_anova handles fewer lineages than num_lineages_bottom (D10)", {
  seurat_object <- .plot_seurat_object(num_lineages = 4)
  score_vec <- stats::setNames(stats::rnorm(ncol(seurat_object)),
                               colnames(seurat_object))
  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = score_vec,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = "day7")

  expect_true(inherits(plot1, "ggplot"))
})