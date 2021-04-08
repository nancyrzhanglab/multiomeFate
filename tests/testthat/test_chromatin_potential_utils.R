context("Test chromatin protential utility functions")

## .chrom_options is correct

test_that(".chrom_options works", {
  res <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                        cand_method = "nn_xonly", rec_method = "nn_yonly",
                        options = list())
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("form_options", "est_options", "cand_options", "rec_options"))))
})

test_that(".chrom_options works with custom options", {
  res <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                        cand_method = "nn_xonly", rec_method = "nn_yonly",
                        options = list(est_family = "bernoulli", est_switch = FALSE, 
                                       cand_nn = 20))
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("form_options", "est_options", "cand_options", "rec_options"))))
  expect_true(res$est_options$family == "bernoulli")
  expect_true(res$est_options$switch == F)
  expect_true(res$cand_options$nn == 20)
})

test_that(".chrom_options can throw warnings", {
  expect_warning(.chrom_options(form_method = "literal", est_method = "glmnet", 
                                cand_method = "nn_xonly", rec_method = "nn_yonly",
                        options = list(est_family = "poisson", asdf = F, random = 10)))
  
  expect_warning(.chrom_options(form_method = "literal", est_method = "glmnet", 
                                cand_method = "nn_xonly", rec_method = "nn_yonly",
                                options = list(est_family = "poisson", asdf = F, random = 10, est_asdf = 50)))
  
  expect_warning(.chrom_options(form_method = "literal", est_method = "glmnet", 
                                cand_method = "nn_xonly", rec_method = "nn_yonly",
                                options = list(est_asdf = 50)))
})

#########################

## .gene_peak_map is correct

test_that(".gene_peak_map works", {
  set.seed(10)
  p1 <- 20; p2 <- 5; genome_length <- 1000
  df <- generate_df_simple(p1 = p1, p2 = p2, genome_length = genome_length, window = 10)
  options <- .chrom_options(form_method = "literal", est_method = "glmnet", 
                            cand_method = "nn_xonly", rec_method = "nn_yonly",
                        options = list())
  
  res <- .gene_peak_map(df$df_x, df$df_y, options$est_options)
  
  expect_true(all(names(options$est_options) %in% names(res)))
  for(i in names(options$est_options)){
    expect_equal(options$est_options[[i]], res[[i]])
  }
  expect_true(class(res$ht_map) == "hash")
  expect_true(length(res$ht_map) == p2)
})
