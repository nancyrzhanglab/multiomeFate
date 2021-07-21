context("Test chromatin protential utility functions")

## .chrom_options is correct

test_that(".chrom_options works", {
  res <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                        form_method = "average", est_method = "glmnet", 
                        cand_method = "nn_any", rec_method = "distant_cor",
                        options = list())
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dim_options", "nn_options",
                                             "form_options", "est_options", 
                                             "cand_options", "rec_options"))))
})

test_that(".chrom_options works with custom options", {
  res <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                        form_method = "average", est_method = "glmnet", 
                        cand_method = "nn_any", rec_method = "distant_cor",
                        options = list(est_switch = FALSE, 
                                       cand_num_cand = 20))
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dim_options", "nn_options",
                                             "form_options", "est_options", 
                                             "cand_options", "rec_options"))))
  expect_true(res$est_options$switch == F)
  expect_true(res$cand_options$num_cand == 20)
})

test_that(".chrom_options can throw warnings", {
  expect_warning(.chrom_options(dim_method = "pca", nn_method = "annoy",
                                form_method = "average", est_method = "glmnet", 
                                cand_method = "nn_any", rec_method = "distant_cor",
                                options = list(asdf = F, random = 10)))
  
  expect_warning(.chrom_options(dim_method = "pca", nn_method = "annoy",
                                form_method = "average", est_method = "glmnet", 
                                cand_method = "nn_any", rec_method = "distant_cor",
                                options = list(asdf = F, random = 10, est_asdf = 50)))
  
  expect_warning(.chrom_options(dim_method = "pca", nn_method = "annoy",
                                form_method = "average", est_method = "glmnet", 
                                cand_method = "nn_any", rec_method = "distant_cor",
                                options = list(est_asdf = 50)))
})

#########################

## .gene_peak_map is correct

test_that(".gene_peak_map works", {
  set.seed(10)
  g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                   directed = F)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
  idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
  combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                            num_per_wave = num_per_wave, 
                                            distinct_waves = distinct_waves)
  
  df <- simulate_data_input(combn_wave_mat)
  df_x <- df$df_x; df_y <- df$df_y
  
  options <- .chrom_options(dim_method = "pca", nn_method = "annoy",
                            form_method = "average", est_method = "glmnet", 
                            cand_method = "nn_any", rec_method = "distant_cor",
                            options = list())
  
  res <- .gene_peak_map(df_x, df_y, options$est_options)
  
  expect_true(all(names(options$est_options) %in% names(res)))
  for(i in names(options$est_options)){
    expect_equal(options$est_options[[i]], res[[i]])
  }
  expect_true(class(res$ht_map) == "hash")
  expect_true(length(res$ht_map) == nrow(df_y))
})
