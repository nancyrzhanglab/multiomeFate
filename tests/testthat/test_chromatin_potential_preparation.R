context("Test chromatin potential preparation")

## chromatin_potential_prepare is correct

test_that("chromatin_potential_prepare works", {
  set.seed(10)
  g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                   directed = F)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
  g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
  idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
  combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                            num_per_wave = num_per_wave, 
                                            distinct_waves = distinct_waves)
  
  res <- simulate_data_input(combn_wave_mat)
  df_x <- res$df_x; df_y <- res$df_y
  list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise
  df_cell <- simulate_df_cell(100, time_max = max(df_y$time_end_scaffold, na.rm = T),
                              num_branch = 3)
  
  set.seed(10)
  dat <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
  mat_x <- dat$obs_x; mat_y <- dat$obs_y
  
  vec_start <- which(df_cell$time <= 10)
  list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
    intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
  })
  
  set.seed(10)
  res <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, vec_start, list_end)
  
  expect_true(is.list(res))
  expect_true(class(res) == "chromatin_potential_prep")
  expect_true(all(sort(names(res)) == sort(c("mat_x", "mat_y", "df_x", "df_y",
                                             "df_res", "dim_reduc_obj", 
                                             "nn_g", "nn_mat", "nn_obj", 
                                             "list_diagnos", "options"))))
})

##############

## .init_chrom_df is correct

test_that(".init_chrom_df works", {
  res <- .init_chrom_df(50, 1:10, list(11:20), paste0("n", 1:50))
  
  expect_true(is.data.frame(res))
  expect_true(all(sort(colnames(res)) == sort(c("idx", "init_state", "num_cand", "order_rec"))))
  expect_true(all(res$num_cand == 0))
  expect_true(all(res$init_state[1:10] == -1))
  expect_true(all(res$init_state[11:20] == 1))
  expect_true(all(is.na(res$init_state[-(1:20)])))
  expect_true(all(res$order_rec[11:20] == 0))
  expect_true(all(is.na(res$order_rec[-c(11:20)])))
})

########

## .init_chrom_ht is correct

test_that(".init_chrom_ht works", {
  res <- .init_chrom_ht(c(11:20))
  
  expect_true(class(res) == "hash")
  for(i in hash::keys(res)){
    expect_true(res[[i]] == as.numeric(i))
  }
})
