context("Test chromatin potential")

## .update_chrom_df_cand is correct

test_that(".update_chrom_df_cand works", {
  df_res <- .init_chrom_df(50, 1:10, list(11:20), paste0("n", 1:50))
  res <- .update_chrom_df_cand(df_res, c(41:50))
  
  expect_true(all(df_res$num_cand[41:50]+1 == res$num_cand[41:50]))
  expect_true(all(df_res$num_cand[-c(41:50)] == res$num_cand[-c(41:50)]))
})

######################

## .update_chrom_ht is correct

test_that(".update_chrom_ht works", {
  ht_neighbor <- .init_chrom_ht(c(11:20))
  res <- .update_chrom_ht(ht_neighbor, list(list(to = c(11,12), from = c(1:4)),
                                            list(to = c(16:18), from = c(5:8))), F)
  
  expect_true(class(res) == "hash")
  expect_true(all(res[["1"]] == 11:12))
  expect_true(all(res[["2"]] == 11:12))
  expect_true(all(res[["5"]] ==  16:18))
})

########################

## chromatin_potential is correct

test_that("chromatin potential works", {
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
  prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                          vec_start, list_end)
  res <- chromatin_potential(prep_obj, verbose = F)
  
  n <- nrow(dat$obs_x)
  expect_true(class(res) == "chromatin_potential")
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("res_g", "df_res", "ht_neighbor", "options",
                                             "mat_x", "mat_y", "df_x", "df_y", "list_diagnos"))))
  expect_true(is.matrix(res$res_g$mat_g))
  expect_true(all(dim(res$res_g$mat_g) == c(p1,p2)))
  expect_true(length(res$res_g$vec_g) == p2)
  expect_true(nrow(res$df_res) == n)
  expect_true(is.data.frame(res$df_res))
  expect_true(length(res$ht_neighbor) == n)
  expect_true(class(res$ht_neighbor) == "hash")
  expect_true(is.list(res$options))
  expect_true(all(sort(names(res$options)) == sort(c("form_options", "est_options", "rec_options", "cand_options",
                                                     "dim_options", "nn_options"))))
})
