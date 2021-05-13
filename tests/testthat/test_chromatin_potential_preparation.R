context("Test chromatin potential preparation")

## chromatin_potential_prepare is correct

test_that("chromatin_potential_prepare works", {
  set.seed(10)
  p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
  df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
  mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
  timepoints <- 20
  mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
  obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                   list(mat_traj), verbose = F)
  dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)
  
  vec_start <- which(dat$df_info$time <= 0.1)
  list_end <- list(which(dat$df_info$time >= 0.9))
  
  set.seed(10)
  res <- chromatin_potential_prepare(dat$obs_x, dat$obs_y, df$df_x, df$df_y,
                                     vec_start, list_end)
  
  expect_true(is.list(res))
  expect_true(class(res) == "chromatin_potential_prep")
  expect_true(all(sort(names(res)) == sort(c("mat_x", "mat_y", "df_x", "df_y",
                                             "df_res", "dim_reduc_obj", 
                                             "ht_neighbor", "nn_mat", "nn_obj", 
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
  res <- .init_chrom_ht(list(11:20, 21:30))
  
  expect_true(class(res) == "hash")
  for(i in hash::keys(res)){
    expect_true(res[[i]] == as.numeric(i))
  }
})
