context("Test chromatin potential")

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

###############

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
  ht_neighbor <- .init_chrom_ht(list(11:20, 21:30))
  res <- .update_chrom_ht(ht_neighbor, c(31,32), list(c(11:15), c(21:22)))
  
  expect_true(class(res) == "hash")
  expect_true(all(res[["31"]] == 11:15))
  expect_true(all(res[["32"]] == 21:22))
})

##################

## .update_chrom_df_rec is correct

test_that(".update_chrom_df_rec works", {
  df_res <- .init_chrom_df(50, 1:10, list(11:20), paste0("n", 1:50))
  res <- .update_chrom_df_rec(df_res, 21:22, 1)
  
  expect_true(all(res$order_rec[21:22] == 1))
})

########################

## chromatin_potential is correct

test_that("chromatin potential works", {
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
  res <- chromatin_potential(dat$obs_x, dat$obs_y, dat$df_x, dat$df_y,
                             vec_start, list_end, verbose = F)
  
  n <- nrow(dat$obs_x)
  expect_true(class(res) == "chromatin_potential")
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("res_g", "df_res", "ht_neighbor", "options"))))
  expect_true(is.matrix(res$res_g$mat_g))
  expect_true(all(dim(res$res_g$mat_g) == c(p1,p2)))
  expect_true(length(res$res_g$vec_g) == p2)
  expect_true(nrow(res$df_res) == n)
  expect_true(is.data.frame(res$df_res))
  expect_true(length(res$ht_neighbor) == n)
  expect_true(class(res$ht_neighbor) == "hash")
  expect_true(is.list(res$options))
  expect_true(all(sort(names(res$options)) == sort(c("form_options", "est_options", "rec_options", "cand_options"))))
})
