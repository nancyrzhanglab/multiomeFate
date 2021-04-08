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
