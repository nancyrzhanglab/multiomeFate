context("Test .construct_lineage_data (test-fixture correctness)")

## Plan item P1 (additional_context/test-plan_2026-08-06_kevin.md).
##
## `.construct_lineage_data()` must generate `lineage_future_count` from the
## CYFER model itself:
##
##   y_l ~ Poisson( sum_{i in l} exp(beta0 + x_i' beta) )
##
## `exp(coefficient_vec %*% cell_features[i,,drop=F])` conforms the bare
## length-p vector as p x 1 against a 1 x p matrix, so it evaluates the p x p
## OUTER product and then sums exp() over all p^2 entries. That is a different
## generative model. It is correct only when p == 1, where the product happens
## to be 1 x 1 -- which is exactly why the existing p == 1 tests never caught
## it, and why every p == 2 fixture in the suite carries counts that no
## coefficient vector can reproduce.
##
## Shape of the test: drive both feature variances to ~0 so every cell shares
## the feature vector rep(1, p). The Poisson mean is then identical across
## lineages and exactly sum_{i in l} exp(x_i' beta), so the sample mean of
## `lineage_future_count` over many lineages pins the generative model using
## nothing stronger than the law of large numbers. `p = 1` is swept alongside
## `p = 2` and `p = 3` as a control: it must pass both before and after the fix.

.degenerate_lineage_data <- function(coefficient_vec,
                                     L = 4000,
                                     n_each = 5,
                                     seed = 10){
  .construct_lineage_data(coefficient_vec = coefficient_vec,
                          L = L,
                          n_each = n_each,
                          p = length(coefficient_vec),
                          variance_across_lineage = 1e-12,
                          variance_within_lineage = 1e-12,
                          seed = seed)
}

test_that(".construct_lineage_data draws counts from the CYFER model mean", {
  coefficient_list <- list(c(0.8),
                           c(0.8, -0.5),
                           c(0.4, -0.3, 0.2))

  for(coefficient_vec in coefficient_list){
    p <- length(coefficient_vec)
    label <- paste0("p=", p)
    res <- .degenerate_lineage_data(coefficient_vec = coefficient_vec)

    # `res$cell_features` carries the Intercept column and `res$coefficient_vec`
    # carries the matching Intercept = 0 entry, so this product is the model's
    # per-cell rate exactly as `lineage_imputation()` would compute it.
    rate_vec <- exp(as.numeric(res$cell_features %*% res$coefficient_vec))
    mu_vec <- tapply(rate_vec, res$cell_lineage, sum)
    mu_vec <- mu_vec[names(res$lineage_future_count)]

    expect_equal(mean(res$lineage_future_count),
                 mean(mu_vec),
                 tolerance = 0.05,
                 info = label)
  }
})

## Plan item P2. `.lineage_cleanup()` sorts the lineages, so a fixture whose
## sort order tracks its appearance order cannot distinguish "indexed by name"
## from "indexed by position" -- the condition under which the 1.0.2.000 factor
## bug stayed invisible. The stock fixture is only partly protected here: string
## sorting puts "lin:10" between "lin:1" and "lin:2", so its ordering is
## perturbed by exactly one lineage and otherwise ascending.
## `.scramble_lineage_names()` below reverses the ordering completely, which is
## the adversarial case, and is reused by the invariance tests in the other
## _claude test files.

.scramble_lineage_names <- function(res){
  uniq_lineages <- unique(res$cell_lineage)
  # Deliberately reverse-ordered labels: sorting these puts the lineages in the
  # opposite order to their appearance, so any positional indexing misaligns.
  new_names <- paste0("z", sprintf("%03d", rev(seq_along(uniq_lineages))))
  name_map <- stats::setNames(new_names, uniq_lineages)

  res$cell_lineage <- unname(name_map[res$cell_lineage])
  names(res$lineage_future_count) <-
    unname(name_map[names(res$lineage_future_count)])
  res
}

test_that(".scramble_lineage_names breaks the sort-order/appearance-order tie", {
  set.seed(10)
  res <- .construct_lineage_data()
  res_scrambled <- .scramble_lineage_names(res)

  uniq_original <- unique(res$cell_lineage)
  uniq_scrambled <- unique(res_scrambled$cell_lineage)

  # appearance order is the exact reverse of sort order, so no positional index
  # can coincide with a name-based one
  expect_equal(uniq_scrambled, rev(sort(uniq_scrambled)))

  # the counts must still travel with the right lineage
  expect_setequal(names(res_scrambled$lineage_future_count), uniq_scrambled)
  expect_equal(unname(res_scrambled$lineage_future_count[uniq_scrambled]),
               unname(res$lineage_future_count[uniq_original]))
})

