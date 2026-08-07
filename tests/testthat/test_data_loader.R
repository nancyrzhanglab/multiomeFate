context("Test data_loader")

## `data_loader()` itself reads .RData from a hardcoded lab path, so it is not
## testable in CI and is deliberately unexported as of 1.0.2.001. Its input
## validation and its one pure helper are testable without any files, and the
## helper is the one that fixed a real Seurat v5 bug.

## DL1: `.fix_dimreduc_assay()` is why loading a DimReduc without also loading
## the RNA assay stopped erroring with "Cannot find assay 'RNA'". Seurat v5
## validates `@assay.used` against the object, and every DimReduc file in the
## paper's data was computed with `assay.used = "RNA"`.
test_that(".fix_dimreduc_assay repoints a missing assay (DL1)", {
  set.seed(10)
  count_mat <- matrix(stats::rpois(30 * 20, lambda = 3),
                      nrow = 30,
                      dimnames = list(paste0("gene", seq_len(30)),
                                      paste0("cell", seq_len(20))))
  seurat_object <- suppressWarnings(
    Seurat::CreateSeuratObject(counts = count_mat, assay = "ATAC"))

  embedding_mat <- matrix(stats::rnorm(20 * 2), nrow = 20, ncol = 2,
                          dimnames = list(colnames(seurat_object),
                                          paste0("PC_", seq_len(2))))

  # Case 1: the referenced assay is absent, so it must be repointed to the
  # first available non-"Empty" assay.
  dimreduc_missing <- Seurat::CreateDimReducObject(embeddings = embedding_mat,
                                                   key = "PC_",
                                                   assay = "ATAC")
  dimreduc_missing@assay.used <- "RNA"
  res_missing <- .fix_dimreduc_assay(dimreduc_missing, seurat_object)
  expect_true(res_missing@assay.used == "ATAC")
  expect_true(res_missing@assay.used %in% Seurat::Assays(seurat_object))

  # Case 2: the referenced assay is present, so the object is returned untouched.
  dimreduc_present <- Seurat::CreateDimReducObject(embeddings = embedding_mat,
                                                   key = "PC_",
                                                   assay = "ATAC")
  res_present <- .fix_dimreduc_assay(dimreduc_present, seurat_object)
  expect_true(res_present@assay.used == "ATAC")
  expect_true(identical(res_present, dimreduc_present))

  # Case 3: not a DimReduc at all, so it passes straight through.
  expect_true(identical(.fix_dimreduc_assay(1:5, seurat_object), 1:5))
  expect_true(identical(.fix_dimreduc_assay("not a dimreduc", seurat_object),
                        "not a dimreduc"))
  expect_true(is.null(.fix_dimreduc_assay(NULL, seurat_object)))
})

## DL2: the `which_files` check is the only input validation the function has,
## and it fires before any file is touched -- so it is testable even though the
## rest of the function is not.
test_that("data_loader rejects an unknown which_files entry (DL2)", {
  valid_vec <- c("atac", "chromvar", "lineage", "rna", "saver",
                 "fasttopics", "peakvi", "rna_dimred", "wnn")

  expect_error(data_loader(which_files = "not_a_real_file"))
  expect_error(data_loader(which_files = c("rna", "not_a_real_file")))
  expect_error(data_loader(which_files = c("RNA")))

  # Every documented name must be accepted by the check itself. The call then
  # fails on the missing lab path rather than on validation, so assert the
  # error is not the validation one.
  for(file_name in valid_vec){
    label <- paste0("which_files = ", file_name)
    error_message <- tryCatch({
      suppressWarnings(data_loader(which_files = file_name))
      ""
    }, error = function(e) conditionMessage(e))
    expect_true(!grepl("which_files", error_message), info = label)
    # It gets as far as trying to read the lab path, which is the proof that
    # validation passed.
    expect_true(grepl("cannot open|No such file", error_message),
                info = label)
  }
})

## data_loader() is a helper for the package authors, not public API -- it reads
## from a path that exists only on the lab machines. Pin that it stays
## unexported, since exporting it is what would force the hardcoded path into
## CRAN's field of view.
test_that("data_loader is not exported", {
  expect_true(!("data_loader" %in% getNamespaceExports("multiomeFate")))
  expect_true(is.function(multiomeFate:::data_loader))
})
