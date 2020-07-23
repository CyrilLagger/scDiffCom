test_that("prepare_seurat_data returns a data.table", {
  data <- prepare_seurat_data(seurat_random_test,
                              assay = "RNA",
                              slot = "data",
                              return_type = "data.table")
  expect_identical(class(data), c("data.table", "data.frame"))
})
