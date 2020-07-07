test_that("prepare_seurat_data returns a data.frame", {
  data <- prepare_seurat_data(seurat_random_test,
                              assay = "RNA",
                              slot = "data",
                              convert_to_human = FALSE,
                              return_type = "data.frame")$data
  expect_identical(class(data), "data.frame")
})
