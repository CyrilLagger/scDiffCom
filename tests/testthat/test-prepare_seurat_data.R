test_that("prepare_seurat_data returns a data.frame", {
  data <- prepare_seurat_data(seurat_random_test,
                              assay = "RNA",
                              slot = "data",
                              convert_to_human = FALSE,
                              return_type = "data.frame")
  expect_identical(class(data), "data.frame")
})


test_that("prepare_seurat_data returns correct value when converting to orthologs", {
  data <- prepare_seurat_data(seurat_random_test,
                              assay = "RNA",
                              slot = "data",
                              convert_to_human = TRUE,
                              return_type = "sparse")
  expect_identical(seurat_random_test$RNA@data["Chd7", "GCGCGATGTCCAGTGC-1-44-1-0"],
                   data["CHD7", "GCGCGATGTCCAGTGC-1-44-1-0"])
})
