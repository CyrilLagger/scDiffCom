context("LRdb")

test_that("Component functions return", {
    ECM_genes = get_ECM_genes(species = 'mouse')
    expect_equal(ECM_genes, c(
      "Col11a1", "Col11a2", "Col1a2", "Col27a1", "Col2a1", "Col3a1", "Col4a1", "Col4a2", "Col4a3", "Col4a4", "Col4a5", "Col4a6", "Col5a1",
      "Col5a2", "Col5a3", "Matn1", "Mmp12", "Mmp13", "Mmp14", "Mmp1a", "Mmp2", "Mmp24", "Mmp7", "Mmp9", "Tecta" 
    ))
})