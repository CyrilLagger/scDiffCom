test_that("get_human_orthologs returns correct symbols", {
  gene_mouse <- c("Apoe", "Cdkn1a", "a")
  gene_human <- c("APOE", "CDKN1A", "ASIP")

  expect_identical(get_human_orthologs(gene_mouse)$human_symbol, gene_human)
})
