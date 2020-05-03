test_that("get_orthologs returns correct symbols", {
  gene_mouse <- c("Apoe", "Cdkn1a", "a")
  gene_human <- c("APOE", "CDKN1A", "ASIP")

  expect_identical(get_orthologs(gene_mouse,
                                 input_species = "mouse")$human_symbol, gene_human)
})
