#' A database of human ligand-receptor interactions.
#'
#' The DB has been compiled from eight external studies.
#'
#' Can be obtained by running scDiffCom:::build_LRdb(species = "human")
#' @docType data
#'
#' @usage data(LRdb_human)
#'
#' @format A list of 3 data.tables: a non-curated table, a curated table and a table of corresponding GO terms
#'
#' @keywords datasets
#'
#' @references to add
#' (\href{to_add}{PubMed})
#'
"LRdb_human"

#' A database of mouse ligand-receptor interactions.
#'
#' The DB has been compiled from eight external studies.
#'
#' Can be obtained by running scDiffCom:::build_LRdb(species = "mouse")
#' @docType data
#'
#' @usage data(LRdb_mouse)
#'
#' @format A list of 3 data.tables: a non-curated table, a curated table and a table of corresponding GO terms
#'
#' @keywords datasets
#'
#' @references to add
#' (\href{to_add}{PubMed})
#'
"LRdb_mouse"

#' A (sampled) Seurat object from the Tabula Muris Senis liver.
#'
#' This is only a sample from the original object. It only contains
#' two cell-types and roughly 2000 genes. The slot "data" contains
#' log-normalized values with normalization performed before the down-sampling.
#' It is only intended to be used as a toy model to show the functionalities of the package.
#' It does not convey meaningful biological information.
#'
#' @docType data
#'
#' @usage data(seurat_sample_tms_liver)
#'
#' @format A Seurat object.
#'
#' @keywords datasets
#'
#' @references to add
#' (\href{to_add}{PubMed})
#'
"seurat_sample_tms_liver"
