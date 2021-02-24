#' A database of human ligand-receptor interactions.
#'
#' This selection of ligand-receptor interactions (LRIs) is a compilation of
#' eight previous databases (see References). We have distinguished curated
#' LRIs from bioinformatically inferred LRIs (only curated ones are effectively
#' used). We also provide GO terms and KEGG pathways attached to each LRI.
#'
#' The code to build the database is  scDiffCom:::build_LRdb(species = "mouse")
#'
#' @docType data
#'
#' @usage data(LRdb_human)
#'
#' @format A list of data.tables
#'
#' @keywords datasets
#'
#' @references CellChat (\href{https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1}{bioRxiv}),
#'  CellPhoneDB (\href{https://pubmed.ncbi.nlm.nih.gov/32103204/}{PMID: 32103204}),
#'  CellTalkDB (\href{https://pubmed.ncbi.nlm.nih.gov/33147626/}{PMID: 33147626}),
#'  connectomeDB2020 (\href{https://pubmed.ncbi.nlm.nih.gov/33024107/}{PMID: 33024107}),
#'  ICELLNET (\href{https://pubmed.ncbi.nlm.nih.gov/33597528/}{PMID: 33597528}),
#'  NicheNet (\href{https://pubmed.ncbi.nlm.nih.gov/31819264/}{PMID: 31819264}),
#'  SingleCellSignalR (\href{https://pubmed.ncbi.nlm.nih.gov/32196115/}{PMID: 32196115}),
#'  scTensor (\href{https://www.biorxiv.org/content/10.1101/566182v1}{bioRxiv})
#'
#'
"LRdb_human"

#' A database of mouse ligand-receptor interactions.
#'
#' This selection of ligand-receptor interactions (LRIs) is a compilation of
#' eight previous databases (see References). We have distinguished curated
#' LRIs from bioinformatically inferred LRIs (only curated ones are effectively
#' used). We also provide GO terms and KEGG pathways attached to each LRI.
#'
#' The code to build the database is  scDiffCom:::build_LRdb(species = "mouse")
#'
#' @docType data
#'
#' @usage data(LRdb_mouse)
#'
#' @format A list of data.tables
#'
#' @keywords datasets
#'
#' @references CellChat (\href{https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1}{bioRxiv}),
#'  CellPhoneDB (\href{https://pubmed.ncbi.nlm.nih.gov/32103204/}{PMID: 32103204}),
#'  CellTalkDB (\href{https://pubmed.ncbi.nlm.nih.gov/33147626/}{PMID: 33147626}),
#'  connectomeDB2020 (\href{https://pubmed.ncbi.nlm.nih.gov/33024107/}{PMID: 33024107}),
#'  ICELLNET (\href{https://pubmed.ncbi.nlm.nih.gov/33597528/}{PMID: 33597528}),
#'  NicheNet (\href{https://pubmed.ncbi.nlm.nih.gov/31819264/}{PMID: 31819264}),
#'  SingleCellSignalR (\href{https://pubmed.ncbi.nlm.nih.gov/32196115/}{PMID: 32196115}),
#'  scTensor (\href{https://www.biorxiv.org/content/10.1101/566182v1}{bioRxiv})
#'
#'
"LRdb_mouse"

#' A down-sampled Seurat object used for testing and benchmarking
#'
#' This object of class Seurat has been down-sampled from the original
#' Tabula Muris Senis Liver object. Pre-processing and normalization has
#' been performed before down-sampling. It contains 726 features (genes) and
#' 468 samples (cells). It is only intended to be used for testing and
#' benchmarking and does not contain meaningful biological information.
#'
#' @docType data
#'
#' @usage data(seurat_sample_tms_liver)
#'
#' @format An S4 object of class Seurat.
#'
#' @keywords datasets
#'
#' @references \emph{A single-cell transcriptomic atlas characterizes
#'  ageing tissues in the mouse}, Tabula Muris Consortium (2020)
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32669714/}{PMID: 32669714})
#'
"seurat_sample_tms_liver"
