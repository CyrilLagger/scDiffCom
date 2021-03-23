#' A database of human ligand-receptor interactions.
#'
#' A list that contains human ligand-receptor interactions (LRIs)
#'  with GO-term and KEGG-pathway annotations. The LRIs have been retrieved
#'  from eight databases (see References). Curated and non-curated LRIs
#'  are distinguished (note that only curated LRIs are used in scDiffCom
#'  interaction analysis, non-curated LRIs are given for reference).
#'
#' The dataset is built internally in scDiffCom according to
#'  scDiffCom:::build_LRI(species = "human").
#'
#' @docType data
#'
#' @usage data(LRI_human)
#'
#' @format A list with seven items:
#' \enumerate{
#'   \item LRI_not_curated: a data.table of curated and non-curated LRIs
#'   \item LRI_curated: a data.table of curated LRIs
#'   \item LRI_curated_GO: a data.table with GO terms attached to
#'    curated LRIs
#'   \item LRI_curated_KEGG: a data.table with KEGG pathways attached to
#'    curated LRIs
#'   \item LRI_retrieved_dates: dates at which data have been retrieved
#'    from the eight external databases
#'   \item LRI_retrieved_from: paths or packages from where data have
#'    been retrieved
#'   \item LRI_biomart_ensembl_version: version of ensembl used for
#'    GO annotation
#' }
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
"LRI_human"

#' A database of mouse ligand-receptor interactions.
#'
#' A list that contains mouse ligand-receptor interactions (LRIs)
#'  with GO-term and KEGG-pathway annotations. The LRIs have been retrieved
#'  from eight databases (see References). Curated and non-curated LRIs
#'  are distinguished (note that only curated LRIs are used in scDiffCom
#'  interaction analysis, non-curated LRIs are given for reference).
#'
#' The code to build the database is  scDiffCom:::build_LRI(species = "mouse")
#'
#' @docType data
#'
#' @usage data(LRI_mouse)
#'
#' @format A list with seven items:
#' \enumerate{
#'   \item LRI_not_curated: a data.table of curated and non-curated LRIs
#'   \item LRI_curated: a data.table of curated LRIs
#'   \item LRI_curated_GO: a data.table with GO terms attached to
#'    curated LRI
#'   \item LRI_curated_KEGG: a data.table with KEGG pathways attached to
#'    curated LRIs
#'    \item LRI_retrieved_dates: dates at which data have been retrieved
#'    from the eight external databases
#'   \item LRI_retrieved_from: paths or packages from where data have
#'    been retrieved
#'   \item LRI_biomart_ensembl_version: version of ensembl used for
#'    GO annotation and orthology conversion
#' }
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
"LRI_mouse"

#' A down-sampled Seurat object used for testing and benchmarking
#'
#' This Seurat object has been down-sampled from the original
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
