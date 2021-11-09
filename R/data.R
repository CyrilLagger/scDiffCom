#' A collection of human ligand-receptor interactions.
#'
#' This dataset contains a data.table of curated human ligand-receptor
#'  interactions as well as related annotations (GO Terms, KEGG Pathways) and
#'  metadata.
#'
#' The dataset has been built internally in scDiffCom according to
#'  \code{scDiffCom:::build_LRI(species = "human")}. The LRIs have been retrieved
#'  from seven databases (see References). Note that only curated LRIs
#'  have been kept.
#'
#' @docType data
#'
#' @usage data(LRI_human)
#'
#' @format A list with the following items:
#' \enumerate{
#'   \item \code{LRI_curated}: a data.table of curated LRIs
#'   \item \code{LRI_curated_GO}: a data.table with GO terms attached to
#'    curated LRIs
#'   \item \code{LRI_curated_KEGG}: a data.table with KEGG pathways attached to
#'    curated LRIs
#'   \item \code{LRI_retrieved_dates}: dates at which data have been retrieved
#'    from the seven external databases
#'   \item \code{LRI_retrieved_from}: paths or packages from where data have
#'    been retrieved
#'   \item \code{LRI_biomart_ensembl_version}: version of ensembl used for
#'    GO annotation
#' }
#'
#' @keywords datasets
#'
#' @references
#'  CellChat (\href{https://pubmed.ncbi.nlm.nih.gov/33597522/}{PMID: 33597522}),
#'  CellPhoneDB (\href{https://pubmed.ncbi.nlm.nih.gov/32103204/}{PMID: 32103204}),
#'  CellTalkDB (\href{https://pubmed.ncbi.nlm.nih.gov/33147626/}{PMID: 33147626}),
#'  connectomeDB2020 (\href{https://pubmed.ncbi.nlm.nih.gov/33024107/}{PMID: 33024107}),
#'  ICELLNET (\href{https://pubmed.ncbi.nlm.nih.gov/33597528/}{PMID: 33597528}),
#'  NicheNet (\href{https://pubmed.ncbi.nlm.nih.gov/31819264/}{PMID: 31819264}),
#'  SingleCellSignalR (\href{https://pubmed.ncbi.nlm.nih.gov/32196115/}{PMID: 32196115})
#'
#'
"LRI_human"

#' A collection of mouse ligand-receptor interactions.
#'
#' This dataset contains a data.table of curated mouse ligand-receptor
#'  interactions as well as related annotations (GO Terms, KEGG Pathways) and
#'  metadata.
#'
#' The dataset has been built internally in scDiffCom according to
#'  \code{scDiffCom:::build_LRI(species = "mouse")}. The LRIs have been retrieved
#'  from seven databases (see References). Note that only curated LRIs
#'  have been kept.
#'
#' @docType data
#'
#' @usage data(LRI_mouse)
#'
#' @format A list with the following items:
#' \enumerate{
#'   \item \code{LRI_curated}: a data.table of curated LRIs
#'   \item \code{LRI_curated_GO}: a data.table with GO terms attached to
#'    curated LRI
#'   \item \code{LRI_curated_KEGG}: a data.table with KEGG pathways attached to
#'    curated LRIs
#'    \item \code{LRI_retrieved_dates}: dates at which data have been retrieved
#'    from the seven external databases
#'   \item \code{LRI_retrieved_from}: paths or packages from where data have
#'    been retrieved
#'   \item \code{LRI_biomart_ensembl_version}: version of ensembl used for
#'    GO annotation and orthology conversion
#' }
#'
#' @keywords datasets
#'
#' @references
#'  CellChat (\href{https://pubmed.ncbi.nlm.nih.gov/33597522/}{PMID: 33597522}),
#'  CellPhoneDB (\href{https://pubmed.ncbi.nlm.nih.gov/32103204/}{PMID: 32103204}),
#'  CellTalkDB (\href{https://pubmed.ncbi.nlm.nih.gov/33147626/}{PMID: 33147626}),
#'  connectomeDB2020 (\href{https://pubmed.ncbi.nlm.nih.gov/33024107/}{PMID: 33024107}),
#'  ICELLNET (\href{https://pubmed.ncbi.nlm.nih.gov/33597528/}{PMID: 33597528}),
#'  NicheNet (\href{https://pubmed.ncbi.nlm.nih.gov/31819264/}{PMID: 31819264}),
#'  SingleCellSignalR (\href{https://pubmed.ncbi.nlm.nih.gov/32196115/}{PMID: 32196115})
#'
#'
"LRI_mouse"

#' A collection of rat ligand-receptor interactions.
#'
#' This dataset contains a data.table of curated rat ligand-receptor
#'  interactions as well as related annotations (GO Terms, KEGG Pathways) and
#'  metadata.
#'
#' The dataset has been built internally in scDiffCom according to
#'  \code{scDiffCom:::build_LRI(species = "rat")}. The LRIs have been retrieved
#'  from seven databases (see References). Note that only curated LRIs
#'  have been kept.
#'
#' @docType data
#'
#' @usage data(LRI_rat)
#'
#' @format A list with the following items:
#' \enumerate{
#'   \item \code{LRI_curated}: a data.table of curated LRIs
#'   \item \code{LRI_curated_GO}: a data.table with GO terms attached to
#'    curated LRI
#'   \item \code{LRI_curated_KEGG}: a data.table with KEGG pathways attached to
#'    curated LRIs
#'    \item \code{LRI_retrieved_dates}: dates at which data have been retrieved
#'    from the seven external databases
#'   \item \code{LRI_retrieved_from}: paths or packages from where data have
#'    been retrieved
#'   \item \code{LRI_biomart_ensembl_version}: version of ensembl used for
#'    GO annotation and orthology conversion
#' }
#'
#' @keywords datasets
#'
#' @references
#'  CellChat (\href{https://pubmed.ncbi.nlm.nih.gov/33597522/}{PMID: 33597522}),
#'  CellPhoneDB (\href{https://pubmed.ncbi.nlm.nih.gov/32103204/}{PMID: 32103204}),
#'  CellTalkDB (\href{https://pubmed.ncbi.nlm.nih.gov/33147626/}{PMID: 33147626}),
#'  connectomeDB2020 (\href{https://pubmed.ncbi.nlm.nih.gov/33024107/}{PMID: 33024107}),
#'  ICELLNET (\href{https://pubmed.ncbi.nlm.nih.gov/33597528/}{PMID: 33597528}),
#'  NicheNet (\href{https://pubmed.ncbi.nlm.nih.gov/31819264/}{PMID: 31819264}),
#'  SingleCellSignalR (\href{https://pubmed.ncbi.nlm.nih.gov/32196115/}{PMID: 32196115})
#'
#'
"LRI_rat"

#' A down-sampled Seurat object to use for testing and benchmarking
#'
#' This Seurat object has been down-sampled from the original
#' Tabula Muris Senis liver object. Pre-processing and normalization has
#' been performed before down-sampling. It contains 726 features (genes) and
#' 468 samples (cells). It is only intended to be used for testing and
#' benchmarking and does not contain meaningful biological information.
#'
#' @docType data
#'
#' @usage data(seurat_sample_tms_liver)
#'
#' @format An object of class Seurat.
#'
#' @keywords datasets
#'
#' @references \emph{A single-cell transcriptomic atlas characterizes
#'  ageing tissues in the mouse}, Tabula Muris Consortium (2020)
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32669714/}{PMID: 32669714})
#'
"seurat_sample_tms_liver"

#' All gene ontology terms annotated with their levels
#'
#' This data.table contains all GO terms retrieved from the package
#' ontoProc. Each term is annotated with its number of ancestors,
#' parents and children, as well as with its level (i.e. depth) in
#' the gene ontology graph. Levels are computed by scDiffCom according to
#' \code{scDiffCom:::get_GO_LEVELS()}.
#'
#' @docType data
#'
#' @usage data(gene_ontology_level)
#'
#' @format A data.table
#'
#' @keywords datasets
#'
#' @references \href{https://doi.org/doi:10.18129/B9.bioc.ontoProc}{ontoProc}
#'
"gene_ontology_level"


