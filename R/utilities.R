#' Convert mouse gene Symbols to human ortholog Symbols.
#'
#' Only keeps genes with orthologogy_confidence == 1
#' and "ortholog_one2one" mapping. The function deals
#' with none 1:1 mapping as (to explain). Is it to
#' stringent? E.g. CDKN2A is not returned...
#'
#' @param genes A character vector of mouse gene Symbols.
#'
#' @return A dataframe with mouse and human Symbols as columns.
#' @export
#'
#' @examples
#' get_human_orthologs(c("Apoe", "Cdkn2a"))
get_human_orthologs <- function(genes) {
  mart_mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  ensembl_mouse <- biomaRt::getBM(attributes = c("mgi_symbol", "ensembl_gene_id"),
                                  filters = "mgi_symbol",
                                  mart = mart_mouse,
                                  values = genes)
  ensembl_conv <-  biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                 "hsapiens_homolog_associated_gene_name",
                                                 "hsapiens_homolog_orthology_confidence",
                                                 "hsapiens_homolog_perc_id_r1",
                                                 "hsapiens_homolog_orthology_type"),
                                  filters = "ensembl_gene_id",
                                  mart = mart_mouse,
                                  value = ensembl_mouse$ensembl_gene_id)

  ensembl_all <- merge(ensembl_mouse, ensembl_conv)
  ensembl_all <- ensembl_all[ensembl_all$hsapiens_homolog_orthology_type == 'ortholog_one2one' &
                               ensembl_all$hsapiens_homolog_orthology_confidence == 1, ]
  ensembl_all <- ensembl_all[,c('mgi_symbol','hsapiens_homolog_associated_gene_name')]
  ensembl_all <- dplyr::distinct(ensembl_all)
  colnames(ensembl_all) <- c('mouse_symbol', 'human_symbol')
  #check for remaining duplicate
  ensembl_all$ml <- tolower(ensembl_all$mouse_symbol)
  ensembl_all$hl <- tolower(ensembl_all$human_symbol)
  duplicate_mouse <- ensembl_all$mouse_symbol[duplicated(ensembl_all$mouse_symbol)]
  if(length(duplicate_mouse > 0)) {
    for(i in 1:length(duplicate_mouse)) {
      id <- which(ensembl_all$mouse_symbol == duplicate_mouse[[i]])
      id_keep <- which((ensembl_all$mouse_symbol == duplicate_mouse[[i]]) & (ensembl_all$ml == ensembl_all$hl))
      if(length(id_keep) > 0) { id_keep == id_keep[[1]]} else {id_keep == id[[1]]}
      id_remove <- id[id!=id_keep]
      ensembl_all <- ensembl_all[-id_remove,]
    }
  }
  ensembl_all <- ensembl_all[,c(1,2)]
  if(!(length(unique(ensembl_all[,1])) == length(unique(ensembl_all[,2])) &
       length(unique(ensembl_all[,1])) == dim(ensembl_all)[[1]])) {
    stop("Problem of ortholog conversion.")
  } else {
    message(paste0("Percentage of orthologs returned: ", nrow(ensembl_all)/length(genes)*100))
    return(ensembl_all)
  }
}
