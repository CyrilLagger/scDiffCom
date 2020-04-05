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
    message(paste0("Percentage of orthologs returned: ", nrow(ensembl_all)/length(genes)*100), "%")
    return(ensembl_all)
  }
}


#' Prepare Seurat data matrix for downstream analysis.
#'
#'
#'
#' @param seurat_obj a Seurat object
#' @param assay character indicating the assay to choose from
#' @param slot character indicating the slot to choose from
#' @param log_scale logical indicating if returning log-normalized value (only for slot data)
#' @param convert_to_human logical indicating if gene names have to be converted to human orthologs
#' @param return_type character indicating the class of the return data (sparse, dense or data.frame)
#'
#' @return a dgCMatrix, a matrix or a data.frame
#' @export
#'
#' @examples
#' prepare_seurat_data(seurat_random_test)
prepare_seurat_data <- function(seurat_obj,
                            assay = "RNA",
                            slot = "data",
                            log_scale = TRUE,
                            convert_to_human = FALSE,
                            return_type = "dense") {
  data <- Seurat::GetAssayData(object = seurat_obj,
                               slot = slot,
                               assay = assay)
  if(slot == "data" & !log_scale) {
    data <- expm1(data)
  } else if(!log_scale) {
    message("There is no log option for slot counts and scale.data.")
  }
  if(convert_to_human) {
    message("Converting mouse genes to human orthologs.")
    gene_mapping <- get_human_orthologs(rownames(data))
    ng <- nrow(data)
    data <- data[rownames(data) %in% gene_mapping$mouse_symbol, ]
    message(paste0("Removing ", ng - nrow(data), " genes with no orthologs."))
    gene_mapping <- gene_mapping[order(match(gene_mapping$mouse_symbol, rownames(data))),]
    if(identical(rownames(data), gene_mapping$mouse_symbol)) {
      rownames(data) <- gene_mapping$human_symbol
    } else {
      stop("Problem in ordering of the rows (genes). To solve later on.")
    }
  }
  if(class(data) == "dgCMatrix") {
    if(return_type == "sparse") {
      message("Initial class (sparse) dgCMatrix, returning dgCMatrix.")
      return(data)
    } else if(return_type == "dense") {
      message("Initial class (sparse) dgCMatrix, returning (dense) matrix.")
      return(as.matrix(data))
    } else {
      message("Initial class (sparse) dgCMatrix, returning (dense) data.frame.")
      return(as.data.frame(as.matrix(data)))
    }
  } else if(class(data) == "matrix") {
    if((return_type == "dense") | (return_type == "sparse")) {
      message("Initial class (dense) Matrix, returning (dense) Matrix.")
      return(data)
    } else {
      message("Initial class (dense) Matrix, returning (dense) data.frame.")
      return(as.data.frame(data))
    }
  } else {
    stop(paste0("Class ", class(data), " is not recognized."))
  }
}

#' Prepare Seurat metadata for downstream analysis
#'
#' @param seurat_obj a Seurat object
#' @param seurat_cell_type_id a character indicating the column of cell-types in the Seurat object
#' @param condition_id a character indicating a column with some condition on the cells
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' prepare_seurat_metadata(seurat_random_test)
prepare_seurat_metadata <- function(seurat_obj,
                                    seurat_cell_type_id,
                                    condition_id = NULL
) {
  if(is.null(condition_id)) {
    return(data.frame(cell_id = rownames(seurat_obj@meta.data),
                      cell_type = seurat_obj@meta.data[, seurat_cell_type_id]))
  } else {
    return(data.frame(cell_id = rownames(seurat_obj@meta.data),
                      cell_type = seurat_obj@meta.data[, seurat_cell_type_id],
                      condition = seurat_obj@meta.data[, condition_id]))
  }
}
