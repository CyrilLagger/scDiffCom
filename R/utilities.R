#' Conversion between mouse and human genes using orthologs.
#'
#' Conversion is based on orthologogy_confidence == 1
#' and "ortholog_one2one" mapping.
#'
#' @param genes character vector of the symbols of the genes (mgi_symbol or hgnc_symbol)
#' @param input_species character indicating the species of the input genes; either "mouse" or "human".
#' @param one2one, logical indicating if using one2one orthology relationship
#'
#' @return 2-column data-frame mouse and human orthologs
#' @export
get_orthologs <- function(
  genes,
  input_species,
  one2one = TRUE
) {
  if(input_species == "mouse") {
    id_in <- "mmusculus"
    id_out <- "hsapiens"
    id_gene <- "mgi_symbol"
    name_in <- "mouse"
    name_out <- "human"
  } else if(input_species == "human") {
    id_in <- "hsapiens"
    id_out <- "mmusculus"
    id_gene <- "hgnc_symbol"
    name_in <- "human"
    name_out <- "mouse"
  } else {
    stop("Species not supported in function get_orthologs.")
  }
  mart <- biomaRt::useMart(
    "ensembl",
    dataset = paste0(id_in, "_gene_ensembl")
  )
  ensembl <- biomaRt::getBM(
    attributes = c(id_gene,
                   "ensembl_gene_id"),
    filters = id_gene,
    mart = mart,
    values = genes
  )
  ensembl_conv <- biomaRt::getBM(
    attributes = c("ensembl_gene_id",
                   paste0(id_out, "_homolog_associated_gene_name"),
                   paste0(id_out, "_homolog_orthology_confidence"),
                   paste0(id_out, "_homolog_perc_id_r1"),
                   paste0(id_out, "_homolog_orthology_type")
    ),
    filters = "ensembl_gene_id",
    mart = mart,
    value = ensembl$ensembl_gene_id)
  ensembl_all <- merge(ensembl, ensembl_conv)
  if(one2one) {
    ensembl_all <- ensembl_all[ensembl_all[[paste0(id_out, "_homolog_orthology_type")]] == 'ortholog_one2one' &
                                 ensembl_all[[paste0(id_out, "_homolog_orthology_confidence")]] == 1, ]
  } else {
    ensembl_all <- ensembl_all[ensembl_all[[paste0(id_out, "_homolog_orthology_confidence")]] == 1, ]
  }
  ensembl_all <- ensembl_all[,c(id_gene, paste0(id_out, "_homolog_associated_gene_name"))]
  ensembl_all <- dplyr::distinct(ensembl_all)
  colnames(ensembl_all) <- c("input", "output")
  if(one2one) {
    #check for remaining duplicate
    ensembl_all$inl <- tolower(ensembl_all$input)
    ensembl_all$outl <- tolower(ensembl_all$output)
    duplicate_genes <- ensembl_all$input[duplicated(ensembl_all$input)]
    if(length(duplicate_genes > 0)) {
      for(i in 1:length(duplicate_genes)) {
        id <- which(ensembl_all$input == duplicate_genes[[i]])
        id_keep <- which((ensembl_all$input == duplicate_genes[[i]]) & (ensembl_all$inl == ensembl_all$outl))
        if(length(id_keep) > 0) { id_keep == id_keep[[1]]} else {id_keep == id[[1]]}
        id_remove <- id[id!=id_keep]
        ensembl_all <- ensembl_all[-id_remove,]
      }
    }
    ensembl_all <- ensembl_all[,c(1,2)]
    colnames(ensembl_all) <- c(paste0(name_in, "_symbol"), paste0(name_out, "_symbol"))
    if(!(length(unique(ensembl_all[,1])) == length(unique(ensembl_all[,2])) &
         length(unique(ensembl_all[,1])) == dim(ensembl_all)[[1]])) {
      stop("Problem of ortholog conversion.")
    } else {
      message(paste0("Percentage of orthologs returned: ", nrow(ensembl_all)/length(genes)*100), "%")
      return(ensembl_all)
    }
  } else {
    colnames(ensembl_all) <- c(paste0(name_in, "_symbol"), paste0(name_out, "_symbol"))
    return(ensembl_all)
  }
}

#' Return data and metadata of a Seurat object for downstream analysis.
#'
#' @param seurat_object Seurat object
#' @param assay character indicating the assay of the Seurat object where to pull the data from; default is "RNA".
#' @param slot character indicating the slot of the Seurat object where to pull the data from; default is "data".
#' @param log_scale logical indicating if using log-normalized data (TRUE) or normalized data (FALSE); default is "TRUE".
#' Only considered if slot == "data".
#' @param convert_to_human logical indicating if converting mouse gene names to human gene names.
#' @param return_type character indicating the class of the returned data (sparse, dense or data.frame)
#' @param seurat_cell_type_id character indicating the column specifying the cell-types of the cells
#' @param condition_id character indicating the column specifying the conditions on the cells. Set to NULL to run the analysis
#' on no condition; default is NULL.
#' @param min_cells numeric indicating the minimal number of cells each cluster need to contain to be considered in the analysis.
#'
#' @return A list containing data, metadata, cell-types and possibly orthology mapping.
preprocess_seurat <- function(
  seurat_object,
  assay,
  slot,
  log_scale,
  convert_to_human,
  return_type,
  seurat_cell_type_id,
  condition_id,
  min_cells
) {
  prep_data <- prepare_seurat_data(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    convert_to_human = convert_to_human,
    return_type = return_type
  )
  prep_meta <- prepare_seurat_metadata(
    seurat_object = seurat_object,
    seurat_cell_type_id = seurat_cell_type_id,
    condition_id = condition_id
  )
  cell_types_filtered <- filter_cell_types(
    metadata = prep_meta,
    min_cells = min_cells
  )
  metadata <- prep_meta[prep_meta$cell_type %in% cell_types_filtered, ]
  data <- prep_data$data[, colnames(prep_data$data) %in% metadata$cell_id]
  return(list(
    data = data,
    metadata = metadata,
    cell_types = cell_types_filtered,
    gene_mapping = prep_data$gene_mapping)
  )
}

#' Prepare Seurat data matrix for downstream analysis.
#'
#' @param seurat_object a Seurat object
#' @param assay character indicating the assay of the Seurat object where to pull the data from; default is "RNA".
#' @param slot character indicating the slot of the Seurat object where to pull the data from; default is "data".
#' @param log_scale logical indicating if using log-normalized data (TRUE) or normalized data (FALSE); default is "TRUE".
#' Only considered if slot == "data".
#' @param convert_to_human logical indicating if converting mouse gene names to human gene names.
#' @param return_type character indicating the class of the return data (sparse, dense or data.frame)
#'
#' @return A list with data as first argument (a dgCMatrix, a matrix or a data.frame) and the gene mapping if converstion to orthologs
prepare_seurat_data <- function(
  seurat_object,
  assay = "RNA",
  slot = "data",
  log_scale = TRUE,
  convert_to_human = FALSE,
  return_type = "dense"
) {
  data <- Seurat::GetAssayData(
    object = seurat_object,
    slot = slot,
    assay = assay
    )
  if(slot == "data" & !log_scale) {
    data <- expm1(data)
  } else if(!log_scale) {
    message("There is no log option for slot counts and scale.data.")
  }
  gene_mapping <- NULL
  if(convert_to_human) {
    message("Converting mouse genes to human orthologs.")
    gene_mapping <- get_orthologs(rownames(data),
                                  input_species = "mouse")
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
      #message("Initial class (sparse) dgCMatrix, returning dgCMatrix.")
      return(list(data = data, gene_mapping = gene_mapping))
    } else if(return_type == "dense") {
      #message("Initial class (sparse) dgCMatrix, returning (dense) matrix.")
      return(list(data = as.matrix(data), gene_mapping = gene_mapping))
    } else {
      #message("Initial class (sparse) dgCMatrix, returning (dense) data.frame.")
      return(list(data = as.data.frame(as.matrix(data)), gene_mapping = gene_mapping))
    }
  } else if(class(data) == "matrix") {
    if((return_type == "dense") | (return_type == "sparse")) {
      #message("Initial class (dense) Matrix, returning (dense) Matrix.")
      return(list(data = data, gene_mapping = gene_mapping))
    } else {
      #message("Initial class (dense) Matrix, returning (dense) data.frame.")
      return(list(data = as.data.frame(data), gene_mapping = gene_mapping))
    }
  } else {
    stop(paste0("Class ", class(data), " is not recognized."))
  }
}

#' Prepare Seurat metadata for downstream analysis
#'
#' @param seurat_object a Seurat object
#' @param seurat_cell_type_id a character indicating the column of cell-types in the Seurat object
#' @param condition_id a character indicating a column with some condition on the cells
#'
#' @return a dataframe
prepare_seurat_metadata <- function(
  seurat_object,
  seurat_cell_type_id,
  condition_id = NULL
) {
  if(is.null(condition_id)) {
    return(data.frame(
      cell_id = rownames(seurat_object@meta.data),
                      cell_type = as.character(seurat_object@meta.data[, seurat_cell_type_id]),
                      stringsAsFactors = FALSE)
      )
  } else {
    cond <- seurat_object@meta.data[, condition_id]
    if(length(unique(cond)) != 2)
    {
      stop("Comparison is only possible between two conditions.")
    } else {
      return(data.frame(cell_id = rownames(seurat_object@meta.data),
                        cell_type = as.character(seurat_object@meta.data[, seurat_cell_type_id]),
                        condition = cond,
                        stringsAsFactors = FALSE))
    }
  }
}

#' Return cell-types containing a minimum of cells in each condition
#'
#' @param metadata Result from prepare_seurat_metadata
#' @param min_cells Numeric indicating minimal number of cells in cell-types
#'
#' @return A vector of cell-types with at least min_cells cells
filter_cell_types <- function(
  metadata,
  min_cells

) {
  if("condition" %in% colnames(metadata)) {
    filt <- apply(table(metadata$cell_type,
                        metadata$condition) >= min_cells,
                  MARGIN = 1,
                  FUN = all)
  } else {
    filt <- table(metadata$cell_type) >= min_cells
  }
  cell_type_filt <- names(filt[filt])
  return(cell_type_filt)
}

#' Correspondance between dataset genes and LR-pairs
#'
#' @param data matrix or data.frame with genes as rows and cells as columns.
#' @param LR_data 3-column data.frame of Ligand-Receptor pairs,
#'  satisfying colnames(LR_data) == c("GENESYMB_L", "GENESYMB_R", "SYMB_LR").
#'
#' @return A list containing the subsetted data and the subsetted LR-pairs.
preprocess_LR <- function(
  data,
  LR_data
) {
  if(!identical(colnames(LR_data), c("GENESYMB_L", "GENESYMB_R", "SYMB_LR"))) {
    stop("Wrong formating of LR_data.")
  }
  message(paste0("Number of considered LR pairs: ", length(unique(LR_data$SYMB_LR)), "."))
  LR_keep <- LR_data[LR_data$GENESYMB_L %in% rownames(data) &
                       LR_data$GENESYMB_R %in% rownames(data), ]
  LR_genes <- unique(c(unique(LR_keep$GENESYMB_L), unique(LR_keep$GENESYMB_R)))
  data_keep <- data[rownames(data) %in% LR_genes, ]
  message(paste0("Number of LR pairs in the dataset: ", length(unique(LR_keep$SYMB_LR)), "."))
  LR_keep$L_GENE <- LR_keep$GENESYMB_L
  LR_keep$R_GENE <- LR_keep$GENESYMB_R
  LR_keep$LR_GENES <- LR_keep$SYMB_LR
  LR_keep$SYMB_LR <- NULL
  LR_keep$GENESYMB_L <- NULL
  LR_keep$GENESYMB_R <- NULL
  return(list(data = data_keep, LR_df = LR_keep))
}

#' Determine if two values are bigger than a threshold at the same time.
#'
#' @param x x
#' @param y x
#' @param thr x
#'
#' @return x
is_detected <- Vectorize(function(
  x,
  y,
  thr
) {
  if (x >= thr & y >= thr) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

