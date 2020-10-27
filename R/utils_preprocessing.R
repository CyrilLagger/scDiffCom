preprocess_condition <- function(
  condition_col_id,
  cond1_name,
  cond2_name,
  metadata,
  verbose
) {
  if(is.null(condition_col_id)) {
    if(!is.null(cond1_name) | !is.null(cond1_name)) {
      warning("You specified at least one condition name, but no corresponding column name. Performing the analysis without condition.")
    }
    cond_info <- list(
      is_cond = FALSE
    )
    if(verbose) message("Starting CCI analysis without condition.")
  } else {
    conds <- unique(metadata$condition)
    if(cond1_name == conds[[1]] & cond2_name == conds[[2]]) {
      cond_info <- list(
        is_cond = TRUE,
        cond1 = conds[[1]],
        cond2 = conds[[2]]
      )
      if(verbose) message(paste0("Starting CCI analysis with conditions ", conds[[1]], " and ", conds[[2]], "."))
    } else if(cond1_name == conds[[2]] & cond2_name == conds[[1]]) {
      cond_info <- list(
        is_cond = TRUE,
        cond1 = conds[[2]],
        cond2 = conds[[1]]
      )
      if(verbose) message(paste0("Starting CCI analysis with conditions ", conds[[2]], " and ", conds[[1]], "."))
    } else {
      stop(paste0(
        "The names of the conditions do not match with the Seurat metadata: ",
        conds[[1]],
        " and ",
        conds[[2]],
        "."))
    }
  }
  return(cond_info)
}

preprocess_LR <- function(
  data,
  LR_object,
  verbose
) {
  if (verbose) message("Preprocessing LR interactions.")
  LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  cols_compulsory <- c("LR_SORTED", "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")
  LR_keep <- LR_object[, cols_compulsory, with = FALSE]
  LR_keep <- unique(LR_keep)
  if (verbose) message(paste0("Number of interactions in the LR database: ", length(unique(LR_keep$LR_SORTED)), "."))
  LR_keep <- LR_keep[LIGAND_1 %in% rownames(data) &
                       RECEPTOR_1 %in% rownames(data) &
                       LIGAND_2 %in% c(rownames(data), NA) &
                       RECEPTOR_2 %in% c(rownames(data), NA) &
                       RECEPTOR_3 %in% c(rownames(data), NA), ]
  LR_genes <- unique(c(unique(LR_keep$LIGAND_1),
                       unique(LR_keep$LIGAND_2),
                       unique(LR_keep$RECEPTOR_1),
                       unique(LR_keep$RECEPTOR_2),
                       unique(LR_keep$RECEPTOR_3)))
  data_keep <- data[rownames(data) %in% LR_genes, ]
  n_ID <- length(unique(LR_keep$LR_SORTED))
  if(n_ID == 0) {
    stop("There are no LR interactions in the Seurat dataset.")
  } else {
    if (verbose) message(paste0("Number of LR interactions in the Seurat dataset: ", n_ID, "."))
  }
  if(all(is.na(LR_keep$RECEPTOR_1)) | all(is.na(LR_keep$LIGAND_1))) {
    stop("Error of formatting in the LR database: only NA's in LIGAND_1 or RECEPTOR_1.")
  } else {
    if(all(is.na(LR_keep$LIGAND_2))) {
      max_nL <- 1
      LR_keep <- base::subset(LR_keep, select = -c(LIGAND_2))
    } else {
      max_nL <- 2
    }
    if(all(is.na(LR_keep$RECEPTOR_2)) & all(is.na(LR_keep$RECEPTOR_3))) {
      max_nR <- 1
      LR_keep <- base::subset(LR_keep, select = -c(RECEPTOR_2, RECEPTOR_3))
    } else if(all(is.na(LR_keep$RECEPTOR_2))) {
      stop("Error of formatting in the LR database: only NA's in RECEPTOR_2 but not in RECEPTOR_3.")
    } else if(all(is.na(LR_keep$RECEPTOR_3))) {
      max_nR <- 2
      LR_keep <- base::subset(LR_keep, select = -c(RECEPTOR_3))
    } else {
      max_nR <- 3
    }
  }
  return(list(data = data_keep, LR_db = LR_keep, max_nL = max_nL, max_nR = max_nR))
}

preprocess_seurat <- function(
  seurat_object,
  celltype_col_id,
  condition_col_id,
  assay,
  slot,
  log_scale,
  return_type,
  min_cells,
  verbose
) {
  cell_type <- NULL
  if (verbose) message("Preprocessing Seurat object.")
  prep_data <- extract_seurat_data(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    return_type = return_type,
    verbose = verbose
  )
  prep_meta <- extract_seurat_metadata(
    seurat_object = seurat_object,
    celltype_col_id = celltype_col_id,
    condition_col_id = condition_col_id
  )
  celltypes_filtered <- filter_celltypes(
    metadata = prep_meta,
    min_cells = min_cells
  )
  metadata <- prep_meta[cell_type %in% celltypes_filtered, ]
  data <- prep_data[, colnames(prep_data) %in% metadata$cell_id]
  return(
    list(
      data = data,
      metadata = metadata,
      cell_types = celltypes_filtered
    )
  )
}

extract_seurat_data <- function(
  seurat_object,
  assay,
  slot,
  log_scale,
  return_type,
  verbose
) {
  if(slot == "data") {
    if(verbose) message(
      paste0(
        "Extracting data from assay ",
        assay,
        " and slot 'data' (assuming normalized log-transformed values)."
      )
    )
  } else if(slot == "counts") {
    if(verbose) message(
      paste0("Extracting data from assay ",
             assay,
             " and slot 'counts' (assuming normalized non-log-transformed values)."
      )
    )
  } else {
    stop("The slot should be 'counts' or 'data'.")
  }
  temp_data <- Seurat::GetAssayData(
    object = seurat_object,
    slot = slot,
    assay = assay
  )
  if(slot == "data" & !log_scale) {
    if(verbose) message("Converting values from log-transformed to non-log-transformed.")
    temp_data <- expm1(temp_data)
  }
  if(slot == "counts" & log_scale) {
    if(verbose) message("Converting values from non-log-transformed to log-transformed.")
    temp_data <- log1p(temp_data)
  }
  if(class(temp_data) == "dgCMatrix") {
    if(return_type == "sparse") {
      if(verbose) message("Returning a sparse data matrix from the Seurat object.")
      return(temp_data)
    } else if(return_type == "dense") {
      if(verbose) message("Returning a dense data matrix from the Seurat object.")
      return(as.matrix(temp_data))
    } else {
      if(verbose) message("Returning a data.table from the Seurat object.")
      return(data.table::as.data.table(as.matrix(temp_data), keep.rownames = TRUE))
    }
  } else if(class(temp_data) == "matrix") {
    if((return_type == "dense") | (return_type == "sparse")) {
      if(verbose) message("Returning a dense data matrix from the Seurat object.")
      return(temp_data)
    } else {
      if(verbose) message("Return a data.table from the Seurat object.")
      return(data.table::as.data.table(temp_data, keep.rownames = TRUE))
    }
  } else {
    stop(paste0("Class ", class(temp_data), " is not recognized."))
  }
}

extract_seurat_metadata <- function(
  seurat_object,
  celltype_col_id,
  condition_col_id
) {
  cols_to_keep <- c("cell_id", celltype_col_id, condition_col_id)
  temp_md <- data.table::copy(x = seurat_object[[]])
  temp_md <- data.table::setDT(
    x = temp_md,
    keep.rownames = "cell_id"
  )
  temp_md <- temp_md[, cols_to_keep, with = FALSE]
  temp_md[, names(temp_md) := lapply(.SD, as.character)]
  if(!is.null(condition_col_id)) {
    if(length(unique(temp_md[[condition_col_id]])) != 2) {
      stop(paste0("The column ", condition_col_id, " of the Seurat object does not contain exactly 2 conditions."))
    } else {
      new_colnames <- c("cell_id", "cell_type", "condition")
    }
  } else {
    new_colnames <- c("cell_id", "cell_type")
  }
  data.table::setnames(
    x = temp_md,
    old = cols_to_keep,
    new = new_colnames
  )
  return(temp_md)
}

filter_celltypes <- function(
  metadata,
  min_cells
) {
  if("condition" %in% colnames(metadata)) {
    filt <- apply(
      X = table(metadata$cell_type,
                metadata$condition) >= min_cells,
      MARGIN = 1,
      FUN = all
    )
  } else {
    filt <- table(metadata$cell_type) >= min_cells
  }
  cell_type_filt <- names(filt[filt])
  return(cell_type_filt)
}
