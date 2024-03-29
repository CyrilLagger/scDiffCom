extract_analysis_inputs <- function(
  seurat_object,
  celltype_column_id,
  sample_column_id,
  condition_column_id,
  cond1_name,
  cond2_name,
  assay,
  slot,
  log_scale,
  threshold_min_cells,
  LRI_table,
  LRI_species,
  verbose
) {
  seurat_inputs <- extract_seurat_inputs(
    seurat_object = seurat_object,
    celltype_column_id = celltype_column_id,
    sample_column_id = sample_column_id,
    condition_column_id = condition_column_id,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    threshold_min_cells = threshold_min_cells,
    verbose = verbose
  )
  LRI_inputs <- extract_LRI_inputs(
    data = seurat_inputs$data,
    LRI_table = LRI_table,
    LRI_species = LRI_species,
    verbose = verbose
  )
  condition_inputs <- extract_condition_inputs(
    sample_column_id = sample_column_id,
    condition_column_id = condition_column_id,
    cond1_name = cond1_name,
    cond2_name = cond2_name,
    metadata = seurat_inputs$metadata,
    verbose = verbose
  )
  list(
    data_tr = DelayedArray::t(LRI_inputs$data),
    metadata = seurat_inputs$metadata,
    cell_types = seurat_inputs$cell_types,
    LRI = LRI_inputs$LRI,
    max_nL = LRI_inputs$max_nL,
    max_nR = LRI_inputs$max_nR,
    condition = condition_inputs
  )
}

extract_seurat_inputs <- function(
  seurat_object,
  celltype_column_id,
  sample_column_id,
  condition_column_id,
  assay,
  slot,
  log_scale,
  threshold_min_cells,
  verbose
) {
  cell_type <- NULL
  mes <- "Extracting data from assay '"
  if (slot == "data") {
    mes <- paste0(
      mes,
      assay,
      "' and slot 'data' (assuming normalized log1p-transformed data)."
    )
  } else if (slot == "counts") {
    mes <- paste0(
      mes,
      assay,
      " and slot 'counts' (assuming normalized non-log1p-transformed data)."
    )
  }
  if (verbose) message(mes)
  temp_data <- Seurat::GetAssayData(
    object = seurat_object,
    slot = slot,
    assay = assay
  )
  if(!methods::is(temp_data, "dgCMatrix")) {
    stop(
      paste0(
        "slot ",
        slot,
        " of 'seurat_object' must be of class 'dgCMatrix'"
      )
    )
  }
  if (slot == "data" & !log_scale) {
    if (verbose) {
      message(
        paste0(
          "Converting normalized data from log1p-transformed ",
          "to non-log1p-transformed."
        )
      )
    }
    temp_data <- expm1(temp_data)
  }
  if (slot == "counts" & log_scale) {
    if (verbose) {
      message(
        "Converting data from non-log1p-transformed to log1p-transformed."
      )
    }
    temp_data <- log1p(temp_data)
  }
  if (
    !identical(
      rownames(seurat_object[[]]),
      colnames(seurat_object)
    )
  ) stop(
    "Column names of 'seurat_object' must be identical to ",
    "row names of 'seurat_object[[]]'"
  )
  temp_md <- copy(x = seurat_object[[]])
  if (!(celltype_column_id %in% names(temp_md))) {
    stop(
      paste0(
        "Can't find column '",
        celltype_column_id,
        "' in the meta.data of 'seurat_object'")
    )
  }
  if (celltype_column_id == "cell_id") {
    colnames(temp_md)[colnames(temp_md) == "cell_id"] <- "cell_type"
    celltype_column_id <- "cell_type"
  }
  if (!is.null(sample_column_id)) {
    if (!(sample_column_id %in% names(temp_md))) {
      stop(
        paste0(
          "Can't find column '",
          sample_column_id,
          "' in the meta.data of 'seurat_object'"
        )
      )
    }
  }
  if (!is.null(condition_column_id)) {
    if (!(condition_column_id %in% names(temp_md))) {
      stop(
        paste0(
          "Can't find column '",
          condition_column_id,
          "' in the meta.data of 'seurat_object'"
        )
      )
    }
  }
  cols_to_keep <- c(
    celltype_column_id,
    sample_column_id,
    condition_column_id
  )
  temp_md <- temp_md[, cols_to_keep, drop = FALSE]
  temp_md <- setDT(
    x = temp_md,
    keep.rownames = "cell_id"
  )
  temp_md[, names(temp_md) := lapply(.SD, as.character)]
  if (!is.null(condition_column_id)) {
    temp_cond <- unique(temp_md[[condition_column_id]])
    if (length(temp_cond) != 2) {
      stop(
        paste0(
          "meta.data ",
          condition_column_id,
          " of 'seurat_object' must contain exactly two groups (",
          length(temp_cond),
          " supplied)."
        )
      )
    }
    if(!is.null(sample_column_id)) {
      temp_cols <- c(sample_column_id, condition_column_id)
      temp_md_sample <- unique(temp_md[, temp_cols, with = FALSE])
      temp_samples <- unique(temp_md[[sample_column_id]])
      if (length(temp_samples) != nrow(temp_md_sample)) {
        stop(
          paste0(
            "Column ",
            sample_column_id,
            " of 'seurat_object' must match to column",
            condition_column_id
          )
        )
      }
      new_colnames <- c(
        "cell_id",
        "cell_type",
        "sample_id",
        "condition"
      )
    } else {
      new_colnames <- c("cell_id", "cell_type", "condition")
    }
  } else {
    if(!is.null(sample_column_id)) {
      stop(
        paste0(
          "Parameter 'seurat_column_id' must be supplied ",
          "when parameter 'seurat_sample_id' is not NULL"
        )
      )
    } else {
      new_colnames <- c("cell_id", "cell_type")
    }
  }
  setnames(
    x = temp_md,
    old = c("cell_id", cols_to_keep),
    new = new_colnames
  )
  if(any(grepl("_", temp_md[["cell_type"]], fixed = TRUE))) {
    warning(
      "Underscores ('_') are not allowed in cell-type names: replacing with '-'"
    )
    temp_md[, cell_type := gsub(
      pattern = "_",
      replacement = "-",
      x = cell_type,
      fixed = TRUE)]
  }
  celltypes_filtered <- filter_celltypes(
    metadata = temp_md,
    threshold_min_cells = threshold_min_cells
  )
  metadata <- temp_md[cell_type %in% celltypes_filtered, ]
  data <- temp_data[, colnames(temp_data) %in% metadata$cell_id]
  mes <- paste0(
    "Input data: ",
    nrow(data),
    " genes, ",
    ncol(data),
    " cells and ",
    length(unique(metadata$cell_type)),
    " cell-types."
  )
  if (verbose) message(mes)
  list(
    data = data,
    metadata = metadata,
    cell_types = celltypes_filtered
  )
}

filter_celltypes <- function(
  metadata,
  threshold_min_cells
) {
  if ("condition" %in% colnames(metadata)) {
    filt <- apply(
      X = table(
        metadata$cell_type,
        metadata$condition
      ) >= threshold_min_cells,
      MARGIN = 1,
      FUN = all
    )
  } else {
    filt <- table(metadata$cell_type) >= threshold_min_cells
  }
  res <- names(filt[filt])
  if(length(res) < 2) {
    stop(
      paste0(
        "Inputs must contain at least 2 cell-types with at least ",
        threshold_min_cells,
        " cells (in each condition)."
      )
    )
  }
  res
}

extract_LRI_inputs <- function(
  data,
  LRI_table,
  LRI_species,
  verbose
) {
  LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  cols_compulsory <- c(
    "LRI",
    "LIGAND_1",
    "LIGAND_2",
    "RECEPTOR_1",
    "RECEPTOR_2",
    "RECEPTOR_3"
  )
  if (!all(cols_compulsory %in% names(LRI_table))) {
    stop(
      paste0(
        "'LRI_table' must contain the columns ",
        paste0(
          cols_compulsory,
          collapse = ", "
        )
      )
    )
  }
  LRI_keep <- LRI_table[, cols_compulsory, with = FALSE]
  LRI_keep <- unique(LRI_keep)
  LRI_keep <- LRI_keep[
    LIGAND_1 %in% rownames(data) &
      RECEPTOR_1 %in% rownames(data) &
      LIGAND_2 %in% c(rownames(data), NA) &
      RECEPTOR_2 %in% c(rownames(data), NA) &
      RECEPTOR_3 %in% c(rownames(data), NA),
  ]
  LRI_genes <- unique(c(
    unique(LRI_keep$LIGAND_1),
    unique(LRI_keep$LIGAND_2),
    unique(LRI_keep$RECEPTOR_1),
    unique(LRI_keep$RECEPTOR_2),
    unique(LRI_keep$RECEPTOR_3)
  ))
  data_keep <- data[rownames(data) %in% LRI_genes, ]
  n_ID <- length(unique(LRI_keep$LRI))
  if (n_ID == 0) {
    stop(
      paste0(
        "There are no genes from 'LRI_table' in 'seurat_object'.",
        " Have you supplied 'seurat_object' with correctly formatted gene ",
        " symbols for your species (namely HGNC or MGI)? "
      )
    )
  }
  if (n_ID <= 10) {
    warning(
      paste0(
        "Only ",
        n_ID,
        " ligand-receptor interactions found in the dataset.",
        " Have you supplied 'seurat_object' with correctly formatted gene ",
        " symbols for your species (namely HGNC or MGI)? "
      )
    )
  }
  if (all(is.na(LRI_keep$RECEPTOR_1)) | all(is.na(LRI_keep$LIGAND_1))) {
    stop(
      paste0(
        "'LRI_table' must not contain only NA in columns ",
        "'LIGAND_1'' or 'RECEPTOR_1'."
      )
    )
  } else {
    if (all(is.na(LRI_keep$LIGAND_2))) {
      max_nL <- 1
      LRI_keep <- base::subset(LRI_keep, select = -c(LIGAND_2))
    } else {
      max_nL <- 2
    }
    if (all(is.na(LRI_keep$RECEPTOR_2)) & all(is.na(LRI_keep$RECEPTOR_3))) {
      max_nR <- 1
      LRI_keep <- base::subset(
        LRI_keep,
        select = -c(RECEPTOR_2, RECEPTOR_3)
      )
    } else if (all(is.na(LRI_keep$RECEPTOR_2))) {
      stop(
        paste0(
          "'LRI_table' must not contain only NA in column ",
          "'RECEPTOR_2'' and non-NA in column 'RECEPTOR_3'."
        )
      )
    } else if (all(is.na(LRI_keep$RECEPTOR_3))) {
      max_nR <- 2
      LRI_keep <- base::subset(
        LRI_keep,
        select = -c(RECEPTOR_3)
      )
    } else {
      max_nR <- 3
    }
  }
  mes <- paste0(
    "Input ligand-receptor database: ",
    length(unique(LRI_table$LRI)),
    " ",
    LRI_species,
    " interactions.\n",
    "Number of LRIs that match to genes present in the dataset: ",
    n_ID,
    "."
  )
  if (verbose) message(mes)
  list(
    data = data_keep,
    LRI = LRI_keep,
    max_nL = max_nL,
    max_nR = max_nR
  )
}

extract_condition_inputs <- function(
  sample_column_id,
  condition_column_id,
  cond1_name,
  cond2_name,
  metadata,
  verbose
) {
  if (is.null(condition_column_id)) {
    cond_info <- list(
      is_cond = FALSE,
      is_samp = FALSE
    )
    if (!is.null(cond1_name) | !is.null(cond2_name)) {
      warning(
        paste0(
          "'condition_column_id' is NULL but either 'cond1_name' or ",
          "'cond2_name' is not NULL."
        )
      )
    }
  } else {
    if (is.null(sample_column_id)) {
      cond_info <- list(
        is_cond = TRUE,
        is_samp = FALSE
      )
    } else {
      cond_info <- list(
        is_cond = TRUE,
        is_samp = TRUE
      )
    }
    conds <- unique(metadata$condition)
    if (cond1_name == conds[[1]] & cond2_name == conds[[2]]) {
      cond_info$cond1 <- conds[[1]]
      cond_info$cond2 <- conds[[2]]
    } else if (cond1_name == conds[[2]] & cond2_name == conds[[1]]) {
      cond_info$cond1 <- conds[[2]]
      cond_info$cond2 <- conds[[1]]
    } else {
      stop(
        paste0(
          "Either 'cond1_name' or 'cond2_name' does not match ",
          "with the content of 'condition_column_id':",
          conds[[1]],
          " and ",
          conds[[2]],
          "."
        )
      )
    }
  }
  mes <- "Type of analysis to be performed:"
  if (!cond_info$is_cond) {
    mes <- paste0(
      mes,
      " detection analysis without conditions."
    )
  } else {
    mes <- paste0(
      mes,
      " differential analysis between ",
      cond_info$cond1,
      " and ",
      cond_info$cond2,
      " cells"
    )
    if(cond_info$is_samp) {
      mes <- paste0(
        mes,
        " (based on samples resampling and cells permutation)."
      )
    } else {
      mes <- paste0(
        mes,
        "."
      )
    }
  }
  if (verbose) message(mes)
  cond_info
}
