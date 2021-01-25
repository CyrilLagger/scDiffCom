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
  LRdb_table,
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
  LRdb_inputs <- extract_LRdb_inputs(
    data = seurat_inputs$data,
    LRdb_table = LRdb_table,
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
    data_tr = DelayedArray::t(LRdb_inputs$data),
    metadata = seurat_inputs$metadata,
    cell_types = seurat_inputs$cell_types,
    LRdb = LRdb_inputs$LRdb,
    max_nL = LRdb_inputs$max_nL,
    max_nR = LRdb_inputs$max_nR,
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
  mes <- "Extracting data from assay `"
  if (slot == "data") {
    mes <- paste0(
      mes,
      assay,
      " and slot 'data' (assuming normalized log1p-transformed data)."
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
    stop("`slot` of `seurat_object` must be of class `dgCMatrix`")
  }
  if (slot == "data" & !log_scale) {
    if (verbose) message("Converting data from log1p-transformed to non-log1p-transformed.")
    temp_data <- expm1(temp_data)
  }
  if (slot == "counts" & log_scale) {
    if (verbose) message("Converting data from non-log1p-transformed to log1p-transformed.")
    temp_data <- log1p(temp_data)
  }
  temp_md <- copy(x = seurat_object[[]])
  temp_md <- setDT(
    x = temp_md,
    keep.rownames = "cell_id"
  )
  if (!(celltype_column_id %in% names(temp_md))) {
    stop(paste0("Can't find column `", celltype_column_id, "` in the meta.data of `seurat_object`"))
  }
  if (!is.null(sample_column_id)) {
    if (!(sample_column_id %in% names(temp_md))) {
      stop(paste0("Can't find column `", sample_column_id, "` in the meta.data of `seurat_object`"))
    }
  }
  if (!is.null(condition_column_id)) {
    if (!(condition_column_id %in% names(temp_md))) {
      stop(paste0("Can't find column `", condition_column_id, "` in the meta.data of `seurat_object`"))
    }
  }
  cols_to_keep <- c("cell_id", celltype_column_id, sample_column_id, condition_column_id)
  temp_md <- temp_md[, cols_to_keep, with = FALSE]
  temp_md[, names(temp_md) := lapply(.SD, as.character)]
  if (!is.null(condition_column_id)) {
    temp_cond <- unique(temp_md[[condition_column_id]])
    if (length(temp_cond) != 2) {
      stop(paste0("Column ", condition_column_id, " of `seurat_object` must contain exactly two groups:
                  You have supplied ", length(temp_cond), " groups."))
    }
    if(!is.null(sample_column_id)) {
      temp_cols <- c(sample_column_id, condition_column_id)
      temp_md_sample <- unique(temp_md[, temp_cols, with = FALSE])
      temp_samples <- unique(temp_md[[sample_column_id]])
      if (length(temp_samples) != nrow(temp_md_sample)) {
        stop(paste0("Column ", sample_column_id, " of 'seurat_object' must match to column", condition_column_id))
      }
      new_colnames <- c("cell_id", "cell_type", "sample_id", "condition")
    } else {
      new_colnames <- c("cell_id", "cell_type", "condition")
    }
  } else {
    if(!is.null(sample_column_id)) {
      stop("Parameter `seurat_column_id` must be supplied when parameter `seurat_sample_id` is not NULL")
    } else {

      new_colnames <- c("cell_id", "cell_type")
    }
  }
  setnames(
    x = temp_md,
    old = cols_to_keep,
    new = new_colnames
  )
  celltypes_filtered <- filter_celltypes(
    metadata = temp_md,
    threshold_min_cells = threshold_min_cells
  )
  metadata <- temp_md[cell_type %in% celltypes_filtered, ]
  data <- temp_data[, colnames(temp_data) %in% metadata$cell_id]
  mes <- paste0(
    "Input data for analysis: ",
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
  names(filt[filt])
}

extract_LRdb_inputs <- function(
  data,
  LRdb_table,
  verbose
) {
  LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  cols_compulsory <- c("LR_GENES", "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")
  if (!all(cols_compulsory %in% names(LRdb_table))) {
    stop(paste0("`LRdb_table` must contain the columns ", paste0(cols_compulsory, collapse = ", ")))
  }
  LRdb_keep <- LRdb_table[, cols_compulsory, with = FALSE]
  LRdb_keep <- unique(LRdb_keep)
  LRdb_keep <- LRdb_keep[
    LIGAND_1 %in% rownames(data) &
      RECEPTOR_1 %in% rownames(data) &
      LIGAND_2 %in% c(rownames(data), NA) &
      RECEPTOR_2 %in% c(rownames(data), NA) &
      RECEPTOR_3 %in% c(rownames(data), NA),
    ]
  LRdb_genes <- unique(c(
    unique(LRdb_keep$LIGAND_1),
    unique(LRdb_keep$LIGAND_2),
    unique(LRdb_keep$RECEPTOR_1),
    unique(LRdb_keep$RECEPTOR_2),
    unique(LRdb_keep$RECEPTOR_3)
  ))
  data_keep <- data[rownames(data) %in% LRdb_genes, ]
  n_ID <- length(unique(LRdb_keep$LR_GENES))
  if (n_ID == 0) {
    stop("There are no genes from `LRdb_table` in `seurat_object`")
  }
  if (all(is.na(LRdb_keep$RECEPTOR_1)) | all(is.na(LRdb_keep$LIGAND_1))) {
    stop("`LRdb_table` must not contain only NA in columns `LIGAND_1`` or `RECEPTOR_1`.")
  } else {
    if (all(is.na(LRdb_keep$LIGAND_2))) {
      max_nL <- 1
      LRdb_keep <- base::subset(LRdb_keep, select = -c(LIGAND_2))
    } else {
      max_nL <- 2
    }
    if (all(is.na(LRdb_keep$RECEPTOR_2)) & all(is.na(LRdb_keep$RECEPTOR_3))) {
      max_nR <- 1
      LRdb_keep <- base::subset(LRdb_keep, select = -c(RECEPTOR_2, RECEPTOR_3))
    } else if (all(is.na(LRdb_keep$RECEPTOR_2))) {
      stop("`LRdb_table` must not have only NA in column `RECEPTOR_2`` and non-NA in column `RECEPTOR_3`.")
    } else if (all(is.na(LRdb_keep$RECEPTOR_3))) {
      max_nR <- 2
      LRdb_keep <- base::subset(LRdb_keep, select = -c(RECEPTOR_3))
    } else {
      max_nR <- 3
    }
  }
  mes <- paste0(
    "Using a ligand-receptor database of ",
    length(unique(LRdb_table$LR_GENES)),
    " interactions (",
    n_ID,
    " present in the data).\n",
    "Subsetting data (keeping ",
    nrow(data_keep),
    " genes)."
  )
  if (verbose) message(mes)
  list(
    data = data_keep,
    LRdb = LRdb_keep,
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
      is_cond = FALSE
    )
    if (!is.null(cond1_name) | !is.null(cond2_name)) {
      warning("`condition_column_id` is NULL but either `cond1_name` or `cond2_name` is not NULL.")
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
      stop(paste0(
        "Either `cond1_name` or `cond2_name` does not match with the content of `condition_column_id`:",
        conds[[1]],
        " and ",
        conds[[2]],
        "."
      ))
    }
  }
  mes <- "Interaction analysis will be performed"
  if (!cond_info$is_cond) {
    mes <- paste0(
      mes,
      " on the full dataset at once (no differential expression analysis)."
    )
  } else {
    mes <- paste0(
      mes,
      " with differential expression between the groups ",
      cond_info$cond1,
      " and ",
      cond_info$cond2
    )
    if(cond_info$is_samp) {
      mes <- paste0(
        mes,
        " (based on samples resampling and cells permutation)."
      )
    } else {
      mes <- paste0(
        mes,
        " (based on cells permutation)."
      )
    }
  }
  if (verbose) message(mes)
  return(cond_info)
}
