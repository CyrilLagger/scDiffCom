#' Run the full differential analysis from a Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param LR_data A ligand-receptor data.frame (a specific format is expected)
#' @param assay The Seurat assay to pull data from; default is "RNA"
#' @param slot The Seurat slot to pull data from; default is "data"
#' @param log_scale Whether to return log-normalized or normalized data (only relevant when slot = "data"); default is TRUE
#' @param seurat_cell_type_id Name of the column of the metadata data.frame containing the cell-type ids
#' @param condition_id Name of the column of the metadata data.frame containing the the condition on the cells. Set to NULL for no conditions
#' @param min_cells Minimum number of cells (per condition if relevant) required to keep a cell-type
#' @param threshold Minimum percentage of cells that have to express a gene to consider this gene in the analysis
#' @param convert_to_human Whether to convert LR genes from mouse to human; default is FALSE
#' @param statistical_analysis Whether to perform the permutation analysis
#' @param iterations Number of permutations when runnning the statistical analysis
#' @param return_type Wheter to return an array or a data.table
#'
#' @return Depending on return_type, either return an array or a data.table
#' @export
#'
#' @examples
run_diffcom_from_seurat <- function(seurat_obj,
                                    LR_data,
                                    assay = "RNA",
                                    slot = "data",
                                    log_scale = FALSE,
                                    seurat_cell_type_id = "cell_ontology_class",
                                    condition_id = NULL,
                                    min_cells = 5,
                                    threshold = 0.1,
                                    convert_to_human = FALSE,
                                    statistical_analysis = TRUE,
                                    iterations = 1000,
                                    return_type = "data.table"
) {
  prep <- prepare_seurat_data(seurat_obj = seurat_obj,
                              assay = assay,
                              slot = slot,
                              log_scale = log_scale,
                              convert_to_human = convert_to_human,
                              return_type = "dense")
  data <- prep$data
  metadata <- prepare_seurat_metadata(seurat_obj = seurat_obj,
                                      seurat_cell_type_id = seurat_cell_type_id,
                                      condition_id = condition_id)
  cell_type_filt <- filter_cell_types(metadata = metadata,
                                      min_cells = min_cells)
  metadata <- metadata[metadata$cell_type %in% cell_type_filt, ]
  data <- data[, colnames(data) %in% metadata$cell_id]
  if(convert_to_human) {
    gene_mapping <- prep$gene_mapping
    stop("Conversion to human orthologs not supported yet in this function, stay tuned :)")
  } else {
    sub <- subset_by_LR(data = data,
                        LR_data = LR_data)
    data <- sub$data
    LR_keep <- sub$LR_data
  }
  data_t <- t(data)
  message("Computing scores.")
  array_res <- get_cci_score(data_t = data_t,
                                metadata = metadata,
                                LR_data = LR_keep,
                                use_thr = TRUE,
                                detection_thr = threshold,
                                condition_id = condition_id)
  if(statistical_analysis) {
    message("Starting statistical analysis.")
    diff_noPerm <- get_diff_score_fast(data_t = data_t,
                                       metadata = metadata,
                                       LR_data = LR_keep,
                                       permutation_test = FALSE)
    distr <- run_permutation(data_t = data_t,
                             metadata = metadata,
                             LR_data = LR_keep,
                             iterations = iterations)
    distr <- cbind(distr, diff_noPerm)
    #return(distr)
    pvals <- rowSums(abs(distr[,1:iterations]) >= abs(distr[,(iterations+1)]))/iterations
    array_res[,,,7,1] <- array(pvals, dim = dim(array_res[,,,1,1]) )
  } else {
    stop("To do")
  }
  if(return_type == "array") {
    return(array_res)
  } else if(return_type == "data.table") {
    dt <- data.table::dcast.data.table(data.table::as.data.table(array_res),
                                       formula = V1 + V2 + V3  ~ V4 + V5 ,
                                       value.var = "value")
    names(dt)[names(dt) == "V1"] <- "LR_pair"
    names(dt)[names(dt) == "V2"] <- "Ligand_cell_type"
    names(dt)[names(dt) == "V3"] <- "Receptor_cell_type"
    return(dt)
  } else {
    stop("Return type not supported.")
  }

}

#' Title
#'
#' @param data_t x
#' @param metadata x
#' @param LR_data x
#' @param use_thr x
#' @param detection_thr x
#' @param condition_id x
#'
#' @return x
#' @export
#'
#' @examples
get_cci_score <- function(data_t,
                          metadata,
                          LR_data,
                          use_thr,
                          detection_thr,
                          condition_id
) {
  if(is.null(condition_id)) {
    stop("To do later, only one condition")
  } else {
    conds <- unique(metadata$condition)
    if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
    ct_sorted <- sort(unique(metadata$cell_type))
    full_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 7, 2))
    dimnames(full_array) <- list(LR_data$SYMB_LR,
                                 ct_sorted,
                                 ct_sorted,
                                 c('LR_score',
                                   'LR_detection',
                                   'Ligand_expr',
                                   'Receptor_expr',
                                   'Ligand_detection',
                                   'Receptor_detection',
                                   'Pvalue'),
                                 conds)
    for(cond in conds) {
      meta <- metadata[metadata$condition == cond, ]
      data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
      data_average <- aggregate_means(data = data_keep,
                                      group = meta$cell_type)
      data_average <- data_average[ct_sorted,]
      detec_rate <- aggregate_means(data = 1*(data_keep > 0),
                                    group = meta$cell_type)
      detec_rate <- detec_rate[ct_sorted,]
      for(i in 1:nrow(LR_data)) {
        LR_scores <- outer(data_average[, LR_data[i,'GENESYMB_L']],
                           data_average[, LR_data[i,'GENESYMB_R']],
                           FUN = "+")/2
        LR_detec <- outer(detec_rate[, LR_data[i,'GENESYMB_L']],
                          detec_rate[, LR_data[i,'GENESYMB_R']],
                          FUN = is_detected,
                          detection_thr)
        if(use_thr) {
          LR_scores[LR_detec == FALSE] <- 0
        }
        full_array[i,,,1, cond] <- LR_scores
        full_array[i,,,2, cond] <- LR_detec
        full_array[i,,,3, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
        full_array[i,,,4, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
        full_array[i,,,5, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
        full_array[i,,,6, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
      }
    }
  }
  return(full_array)
}


#' Title
#'
#' @param data_t x
#' @param metadata x
#' @param LR_data x
#' @param permutation_test x
#'
#' @return x
#' @export
#'
#' @examples
get_diff_score_fast <- function(
  data_t,
  metadata,
  LR_data,
  permutation_test = TRUE
) {
  if(permutation_test) {
    #metadata$cell_id <- sample(metadata$cell_id)
    for(x in unique(metadata$cell_type)) {
      metadata$condition[metadata$cell_type == x] <- sample(metadata$condition[metadata$cell_type == x])
    }
  }
  conds <- unique(metadata$condition)
  ct_sorted <- sort(unique(metadata$cell_type))
  temp_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 2))
  dimnames(temp_array) <- list(NULL, NULL, NULL, conds)
  for(cond in conds) {
    meta <- metadata[metadata$condition == cond, ]
    data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
    data_average <- aggregate_means(data = data_keep,
                                    group = meta$cell_type)
    data_average <- data_average[ct_sorted,]
    for(i in 1:nrow(LR_data)) {
      #be careful Rfast::Outer returns the transpose of base::outer!!
      temp_array[i,,,cond] <- Rfast::transpose(Rfast::Outer(data_average[, LR_data[i,'GENESYMB_L']],
                                      data_average[, LR_data[i,'GENESYMB_R']],
                                      oper = "+")/2)
    }
  }
  return(as.vector(temp_array[,,,2]-temp_array[,,,1]))
}

#' Title
#'
#' @param data_t x
#' @param metadata x
#' @param LR_data x
#' @param iterations x
#'
#' @return x
#' @export
#'
#' @examples
run_permutation <- function(data_t,
                            metadata,
                            LR_data,
                            iterations
) {
  replicate(iterations, get_diff_score_fast(data_t = data_t,
                                           metadata = metadata,
                                           LR_data = LR_data,
                                           permutation_test = TRUE))
}


#' #' Compute arithmetic mean
#' #'
#' #' @param x
#' #' @param y
#' #' @param scale_factor
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cci_amean <- Vectorize(function(
#'   x,
#'   y,
#'   scale_factor = 0
#' ) {
#'   m <- (x+y)/2
#'   if(!(scale_factor == 0)) {
#'     return(m/(m+scale_factor))
#'   } else {
#'     return(m)
#'   }
#' })


#' Title
#'
#' @param x x
#' @param y x
#' @param thr x
#'
#' @return xxx
#' @export
#'
#' @examples
is_detected <- Vectorize(function(x,
                                  y,
                                  thr
) {
  if(x > thr & y > thr) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

create_diffcom_cci <- function(array_score)
{

}


