

run_diffcom_from_seurat <- function(seurat_obj,
                                    LR_data,
                                    assay = "RNA",
                                    slot = "data",
                                    log_scale = TRUE,
                                    seurat_cell_type_id = "cell_ontology_class",
                                    condition_id = NULL,
                                    min_cells = 5,
                                    threshold = 0.1,
                                    convert_to_human = FALSE,
                                    statistical_analysis = TRUE,
                                    iterations = 1000
) {
  pp_seurat <- preprocess_seurat(seurat_obj = seurat_obj,
                                 assay = assay,
                                 slot = slot,
                                 log_scale = log_scale,
                                 convert_to_human = convert_to_human,
                                 return_type = "dense",
                                 seurat_cell_type_id = seurat_cell_type_id,
                                 condition_id = condition_id,
                                 min_cells = min_cells)
  if(convert_to_human) {
    stop("Conversion to human orthologs not supported yet in this function, stay tuned :)")
  }
  pp_LR <- preprocess_LR(data = pp_seurat$data,
                         LR_data = LR_data)
  cci_analysis <- run_cci_analysis(expr_tr = t(pp_LR$data),
                                   metadata = pp_seurat$metadata,
                                   cell_types = pp_seurat$cell_types,
                                   LR_df = pp_LR$LR_df,
                                   threshold = threshold,
                                   condition_id = condition_id,
                                   statistical_analysis = statistical_analysis,
                                   iterations = iterations)
}

run_cci_analysis <- function(expr_tr,
                             metadata,
                             cell_types,
                             LR_df,
                             threshold,
                             condition_id,
                             statistical_analysis,
                             iterations
) {
  if(is.null(condition_id)) {
    array_noCond <- run_simple_analysis(expr_tr = expr_tr,
                                        metadata = metadata,
                                        cell_types = cell_types,
                                        cond = NULL,
                                        LR_df = LR_df,
                                        threshold = threshold,
                                        compute_fast = FALSE)
    if(!statistical_analysis) {
      cci_dt <- data.table::dcast.data.table(data.table::as.data.table(array_noCond),
                                             formula = V1 + V2 + V3  ~ V4 ,
                                             value.var = "value")
    } else {
      stat_res <- run_stat_analysis(expr_tr = expr_tr,
                                       metadata = metadata,
                                       cell_types = cell_types,
                                       LR_df = LR_df,
                                       cond1 = NULL,
                                       cond2 = NULL,
                                       iterations = iterations)
      big_array <- array(data = 0,
                         dim = c(dim(array_noCond)[1:3], dim(array_noCond)[4] + 1))
      dnames <- dimnames(array_noCond)
      dimnames(big_array) <- list(dnames[[1]], dnames[[2]], dnames[[3]], c(dnames[[4]], "pval_spec"))
      big_array[,,,1:6] <- array_noCond
      big_array[,,,7] <- array(data = stat_res, dim = dim(array_noCond)[1:3])
      cci_dt <- data.table::dcast(data.table::as.data.table(big_array),
                                  formula = V1 + V2 + V3 ~ V4,
                                  value.var = "value")
    }

    names(cci_dt)[names(cci_dt) == "V1"] <- "LR_pair"
    names(cci_dt)[names(cci_dt) == "V2"] <- "Ligand_cell_type"
    names(cci_dt)[names(cci_dt) == "V3"] <- "Receptor_cell_type"
  } else{
    conds <- unique(metadata$condition)
    if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
    cond1 <- conds[[1]]
    cond2 <- conds[[2]]
    array_cond1 <- run_simple_analysis(expr_tr = expr_tr,
                                       metadata = metadata,
                                       cell_types = cell_types,
                                       cond = cond1,
                                       LR_df = LR_df,
                                       threshold = threshold,
                                       compute_fast = FALSE)
    array_cond2 <- run_simple_analysis(expr_tr = expr_tr,
                                       metadata = metadata,
                                       cell_types = cell_types,
                                       cond = cond2,
                                       LR_df = LR_df,
                                       threshold = threshold,
                                       compute_fast = FALSE)
    if(!statistical_analysis) {
      cci_dt1 <- data.table::dcast.data.table(data.table::as.data.table(array_cond1, sorted = FALSE),
                                              formula = V1 + V2 + V3  ~ V4 ,
                                              value.var = "value")
      cci_dt2 <- data.table::dcast.data.table(data.table::as.data.table(array_cond2),
                                              formula = V1 + V2 + V3  ~ V4 ,
                                              value.var = "value")
      cci_dt <- data.table::merge.data.table(x = cci_dt1,
                                             y = cci_dt2,
                                             by = c("V1", "V2", "V3"),
                                             sort = FALSE,
                                             suffixes = c(paste0("_",cond1), paste0("_", cond2)))
    } else {
      stat_res <- run_stat_analysis(expr_tr = expr_tr,
                                    metadata = metadata,
                                    cell_types = cell_types,
                                    LR_df = LR_df,
                                    cond1 = cond1,
                                    cond2 = cond2,
                                    iterations = iterations)
      #to keep the order of each cci the same
      big_array <- array(data = 0,
                         dim = c(dim(array_cond1)[1:3], dim(array_cond1)[4] + 2, 2))
      dnames <- dimnames(array_cond1)
      dimnames(big_array) <- list(dnames[[1]], dnames[[2]], dnames[[3]], c(dnames[[4]], "pval_spec", "pval_diff"), c(cond1, cond2))
      big_array[,,,1:6,1] <- array_cond1
      big_array[,,,7,1] <- array(data = stat_res$pvals_cond1, dim = dim(array_cond1)[1:3])
      big_array[,,,8,1] <- array(data = stat_res$pvals_diff, dim = dim(array_cond1)[1:3])
      big_array[,,,1:6,2] <- array_cond2
      big_array[,,,7,2] <- array(data = stat_res$pvals_cond2, dim = dim(array_cond1)[1:3])
      big_array[,,,8,2] <- array(data = stat_res$pvals_diff, dim = dim(array_cond1)[1:3])

      cci_dt <- data.table::dcast(data.table::as.data.table(big_array),
                                  formula = V1 + V2 + V3 ~ V4 + V5,
                                  value.var = "value")

    }
    names(cci_dt)[names(cci_dt) == "V1"] <- "LR_pair"
    names(cci_dt)[names(cci_dt) == "V2"] <- "Ligand_cell_type"
    names(cci_dt)[names(cci_dt) == "V3"] <- "Receptor_cell_type"

  }
  return(cci_dt)
}

run_stat_analysis <- function(expr_tr,
                              metadata,
                              cell_types,
                              LR_df,
                              cond1,
                              cond2,
                              iterations
) {
  if(is.null(cond1) | is.null(cond2)) {
    cci_noCond <- as.vector(run_simple_analysis(expr_tr = expr_tr,
                                                metadata = metadata,
                                                cell_types = cell_types,
                                                cond = NULL,
                                                LR_df = LR_df,
                                                threshold = 0,
                                                compute_fast = TRUE))
    cci_perm <- replicate(iterations,
                          run_stat_iteration(expr_tr = expr_tr,
                                             metadata = metadata,
                                             cell_types = cell_types,
                                             LR_df = LR_df,
                                             cond1 = NULL,
                                             cond2 = NULL),
                          simplify = "array")
    distr <- cbind(cci_perm, cci_noCond)
    pvals <- rowSums(distr[,1:iterations] >= distr[,(iterations+1)])/iterations
    return(pvals)
  } else {
    cci_cond1 <- as.vector(run_simple_analysis(expr_tr = expr_tr,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond1,
                                               LR_df = LR_df,
                                               threshold = 0,
                                               compute_fast = TRUE))
    cci_cond2 <- as.vector(run_simple_analysis(expr_tr = expr_tr,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond2,
                                               LR_df = LR_df,
                                               threshold = 0,
                                               compute_fast = TRUE))
    array_perm <- replicate(iterations,
                            run_stat_iteration(expr_tr = expr_tr,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               LR_df = LR_df,
                                               cond1 = cond1,
                                               cond2 = cond2),
                            simplify = "array")
    distr_diff <- cbind(array_perm[,1,], cci_cond2 - cci_cond1)
    pvals_diff <- rowSums(abs(distr_diff[,1:iterations]) >= abs(distr_diff[,(iterations+1)]))/iterations
    distr_cond1 <- cbind(array_perm[,2,], cci_cond1)
    pvals_cond1 <- rowSums(distr_cond1[,1:iterations] >= distr_cond1[,(iterations+1)])/iterations
    distr_cond2 <- cbind(array_perm[,3,], cci_cond2)
    pvals_cond2 <- rowSums(distr_cond2[,1:iterations] >= distr_cond2[,(iterations+1)])/iterations
    return(list(pvals_diff = pvals_diff,
                pvals_cond1 = pvals_cond1,
                pvals_cond2 = pvals_cond2))
  }
}

run_stat_iteration <- function(expr_tr,
                               metadata,
                               cell_types,
                               LR_df,
                               cond1,
                               cond2
) {
  if(is.null(cond1) | is.null(cond2)) {
    meta_ct <- metadata
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    cci_permct <- run_simple_analysis(expr_tr = expr_tr,
                                      metadata = meta_ct,
                                      cell_types = cell_types,
                                      cond = NULL,
                                      LR_df = LR_df,
                                      threshold = 0,
                                      compute_fast = TRUE)
    return(as.vector(cci_permct))
  } else {
    meta_cond <- metadata
    for(x in cell_types) {
      meta_cond$condition[meta_cond$cell_type == x] <- sample(meta_cond$condition[meta_cond$cell_type == x])
    }
    cci_permcond_cond1 <- run_simple_analysis(expr_tr = expr_tr,
                                              metadata = meta_cond,
                                              cell_types = cell_types,
                                              cond = cond1,
                                              LR_df = LR_df,
                                              threshold = 0,
                                              compute_fast = TRUE)
    cci_permcond_cond2 <- run_simple_analysis(expr_tr = expr_tr,
                                              metadata = meta_cond,
                                              cell_types = cell_types,
                                              cond = cond2,
                                              LR_df = LR_df,
                                              threshold = 0,
                                              compute_fast = TRUE)
    diff_permcond <- as.vector(cci_permcond_cond2 - cci_permcond_cond1)

    meta_ct <- metadata
    meta_ct$cell_type[meta_ct$condition == cond1] <- sample(meta_ct$cell_type[meta_ct$condition == cond1])
    meta_ct$cell_type[meta_ct$condition == cond2] <- sample(meta_ct$cell_type[meta_ct$condition == cond2])
    cci_permct_cond1 <- as.vector(run_simple_analysis(expr_tr = expr_tr,
                                                      metadata = meta_ct,
                                                      cell_types = cell_types,
                                                      cond = cond1,
                                                      LR_df = LR_df,
                                                      threshold = 0,
                                                      compute_fast = TRUE))
    cci_permct_cond2 <- as.vector(run_simple_analysis(expr_tr = expr_tr,
                                                      metadata = meta_ct,
                                                      cell_types = cell_types,
                                                      cond = cond2,
                                                      LR_df = LR_df,
                                                      threshold = 0,
                                                      compute_fast = TRUE))
    return(cbind(diff_permcond, cci_permcond_cond1, cci_permct_cond2))
    #return(list(diff_permcond = diff_permcond,
    #            permct_cond1 = as.vector(cci_permct_cond1),
    #            permct_cond2 = as.vector(cci_permct_cond2) ))
  }
}


run_simple_analysis <- function(expr_tr,
                                metadata,
                                cell_types,
                                cond,
                                LR_df,
                                threshold,
                                compute_fast
) {
  averaged_expr <- aggregate_cells(expr_tr = expr_tr,
                                   metadata = metadata,
                                   cell_types = cell_types,
                                   cond = cond)
  cci <- build_cci_array(averaged_expr = averaged_expr,
                         LR_df = LR_df)
  if(compute_fast) {
    return(cci)
  } else {
    detection_rate <- aggregate_cells(expr_tr = 1*(expr_tr > 0),
                                      metadata = metadata,
                                      cell_types = cell_types,
                                      cond = cond)
    drate <- build_drate_array(detection_rate = detection_rate,
                               detection_thr = threshold,
                               LR_df = LR_df )
    full_array <- array(data = 0,
                        dim = c(dim(cci), 6))
    dimnames(full_array) <- list(paste(LR_df$ligand, LR_df$receptor, sep = "_"),
                                 cell_types,
                                 cell_types,
                                 c('LR_score',
                                   'LR_detection',
                                   'Ligand_expr',
                                   'Receptor_expr',
                                   'Ligand_detec',
                                   'Receptor_detec'
                                 )
    )
    full_array[,,,1] <- cci
    full_array[,,,2] <- drate
    for(i in 1:nrow(LR_df)) {
      full_array[i,,,3] <- matrix(averaged_expr[, LR_df[i,'ligand']],
                                  ncol = nrow(averaged_expr), nrow = nrow(averaged_expr), byrow = FALSE)
      full_array[i,,,4] <- matrix(averaged_expr[, LR_df[i,'receptor']],
                                  ncol = nrow(averaged_expr), nrow = nrow(averaged_expr), byrow = TRUE)
      full_array[i,,,5] <- matrix(detection_rate[, LR_df[i,'ligand']],
                                  ncol = nrow(detection_rate), nrow = nrow(detection_rate), byrow = FALSE)
      full_array[i,,,6] <- matrix(detection_rate[, LR_df[i,'receptor']],
                                  ncol = nrow(detection_rate), nrow = nrow(detection_rate), byrow = TRUE)
    }
    return(full_array)
  }
}

#' Aggregate data per cell-types and conditions
#'
#' @param data Data matrix
#' @param group Vector indicating groups
#'
#' @return A matrix with the mean expression of each gene for each cell-type
#' @export
#'
#' @examples
aggregate_cells <- function(expr_tr,
                            metadata,
                            cell_types,
                            cond
) {
  if(is.null(cond)) {
    meta <- metadata
    expr <- expr_tr
  } else {
    meta <- metadata[metadata$condition == cond, ]
    expr <- expr_tr[rownames(expr_tr) %in% meta$cell_id, ]
  }
  sums <- rowsum(x = expr,
                 group = meta$cell_type)
  aggr <- sums/as.vector(table(meta$cell_type))
  return(aggr[cell_types, ])
}

#' Title
#'
#' @param averaged_expr
#' @param LR_df
#'
#' @return
#' @export
#'
#' @examples
build_cci_array <- function(averaged_expr,
                            LR_df
) {
  cci_array <- array(data = 0,
                     dim = c(nrow(LR_df), nrow(averaged_expr), nrow(averaged_expr)))
  for(i in 1:nrow(LR_df)) {
    cci_array[i, , ] <- compute_LR_score(ligand = LR_df[i, "ligand"],
                                         receptor = LR_df[i, "receptor"],
                                         averaged_expr= averaged_expr)
  }
  return(cci_array)
}

#' Title
#'
#' @param averaged_expr
#' @param LR_df
#'
#' @return
#' @export
#'
#' @examples
build_drate_array <- function(detection_rate,
                              detection_thr,
                              LR_df
) {
  drate_array <- array(data = 0,
                       dim = c(nrow(LR_df), nrow(detection_rate), nrow(detection_rate)))
  for(i in 1:nrow(LR_df)) {
    drate_array[i, , ] <- compute_LR_drate(ligand = LR_df[i, "ligand"] ,
                                           receptor = LR_df[i, "receptor"],
                                           detection_rate = detection_rate,
                                           detection_thr = detection_thr)
  }
  return(drate_array)
}

#' Compute LR scores for one LR-pair over all cell-types
#'
#' @param ligand Name of the ligand
#' @param receptor Name of the receptor
#' @param average_expression A matrix with the mean expression of a each gene per cell-type
#'
#' @return A matrix with the score of LR for each cell-type pair combination
#' @export
#'
#' @examples
compute_LR_score <- function(ligand,
                             receptor,
                             averaged_expr
) {
  #be careful Rfast::Outer returns the transpose of base::outer!
  Rfast::transpose(Rfast::Outer(averaged_expr[, ligand],
                                averaged_expr[, receptor],
                                oper = "+")/2)
}

#' Compute LR detection rate for one LR-pair over all cell-types
#'
#' @param ligand Name of the ligand
#' @param receptor Name of the receptor
#' @param average_expression A matrix with the detection rate of a each gene per cell-type
#'
#' @return A matrix with the detection rate of LR for each cell-type pair combination
#' @export
#'
#' @examples
compute_LR_drate <- function(ligand,
                             receptor,
                             detection_rate,
                             detection_thr
) {
  outer(detection_rate[, ligand],
        detection_rate[, receptor],
        FUN = is_detected,
        detection_thr)
}

#' Determine if two values are bigger than a threshold at the same time.
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
  if(x >= thr & y >= thr) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

###############################

###############


#'
#'
#'
#'
#' #' Run the full differential analysis from a Seurat object.
#' #'
#' #' @param seurat_obj A Seurat object
#' #' @param LR_data A ligand-receptor data.frame (a specific format is expected)
#' #' @param assay The Seurat assay to pull data from; default is "RNA"
#' #' @param slot The Seurat slot to pull data from; default is "data"
#' #' @param log_scale Whether to return log-normalized or normalized data (only relevant when slot = "data"); default is TRUE
#' #' @param seurat_cell_type_id Name of the column of the metadata data.frame containing the cell-type ids
#' #' @param condition_id Name of the column of the metadata data.frame containing the the condition on the cells. Set to NULL for no conditions
#' #' @param min_cells Minimum number of cells (per condition if relevant) required to keep a cell-type
#' #' @param threshold Minimum percentage of cells that have to express a gene to consider this gene in the analysis
#' #' @param convert_to_human Whether to convert LR genes from mouse to human; default is FALSE
#' #' @param statistical_analysis Whether to perform the permutation analysis
#' #' @param iterations Number of permutations when runnning the statistical analysis
#' #' @param return_type Wheter to return an array or a data.table
#' #'
#' #' @return Depending on return_type, either return an array or a data.table
#' #' @export
#' #'
#' #' @examples
#' run_diffcom_from_seurat <- function(seurat_obj,
#'                                     LR_data,
#'                                     assay = "RNA",
#'                                     slot = "data",
#'                                     log_scale = FALSE,
#'                                     seurat_cell_type_id = "cell_ontology_class",
#'                                     condition_id = NULL,
#'                                     min_cells = 5,
#'                                     threshold = 0.1,
#'                                     convert_to_human = FALSE,
#'                                     statistical_analysis = TRUE,
#'                                     iterations = 1000,
#'                                     return_type = "data.table"
#' ) {
#'   prep <- prepare_seurat_data(seurat_obj = seurat_obj,
#'                               assay = assay,
#'                               slot = slot,
#'                               log_scale = log_scale,
#'                               convert_to_human = convert_to_human,
#'                               return_type = "dense")
#'   data <- prep$data
#'   metadata <- prepare_seurat_metadata(seurat_obj = seurat_obj,
#'                                       seurat_cell_type_id = seurat_cell_type_id,
#'                                       condition_id = condition_id)
#'   cell_type_filt <- filter_cell_types(metadata = metadata,
#'                                       min_cells = min_cells)
#'   metadata <- metadata[metadata$cell_type %in% cell_type_filt, ]
#'   data <- data[, colnames(data) %in% metadata$cell_id]
#'   if(convert_to_human) {
#'     gene_mapping <- prep$gene_mapping
#'     stop("Conversion to human orthologs not supported yet in this function, stay tuned :)")
#'   } else {
#'     sub <- subset_by_LR(data = data,
#'                         LR_data = LR_data)
#'     data <- sub$data
#'     LR_keep <- sub$LR_data
#'   }
#'   data_t <- t(data)
#'   message("Computing scores.")
#'   array_res <- get_cci_score(data_t = data_t,
#'                              metadata = metadata,
#'                              LR_data = LR_keep,
#'                              use_thr = TRUE,
#'                              detection_thr = threshold,
#'                              condition_id = condition_id)
#'   if(statistical_analysis) {
#'     message("Starting statistical analysis.")
#'     diff_noPerm <- get_diff_score_fast(data_t = data_t,
#'                                        metadata = metadata,
#'                                        LR_data = LR_keep,
#'                                        permutation_test = FALSE)
#'     distr <- run_permutation(data_t = data_t,
#'                              metadata = metadata,
#'                              LR_data = LR_keep,
#'                              iterations = iterations)
#'     distr <- cbind(distr, diff_noPerm)
#'     #return(distr)
#'     pvals <- rowSums(abs(distr[,1:iterations]) >= abs(distr[,(iterations+1)]))/iterations
#'     array_res[,,,7,1] <- array(pvals, dim = dim(array_res[,,,1,1]) )
#'   } else {
#'     stop("To do")
#'   }
#'   if(return_type == "array") {
#'     return(array_res)
#'   } else if(return_type == "data.table") {
#'     dt <- data.table::dcast.data.table(data.table::as.data.table(array_res),
#'                                        formula = V1 + V2 + V3  ~ V4 + V5 ,
#'                                        value.var = "value")
#'     names(dt)[names(dt) == "V1"] <- "LR_pair"
#'     names(dt)[names(dt) == "V2"] <- "Ligand_cell_type"
#'     names(dt)[names(dt) == "V3"] <- "Receptor_cell_type"
#'     return(dt)
#'   } else {
#'     stop("Return type not supported.")
#'   }
#'
#' }
#'
#' #' Title
#' #'
#' #' @param data_t x
#' #' @param metadata x
#' #' @param LR_data x
#' #' @param use_thr x
#' #' @param detection_thr x
#' #' @param condition_id x
#' #'
#' #' @return x
#' #' @export
#' #'
#' #' @examples
#' get_cci_score <- function(data_t,
#'                           metadata,
#'                           LR_data,
#'                           use_thr,
#'                           detection_thr,
#'                           condition_id
#' ) {
#'   if(is.null(condition_id)) {
#'     stop("To do later, only one condition")
#'   } else {
#'     conds <- unique(metadata$condition)
#'     if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
#'     ct_sorted <- sort(unique(metadata$cell_type))
#'     full_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 7, 2))
#'     dimnames(full_array) <- list(LR_data$SYMB_LR,
#'                                  ct_sorted,
#'                                  ct_sorted,
#'                                  c('LR_score',
#'                                    'LR_detection',
#'                                    'Ligand_expr',
#'                                    'Receptor_expr',
#'                                    'Ligand_detection',
#'                                    'Receptor_detection',
#'                                    'Pvalue'),
#'                                  conds)
#'     for(cond in conds) {
#'       meta <- metadata[metadata$condition == cond, ]
#'       data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
#'       data_average <- aggregate_cells(data = data_keep,
#'                                       group = meta$cell_type)
#'       data_average <- data_average[ct_sorted,]
#'       detec_rate <- aggregate_cells(data = 1*(data_keep > 0),
#'                                     group = meta$cell_type)
#'       detec_rate <- detec_rate[ct_sorted,]
#'       for(i in 1:nrow(LR_data)) {
#'         LR_scores <- outer(data_average[, LR_data[i,'GENESYMB_L']],
#'                            data_average[, LR_data[i,'GENESYMB_R']],
#'                            FUN = "+")/2
#'         LR_detec <- outer(detec_rate[, LR_data[i,'GENESYMB_L']],
#'                           detec_rate[, LR_data[i,'GENESYMB_R']],
#'                           FUN = is_detected,
#'                           detection_thr)
#'         if(use_thr) {
#'           LR_scores[LR_detec == FALSE] <- 0
#'         }
#'         full_array[i,,,1, cond] <- LR_scores
#'         full_array[i,,,2, cond] <- LR_detec
#'         full_array[i,,,3, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
#'         full_array[i,,,4, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
#'         full_array[i,,,5, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
#'         full_array[i,,,6, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
#'       }
#'     }
#'   }
#'   return(full_array)
#' }
#'
#'
#' #' Title
#' #'
#' #' @param data_t x
#' #' @param metadata x
#' #' @param LR_data x
#' #' @param permutation_test x
#' #'
#' #' @return x
#' #' @export
#' #'
#' #' @examples
#' get_diff_score_fast <- function(
#'   data_t,
#'   metadata,
#'   LR_data,
#'   permutation_test = TRUE
#' ) {
#'   if(permutation_test) {
#'     #metadata$cell_id <- sample(metadata$cell_id)
#'     for(x in unique(metadata$cell_type)) {
#'       metadata$condition[metadata$cell_type == x] <- sample(metadata$condition[metadata$cell_type == x])
#'     }
#'   }
#'   conds <- unique(metadata$condition)
#'   ct_sorted <- sort(unique(metadata$cell_type))
#'   temp_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 2))
#'   dimnames(temp_array) <- list(NULL, NULL, NULL, conds)
#'   for(cond in conds) {
#'     meta <- metadata[metadata$condition == cond, ]
#'     data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
#'     data_average <- aggregate_cells(data = data_keep,
#'                                     group = meta$cell_type)
#'     data_average <- data_average[ct_sorted,]
#'     for(i in 1:nrow(LR_data)) {
#'       #be careful Rfast::Outer returns the transpose of base::outer!!
#'       temp_array[i,,,cond] <- Rfast::transpose(Rfast::Outer(data_average[, LR_data[i,'GENESYMB_L']],
#'                                                             data_average[, LR_data[i,'GENESYMB_R']],
#'                                                             oper = "+")/2)
#'     }
#'   }
#'   return(as.vector(temp_array[,,,2]-temp_array[,,,1]))
#' }
#'
#' #' Title
#' #'
#' #' @param data_t x
#' #' @param metadata x
#' #' @param LR_data x
#' #' @param iterations x
#' #'
#' #' @return x
#' #' @export
#' #'
#' #' @examples
#' run_permutation <- function(data_t,
#'                             metadata,
#'                             LR_data,
#'                             iterations
#' ) {
#'   replicate(iterations, get_diff_score_fast(data_t = data_t,
#'                                             metadata = metadata,
#'                                             LR_data = LR_data,
#'                                             permutation_test = TRUE))
#' }
#'
#'
#' get_score_fast <- function(
#'   data_t,
#'   metadata,
#'   LR_data,
#'   permutation_test = TRUE
#' ) {
#'   if(permutation_test) {
#'     #metadata$cell_id <- sample(metadata$cell_id)
#'     for(x in unique(metadata$cell_type)) {
#'       metadata$condition[metadata$cell_type == x] <- sample(metadata$condition[metadata$cell_type == x])
#'     }
#'   }
#'   conds <- unique(metadata$condition)
#'   ct_sorted <- sort(unique(metadata$cell_type))
#'   temp_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 2))
#'   dimnames(temp_array) <- list(NULL, NULL, NULL, conds)
#'   for(cond in conds) {
#'     meta <- metadata[metadata$condition == cond, ]
#'     data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
#'     data_average <- aggregate_cells(data = data_keep,
#'                                     group = meta$cell_type)
#'     data_average <- data_average[ct_sorted,]
#'     for(i in 1:nrow(LR_data)) {
#'       temp_array[i,,,cond] <- compute_single_LR(ligand = LR_data[i,'GENESYMB_L'],
#'                                                 recetpor = LR_data[i,'GENESYMB_R'],
#'                                                 average_expression = data_average)
#'     }
#'   }
#'   return(as.vector(temp_array[,,,2]-temp_array[,,,1]))
#' }
#'
#' get_array_score_fast <- function(data_t,
#'                                  metadata,
#'                                  LR_data,
#'                                  is_condition
#' ) {
#'   conds <- unique(metadata$condition)
#'   ct_sorted <- sort(unique(metadata$cell_type))
#'   temp_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 2))
#'   dimnames(temp_array) <- list(NULL, NULL, NULL, conds)
#'   for(cond in conds) {
#'     meta <- metadata[metadata$condition == cond, ]
#'     data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
#'     data_average <- aggregate_cells(data = data_keep,
#'                                     group = meta$cell_type)
#'     data_average <- data_average[ct_sorted,]
#'     for(i in 1:nrow(LR_data)) {
#'       temp_array[i,,,cond] <- compute_LR(ligand = LR_data[i,'GENESYMB_L'],
#'                                          recetpor = LR_data[i,'GENESYMB_R'],
#'                                          average_expression = data_average)
#'     }
#'   }
#'   return(temp_array)
#' }
#'
#' #################
#'
#'
#' get_cci_score <- function(data_t,
#'                           metadata,
#'                           LR_data,
#'                           use_thr,
#'                           detection_thr,
#'                           condition_id
#' ) {
#'   if(is.null(condition_id)) {
#'     stop("To do later, only one condition")
#'   } else {
#'     conds <- unique(metadata$condition)
#'     if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
#'     ct_sorted <- sort(unique(metadata$cell_type))
#'     full_array <- array(data = 0, dim = c(nrow(LR_data), length(ct_sorted), length(ct_sorted), 7, 2))
#'     dimnames(full_array) <- list(LR_data$SYMB_LR,
#'                                  ct_sorted,
#'                                  ct_sorted,
#'                                  c('LR_score',
#'                                    'LR_detection',
#'                                    'Ligand_expr',
#'                                    'Receptor_expr',
#'                                    'Ligand_detection',
#'                                    'Receptor_detection',
#'                                    'Pvalue'),
#'                                  conds)
#'     for(cond in conds) {
#'       meta <- metadata[metadata$condition == cond, ]
#'       data_keep <- data_t[rownames(data_t) %in% meta$cell_id, ]
#'       data_average <- aggregate_cells(data = data_keep,
#'                                       group = meta$cell_type)
#'       data_average <- data_average[ct_sorted,]
#'       detec_rate <- aggregate_cells(data = 1*(data_keep > 0),
#'                                     group = meta$cell_type)
#'       detec_rate <- detec_rate[ct_sorted,]
#'       for(i in 1:nrow(LR_data)) {
#'         LR_scores <- outer(data_average[, LR_data[i,'GENESYMB_L']],
#'                            data_average[, LR_data[i,'GENESYMB_R']],
#'                            FUN = "+")/2
#'         LR_detec <- outer(detec_rate[, LR_data[i,'GENESYMB_L']],
#'                           detec_rate[, LR_data[i,'GENESYMB_R']],
#'                           FUN = is_detected,
#'                           detection_thr)
#'         if(use_thr) {
#'           LR_scores[LR_detec == FALSE] <- 0
#'         }
#'         full_array[i,,,1, cond] <- LR_scores
#'         full_array[i,,,2, cond] <- LR_detec
#'         full_array[i,,,3, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
#'         full_array[i,,,4, cond] <- matrix(data_average[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
#'         full_array[i,,,5, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_L']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = FALSE)
#'         full_array[i,,,6, cond] <- matrix(detec_rate[, LR_data[i,'GENESYMB_R']], ncol = nrow(data_average), nrow = nrow(data_average), byrow = TRUE)
#'       }
#'     }
#'   }
#'   return(full_array)
#' }
#'
#'
