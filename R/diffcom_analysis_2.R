#' Determine (differential) intercellular communication scores (with significance)
#'
#' Read a Seurat object and a data.frame of Ligand-Receptor interactions.
#' Compute the score of each cell-cell interaction (CCI) for one or two conditions on the cells.
#' There is the possibility to return a specicifity p-value for each CCI over the conditions.
#' When two conditions are considered, it computes the changes of each CCI with the option to return a p-value to assess
#' the significance of each change.
#' Statistical assessements are based on a permutation test (either permuting the condition or cluster labels.)
#'
#' @param seurat_obj Seurat object of several cell-types and human/mouse genes.
#' @param LR_data 3-column data.frame of Ligand-Receptor pairs,
#'  satisfying colnames(LR_data) == c("GENESYMB_L", "GENESYMB_R", "SYMB_LR").
#' @param seurat_cell_type_id character indicating the column specifying the cell-types of the cells
#' @param condition_id character indicating the column specifying the conditions on the cells. Set to NULL to run the analysis
#' on no condition; default is NULL.
#' @param assay character indicating the assay of the Seurat object where to pull the data from; default is "RNA".
#' @param slot character indicating the slot of the Seurat object where to pull the data from; default is "data".
#' @param log_scale logical indicating if using log-normalized data (TRUE) or normalized data (FALSE); default is "TRUE".
#' Only considered if slot == "data".
#' @param min_cells numeric indicating the minimal number of cells each cluster need to contain to be considered in the analysis.
#' @param threshold numeric indicating the percentage of cells that need to express a gene in a cluster for the gene to be
#' considered detected.
#' @param differential_analysis logical indicating if performing the permutation test for the significance of differential
#' expression between the conditions. Only considered when condidition_id is not NULL; default is TRUE.
#' @param specificity_analysis logical indicating if performing the permutation test for the specificity of each CCI on each condition;
#' default is FALSE.
#' @param iterations integer indicating the number of iterations during the permutation test.
#'
#' @return A data.table where each row is CCI. The columns vary in functions of the parameters used when calling the function.
#' It includes the CCI information and for each condition the scores, detection rates and possibly p-values.
#' @export
run_diffcom_2 <- function(seurat_obj,
                          LR_data,
                          seurat_cell_type_id = "cell_ontology_class",
                          condition_id = NULL,
                          assay = "RNA",
                          slot = "data",
                          log_scale = TRUE,
                          min_cells = 5,
                          threshold = 0.1,
                          differential_analysis = TRUE,
                          specificity_analysis = FALSE,
                          iterations = 1000
) {
  message("Preprocessing Seurat object.")
  pp_seurat <- preprocess_seurat(seurat_obj = seurat_obj,
                                 assay = assay,
                                 slot = slot,
                                 log_scale = log_scale,
                                 convert_to_human = FALSE,
                                 return_type = "dense",
                                 seurat_cell_type_id = seurat_cell_type_id,
                                 condition_id = condition_id,
                                 min_cells = min_cells)
  message("Preprocessing ligand-receptor pairs.")
  pp_LR <- preprocess_LR(data = pp_seurat$data,
                         LR_data = LR_data)
  message("Start CCI analysis.")
  LR_df <- pp_LR$LR_df
  LR_df$keep <- TRUE
  cci_analysis <- run_cci_analysis_2(expr_tr = t(pp_LR$data),
                                     metadata = pp_seurat$metadata,
                                     cell_types = pp_seurat$cell_types,
                                     LR_df = LR_df,
                                     threshold = threshold,
                                     condition_id = condition_id,
                                     differential_analysis = differential_analysis,
                                     specificity_analysis = specificity_analysis,
                                     iterations = iterations)
  #formatting output data.table
  if(is.null(condition_id)) {
    if(!specificity_analysis) {
      cci_analysis <- cci_analysis[, c(1,2,3,5,4,7,6,9,8)]
    } else {
      cci_analysis <- cci_analysis[, c(1,2,3,5,4,10,7,6,9,8)]
    }

  } else {
    if(!specificity_analysis) {
      if(!differential_analysis) {
        cci_analysis <- cci_analysis[, c(1,2,3,5,11,4,10,7,13,6,12,9,15,8,14)]
      } else {
        cci_analysis <- cci_analysis[, c(1,2,3,16,6,7,4,5,10,11,8,9,14,15,12,13)]
        colnames(cci_analysis)[[4]] <- "pval_diff"
      }
    } else {
      if(!differential_analysis) {
        cci_analysis <- cci_analysis[, c(1,2,3,6,7,4,5,16,17,10,11,8,9,14,15,12,13)]
      } else {
        cci_analysis <- cci_analysis[, c(1,2,3,16,6,7,4,5,18,19,10,11,8,9,14,15,12,13)]
        colnames(cci_analysis)[[4]] <- "pval_diff"
      }
    }
  }
  return(cci_analysis)
}

#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param cell_types x
#' @param LR_df x
#' @param threshold x
#' @param condition_id x
#' @param differential_analysis x
#' @param specificity_analysis x
#' @param iterations x
#'
#' @return x
run_cci_analysis_2 <- function(expr_tr,
                               metadata,
                               cell_types,
                               LR_df,
                               threshold,
                               condition_id,
                               differential_analysis,
                               specificity_analysis,
                               iterations
) {
  if(is.null(condition_id)) {
    message("Performing simple analysis without condition.")
    array_noCond <- run_simple_analysis_2(expr_tr = expr_tr,
                                          metadata = metadata,
                                          cell_types = cell_types,
                                          cond = NULL,
                                          LR_df = LR_df,
                                          threshold = threshold,
                                          compute_fast = FALSE)
    if(!specificity_analysis) {
      cci_dt <- data.table::dcast.data.table(data.table::as.data.table(array_noCond),
                                             formula = V1 + V2 + V3 ~ V4 ,
                                             value.var = "value")
    } else {
      message("Performing specificity analysis without condition.")
      LR_df$keep <- filter_detected_LR_2(array_dr1 = array_noCond[,,,"LR_detection"],
                                         array_dr2 = NULL,
                                         LR_df = LR_df)
      stat_res <- run_stat_analysis_2(expr_tr = expr_tr,
                                      metadata = metadata,
                                      cell_types = cell_types,
                                      LR_df = LR_df,
                                      cond1 = NULL,
                                      cond2 = NULL,
                                      iterations = iterations,
                                      use_case = "no_cond")
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
    message("Performing simple analysis on the two conditions.")
    array_cond1 <- run_simple_analysis_2(expr_tr = expr_tr,
                                         metadata = metadata,
                                         cell_types = cell_types,
                                         cond = cond1,
                                         LR_df = LR_df,
                                         threshold = threshold,
                                         compute_fast = FALSE)
    array_cond2 <- run_simple_analysis_2(expr_tr = expr_tr,
                                         metadata = metadata,
                                         cell_types = cell_types,
                                         cond = cond2,
                                         LR_df = LR_df,
                                         threshold = threshold,
                                         compute_fast = FALSE)
    LR_df$keep <- filter_detected_LR_2(array_dr1 = array_cond1[,,,"LR_detection"],
                                       array_dr2 = array_cond2[,,,"LR_detection"],
                                       LR_df = LR_df
    )
    if(!differential_analysis) {
      if(!specificity_analysis) {
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
        message("Performing specificity analysis on the two conditions.")
        stat_res <- run_stat_analysis_2(expr_tr = expr_tr,
                                        metadata = metadata,
                                        cell_types = cell_types,
                                        LR_df = LR_df,
                                        cond1 = cond1,
                                        cond2 = cond2,
                                        iterations = iterations,
                                        use_case = "cond_spec")
        #to keep the order of each cci the same
        big_array <- array(data = 0,
                           dim = c(dim(array_cond1)[1:3], dim(array_cond1)[4] + 1, 2))
        dnames <- dimnames(array_cond1)
        dimnames(big_array) <- list(dnames[[1]], dnames[[2]], dnames[[3]], c(dnames[[4]], "pval_spec"), c(cond1, cond2))
        big_array[,,,1:6,1] <- array_cond1
        big_array[,,,7,1] <- array(data = stat_res$pvals_cond1, dim = dim(array_cond1)[1:3])
        big_array[,,,1:6,2] <- array_cond2
        big_array[,,,7,2] <- array(data = stat_res$pvals_cond2, dim = dim(array_cond1)[1:3])
        cci_dt <- data.table::dcast(data.table::as.data.table(big_array),
                                    formula = V1 + V2 + V3 ~ V4 + V5,
                                    value.var = "value")
      }

    } else {
      if(!specificity_analysis) {
        message("Performing differential analysis between the two conditions.")
        stat_res <- run_stat_analysis_2(expr_tr = expr_tr,
                                        metadata = metadata,
                                        cell_types = cell_types,
                                        LR_df = LR_df,
                                        cond1 = cond1,
                                        cond2 = cond2,
                                        iterations = iterations,
                                        use_case = "cond_stat")
        #to keep the order of each cci the same
        big_array <- array(data = 0,
                           dim = c(dim(array_cond1)[1:3], dim(array_cond1)[4] + 1, 2))
        dnames <- dimnames(array_cond1)
        dimnames(big_array) <- list(dnames[[1]], dnames[[2]], dnames[[3]], c(dnames[[4]], "pval_diff"), c(cond1, cond2))
        big_array[,,,1:6,1] <- array_cond1
        big_array[,,,7,1] <- array(data = stat_res, dim = dim(array_cond1)[1:3])
        big_array[,,,1:6,2] <- array_cond2
        big_array[,,,7,2] <- array(data = stat_res, dim = dim(array_cond1)[1:3])
        cci_dt <- data.table::dcast(data.table::as.data.table(big_array),
                                    formula = V1 + V2 + V3 ~ V4 + V5,
                                    value.var = "value")

      } else {
        message("Performing differential analysis between the two conditions and specificity analysis.")
        stat_res <- run_stat_analysis_2(expr_tr = expr_tr,
                                        metadata = metadata,
                                        cell_types = cell_types,
                                        LR_df = LR_df,
                                        cond1 = cond1,
                                        cond2 = cond2,
                                        iterations = iterations,
                                        use_case = "cond_stat_spec")
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
    }
    names(cci_dt)[names(cci_dt) == "V1"] <- "LR_pair"
    names(cci_dt)[names(cci_dt) == "V2"] <- "Ligand_cell_type"
    names(cci_dt)[names(cci_dt) == "V3"] <- "Receptor_cell_type"

  }
  return(cci_dt)
}

#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param cell_types x
#' @param cond x
#' @param LR_df x
#' @param threshold x
#' @param compute_fast x
#'
#' @return x
run_simple_analysis_2 <- function(expr_tr,
                                  metadata,
                                  cell_types,
                                  cond,
                                  LR_df,
                                  threshold,
                                  compute_fast
) {
  averaged_expr <- aggregate_cells_slow(expr_tr = expr_tr,
                                   metadata = metadata,
                                   cell_types = cell_types,
                                   cond = cond)
  cci <- build_cci_array_2(averaged_expr = averaged_expr,
                           LR_df = LR_df)
  if(compute_fast) {
    return(cci)
  } else {
    detection_rate <- aggregate_cells_slow(expr_tr = 1*(expr_tr > 0),
                                      metadata = metadata,
                                      cell_types = cell_types,
                                      cond = cond)
    drate <- build_drate_array_2(detection_rate = detection_rate,
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

#' Indicates which LR has been detected in at least one CCI
#'
#' @param array_dr1 3D array for noCond or cond1
#' @param array_dr2 3D array for cond2 (NULL if noCond)
#' @param LR_df LR data.frame
#'
#' @return logical vector
filter_detected_LR_2 <- function(array_dr1,
                                 array_dr2,
                                 LR_df
) {
  if(is.null(array_dr2)) {
    sums <- sapply(1:nrow(LR_df), function(i) {
      sum(array_dr1[i,,])
    })
    res <- sums > 0
  } else {
    sums_cond1 <- sapply(1:nrow(LR_df), function(i) {
      sum(array_dr1[i,,])
    })
    sums_cond2 <- sapply(1:nrow(LR_df), function(i) {
      sum(array_dr2[i,,])
    })
    res <- (sums_cond1 + sums_cond2) > 0
  }
  message(paste0("Number of detected LR pairs: ", sum(res)))
  return(res)
}

#' Title
#'
#' @param averaged_expr x
#' @param LR_df x
#'
#' @return x
build_cci_array_2 <- function(averaged_expr,
                              LR_df
) {
  cci_array <- array(data = 0,
                     dim = c(nrow(LR_df), nrow(averaged_expr), nrow(averaged_expr)))
  for(i in 1:nrow(LR_df)) {
    if(LR_df[i, "keep"]) {
      cci_array[i, , ] <- compute_LR_score_2(ligand = LR_df[i, "ligand"],
                                             receptor = LR_df[i, "receptor"],
                                             averaged_expr= averaged_expr)
    } else {
      cci_array[i, , ] <- 0
    }

  }
  return(cci_array)
}

#' Compute LR scores for one LR-pair over all cell-types
#'
#' @param ligand character indicating the name of the ligand
#' @param receptor character indicating the name of the receptor
#' @param averaged_expr matrix with the mean expression of a each gene per cell-type
#'
#' @return matrix with the score of LR for each cell-type pair combination
compute_LR_score_2 <- function(ligand,
                               receptor,
                               averaged_expr
) {
  #be careful Rfast::Outer returns the transpose of base::outer!
  Rfast::transpose(Rfast::Outer(averaged_expr[, ligand],
                                averaged_expr[, receptor],
                                oper = "+")/2)
}

#' Title
#'
#' @param detection_rate x
#' @param detection_thr x
#' @param LR_df x
#'
#' @return x
build_drate_array_2 <- function(detection_rate,
                                detection_thr,
                                LR_df
) {
  drate_array <- array(data = 0,
                       dim = c(nrow(LR_df), nrow(detection_rate), nrow(detection_rate)))
  for(i in 1:nrow(LR_df)) {
    drate_array[i, , ] <- compute_LR_drate_2(ligand = LR_df[i, "ligand"] ,
                                             receptor = LR_df[i, "receptor"],
                                             detection_rate = detection_rate,
                                             detection_thr = detection_thr)
  }
  return(drate_array)
}

#' Compute LR detection rate for one LR-pair over all cell-types
#'
#' @param ligand character indicating the name of the ligand
#' @param receptor character indicating the name of the receptor
#' @param detection_rate x
#' @param detection_thr x
#'
#' @return matrix with the detection rate of LR for each cell-type pair combination
compute_LR_drate_2 <- function(ligand,
                               receptor,
                               detection_rate,
                               detection_thr
) {
  outer(detection_rate[, ligand],
        detection_rate[, receptor],
        FUN = is_detected_2,
        detection_thr)
}

#' Determine if two values are bigger than a threshold at the same time.
#'
#' @param x x
#' @param y x
#' @param thr x
#'
#' @return x
is_detected_2 <- Vectorize(function(x,
                                  y,
                                  thr
) {
  if(x >= thr & y >= thr) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

compute_pvalue_2 <- function(expr_tr,
                           metadata,
                           cell_types,
                           LR_df,
                           iterations,
                           sample_matrix_spec,
                           sample_matrix_diff,
                           cond1,
                           cond2,
                           use_case
) {
  n_ct <- length(cell_types)
  if(use_case == "no_cond") {
    pvals <- vector()
    for(i in 1:nrow(LR_df)) {
      if(LR_df[i, "keep"]) {
        ligand <- LR_df[i, "ligand"]
        receptor <- LR_df[i, "receptor"]
        expr_tr_sub <- expr_tr[, c(ligand, receptor)]
        distr <- array(data = 0,
                       dim = c(n_ct*n_ct, iterations + 1))
        averaged_expr <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                         metadata = metadata,
                                         cell_types = cell_types,
                                         cond = NULL)
        distr[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                              receptor = receptor,
                                                              averaged_expr = averaged_expr))
        for(i in 1:iterations) {
          md <- metadata
          md$cell_type <- sample_matrix_spec[i, ]
          averaged_expr <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                           metadata = md,
                                           cell_types = cell_types,
                                           cond = NULL)
          distr[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                     receptor = receptor,
                                                     averaged_expr = averaged_expr))
        }
        pvals <- c(pvals, rowSums(distr[,1:iterations] >= distr[,(iterations+1)])/iterations)
      } else {
        pvals <- c(pvals, rep(1, n_ct*n_ct))
      }
    }
    pvals <- matrix(data = pvals, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    return(as.vector(t(pvals)))
  } else if (use_case == "cond_spec") {
    pvals_cond1 <- vector()
    pvals_cond2 <- vector()
    for(i in 1:nrow(LR_df)) {
      if(LR_df[i, "keep"]) {
        ligand <- LR_df[i, "ligand"]
        receptor <- LR_df[i, "receptor"]
        expr_tr_sub <- expr_tr[, c(ligand, receptor)]
        distr_cond1 <- array(data = 0,
                             dim = c(n_ct*n_ct, iterations + 1))
        distr_cond2 <- array(data = 0,
                             dim = c(n_ct*n_ct, iterations + 1))
        averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond1)
        averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond2)
        distr_cond1[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                    receptor = receptor,
                                                                    averaged_expr = averaged_expr_cond1))
        distr_cond2[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                    receptor = receptor,
                                                                    averaged_expr = averaged_expr_cond2))
        for(i in 1:iterations) {
          md <- metadata
          md$cell_type <- sample_matrix_spec[i, ]
          averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                 metadata = md,
                                                 cell_types = cell_types,
                                                 cond = cond1)
          averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                 metadata = md,
                                                 cell_types = cell_types,
                                                 cond = cond2)
          distr_cond1[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                           receptor = receptor,
                                                           averaged_expr = averaged_expr_cond1))
          distr_cond2[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                           receptor = receptor,
                                                           averaged_expr = averaged_expr_cond2))
        }
        pvals_cond1 <- c(pvals_cond1, rowSums(distr_cond1[,1:iterations] >= distr_cond1[,(iterations+1)])/iterations)
        pvals_cond2 <- c(pvals_cond2, rowSums(distr_cond2[,1:iterations] >= distr_cond2[,(iterations+1)])/iterations)
      } else {
        pvals_cond1 <- c(pvals_cond1, rep(1, n_ct*n_ct))
        pvals_cond2 <- c(pvals_cond2, rep(1, n_ct*n_ct))
      }
    }

    pvals_cond1 <- matrix(data = pvals_cond1, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    pvals_cond2 <- matrix(data = pvals_cond2, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    return(list(pvals_cond1 = as.vector(t(pvals_cond1)),
                pvals_cond2 = as.vector(t(pvals_cond2))))
  } else if(use_case == "cond_stat") {
    pvals <- vector()
    for(i in 1:nrow(LR_df)) {
      if(LR_df[i, "keep"]) {
        ligand <- LR_df[i, "ligand"]
        receptor <- LR_df[i, "receptor"]
        expr_tr_sub <- expr_tr[, c(ligand, receptor)]
        distr <- array(data = 0,
                       dim = c(n_ct*n_ct, iterations + 1))
        averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond1)
        averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond2)
        distr[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                              receptor = receptor,
                                                              averaged_expr = averaged_expr_cond2) -
                                             compute_LR_score_2(ligand = ligand,
                                                                receptor = receptor,
                                                                averaged_expr = averaged_expr_cond1))
        for(i in 1:iterations) {
          md <- metadata
          md$condition <- sample_matrix_diff[i, ]
          averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                           metadata = md,
                                           cell_types = cell_types,
                                           cond = cond1)
          averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                 metadata = md,
                                                 cell_types = cell_types,
                                                 cond = cond2)
          distr[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                     receptor = receptor,
                                                     averaged_expr = averaged_expr_cond2) -
                                    compute_LR_score_2(ligand = ligand,
                                                       receptor = receptor,
                                                       averaged_expr = averaged_expr_cond1))
        }
        pvals <- c(pvals, rowSums(abs(distr[,1:iterations]) >= abs(distr[,(iterations+1)]))/iterations)
      } else {
        pvals <- c(pvals, rep(1, n_ct*n_ct))
      }
    }
    pvals <- matrix(data = pvals, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    return(as.vector(t(pvals)))
  } else if (use_case == "cond_stat_spec") {
    pvals_diff <- vector()
    pvals_cond1 <- vector()
    pvals_cond2 <- vector()
    for(i in 1:nrow(LR_df)) {
      if(LR_df[i, "keep"]) {
        ligand <- LR_df[i, "ligand"]
        receptor <- LR_df[i, "receptor"]
        expr_tr_sub <- expr_tr[, c(ligand, receptor)]
        distr_diff <- array(data = 0,
                       dim = c(n_ct*n_ct, iterations + 1))
        distr_cond1 <- array(data = 0,
                             dim = c(n_ct*n_ct, iterations + 1))
        distr_cond2 <- array(data = 0,
                             dim = c(n_ct*n_ct, iterations + 1))
        averaged_expr_diff_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond1)
        averaged_expr_diff_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond2)
        averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond1)
        averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                               metadata = metadata,
                                               cell_types = cell_types,
                                               cond = cond2)

        distr_diff[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                              receptor = receptor,
                                                              averaged_expr = averaged_expr_diff_cond2) -
                                             compute_LR_score_2(ligand = ligand,
                                                                receptor = receptor,
                                                                averaged_expr = averaged_expr_diff_cond1))
        distr_cond1[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                    receptor = receptor,
                                                                    averaged_expr = averaged_expr_cond1))
        distr_cond2[, (iterations + 1)] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                    receptor = receptor,
                                                                    averaged_expr = averaged_expr_cond2))
        for(i in 1:iterations) {
          md_diff <- metadata
          md_diff$condition <- sample_matrix_diff[i, ]
          md_cond <- metadata
          md_cond$cell_type <- sample_matrix_spec[i, ]
          averaged_expr_diff_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                      metadata = md_diff,
                                                      cell_types = cell_types,
                                                      cond = cond1)
          averaged_expr_diff_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                      metadata = md_diff,
                                                      cell_types = cell_types,
                                                      cond = cond2)
          averaged_expr_cond1 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                 metadata = md_cond,
                                                 cell_types = cell_types,
                                                 cond = cond1)
          averaged_expr_cond2 <- aggregate_cells_slow(expr_tr = expr_tr_sub,
                                                 metadata = md_cond,
                                                 cell_types = cell_types,
                                                 cond = cond2)
          distr_diff[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                     receptor = receptor,
                                                                     averaged_expr = averaged_expr_diff_cond2) -
                                                    compute_LR_score_2(ligand = ligand,
                                                                       receptor = receptor,
                                                                       averaged_expr = averaged_expr_diff_cond1))
          distr_cond1[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                      receptor = receptor,
                                                                      averaged_expr = averaged_expr_cond1))
          distr_cond2[, i] <- as.vector(compute_LR_score_2(ligand = ligand,
                                                                      receptor = receptor,
                                                                      averaged_expr = averaged_expr_cond2))
        }
        pvals_diff <- c(pvals_diff, rowSums(abs(distr_diff[,1:iterations]) >= abs(distr_diff[,(iterations+1)]))/iterations)
        pvals_cond1 <- c(pvals_cond1, rowSums(distr_cond1[,1:iterations] >= distr_cond1[,(iterations+1)])/iterations)
        pvals_cond2 <- c(pvals_cond2, rowSums(distr_cond2[,1:iterations] >= distr_cond2[,(iterations+1)])/iterations)
      } else {
        pvals_diff <- c(pvals_diff, rep(1, n_ct*n_ct))
        pvals_cond1 <- c(pvals_cond1, rep(1, n_ct*n_ct))
        pvals_cond2 <- c(pvals_cond2, rep(1, n_ct*n_ct))
      }
    }
    pvals_diff <- matrix(data = pvals_diff, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    pvals_cond1 <- matrix(data = pvals_cond1, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    pvals_cond2 <- matrix(data = pvals_cond2, nrow = n_ct*n_ct, ncol = nrow(LR_df), byrow = FALSE)
    return(list(pvals_diff = as.vector(t(pvals_diff)),
                pvals_cond1 = as.vector(t(pvals_cond1)),
                pvals_cond2 = as.vector(t(pvals_cond2))))

  } else {
    stop("Case not supported.")
  }
}


#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param cell_types x
#' @param LR_df x
#' @param cond1 x
#' @param cond2 x
#' @param iterations x
#' @param use_case x
#'
#' @return x
run_stat_analysis_2 <- function(expr_tr,
                                metadata,
                                cell_types,
                                LR_df,
                                cond1,
                                cond2,
                                iterations,
                                use_case
) {
  if(use_case == "no_cond") {
    ct <- metadata$cell_type
    sample_matrix <- permute::shuffleSet(seq_along(ct),
                                         nset = iterations)
    sample_matrix[] <- ct[sample_matrix]
    pvals <- compute_pvalue_2(expr_tr = expr_tr,
                            metadata = metadata,
                            cell_types = cell_types,
                            LR_df = LR_df,
                            iterations = iterations,
                            sample_matrix_spec = sample_matrix,
                            sample_matrix_diff = NULL,
                            cond1 = NULL,
                            cond2 = NULL,
                            use_case = "no_cond")
    return(pvals)
  } else if (use_case == "cond_spec") {
    ct_cond <- metadata$cell_type
    ctrl <- permute::how(blocks = as.factor(metadata$condition))
    sample_matrix <- permute::shuffleSet(seq_along(ct_cond),
                                         nset = iterations,
                                         control = ctrl)
    sample_matrix[] <- ct_cond[sample_matrix]
    ##
    #sample_matrix <- array(data = 0, dim = c(iterations, length(ct_cond)))
    #for(i in 1:iterations) {
    #  meta_ct <- metadata
    #  meta_ct$cell_type[meta_ct$condition == cond1] <- sample(meta_ct$cell_type[meta_ct$condition == cond1])
    #  meta_ct$cell_type[meta_ct$condition == cond2] <- sample(meta_ct$cell_type[meta_ct$condition == cond2])
    #  sample_matrix[i,] <- meta_ct$cell_type
    #}
    ##
    pvals <- compute_pvalue_2(expr_tr = expr_tr,
                            metadata = metadata,
                            cell_types = cell_types,
                            LR_df = LR_df,
                            iterations = iterations,
                            sample_matrix_spec = sample_matrix,
                            sample_matrix_diff = NULL,
                            cond1 = cond1,
                            cond2 = cond2,
                            use_case = "cond_spec")
    return(list(pvals_cond1 = pvals$pvals_cond1,
                pvals_cond2 = pvals$pvals_cond2))
  } else if(use_case == "cond_stat") {
    cond_diff <- metadata$condition
    ctrl <- permute::how(blocks = as.factor(metadata$cell_type))
    sample_matrix <- permute::shuffleSet(seq_along(cond_diff),
                                         nset = iterations,
                                         control = ctrl)
    sample_matrix[] <- cond_diff[sample_matrix]
    pvals <- compute_pvalue_2(expr_tr = expr_tr,
                            metadata = metadata,
                            cell_types = cell_types,
                            LR_df = LR_df,
                            iterations = iterations,
                            sample_matrix_spec = NULL,
                            sample_matrix_diff = sample_matrix,
                            cond1 = cond1,
                            cond2 = cond2,
                            use_case = "cond_stat")
    return(pvals)
  } else if(use_case == "cond_stat_spec") {
    cond_diff <- metadata$condition
    ctrl_diff <- permute::how(blocks = as.factor(metadata$cell_type))
    sample_matrix_diff <- permute::shuffleSet(seq_along(cond_diff),
                                              nset = iterations,
                                              control = ctrl_diff)
    sample_matrix_diff[] <- cond_diff[sample_matrix_diff]
    ct_cond <- metadata$cell_type
    ctrl_cond <- permute::how(blocks = as.factor(metadata$condition))
    sample_matrix_cond <- permute::shuffleSet(seq_along(ct_cond),
                                              nset = iterations,
                                              control = ctrl_cond)
    sample_matrix_cond[] <- ct_cond[sample_matrix_cond]
    pvals <- compute_pvalue_2(expr_tr = expr_tr,
                            metadata = metadata,
                            cell_types = cell_types,
                            LR_df = LR_df,
                            iterations = iterations,
                            sample_matrix_spec = sample_matrix_cond,
                            sample_matrix_diff = sample_matrix_diff,
                            cond1 = cond1,
                            cond2 = cond2,
                            use_case = "cond_stat_spec")

    return(list(pvals_diff = pvals$pvals_diff,
                pvals_cond1 = pvals$pvals_cond1,
                pvals_cond2 = pvals$pvals_cond2))

  } else {
    stop("Case not supported in function run_stat_analysis_2.")
  }
}
