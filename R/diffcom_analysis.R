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
#' @param statistical_analysis logical indicating if performing the permutation test for the significance of differential
#' expression between the conditions. Only considered when condidition_id is not NULL; default is TRUE.
#' @param specificity_analysis logical indicating if performing the permutation test for the specificity of each CCI on each condition;
#' default is FALSE.
#' @param iterations integer indicating the number of iterations during the permutation test.
#'
#' @return A data.table where each row is CCI. The columns vary in functions of the parameters used when calling the function.
#' It includes the CCI information and for each condition the scores, detection rates and possibly p-values.
#' @export
run_diffcom <- function(seurat_obj,
                        LR_data,
                        seurat_cell_type_id = "cell_ontology_class",
                        condition_id = NULL,
                        assay = "RNA",
                        slot = "data",
                        log_scale = TRUE,
                        min_cells = 5,
                        threshold = 0.1,
                        statistical_analysis = TRUE,
                        specificity_analysis = FALSE,
                        iterations = 1000
) {
  pp_seurat <- preprocess_seurat(seurat_obj = seurat_obj,
                                 assay = assay,
                                 slot = slot,
                                 log_scale = log_scale,
                                 convert_to_human = FALSE,
                                 return_type = "dense",
                                 seurat_cell_type_id = seurat_cell_type_id,
                                 condition_id = condition_id,
                                 min_cells = min_cells)
  pp_LR <- preprocess_LR(data = pp_seurat$data,
                         LR_data = LR_data)
  cci_analysis <- run_cci_analysis(expr_tr = t(pp_LR$data),
                                   metadata = pp_seurat$metadata,
                                   cell_types = pp_seurat$cell_types,
                                   LR_df = pp_LR$LR_df,
                                   threshold = threshold,
                                   condition_id = condition_id,
                                   statistical_analysis = statistical_analysis,
                                   specificity_analysis = specificity_analysis,
                                   iterations = iterations)
}

#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param cell_types x
#' @param LR_df x
#' @param threshold x
#' @param condition_id x
#' @param statistical_analysis x
#' @param specificity_analysis x
#' @param iterations x
#'
#' @return x
run_cci_analysis <- function(expr_tr,
                             metadata,
                             cell_types,
                             LR_df,
                             threshold,
                             condition_id,
                             statistical_analysis,
                             specificity_analysis,
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
    if(!specificity_analysis) {
      cci_dt <- data.table::dcast.data.table(data.table::as.data.table(array_noCond),
                                             formula = V1 + V2 + V3 ~ V4 ,
                                             value.var = "value")
    } else {
      stat_res <- run_stat_analysis(expr_tr = expr_tr,
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
        stat_res <- run_stat_analysis(expr_tr = expr_tr,
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
        stat_res <- run_stat_analysis(expr_tr = expr_tr,
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
        stat_res <- run_stat_analysis(expr_tr = expr_tr,
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
        #cci_dt$BH <- stats::p.adjust(cci_dt[[paste0("pval_diff_", cond1)]])
        #cci_dt$log2fc <- log2(cci_dt[[paste0("LR_score_", cond2)]]/cci_dt[[paste0("LR_score_", cond1)]])

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

#' Title
#'
#' @param averaged_expr x
#' @param LR_df x
#'
#' @return x
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

#' Compute LR scores for one LR-pair over all cell-types
#'
#' @param ligand character indicating the name of the ligand
#' @param receptor character indicating the name of the receptor
#' @param averaged_expr matrix with the mean expression of a each gene per cell-type
#'
#' @return matrix with the score of LR for each cell-type pair combination
compute_LR_score <- function(ligand,
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

#' Compute LR detection rate for one LR-pair over all cell-types
#'
#' @param ligand character indicating the name of the ligand
#' @param receptor character indicating the name of the receptor
#' @param detection_rate x
#' @param detection_thr x
#'
#' @return matrix with the detection rate of LR for each cell-type pair combination
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
#' @return x
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
run_stat_analysis <- function(expr_tr,
                              metadata,
                              cell_types,
                              LR_df,
                              cond1,
                              cond2,
                              iterations,
                              use_case
) {
  if(use_case == "no_cond") {
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
                                             cond2 = NULL,
                                             use_case = use_case),
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
                                               cond2 = cond2,
                                               use_case = use_case),
                            simplify = "array")
    if(use_case == "cond_spec") {
      distr_cond1 <- cbind(array_perm[,1,], cci_cond1)
      pvals_cond1 <- rowSums(distr_cond1[,1:iterations] >= distr_cond1[,(iterations+1)])/iterations
      distr_cond2 <- cbind(array_perm[,2,], cci_cond2)
      pvals_cond2 <- rowSums(distr_cond2[,1:iterations] >= distr_cond2[,(iterations+1)])/iterations
      return(list(pvals_cond1 = pvals_cond1,
                  pvals_cond2 = pvals_cond2))
    } else if(use_case == "cond_stat") {
      distr_diff <- cbind(array_perm, cci_cond2 - cci_cond1)
      pvals_diff <- rowSums(abs(distr_diff[,1:iterations]) >= abs(distr_diff[,(iterations+1)]))/iterations
      return(pvals_diff)
    } else if(use_case == "cond_stat_spec") {
      distr_diff <- cbind(array_perm[,1,], cci_cond2 - cci_cond1)
      pvals_diff <- rowSums(abs(distr_diff[,1:iterations]) >= abs(distr_diff[,(iterations+1)]))/iterations
      distr_cond1 <- cbind(array_perm[,2,], cci_cond1)
      pvals_cond1 <- rowSums(distr_cond1[,1:iterations] >= distr_cond1[,(iterations+1)])/iterations
      distr_cond2 <- cbind(array_perm[,3,], cci_cond2)
      pvals_cond2 <- rowSums(distr_cond2[,1:iterations] >= distr_cond2[,(iterations+1)])/iterations
      return(list(pvals_diff = pvals_diff,
                  pvals_cond1 = pvals_cond1,
                  pvals_cond2 = pvals_cond2))

    } else {
      stop("Case not supported in function run_stat_analysis.")
    }
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
#' @param use_case x
#'
#' @return x
run_stat_iteration <- function(expr_tr,
                               metadata,
                               cell_types,
                               LR_df,
                               cond1,
                               cond2,
                               use_case
) {
  if(use_case == "no_cond") {
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
  } else if(use_case == "cond_spec") {
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
    return(cbind(cci_permcond_cond1, cci_permct_cond2))
  } else if(use_case == "cond_stat") {
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
    return(diff_permcond)
  } else if(use_case == "cond_stat_spec") {
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
  } else {
    stop("Case  not supported in function run stat iteration")
  }
}
