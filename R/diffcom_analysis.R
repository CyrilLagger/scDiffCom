#' Determine (differential) intercellular communication scores (with significance)
#'
#' Read a Seurat object and a data.frame of Ligand-Receptor interactions.
#' Compute the score of each cell-cell interaction (CCI) for one or two conditions on the cells.
#' There is the possibility to return a specicifity p-value for each CCI over the conditions.
#' When two conditions are considered, it computes the changes of each CCI with the option to return a p-value to assess
#' the significance of each change.
#' Statistical assessements are based on a permutation test (either permuting the condition or cluster labels.)
#'
#' @param seurat_object Seurat object of several cell-types and human/mouse genes.
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
#' @param permutation_analysis logical indicating if performing permutation test; default is TRUE.
#' @param iterations integer indicating the number of iterations during the permutation test.
#' @param one_sided logical indicating if computing differential p-values from a one-sided or two-sided test.
#' @param return_distr logical indicating if returning the matrix with the distributions obtained from the permutation test.
#'
#' @return A data.table where each row is CCI. The columns vary in functions of the parameters used when calling the function.
#' It includes the CCI information and for each condition the scores, detection rates and possibly p-values.
#' @export
run_diffcom <- function(
  seurat_object,
  LR_data,
  seurat_cell_type_id = "cell_ontology_class",
  condition_id = NULL,
  assay = "RNA",
  slot = "data",
  log_scale = TRUE,
  min_cells = 5,
  threshold = 0.1,
  permutation_analysis = TRUE,
  one_sided = FALSE,
  iterations = 1000,
  return_distr = FALSE
) {
  message("Preprocessing Seurat object.")
  pp_seurat <- preprocess_seurat(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    convert_to_human = FALSE,
    return_type = "dense",
    seurat_cell_type_id = seurat_cell_type_id,
    condition_id = condition_id,
    min_cells = min_cells
  )
  message("Preprocessing ligand-receptor pairs.")
  pp_LR <- preprocess_LR(
    data = pp_seurat$data,
    LR_data = LR_data
  )
  message("Start CCI analysis.")
  template_cci_dt <- prepare_template_cci(
    LR_df = pp_LR$LR_df,
    cell_types = pp_seurat$cell_types
  )
  expr_tr <- t(pp_LR$data)
  metadata <- data.table::setDT(pp_seurat$metadata)
  if(is.null(condition_id)) {
    cond1 <- NULL
    cond2 <- NULL
    message("Performing simple analysis without condition.")
  } else {
    conds <- unique(metadata$condition)
    if (length(conds) != 2)
      stop("Wrong number of groups in cell-type conditions (expected 2).")
    cond1 <- conds[[1]]
    cond2 <- conds[[2]]
    message("Performing simple analysis on the two conditions.")
  }
  cci_dt_simple <- run_simple_analysis(
    expr_tr = expr_tr,
    metadata = metadata,
    template_cci_dt = template_cci_dt,
    cond1 = cond1,
    cond2 = cond2,
    threshold = threshold,
    compute_fast = FALSE
  )
  if (!permutation_analysis) {
    if(return_distr) {
      stop("No permutation distribution to return!")
    }
    cci_analysis <- cci_dt_simple
  } else {
    cci_analysis <- run_stat_analysis(
      cci_dt_simple = cci_dt_simple,
      expr_tr = expr_tr,
      metadata = metadata,
      cond1 = cond1,
      cond2 = cond2,
      iterations = iterations,
      one_sided = one_sided,
      return_distr = return_distr
    )
  }
  if(return_distr) {
    return(cci_analysis)
  }
  cci_analysis <- add_cell_number(
    cci_dt = cci_analysis,
    cond1 = cond1,
    cond2 = cond2,
    metadata = metadata
  )
  cci_analysis <- clean_colnames(
    cci_dt = cci_analysis,
    cond1 = cond1,
    cond2 = cond2,
    permutation_analysis = permutation_analysis
  )
  return(cci_analysis)
}

#' Create a template data.table with all the CCIs
#'
#' @param LR_df data.frame of LR CCI
#' @param cell_types character of the cell_types present in the Seurat object
#'
#' @return a data.table
prepare_template_cci <- function(
  LR_df,
  cell_types
) {
  setDT(LR_df)
  template <- data.table::CJ(
    L_CELLTYPE = cell_types,
    R_CELLTYPE = cell_types,
    LR_GENES = LR_df$LR_GENES
  )
  template <- data.table::merge.data.table(
    x = template,
    y = LR_df,
    by.x = "LR_GENES",
    by.y = "LR_GENES",
    all.x = TRUE,
    sort = FALSE
  )
}

#' Compute the score and detection rate of each CCI
#'
#' @param expr_tr matrix with the (transposed) expression data from the Seurat object
#' @param metadata data.table with the relevant metadata from the Seurat object
#' @param template_cci_dt data.table with a template of all CCIs
#' @param cond1 character indicating the name of the first condition, or NULL if no condition
#' @param cond2 character indicating the name of the first condition, or NULL if no condition
#' @param threshold numeric indicating the percentage of cells that need to express a gene in a cluster for the gene to be
#' considered detected.
#' @param compute_fast logical indicating if doing a fast computation for the permutation test or a full computation with
#' detection rate.
#'
#' @return a data.table or a (list of) numeric vector(s)
run_simple_analysis <- function(
  expr_tr,
  metadata,
  template_cci_dt,
  cond1,
  cond2,
  threshold,
  compute_fast
) {
  if(is.null(cond1) | is.null(cond2)) {
    is_cond = FALSE
  } else {
    is_cond = TRUE
  }
  averaged_expr <- aggregate_cells(
    expr_tr = expr_tr,
    metadata = metadata,
    is_cond = is_cond
  )
  cci_dt <- build_cci_drate_dt(
    averaged_expr = averaged_expr,
    template_cci_dt = template_cci_dt,
    cond1 = cond1,
    cond2 = cond2,
    detection_thr = threshold,
    cci_or_drate = "cci"
  )
  if(compute_fast) {
    if(!is_cond) {
      return(cci_dt[["LR_SCORE"]])
    } else {
      return(list(
        cond1 = cci_dt[[paste0("LR_SCORE_", cond1)]],
        cond2 = cci_dt[[paste0("LR_SCORE_", cond2)]]
      )
      )
    }
  }
  detection_rate <- aggregate_cells(
    expr_tr = 1 * (expr_tr > 0),
    metadata = metadata,
    is_cond = is_cond
  )
  drate_dt <- build_cci_drate_dt(
    averaged_expr = detection_rate,
    template_cci_dt = template_cci_dt,
    cond1 = cond1,
    cond2 = cond2,
    detection_thr = threshold,
    cci_or_drate = "drate"
  )
  dt <- data.table::merge.data.table(
    x = cci_dt,
    y = drate_dt,
    by = c("LR_GENES", "L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE"),
    sort = FALSE
  )
  return(dt)
}

#' Title
#'
#' @param expr_string x
#' @param cond1 x
#' @param cond2 x
#'
#' @return x
split_cond_string <- function(
  expr_string,
  cond1,
  cond2
) {
  cond <- ifelse(
    substr(
      x = expr_string,
      start = 1,
      stop = nchar(cond1)) == cond1,
    cond1,
    cond2)
  cell_type <- ifelse(
    substr(
      x = expr_string,
      start = 1,
      stop = nchar(cond1)) == cond1,
    substr(
      x = expr_string,
      start = nchar(cond1)+2,
      stop = nchar(expr_string)),
    substr(
      x = expr_string,
      start = nchar(cond2)+2,
      stop = nchar(expr_string)))
  return(list(cond = cond, cell_type = cell_type))
}

#' Title
#'
#' @param averaged_expr x
#' @param template_cci_dt x
#' @param cond1 x
#' @param cond2 x
#' @param detection_thr x
#' @param cci_or_drate x
#'
#' @return x
#' @import data.table
build_cci_drate_dt <- function(
  averaged_expr,
  template_cci_dt,
  cond1,
  cond2,
  detection_thr,
  cci_or_drate
) {
  L_GENE <- L_CELLTYPE <- L_EXPRESSION <-  L_DETECTED <-
    R_GENE <- R_CELLTYPE <- R_EXPRESSION <-  R_DETECTED <-
    CONDITION_CELLTYPE <- LR_SCORE <- LR_DETECTED <- GENE <- CELLTYPE <- NULL
  if(is.null(cond1) | is.null(cond2)) {
    dt <- data.table::as.data.table(
      x = averaged_expr,
      keep.rownames = "CELLTYPE",
      sorted = FALSE
    )
    long_dt <- data.table::melt.data.table(
      data = dt,
      id.vars = "CELLTYPE",
      variable.name = "GENE",
      value.name = "EXPRESSION"
    )
    col_long_dt <- c("GENE", "CELLTYPE")
    setkeyv(
      x = long_dt,
      cols = col_long_dt
    )
    col_template_cci_dt <- c("L_GENE", "L_CELLTYPE", "R_GENE", "R_CELLTYPE")
    setkeyv(
      x = template_cci_dt,
      cols = col_template_cci_dt
    )
    full_dt <- long_dt[
      long_dt[
        template_cci_dt,
        on = c("GENE==L_GENE", "CELLTYPE==L_CELLTYPE")],
      on = c("GENE==R_GENE", "CELLTYPE==R_CELLTYPE") ]
    if(cci_or_drate == "cci") {
      data.table::setnames(
        x = full_dt,
        old = c("CELLTYPE", "GENE", "EXPRESSION",
                "i.EXPRESSION", "i.CELLTYPE", "i.GENE"),
        new = c("R_CELLTYPE", "R_GENE", "R_EXPRESSION",
                "L_EXPRESSION", "L_CELLTYPE", "L_GENE")
      )
      full_dt[, LR_SCORE := (L_EXPRESSION + R_EXPRESSION) / 2]
    } else if(cci_or_drate == "drate") {
      data.table::setnames(
        x = full_dt,
        old = c("CELLTYPE", "GENE", "EXPRESSION",
                "i.EXPRESSION", "i.CELLTYPE", "i.GENE"),
        new = c("R_CELLTYPE", "R_GENE", "R_DETECTED",
                "L_DETECTED", "L_CELLTYPE", "L_GENE"))
      full_dt[, LR_DETECTED := is_detected(L_DETECTED, R_DETECTED, detection_thr)]
    } else {
      stop("Error in build_cci_drate_dt.")
    }
  } else {
    dt <- data.table::as.data.table(
      x = averaged_expr,
      keep.rownames = "CONDITION_CELLTYPE",
      sorted = FALSE
    )
    dt[, c("CONDITION", "CELLTYPE") := split_cond_string(CONDITION_CELLTYPE, cond1, cond2)]
    dt[, CONDITION_CELLTYPE := NULL]
    long_dt <- data.table::melt.data.table(
      data = dt,
      id.vars = c("CELLTYPE", "CONDITION"),
      variable.name = "GENE",
      value.name = "EXPRESSION"
    )
    long_dt <- data.table::dcast.data.table(
      data = long_dt,
      formula = CELLTYPE + GENE ~ CONDITION,
      value.var = "EXPRESSION")
    setkey(long_dt, GENE, CELLTYPE)
    setkey(template_cci_dt, L_GENE, L_CELLTYPE, R_GENE, R_CELLTYPE)
    full_dt <- long_dt[
      long_dt[
        template_cci_dt,
        on = c("GENE==L_GENE", "CELLTYPE==L_CELLTYPE")],
      on = c("GENE==R_GENE", "CELLTYPE==R_CELLTYPE") ]
    if(cci_or_drate == "cci") {
      data.table::setnames(
        x = full_dt,
        old = c("CELLTYPE", "GENE", cond1, cond2,
                "i.CELLTYPE", "i.GENE", paste0("i.", cond1), paste0("i.", cond2)),
        new = c("R_CELLTYPE", "R_GENE", paste0("R_EXPRESSION_", cond1), paste0("R_EXPRESSION_", cond2),
                "L_CELLTYPE", "L_GENE", paste0("L_EXPRESSION_", cond1), paste0("L_EXPRESSION_", cond2))
      )
      full_dt[, paste0("LR_SCORE_", c(cond1, cond2)) :=
                .((get(paste0("L_EXPRESSION_", cond1)) + get(paste0("R_EXPRESSION_", cond1))) / 2,
                  (get(paste0("L_EXPRESSION_", cond2)) + get(paste0("R_EXPRESSION_", cond2))) / 2)]
    } else if(cci_or_drate == "drate") {
      data.table::setnames(
        x = full_dt,
        old = c("CELLTYPE", "GENE", cond1, cond2,
                "i.CELLTYPE", "i.GENE", paste0("i.", cond1), paste0("i.", cond2)),
        new = c("R_CELLTYPE", "R_GENE", paste0("R_DETECTED_", cond1), paste0("R_DETECTED_", cond2),
                "L_CELLTYPE", "L_GENE", paste0("L_DETECTED_", cond1), paste0("L_DETECTED_", cond2))
      )
      full_dt[, paste0("LR_DETECTED_", c(cond1, cond2)) :=
                .(is_detected(get(paste0("L_DETECTED_", cond1)), get(paste0("R_DETECTED_", cond1)), detection_thr),
                  is_detected(get(paste0("L_DETECTED_", cond2)), get(paste0("R_DETECTED_", cond2)), detection_thr))]
    } else {
      stop("Error in build_cci_drate_dt.")
    }
  }
  return(full_dt)
}

#' Title
#'
#' @param cci_dt_simple x
#' @param expr_tr x
#' @param metadata x
#' @param cond1 x
#' @param cond2 x
#' @param iterations x
#' @param one_sided logical indicating if computing differential p-values from a one-sided or two-sided test.
#' @param return_distr x
#'
#' @return x
run_stat_analysis <- function(
  cci_dt_simple,
  expr_tr,
  metadata,
  cond1,
  cond2,
  iterations,
  one_sided,
  return_distr = FALSE
) {
  L_GENE <- L_CELLTYPE <- R_GENE <- R_CELLTYPE <-
    LR_GENES <- LR_DETECTED <-
    BH_PVAL <- BH_PVAL_DIFF <- PVAL <- PVAL_DIFF <- NULL
  if(is.null(cond1) | is.null(cond2)) {
    message("Performing permutation analysis without conditions.")
    sub_template_cci_dt <- cci_dt_simple[LR_DETECTED == TRUE]
  } else {
    message("Performing permutation analysis on the two conditions.")
    sub_template_cci_dt <- cci_dt_simple[get(paste0("LR_DETECTED_", cond1)) == TRUE | get(paste0("LR_DETECTED_", cond2)) == TRUE]
  }
  sub_expr_tr <- expr_tr[, colnames(expr_tr) %in%
                           unique(c(sub_template_cci_dt[["L_GENE"]], sub_template_cci_dt[["R_GENE"]]))]
  message(paste0("Initial number of genes: ",
                 ncol(expr_tr),
                 ". Number of detected genes for permutation test: ",
                 ncol(sub_expr_tr)))
  cci_perm <- replicate(
    n = iterations,
    expr = run_stat_iteration(
      expr_tr = sub_expr_tr,
      metadata = metadata,
      template_cci_dt = sub_template_cci_dt[, .(LR_GENES, L_CELLTYPE, R_CELLTYPE, L_GENE, R_GENE)],
      cond1 = cond1,
      cond2 = cond2
    ),
    simplify = "array"
  )
  if(is.null(cond1) | is.null(cond2)) {
    distr <- cbind(cci_perm, sub_template_cci_dt[["LR_SCORE"]])
    if (return_distr) {
      return(distr)
    }
    pvals <- rowSums(distr[, 1:iterations] >= distr[, (iterations + 1)]) / iterations
    sub_template_cci_dt[, PVAL := pvals]
    sub_template_cci_dt[, BH_PVAL := stats::p.adjust(p = pvals, method = "BH")]
    sub_template_cci_dt <- sub_template_cci_dt[, list(LR_GENES, L_CELLTYPE, R_CELLTYPE, L_GENE, R_GENE,
                                                      PVAL, BH_PVAL)]
  } else {
    distr_diff <- cbind(cci_perm[, 1, ], sub_template_cci_dt[[paste0("LR_SCORE_", cond2)]]
                        - sub_template_cci_dt[[paste0("LR_SCORE_", cond1)]])
    distr_cond1 <- cbind(cci_perm[, 2, ], sub_template_cci_dt[[paste0("LR_SCORE_", cond1)]])
    distr_cond2 <- cbind(cci_perm[, 3, ], sub_template_cci_dt[[paste0("LR_SCORE_", cond2)]])
    if (return_distr) {
      return(
        list(
          distr_diff = distr_diff,
          distr_cond1 = distr_cond1,
          distr_cond2 = distr_cond2
        )
      )
    }
    if(one_sided) {
      pvals_diff_cond1 <- rowSums(distr_diff[, 1:iterations] <= distr_diff[, (iterations + 1)]) / iterations
      pvals_diff_cond2 <- rowSums(distr_diff[, 1:iterations] >= distr_diff[, (iterations + 1)]) / iterations
      pvals_diff <- pmin(pvals_diff_cond1, pvals_diff_cond2)
      BH_pvals_diff_cond1 <- stats::p.adjust(p = pvals_diff_cond1, method = "BH")
      BH_pvals_diff_cond2 <- stats::p.adjust(p = pvals_diff_cond2, method = "BH")
      BH_pvals_diff <- pmin(BH_pvals_diff_cond1, BH_pvals_diff_cond2)
    } else {
      pvals_diff <- rowSums(abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations + 1)])) / iterations
      BH_pvals_diff <- stats::p.adjust(p = pvals_diff, method = "BH")
    }
    pvals_cond1 <- rowSums(distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]) / iterations
    pvals_cond2 <- rowSums(distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]) / iterations
    sub_template_cci_dt[, paste0("PVAL_", c(cond1, cond2)) := list(pvals_cond1, pvals_cond2)]
    sub_template_cci_dt[, paste0("BH_PVAL_", c(cond1, cond2)) := list(stats::p.adjust(
      p = pvals_cond1,
      method = "BH"
    ),
    stats::p.adjust(
      p = pvals_cond2,
      method = "BH"
    )
    )]
    sub_template_cci_dt[, PVAL_DIFF := pvals_diff]
    sub_template_cci_dt[, BH_PVAL_DIFF := BH_pvals_diff]
    sub_template_cci_dt <- sub_template_cci_dt[, c("LR_GENES", "L_CELLTYPE", "R_CELLTYPE", "L_GENE", "R_GENE",
                                                   paste0("PVAL_", cond1), paste0("PVAL_", cond2), "PVAL_DIFF",
                                                   paste0("BH_PVAL_", cond1), paste0("BH_PVAL_", cond2), "BH_PVAL_DIFF"), with = FALSE]
  }
  cci_dt <- data.table::merge.data.table(
    x = cci_dt_simple,
    y = sub_template_cci_dt,
    by = colnames(cci_dt_simple)[1:5],
    all.x = TRUE,
    sort = FALSE
  )
  cci_dt[is.na(cci_dt)] <- 1
  return(cci_dt)
}

#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param template_cci_dt x
#' @param cond1 x
#' @param cond2 x
#'
#' @return x
run_stat_iteration <- function(
  expr_tr,
  metadata,
  template_cci_dt,
  cond1,
  cond2
) {
  if(is.null(cond1) | is.null(cond2)) {
    meta_ct <- metadata
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    return(run_simple_analysis(
      expr_tr = expr_tr,
      metadata = meta_ct ,
      template_cci_dt = template_cci_dt ,
      cond1 = cond1,
      cond2 = cond2,
      threshold = 0.1,
      compute_fast = TRUE
    ))
  } else {
    meta_cond <- metadata
    for(x in unique(meta_cond$cell_type)) {
      meta_cond$condition[meta_cond$cell_type == x] <- sample(meta_cond$condition[meta_cond$cell_type == x])
    }
    permcond_dt <- run_simple_analysis(
      expr_tr = expr_tr,
      metadata = meta_cond ,
      template_cci_dt = template_cci_dt ,
      cond1 = cond1,
      cond2 = cond2,
      threshold = 0.1,
      compute_fast = TRUE
    )
    meta_ct <- metadata
    meta_ct$cell_type[meta_ct$condition == cond1] <- sample(meta_ct$cell_type[meta_ct$condition == cond1])
    meta_ct$cell_type[meta_ct$condition == cond2] <- sample(meta_ct$cell_type[meta_ct$condition == cond2])
    permct_dt <- run_simple_analysis(
      expr_tr = expr_tr,
      metadata = meta_ct ,
      template_cci_dt = template_cci_dt ,
      cond1 = cond1,
      cond2 = cond2,
      threshold = 0.1,
      compute_fast = TRUE
    )
    return(cbind(permcond_dt$cond2 - permcond_dt$cond1, permct_dt$cond1, permct_dt$cond2))
  }
}

#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param is_cond x
#'
#' @return x
aggregate_cells <- function(
  expr_tr,
  metadata,
  is_cond
) {
  if(!is_cond) {
    sums <- rowsum(x = expr_tr,
                   group = metadata[["cell_type"]])
    aggr <- sums/as.vector(table(metadata[["cell_type"]]))
  } else {
    group <-paste(metadata[["condition"]], metadata[["cell_type"]], sep = "_")
    sums <- rowsum(x = expr_tr,
                   group = group )
    aggr <- sums/as.vector(table(group))
  }
  return(aggr)
}

#' Title
#'
#' @param cci_dt x
#' @param cond1 x
#' @param cond2 x
#' @param permutation_analysis x
#'
#' @return x
clean_colnames <- function(
  cci_dt,
  cond1,
  cond2,
  permutation_analysis
) {
  first_cols <- c("LR_GENES", "L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE")
  if(is.null(cond1) | is.null(cond2)) {
    last_cols <- c("L_NCELLS", "L_EXPRESSION", "L_DETECTED", "R_NCELLS", "R_EXPRESSION", "R_DETECTED")
    if(!permutation_analysis) {
      ordered_cols <- c(first_cols,
                        "LR_SCORE", "LR_DETECTED",
                        last_cols)
    } else {
      ordered_cols <- c(first_cols,
                        "LR_SCORE", "LR_DETECTED", "PVAL", "BH_PVAL",
                        last_cols)
    }
  } else {
    last_cols <- c(paste0("L_NCELLS_", cond1), paste0("L_EXPRESSION_", cond1), paste0("L_DETECTED_", cond1),
                   paste0("R_NCELLS_", cond1), paste0("R_EXPRESSION_", cond1), paste0("R_DETECTED_", cond1),
                   paste0("L_NCELLS_", cond2), paste0("L_EXPRESSION_", cond2), paste0("L_DETECTED_", cond2),
                   paste0("R_NCELLS_", cond2), paste0("R_EXPRESSION_", cond2), paste0("R_DETECTED_", cond2))
    if(!permutation_analysis) {
      ordered_cols <- c(first_cols,
                        paste0("LR_SCORE_", cond1), paste0("LR_SCORE_", cond2),
                        paste0("LR_DETECTED_", cond1), paste0("LR_DETECTED_", cond2),
                        last_cols)
    } else {
      ordered_cols <- c(first_cols,
                        paste0("LR_SCORE_", cond1), paste0("LR_SCORE_", cond2),
                        paste0("LR_DETECTED_", cond1), paste0("LR_DETECTED_", cond2),
                        "PVAL_DIFF", "BH_PVAL_DIFF",
                        paste0("PVAL_", cond1), paste0("BH_PVAL_", cond1),
                        paste0("PVAL_", cond2), paste0("BH_PVAL_", cond2),
                        last_cols)
    }
  }
  data.table::setcolorder(
    x = cci_dt,
    neworder = ordered_cols
  )
  return(cci_dt)
}

#' Title
#'
#' @param cci_dt x
#' @param cond1 x
#' @param cond2 x
#' @param metadata x
#'
#' @return x
add_cell_number <- function(
  cci_dt,
  cond1,
  cond2,
  metadata
) {
  if(is.null(cond1) | is.null(cond2)) {
    dt_NCELLS <- metadata[, .N, by = "cell_type"]
    cci_dt <- data.table::merge.data.table(
      x = cci_dt,
      y = dt_NCELLS,
      by.x = "L_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE
    )
    cci_dt <- data.table::merge.data.table(
      x = cci_dt,
      y = dt_NCELLS,
      by.x = "R_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("_L", "_R")
    )
    data.table::setnames(
      x = cci_dt,
      old = c("N_L", "N_R"),
      new = c("L_NCELLS", "R_NCELLS")
    )
  } else {
    dt_NCELLS <- metadata[, .N, by = c("cell_type", "condition")]
    dt_NCELLS <- data.table::dcast.data.table(
      data = dt_NCELLS,
      formula = cell_type ~ condition,
      value.var = "N"
    )
    cci_dt <- data.table::merge.data.table(
      x = cci_dt,
      y = dt_NCELLS,
      by.x = "L_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE
    )
    cci_dt <- data.table::merge.data.table(
      x = cci_dt,
      y = dt_NCELLS,
      by.x = "R_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("_L", "_R")
    )
    data.table::setnames(
      x = cci_dt,
      old = c(paste0(cond1, "_L"), paste0(cond2, "_L"), paste0(cond1, "_R"), paste0(cond2, "_R")),
      new = c(paste0("L_NCELLS_", cond1), paste0("L_NCELLS_", cond2), paste0("R_NCELLS_", cond1), paste0("R_NCELLS_", cond2))
    )
  }
  return(cci_dt)
}
