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
run_diffcom <- function(
  seurat_obj,
  LR_data,
  seurat_cell_type_id = "cell_ontology_class",
  condition_id = NULL,
  assay = "RNA",
  slot = "data",
  log_scale = TRUE,
  min_cells = 5,
  threshold = 0.1,
  differential_analysis = TRUE,
  specificity_analysis = TRUE,
  iterations = 1000
) {
  message("Preprocessing Seurat object.")
  pp_seurat <- preprocess_seurat(
    seurat_obj = seurat_obj,
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
  cci_analysis <- run_cci_analysis(
    expr_tr = expr_tr,
    metadata = metadata,
    template_cci_dt = template_cci_dt,
    threshold = threshold,
    condition_id = condition_id,
    differential_analysis = differential_analysis,
    specificity_analysis = specificity_analysis,
    iterations = iterations
  )
  if(is.null(condition_id)) {
    if(!specificity_analysis) {
      data.table::setcolorder(
        x = cci_analysis,
        neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type",
                     "LR_score", "LR_detec", "Ligand_expr", "Ligand_detec_rate", "Receptor_expr", "Receptor_detec_rate")
      )
    } else {
      data.table::setcolorder(
        x = cci_analysis,
        neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type",
                     "LR_score", "LR_detec", "pvals", "Ligand_expr", "Ligand_detec_rate", "Receptor_expr", "Receptor_detec_rate")
      )
    }

  } else {
    conds <- unique(metadata$condition)
    cond1 <- conds[[1]]
    cond2 <- conds[[2]]
    if(!specificity_analysis) {
      if(!differential_analysis) {
        data.table::setcolorder(
          x = cci_analysis,
          neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type",
                       paste0("LR_score_", cond1), paste0("LR_score_", cond2),
                       paste0("LR_detec_", cond1), paste0("LR_detec_", cond2),
                       paste0("Ligand_expr_", cond1), paste0("Ligand_detec_rate_", cond1),
                       paste0("Receptor_expr_", cond1), paste0("Receptor_detec_rate_", cond1),
                       paste0("Ligand_expr_", cond2), paste0("Ligand_detec_rate_", cond2),
                       paste0("Receptor_expr_", cond2), paste0("Receptor_detec_rate_", cond2))
        )
      } else {
        data.table::setcolorder(
          x = cci_analysis,
          neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type", "pvals_diff",
                       paste0("LR_score_", cond1), paste0("LR_score_", cond2),
                       paste0("LR_detec_", cond1), paste0("LR_detec_", cond2),
                       paste0("Ligand_expr_", cond1), paste0("Ligand_detec_rate_", cond1),
                       paste0("Receptor_expr_", cond1), paste0("Receptor_detec_rate_", cond1),
                       paste0("Ligand_expr_", cond2), paste0("Ligand_detec_rate_", cond2),
                       paste0("Receptor_expr_", cond2), paste0("Receptor_detec_rate_", cond2))
        )
      }
    } else {
      if(!differential_analysis) {
        data.table::setcolorder(
          x = cci_analysis,
          neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type",
                       paste0("pvals_", cond1), paste0("pvals_", cond2),
                       paste0("LR_score_", cond1), paste0("LR_score_", cond2),
                       paste0("LR_detec_", cond1), paste0("LR_detec_", cond2),
                       paste0("Ligand_expr_", cond1), paste0("Ligand_detec_rate_", cond1),
                       paste0("Receptor_expr_", cond1), paste0("Receptor_detec_rate_", cond1),
                       paste0("Ligand_expr_", cond2), paste0("Ligand_detec_rate_", cond2),
                       paste0("Receptor_expr_", cond2), paste0("Receptor_detec_rate_", cond2))
        )
      } else {
        data.table::setcolorder(
          x = cci_analysis,
          neworder = c("LR_pair", "ligand", "receptor", "Ligand_cell_type", "Receptor_cell_type", "pvals_diff",
                       paste0("pvals_", cond1), paste0("pvals_", cond2),
                       paste0("LR_score_", cond1), paste0("LR_score_", cond2),
                       paste0("LR_detec_", cond1), paste0("LR_detec_", cond2),
                       paste0("Ligand_expr_", cond1), paste0("Ligand_detec_rate_", cond1),
                       paste0("Receptor_expr_", cond1), paste0("Receptor_detec_rate_", cond1),
                       paste0("Ligand_expr_", cond2), paste0("Ligand_detec_rate_", cond2),
                       paste0("Receptor_expr_", cond2), paste0("Receptor_detec_rate_", cond2))
        )
      }
    }
  }
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
    Ligand_cell_type = cell_types,
    Receptor_cell_type = cell_types,
    LR_pair = LR_df$LR_pair
  )
  template <- data.table::merge.data.table(
    x = template,
    y = LR_df,
    by.x = "LR_pair",
    by.y = "LR_pair",
    all.x = TRUE,
    sort = FALSE
  )
}


#' Main internal routine
#'
#' @param expr_tr matrix with the (transposed) expression data from the Seurat object
#' @param metadata data.table with the relevant metadata from the Seurat object
#' @param template_cci_dt data.table with a template of all CCIs
#' @param threshold numeric indicating the percentage of cells that need to express a gene in a cluster for the gene to be
#' considered detected.
#' @param condition_id character indicating the column specifying the conditions on the cells. Set to NULL to run the analysis
#' on no condition; default is NULL.
#' @param differential_analysis logical indicating if performing the permutation test for the significance of differential
#' expression between the conditions. Only considered when condidition_id is not NULL; default is TRUE.
#' @param specificity_analysis logical indicating if performing the permutation test for the specificity of each CCI on each condition;
#' default is FALSE.
#' @param iterations integer indicating the number of iterations during the permutation test.
#'
#' @return x
run_cci_analysis <- function(
  expr_tr,
  metadata,
  template_cci_dt,
  threshold,
  condition_id,
  differential_analysis,
  specificity_analysis,
  iterations
) {
  if (is.null(condition_id)) {
    message("Performing simple analysis without condition.")
    cci_dt_simple <- run_simple_analysis(
      expr_tr = expr_tr,
      metadata = metadata,
      template_cci_dt = template_cci_dt,
      cond1 = NULL,
      cond2 = NULL,
      threshold = threshold,
      compute_fast = FALSE )
    if (!specificity_analysis) {
      cci_dt <- cci_dt_simple
    } else {
      message("Performing specificity analysis without condition.")
      cci_dt <- run_stat_analysis(
        cci_dt_simple = cci_dt_simple,
        expr_tr = expr_tr,
        metadata = metadata,
        cond1 = NULL,
        cond2 = NULL,
        iterations = iterations,
        use_case = "no_cond_spec",
        return_distr = FALSE
      )
    }
  } else{
    conds <- unique(metadata$condition)
    if (length(conds) != 2)
      stop("Wrong number of groups in cell-type conditions (expected 2).")
    cond1 <- conds[[1]]
    cond2 <- conds[[2]]
    message("Performing simple analysis on the two conditions.")
    cci_dt_simple <- run_simple_analysis(
      expr_tr = expr_tr,
      metadata = metadata,
      template_cci_dt = template_cci_dt,
      cond1 = cond1,
      cond2 = cond2,
      threshold = threshold,
      compute_fast = FALSE )
    if (!differential_analysis) {
      if (!specificity_analysis) {
        cci_dt <- cci_dt_simple
      } else {
        message("Performing specificity analysis on the two conditions.")
        cci_dt <- run_stat_analysis(
          cci_dt_simple = cci_dt_simple,
          expr_tr = expr_tr,
          metadata = metadata,
          cond1 = cond1,
          cond2 = cond2,
          iterations = iterations,
          use_case = "cond_spec",
          return_distr = FALSE
        )
      }
    } else {
      if (!specificity_analysis) {
        message("Performing differential analysis between the two conditions.")
        cci_dt <- run_stat_analysis(
          cci_dt_simple = cci_dt_simple,
          expr_tr = expr_tr,
          metadata = metadata,
          cond1 = cond1,
          cond2 = cond2,
          iterations = iterations,
          use_case = "cond_diff",
          return_distr = FALSE
        )
      } else {
        message("Performing differential analysis between the two conditions and specificity analysis.")
        cci_dt <- run_stat_analysis(
          cci_dt_simple = cci_dt_simple,
          expr_tr = expr_tr,
          metadata = metadata,
          cond1 = cond1,
          cond2 = cond2,
          iterations = iterations,
          use_case = "cond_diff_spec",
          return_distr = FALSE
        )
      }
    }
  }
  return(cci_dt)
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
    averaged_expr <- aggregate_cells(
      expr_tr = expr_tr,
      metadata = metadata,
      is_cond = FALSE
    )
    cci_dt <- build_cci_drate_dt(
      averaged_expr = averaged_expr,
      template_cci_dt = template_cci_dt,
      cond1 = NULL ,
      cond2 = NULL,
      detection_thr = 0,
      cci_or_drate = "cci"
    )
    if(compute_fast) {
      return(cci_dt[["LR_score"]])
    } else {
      detection_rate <- aggregate_cells(
        expr_tr = 1 * (expr_tr > 0),
        metadata = metadata,
        is_cond = FALSE
      )
      drate_dt <- build_cci_drate_dt(
        averaged_expr = detection_rate,
        template_cci_dt = template_cci_dt,
        cond1 = NULL,
        cond2 = NULL,
        detection_thr = threshold,
        cci_or_drate = "drate"
      )
      dt <- data.table::merge.data.table(
        x = cci_dt,
        y = drate_dt,
        by = colnames(cci_dt)[c(1,2,4,5,7)],
        sort = FALSE
      )
      return(dt)
    }
  } else {
    averaged_expr <- aggregate_cells(
      expr_tr = expr_tr,
      metadata = metadata,
      is_cond = TRUE
    )
    cci_dt <- build_cci_drate_dt(
      averaged_expr = averaged_expr,
      template_cci_dt = template_cci_dt,
      cond1 = cond1,
      cond2 = cond2,
      detection_thr = 0,
      cci_or_drate = "cci"
    )
    if(compute_fast) {
      return(list(
        cond1 = cci_dt[[paste0("LR_score_", cond1)]],
        cond2 = cci_dt[[paste0("LR_score_", cond2)]]
      )
      )
    } else {
      detection_rate <- aggregate_cells(
        expr_tr = 1 * (expr_tr > 0),
        metadata = metadata,
        is_cond = TRUE
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
        by = colnames(cci_dt)[c(1, 2, 5, 6, 9)],
        sort = FALSE
      )
      return(dt)
    }
  }
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
  Ligand_expr <- Receptor_expr <- Ligand_detec_rate <- Receptor_detec_rate <- cond_cell_type <-
    LR_score <- LR_detec <- gene <- cell_type <- ligand <- Ligand_cell_type <- receptor <-
    Receptor_cell_type <- NULL

  if(is.null(cond1) | is.null(cond2)) {
    dt <- data.table::as.data.table(
      x = averaged_expr,
      keep.rownames = "cell_type",
      sorted = FALSE
    )
    long_dt <- data.table::melt.data.table(
      data = dt,
      id.vars = "cell_type",
      variable.name = "gene",
      value.name = "avg_expr"
    )
    col_long_dt <- c("gene", "cell_type")
    setkeyv(
      x = long_dt,
      cols = col_long_dt
      )
    col_template_cci_dt <- c("ligand", "Ligand_cell_type", "receptor", "Receptor_cell_type")
    setkeyv(
      x = template_cci_dt,
      cols = col_template_cci_dt
      )
    full_dt <- long_dt[
      long_dt[
        template_cci_dt,
        on = c("gene==ligand", "cell_type==Ligand_cell_type")],
      on = c("gene==receptor", "cell_type==Receptor_cell_type") ]
    if(cci_or_drate == "cci") {
      data.table::setnames(full_dt,
                           c("cell_type", "gene",
                             "avg_expr", "i.avg_expr",
                             "i.cell_type", "i.gene"),
                           c("Receptor_cell_type", "receptor",
                             "Receptor_expr", "Ligand_expr",
                             "Ligand_cell_type", "ligand"))
      full_dt[, LR_score := (Ligand_expr + Receptor_expr) / 2]
    } else if(cci_or_drate == "drate") {
      data.table::setnames(full_dt,
                           c("cell_type", "gene",
                             "avg_expr", "i.avg_expr",
                             "i.cell_type", "i.gene"),
                           c("Receptor_cell_type", "receptor",
                             "Receptor_detec_rate", "Ligand_detec_rate",
                             "Ligand_cell_type", "ligand"))
      full_dt[, LR_detec := is_detected(Ligand_detec_rate, Receptor_detec_rate, detection_thr)]
    } else {
      stop("Error in build_cci_drate_dt.")
    }
  } else {
    dt <- data.table::as.data.table(
      x = averaged_expr,
      keep.rownames = "cond_cell_type",
      sorted = FALSE
    )
    dt[, c("condition", "cell_type") := split_cond_string(cond_cell_type, cond1, cond2)]
    dt[, cond_cell_type := NULL]
    long_dt <- data.table::melt.data.table(
      data = dt,
      id.vars = c("cell_type", "condition"),
      variable.name = "gene",
      value.name = "avg_expr"
    )
    long_dt <- data.table::dcast.data.table(
      data = long_dt,
      formula = cell_type + gene ~ condition,
      value.var = "avg_expr")
    setkey(long_dt, gene, cell_type)
    setkey(template_cci_dt, ligand, Ligand_cell_type, receptor, Receptor_cell_type)
    full_dt <- long_dt[
      long_dt[
        template_cci_dt,
        on = c("gene==ligand", "cell_type==Ligand_cell_type")],
      on = c("gene==receptor", "cell_type==Receptor_cell_type") ]
    if(cci_or_drate == "cci") {
      data.table::setnames(
        full_dt,
        c("cell_type", "gene",
          cond1, cond2,
          "i.cell_type", "i.gene",
          paste0("i.", cond1), paste0("i.", cond2)),
        c("Receptor_cell_type", "receptor",
          paste0("Receptor_expr_", cond1), paste0("Receptor_expr_", cond2),
          "Ligand_cell_type", "ligand",
          paste0("Ligand_expr_", cond1), paste0("Ligand_expr_", cond2)))
      full_dt[, paste0("LR_score_", c(cond1, cond2)) := .((get(paste0("Ligand_expr_", cond1)) + get(paste0("Receptor_expr_", cond1))) / 2,
                                                          (get(paste0("Ligand_expr_", cond2)) + get(paste0("Receptor_expr_", cond2))) / 2)]
    } else if(cci_or_drate == "drate") {
      data.table::setnames(
        full_dt,
        c("cell_type", "gene",
          cond1, cond2,
          "i.cell_type", "i.gene",
          paste0("i.", cond1), paste0("i.", cond2)),
        c("Receptor_cell_type", "receptor",
          paste0("Receptor_detec_rate_", cond1), paste0("Receptor_detec_rate_", cond2),
          "Ligand_cell_type", "ligand",
          paste0("Ligand_detec_rate_", cond1), paste0("Ligand_detec_rate_", cond2)))
      full_dt[, paste0("LR_detec_", c(cond1, cond2)) := .(is_detected(get(paste0("Ligand_detec_rate_", cond1)), get(paste0("Receptor_detec_rate_", cond1)), detection_thr),
                                                          is_detected(get(paste0("Ligand_detec_rate_", cond2)), get(paste0("Receptor_detec_rate_", cond2)), detection_thr))]
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
#' @param use_case x
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
  use_case,
  return_distr = FALSE
) {
  LR_detec <- LR_pair <- Ligand_cell_type <- Receptor_cell_type <- ligand <- receptor <- NULL
  if (use_case == "no_cond_spec") {
    sub_template_cci_dt <- cci_dt_simple[LR_detec == TRUE]
    sub_expr_tr <- expr_tr[,
                           colnames(expr_tr) %in%
                             unique(c(sub_template_cci_dt[["ligand"]], sub_template_cci_dt[["receptor"]]))]
    message(paste0("Number of genes initially: ", ncol(expr_tr), ". Number of genes for permutation test: ", ncol(sub_expr_tr)))
    cci_perm <- replicate(
      n = iterations,
      expr = run_stat_iteration(
        expr_tr = sub_expr_tr,
        metadata = metadata,
        template_cci_dt = sub_template_cci_dt[, .(LR_pair, Ligand_cell_type, Receptor_cell_type, ligand, receptor)],
        cond1 = NULL,
        cond2 = NULL,
        use_case = use_case
      ),
      simplify = "array"
    )
    distr <- cbind(cci_perm, sub_template_cci_dt[["LR_score"]])
    pvals <-
      rowSums(distr[, 1:iterations] >= distr[, (iterations + 1)]) / iterations
    sub_template_cci_dt[, pvals := pvals]
    sub_template_cci_dt <- sub_template_cci_dt[, list(LR_pair, Ligand_cell_type, Receptor_cell_type, ligand, receptor,
                                                      pvals)]
    if (return_distr) {
      return(distr)
    } else {
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
  } else {
    sub_template_cci_dt <- cci_dt_simple[get(paste0("LR_detec_", cond1)) == TRUE | get(paste0("LR_detec_", cond2)) == TRUE]
    sub_expr_tr <- expr_tr[,
                           colnames(expr_tr) %in%
                             unique(c(sub_template_cci_dt[["ligand"]], sub_template_cci_dt[["receptor"]]))]
    message(paste0("Number of genes initially: ", ncol(expr_tr), ". Number of genes for permutation test: ", ncol(sub_expr_tr)))
    cci_perm <- replicate(
      n = iterations,
      expr = run_stat_iteration(
        expr_tr = sub_expr_tr,
        metadata = metadata,
        template_cci_dt = sub_template_cci_dt[, .(LR_pair, Ligand_cell_type, Receptor_cell_type, ligand, receptor)],
        cond1 = cond1,
        cond2 = cond2,
        use_case = use_case
      ),
      simplify = "array"
    )
    if (use_case == "cond_spec") {
      distr_cond1 <- cbind(cci_perm[, 1, ], sub_template_cci_dt[[paste0("LR_score_", cond1)]])
      pvals_cond1 <-
        rowSums(distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]) / iterations
      distr_cond2 <- cbind(cci_perm[, 2, ], sub_template_cci_dt[[paste0("LR_score_", cond2)]])
      pvals_cond2 <-
        rowSums(distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]) / iterations
      sub_template_cci_dt[, paste0("pvals_", c(cond1, cond2)) := list(pvals_cond1, pvals_cond2)]
      sub_template_cci_dt <- sub_template_cci_dt[, c("LR_pair", "Ligand_cell_type", "Receptor_cell_type", "ligand", "receptor",
                                                     paste0("pvals_", cond1), paste0("pvals_", cond2)), with = FALSE]
      if (return_distr) {
        return(list(distr_cond1 = distr_cond1,
                    distr_cond2 = distr_cond2))
      } else {
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
    } else if (use_case == "cond_diff") {
      distr_diff <- cbind(cci_perm, sub_template_cci_dt[[paste0("LR_score_", cond2)]] - sub_template_cci_dt[[paste0("LR_score_", cond1)]])
      pvals_diff <-
        rowSums(abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations +
                                                                       1)])) / iterations
      sub_template_cci_dt[, pvals_diff := pvals_diff]
      sub_template_cci_dt <- sub_template_cci_dt[, list(LR_pair, Ligand_cell_type, Receptor_cell_type, ligand, receptor,
                                                        pvals_diff)]
      if (return_distr) {
        return(distr_diff)
      } else {
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
    } else if (use_case == "cond_diff_spec") {
      distr_diff <- cbind(cci_perm[, 1, ], sub_template_cci_dt[[paste0("LR_score_", cond2)]] - sub_template_cci_dt[[paste0("LR_score_", cond1)]])
      pvals_diff <-
        rowSums(abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations +
                                                                       1)])) / iterations
      distr_cond1 <- cbind(cci_perm[, 2, ], sub_template_cci_dt[[paste0("LR_score_", cond1)]])
      pvals_cond1 <-
        rowSums(distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]) / iterations
      distr_cond2 <- cbind(cci_perm[, 3, ], sub_template_cci_dt[[paste0("LR_score_", cond2)]])
      pvals_cond2 <-
        rowSums(distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]) / iterations
      sub_template_cci_dt[, paste0("pvals_", c(cond1, cond2)) := list(pvals_cond1, pvals_cond2)]
      sub_template_cci_dt[, pvals_diff := pvals_diff]
      sub_template_cci_dt <- sub_template_cci_dt[, c("LR_pair", "Ligand_cell_type", "Receptor_cell_type", "ligand", "receptor",
                                                     paste0("pvals_", cond1), paste0("pvals_", cond2), "pvals_diff"), with = FALSE]
      if (return_distr) {
        return(
          list(
            distr_diff = distr_diff,
            distr_cond1 = distr_cond1,
            distr_cond2 = distr_cond2
          )
        )
      } else {
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
    } else {
      stop("Case not supported in function run_stat_analysis.")
    }
  }
}


#' Title
#'
#' @param expr_tr x
#' @param metadata x
#' @param template_cci_dt x
#' @param cond1 x
#' @param cond2 x
#' @param use_case x
#'
#' @return x
run_stat_iteration <- function(
  expr_tr,
  metadata,
  template_cci_dt,
  cond1,
  cond2,
  use_case
) {
  if(use_case == "no_cond_spec") {
    meta_ct <- metadata
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    return(run_simple_analysis(
      expr_tr = expr_tr,
      metadata = meta_ct ,
      template_cci_dt = template_cci_dt ,
      cond1 = NULL,
      cond2 = NULL ,
      threshold = 0.1,
      compute_fast = TRUE
    ))
  } else if(use_case == "cond_spec") {
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
    return(cbind(permct_dt$cond1, permct_dt$cond2))
  } else if(use_case == "cond_diff") {
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
    return(permcond_dt$cond2 - permcond_dt$cond1)
  } else if(use_case == "cond_diff_spec") {
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
  } else {
    stop("Case  not supported in function run stat iteration")
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

