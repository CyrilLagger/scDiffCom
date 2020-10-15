create_template_cci <- function(
  LR_db,
  cell_types
) {
  template <- data.table::CJ(
    L_CELLTYPE = cell_types,
    R_CELLTYPE = cell_types,
    LR_SORTED = LR_db$LR_SORTED
  )
  template <- data.table::merge.data.table(
    x = template,
    y = LR_db,
    by.x = "LR_SORTED",
    by.y = "LR_SORTED",
    all.x = TRUE,
    sort = FALSE
  )
}

add_cell_number <- function(
  cci_dt,
  condition_info,
  metadata
) {
  if(!condition_info$is_cond) {
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
    new_cols <- c("L_NCELLS", "R_NCELLS")
    data.table::setnames(
      x = cci_dt,
      old = c("N_L", "N_R"),
      new = new_cols
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
    new_cols <- c(
      paste0("L_NCELLS_", condition_info$cond1),
      paste0("L_NCELLS_", condition_info$cond2),
      paste0("R_NCELLS_", condition_info$cond1),
      paste0("R_NCELLS_", condition_info$cond2))
    data.table::setnames(
      x = cci_dt,
      old = c(
        paste0(condition_info$cond1, "_L"),
        paste0(condition_info$cond2, "_L"),
        paste0(condition_info$cond1, "_R"),
        paste0(condition_info$cond2, "_R")),
      new = new_cols
    )
  }
  for(j in new_cols){
    data.table::set(cci_dt, i = which(is.na(cci_dt[[j]])), j = j, value = 0)
  }
  return(cci_dt)
}

run_simple_cci_analysis <- function(
  expr_tr,
  metadata,
  template_cci_dt,
  pp_LR,
  condition_info,
  pct_threshold,
  min_cells = min_cells,
  compute_fast
) {
  averaged_expr <- aggregate_cells(
    expr_tr = expr_tr,
    metadata = metadata,
    is_cond = condition_info$is_cond
  )
  cci_dt <- build_cci_or_drate(
    averaged_expr = averaged_expr,
    template_cci_dt = template_cci_dt,
    pp_LR = pp_LR,
    condition_info = condition_info,
    pct_threshold = pct_threshold,
    min_cells = min_cells,
    cci_or_drate = "cci"
  )
  if(compute_fast) {
    if(!condition_info$is_cond) {
      return(cci_dt[["LR_SCORE"]])
    } else {
      return(list(
        cond1 = cci_dt[[paste0("LR_SCORE_", condition_info$cond1)]],
        cond2 = cci_dt[[paste0("LR_SCORE_", condition_info$cond2)]]
      )
      )
    }
  }
  detection_rate <- aggregate_cells(
    expr_tr = 1 * (expr_tr > 0),
    metadata = metadata,
    is_cond = condition_info$is_cond
  )
  drate_dt <- build_cci_or_drate(
    averaged_expr = detection_rate,
    template_cci_dt = template_cci_dt,
    pp_LR = pp_LR,
    condition_info = condition_info,
    pct_threshold = pct_threshold,
    min_cells = min_cells,
    cci_or_drate = "drate"
  )
  dt <- data.table::merge.data.table(
    x = cci_dt,
    y = drate_dt,
    by = generics::intersect(names(cci_dt), names(drate_dt)),
    sort = FALSE
  )
  return(dt)
}

aggregate_cells <- function(
  expr_tr,
  metadata,
  is_cond
) {
  if(!is_cond) {
    group <- metadata[["cell_type"]]
  } else {
    group <- paste(metadata[["condition"]], metadata[["cell_type"]], sep = "_")
  }
  sums <- DelayedArray::rowsum(
    x = expr_tr,
    group = group,
    reorder = TRUE
  )
  aggr <- sums/as.vector(table(group))
  return(aggr)
}

build_cci_or_drate <- function(
  averaged_expr,
  template_cci_dt,
  pp_LR,
  condition_info,
  pct_threshold,
  min_cells,
  cci_or_drate
) {
  CONDITION_CELLTYPE <- NULL
  full_dt <- data.table::copy(template_cci_dt)
  if(cci_or_drate == "cci") {
    name_tag = "EXPRESSION"
  } else if(cci_or_drate == "drate") {
    name_tag = "DETECTED"
  } else {
    stop("Error in build_cci_drate_dt.")
  }
  if(!condition_info$is_cond) {
    row_id <- "CELLTYPE"
    vars_id <- "CELLTYPE"
    cond1_id <- NULL
    cond2_id <- NULL
    merge_id <- name_tag
    score_id <- "LR_SCORE"
    dr_id <- "LR_DETECTED"
    n_id <- 1
    pmin_id <- NULL
  } else {
    row_id <- "CONDITION_CELLTYPE"
    vars_id <- c("CELLTYPE", "CONDITION")
    cond1_id <- paste0("_", condition_info$cond1)
    cond2_id <- paste0("_", condition_info$cond2)
    merge_id <- c(condition_info$cond1, condition_info$cond2)
    score_id <- paste0("LR_SCORE_", c(condition_info$cond1, condition_info$cond2))
    dr_id <- paste0("LR_DETECTED_", c(condition_info$cond1, condition_info$cond2))
    n_id <- 2
    pmin_id <- c(cond1_id, cond2_id)
  }
  dt <- data.table::as.data.table(
    x = averaged_expr,
    keep.rownames = row_id,
    sorted = FALSE
  )
  if(condition_info$is_cond) {
    dt[, c("CONDITION", "CELLTYPE") := split_cond_string(CONDITION_CELLTYPE, condition_info$cond1, condition_info$cond2)]
    dt[, CONDITION_CELLTYPE := NULL]
  }
  dt <- data.table::melt.data.table(
    data = dt,
    id.vars = vars_id,
    variable.name = "GENE",
    value.name = name_tag
  )
  if(condition_info$is_cond) {
    dt <- data.table::dcast.data.table(
      data = dt,
      formula = CELLTYPE + GENE ~ CONDITION,
      value.var = name_tag)
  }
  dt[is.na(dt)] <- 0
  out_names <- c(
    sapply(
      1:pp_LR$max_nL,
      function(i) {
        paste0("L", i, "_", name_tag, pmin_id)
      }),
    sapply(
      1:pp_LR$max_nR,
      function(i) {
        paste0("R", i, "_", name_tag, pmin_id)
      })
  )
  full_dt[,
          c(out_names) :=
            c(
              sapply(
                1:pp_LR$max_nL,
                function(i) {
                  as.list(
                    dt[.SD,
                       on = c(paste0("GENE==LIGAND_", i), "CELLTYPE==L_CELLTYPE"),
                       mget(paste0("x.", merge_id))
                       ])
                }),
              sapply(
                1:pp_LR$max_nR,
                function(i) {
                  as.list(
                    dt[.SD,
                       on = c(paste0("GENE==RECEPTOR_", i), "CELLTYPE==R_CELLTYPE"),
                       mget(paste0("x.", merge_id))
                       ])
                })
            )
          ]
  if(cci_or_drate == "cci") {
    full_dt[, (score_id) :=
              lapply(1:n_id, function(x) {
                (
                  do.call(pmin,
                          c(lapply(1:pp_LR$max_nL, function(i) {get(paste0("L", i, "_", name_tag, pmin_id[x]))}), na.rm = TRUE))
                  +
                    do.call(pmin,
                            c(lapply(1:pp_LR$max_nR, function(i) {get(paste0("R", i,"_", name_tag, pmin_id[x]))}), na.rm = TRUE))
                )/2
              })
            ]
  } else if(cci_or_drate == "drate") {
    full_dt[, (dr_id) :=
              lapply(1:n_id, function(x) {
                is_detected_full(
                  x_dr = do.call(pmin, c(lapply(1:pp_LR$max_nL, function(i) {get(paste0("L", i, "_", name_tag, pmin_id[x]))}), na.rm = TRUE)),
                  x_ncells = get(paste0("L_NCELLS", pmin_id[x])),
                  y_dr = do.call(pmin, c(lapply(1:pp_LR$max_nR, function(i) {get(paste0("R", i, "_", name_tag, pmin_id[x]))}), na.rm = TRUE)),
                  y_ncells = get(paste0("R_NCELLS", pmin_id[x])),
                  dr_thr = pct_threshold,
                  min_cells = min_cells
                )
              })
            ]
  }
  return(full_dt)
}

clean_colnames <- function(
  cci_dt,
  condition_info,
  max_nL,
  max_nR,
  permutation_analysis
) {
  first_cols <- c("LR_NAME", "LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE")
  LR_names <- c(
    paste0("LIGAND_", 1:max_nL),
    paste0("RECEPTOR_", 1:max_nR )
  )
  first_cols <- c(first_cols, LR_names)
  if(!condition_info$is_cond) {
    last_cols <- c(
      "L_NCELLS",
      paste0("L", 1:max_nL, "_EXPRESSION"),
      paste0("L", 1:max_nL, "_DETECTED"),
      "R_NCELLS",
      paste0("R", 1:max_nR, "_EXPRESSION"),
      paste0("R", 1:max_nR, "_DETECTED")
    )
    if(!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        "LR_SCORE", "LR_DETECTED",
        last_cols
      )
    } else {
      ordered_cols <- c(
        first_cols,
        "LR_SCORE", "LR_DETECTED_AND_SIGNIFICANT", "LR_DETECTED", "PVAL", "BH_PVAL",
        last_cols
      )
    }
  } else {
    last_cols <- c(
      paste0("L_NCELLS_", condition_info$cond1),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_info$cond1),
      paste0("L", 1:max_nL, "_DETECTED_", condition_info$cond1),
      paste0("R_NCELLS_", condition_info$cond1),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_info$cond1),
      paste0("R", 1:max_nR, "_DETECTED_", condition_info$cond1),
      paste0("L_NCELLS_", condition_info$cond2),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_info$cond2),
      paste0("L", 1:max_nL, "_DETECTED_", condition_info$cond2),
      paste0("R_NCELLS_", condition_info$cond2),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_info$cond2),
      paste0("R", 1:max_nR, "_DETECTED_", condition_info$cond2)
    )
    if(!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        "LOGFC", "LOGFC_ABS",
        paste0("LR_SCORE_", condition_info$cond1), paste0("LR_SCORE_", condition_info$cond2),
        paste0("LR_DETECTED_", condition_info$cond1), paste0("LR_DETECTED_", condition_info$cond2),
        last_cols
      )
    } else {
      ordered_cols <- c(
        first_cols,
        "LOGFC", "LOGFC_ABS",
        "DIFFERENTIALLY_EXPRESSED", "DIFFERENTIAL_DIRECTION", "REGULATION", "REGULATION_SIMPLE",
        paste0("LR_SCORE_", condition_info$cond1), paste0("LR_SCORE_", condition_info$cond2),
        paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1), paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2),
        paste0("LR_DETECTED_", condition_info$cond1), paste0("LR_DETECTED_", condition_info$cond2),
        "PVAL_DIFF",
        "BH_PVAL_DIFF",
        paste0("PVAL_", condition_info$cond1),
        paste0("BH_PVAL_", condition_info$cond1),
        paste0("PVAL_", condition_info$cond2),
        paste0("BH_PVAL_", condition_info$cond2),
        last_cols
      )
    }
  }
  data.table::setcolorder(
    x = cci_dt,
    neworder = ordered_cols
  )
  return(cci_dt)
}

is_detected_full <- Vectorize(
  function(
    x_dr,
    x_ncells,
    y_dr,
    y_ncells,
    dr_thr,
    min_cells
  ) {
    if (x_dr >= dr_thr &
        x_dr*x_ncells >= min_cells &
        y_dr >= dr_thr &
        y_dr*y_ncells >= min_cells
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
)

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

