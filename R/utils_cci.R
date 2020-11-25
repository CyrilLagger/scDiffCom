create_cci_template <- function(
  analysis_inputs
) {
  template <- CJ(
    EMITTER_CELLTYPE = analysis_inputs$cell_types,
    RECEIVER_CELLTYPE = analysis_inputs$cell_types,
    LR_SORTED = analysis_inputs$LRdb$LR_SORTED
  )
  template <- merge.data.table(
    x = template,
    y = analysis_inputs$LRdb,
    by.x = "LR_SORTED",
    by.y = "LR_SORTED",
    all.x = TRUE,
    sort = FALSE
  )
  template <- add_cell_number(
    template_table = template,
    condition_inputs = analysis_inputs$condition,
    metadata = analysis_inputs$metadata
  )
}

add_cell_number <- function(
  template_table,
  condition_inputs,
  metadata
) {
  if (!condition_inputs$is_cond) {
    dt_NCELLS <- metadata[, .N, by = "cell_type"]
    template_table <- merge.data.table(
      x = template_table,
      y = dt_NCELLS,
      by.x = "EMITTER_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE
    )
    template_table <- merge.data.table(
      x = template_table,
      y = dt_NCELLS,
      by.x = "RECEIVER_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("_L", "_R")
    )
    new_cols <- c("EMITTER_NCELLS", "RECEIVER_NCELLS")
    setnames(
      x = template_table,
      old = c("N_L", "N_R"),
      new = new_cols
    )
  } else {
    dt_NCELLS <- metadata[, .N, by = c("cell_type", "condition")]
    dt_NCELLS <- dcast.data.table(
      data = dt_NCELLS,
      formula = cell_type ~ condition,
      value.var = "N"
    )
    template_table <- merge.data.table(
      x = template_table,
      y = dt_NCELLS,
      by.x = "EMITTER_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE
    )
    template_table <- merge.data.table(
      x = template_table,
      y = dt_NCELLS,
      by.x = "RECEIVER_CELLTYPE",
      by.y = "cell_type",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("_L", "_R")
    )
    new_cols <- c(
      paste0("EMITTER_NCELLS_", condition_inputs$cond1),
      paste0("EMITTER_NCELLS_", condition_inputs$cond2),
      paste0("RECEIVER_NCELLS_", condition_inputs$cond1),
      paste0("RECEIVER_NCELLS_", condition_inputs$cond2)
    )
    setnames(
      x = template_table,
      old = c(
        paste0(condition_inputs$cond1, "_L"),
        paste0(condition_inputs$cond2, "_L"),
        paste0(condition_inputs$cond1, "_R"),
        paste0(condition_inputs$cond2, "_R")
      ),
      new = new_cols
    )
  }
  for (j in new_cols) {
    set(template_table, i = which(is.na(template_table[[j]])), j = j, value = 0)
  }
  return(template_table)
}

run_simple_cci_analysis <- function(
  analysis_inputs,
  cci_template,
  log_scale,
  threshold_min_cells,
  threshold_pct,
  compute_fast
) {
  LOGFC <- LOGFC_ABS <- NULL
  averaged_expr <- aggregate_cells(
    data_tr = analysis_inputs$data_tr,
    metadata = analysis_inputs$metadata,
    is_cond = analysis_inputs$condition$is_cond
  )
  cci_dt <- build_cci_or_drate(
    averaged_expr = averaged_expr,
    cci_template = cci_template,
    max_nL = analysis_inputs$max_nL,
    max_nR = analysis_inputs$max_nR,
    condition_inputs = analysis_inputs$condition,
    threshold_min_cells = threshold_min_cells,
    threshold_pct = threshold_pct,
    cci_or_drate = "cci"
  )
  if (compute_fast) {
    if (!analysis_inputs$condition$is_cond) {
      return(cci_dt[["CCI_SCORE"]])
    } else {
      return(list(
        cond1 = cci_dt[[paste0("CCI_SCORE_", analysis_inputs$condition$cond1)]],
        cond2 = cci_dt[[paste0("CCI_SCORE_", analysis_inputs$condition$cond2)]]
      ))
    }
  }
  detection_rate <- aggregate_cells(
    data_tr = 1 * (analysis_inputs$data_tr > 0),
    metadata = analysis_inputs$metadata,
    is_cond = analysis_inputs$condition$is_cond
  )
  drate_dt <- build_cci_or_drate(
    averaged_expr = detection_rate,
    cci_template = cci_template,
    max_nL = analysis_inputs$max_nL,
    max_nR = analysis_inputs$max_nR,
    condition_inputs = analysis_inputs$condition,
    threshold_pct = threshold_pct,
    threshold_min_cells = threshold_min_cells,
    cci_or_drate = "drate"
  )
  dt <- merge.data.table(
    x = cci_dt,
    y = drate_dt,
    by = intersect(names(cci_dt), names(drate_dt)),
    sort = FALSE
  )
  if (analysis_inputs$condition$is_cond) {
    if (log_scale) {
      dt[,
         LOGFC := get(paste0("CCI_SCORE_", analysis_inputs$condition$cond2)) -
           get(paste0("CCI_SCORE_", analysis_inputs$condition$cond1))
         ]
    } else {
      dt[,
         LOGFC := log(get(paste0("CCI_SCORE_", analysis_inputs$condition$cond2)) /
                        get(paste0("CCI_SCORE_", analysis_inputs$condition$cond1)))
         ]
      max_logfc <- max(dt[is.finite(LOGFC)][["LOGFC"]])
      min_logfc <- min(dt[is.finite(LOGFC)][["LOGFC"]])
      dt[, LOGFC := ifelse(
        is.finite(LOGFC),
        LOGFC,
        ifelse(
          is.nan(LOGFC),
          0,
          ifelse(
            LOGFC > 0,
            max_logfc,
            min_logfc
          )
        )
      )]
    }
    dt[, LOGFC_ABS := abs(LOGFC)]
  }
  return(dt)
}

aggregate_cells <- function(
  data_tr,
  metadata,
  is_cond
) {
  if (!is_cond) {
    group <- metadata[["cell_type"]]
  } else {
    group <- paste(metadata[["condition"]], metadata[["cell_type"]], sep = "_")
  }
  sums <- DelayedArray::rowsum(
    x = data_tr,
    group = group,
    reorder = TRUE
  )
  aggr <- sums / as.vector(table(group))
  return(aggr)
}

build_cci_or_drate <- function(
  averaged_expr,
  cci_template,
  max_nL,
  max_nR,
  condition_inputs,
  threshold_min_cells,
  threshold_pct,
  cci_or_drate
) {
  CONDITION_CELLTYPE <- NULL
  full_dt <- copy(cci_template)
  if (cci_or_drate == "cci") {
    name_tag <- "EXPRESSION"
  } else if (cci_or_drate == "drate") {
    name_tag <- "DETECTED"
  } else {
    stop("Error in build_cci_drate_dt.")
  }
  if (!condition_inputs$is_cond) {
    row_id <- "CELLTYPE"
    vars_id <- "CELLTYPE"
    cond1_id <- NULL
    cond2_id <- NULL
    merge_id <- name_tag
    score_id <- "CCI_SCORE"
    dr_id <- "CCI_DETECTED"
    n_id <- 1
    pmin_id <- NULL
  } else {
    row_id <- "CONDITION_CELLTYPE"
    vars_id <- c("CELLTYPE", "CONDITION")
    cond1_id <- paste0("_", condition_inputs$cond1)
    cond2_id <- paste0("_", condition_inputs$cond2)
    merge_id <- c(condition_inputs$cond1, condition_inputs$cond2)
    score_id <- paste0("CCI_SCORE_", c(condition_inputs$cond1, condition_inputs$cond2))
    dr_id <- paste0("CCI_DETECTED_", c(condition_inputs$cond1, condition_inputs$cond2))
    n_id <- 2
    pmin_id <- c(cond1_id, cond2_id)
  }
  dt <- as.data.table(
    x = averaged_expr,
    keep.rownames = row_id,
    sorted = FALSE
  )
  if (condition_inputs$is_cond) {
    dt[, c("CONDITION", "CELLTYPE") := split_cond_string(CONDITION_CELLTYPE, condition_inputs$cond1, condition_inputs$cond2)]
    dt[, CONDITION_CELLTYPE := NULL]
  }
  dt <- melt.data.table(
    data = dt,
    id.vars = vars_id,
    variable.name = "GENE",
    value.name = name_tag
  )
  if (condition_inputs$is_cond) {
    dt <- dcast.data.table(
      data = dt,
      formula = CELLTYPE + GENE ~ CONDITION,
      value.var = name_tag
    )
  }
  dt[is.na(dt)] <- 0
  out_names <- c(
    sapply(
      1:max_nL,
      function(i) {
        paste0("L", i, "_", name_tag, pmin_id)
      }
    ),
    sapply(
      1:max_nR,
      function(i) {
        paste0("R", i, "_", name_tag, pmin_id)
      }
    )
  )
  full_dt[
    ,
    c(out_names) :=
      c(
        sapply(
          1:max_nL,
          function(i) {
            as.list(
              dt[.SD,
                 on = c(paste0("GENE==LIGAND_", i), "CELLTYPE==EMITTER_CELLTYPE"),
                 mget(paste0("x.", merge_id))
                 ]
            )
          }
        ),
        sapply(
          1:max_nR,
          function(i) {
            as.list(
              dt[.SD,
                 on = c(paste0("GENE==RECEPTOR_", i), "CELLTYPE==RECEIVER_CELLTYPE"),
                 mget(paste0("x.", merge_id))
                 ]
            )
          }
        )
      )
    ]
  if (cci_or_drate == "cci") {
    full_dt[, (score_id) :=
              lapply(1:n_id, function(x) {
                (
                  do.call(
                    pmin,
                    c(lapply(1:max_nL, function(i) {
                      get(paste0("L", i, "_", name_tag, pmin_id[x]))
                    }), na.rm = TRUE)
                  )
                  +
                    do.call(
                      pmin,
                      c(lapply(1:max_nR, function(i) {
                        get(paste0("R", i, "_", name_tag, pmin_id[x]))
                      }), na.rm = TRUE)
                    )
                ) / 2
              })]
  } else if (cci_or_drate == "drate") {
    full_dt[, (dr_id) :=
              lapply(1:n_id, function(x) {
                is_detected_full(
                  x_dr = do.call(pmin, c(lapply(1:max_nL, function(i) {
                    get(paste0("L", i, "_", name_tag, pmin_id[x]))
                  }), na.rm = TRUE)),
                  x_ncells = get(paste0("EMITTER_NCELLS", pmin_id[x])),
                  y_dr = do.call(pmin, c(lapply(1:max_nR, function(i) {
                    get(paste0("R", i, "_", name_tag, pmin_id[x]))
                  }), na.rm = TRUE)),
                  y_ncells = get(paste0("RECEIVER_NCELLS", pmin_id[x])),
                  dr_thr = threshold_pct,
                  threshold_min_cells = threshold_min_cells
                )
              })]
  }
  return(full_dt)
}

clean_colnames <- function(
  cci_dt,
  condition_inputs,
  max_nL,
  max_nR,
  permutation_analysis
) {
  first_cols <- c("LR_GENES", "LR_SORTED", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "ER_CELLTYPES")
  LR_COLNAMES <- c(
    paste0("LIGAND_", 1:max_nL),
    paste0("RECEPTOR_", 1:max_nR)
  )
  first_cols <- c(first_cols, LR_COLNAMES)
  if (!condition_inputs$is_cond) {
    last_cols <- c(
      "EMITTER_NCELLS",
      paste0("L", 1:max_nL, "_EXPRESSION"),
      paste0("L", 1:max_nL, "_DETECTED"),
      "RECEIVER_NCELLS",
      paste0("R", 1:max_nR, "_EXPRESSION"),
      paste0("R", 1:max_nR, "_DETECTED")
    )
    if (!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        "CCI_SCORE", "CCI_DETECTED",
        last_cols
      )
    } else {
      ordered_cols <- c(
        first_cols,
        "CCI_SCORE", "CCI_DETECTED_AND_SIGNIFICANT", "CCI_DETECTED", "P_VALUE", "BH_P_VALUE",
        last_cols
      )
    }
  } else {
    last_cols <- c(
      paste0("EMITTER_NCELLS_", condition_inputs$cond1),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_inputs$cond1),
      paste0("L", 1:max_nL, "_DETECTED_", condition_inputs$cond1),
      paste0("RECEIVER_NCELLS_", condition_inputs$cond1),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_inputs$cond1),
      paste0("R", 1:max_nR, "_DETECTED_", condition_inputs$cond1),
      paste0("EMITTER_NCELLS_", condition_inputs$cond2),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_inputs$cond2),
      paste0("L", 1:max_nL, "_DETECTED_", condition_inputs$cond2),
      paste0("RECEIVER_NCELLS_", condition_inputs$cond2),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_inputs$cond2),
      paste0("R", 1:max_nR, "_DETECTED_", condition_inputs$cond2)
    )
    if (!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        "LOGFC", "LOGFC_ABS",
        paste0("CCI_SCORE_", condition_inputs$cond1), paste0("CCI_SCORE_", condition_inputs$cond2),
        paste0("CCI_DETECTED_", condition_inputs$cond1), paste0("CCI_DETECTED_", condition_inputs$cond2),
        last_cols
      )
    } else {
      ordered_cols <- c(
        first_cols,
        "LOGFC", "LOGFC_ABS",
        "DIFFERENTIALLY_EXPRESSED", "DIFFERENTIAL_DIRECTION", "REGULATION", "REGULATION_SIMPLE",
        paste0("CCI_SCORE_", condition_inputs$cond1), paste0("CCI_SCORE_", condition_inputs$cond2),
        paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1), paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2),
        paste0("CCI_DETECTED_", condition_inputs$cond1), paste0("CCI_DETECTED_", condition_inputs$cond2),
        "P_VALUE_DE",
        "BH_P_VALUE_DE",
        paste0("P_VALUE_", condition_inputs$cond1),
        paste0("BH_P_VALUE_", condition_inputs$cond1),
        paste0("P_VALUE_", condition_inputs$cond2),
        paste0("BH_P_VALUE_", condition_inputs$cond2),
        last_cols
      )
    }
  }
  setcolorder(
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
    threshold_min_cells
) {
    if (x_dr >= dr_thr &
        x_dr * x_ncells >= threshold_min_cells &
        y_dr >= dr_thr &
        y_dr * y_ncells >= threshold_min_cells
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
      stop = nchar(cond1)
    ) == cond1,
    cond1,
    cond2
  )
  cell_type <- ifelse(
    substr(
      x = expr_string,
      start = 1,
      stop = nchar(cond1)
    ) == cond1,
    substr(
      x = expr_string,
      start = nchar(cond1) + 2,
      stop = nchar(expr_string)
    ),
    substr(
      x = expr_string,
      start = nchar(cond2) + 2,
      stop = nchar(expr_string)
    )
  )
  return(list(cond = cond, cell_type = cell_type))
}
