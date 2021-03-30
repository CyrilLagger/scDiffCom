run_filtering_and_ora <- function(
  object,
  new_threshold_quantile_score,
  new_threshold_p_value_specificity,
  new_threshold_p_value_de,
  new_threshold_logfc,
  skip_ora,
  verbose,
  class_signature
) {
  REGULATION_SIMPLE <- NULL
  temp_param <- object@parameters
  if (!is.null(new_threshold_quantile_score)) {
    temp_param$threshold_quantile_score <-
      new_threshold_quantile_score
  }
  if (!is.null(new_threshold_p_value_specificity)) {
    temp_param$threshold_p_value_specificity <-
      new_threshold_p_value_specificity
  }
  if (!is.null(new_threshold_p_value_de)) {
    temp_param$threshold_p_value_de <-
      new_threshold_p_value_de
  }
  if (!is.null(new_threshold_logfc)) {
    temp_param$threshold_logfc <- new_threshold_logfc
  }
  check_parameters <- validate_parameters(
    params = temp_param,
    from_inputs = FALSE
  )
  if (!is.null(check_parameters$check)) {
    stop(paste0(
      "Invalid parameters: ",
      paste0(check_parameters$check, collapse = " and ")
      ))
  } else {
    temp_param <- check_parameters$params
  }
  object@parameters <- temp_param
  methods::validObject(object)
  condition_inputs <- list(
    is_cond = temp_param$conditional_analysis,
    cond1 = temp_param$seurat_condition_id$cond1_name,
    cond2 = temp_param$seurat_condition_id$cond2_name
  )
  cci_dt <- copy(x = object@cci_table_raw)
  if (verbose) message("Filtering and cleaning `raw` CCIs.")
  cci_dt <- process_cci_raw(
    cci_dt = cci_dt,
    condition_inputs = condition_inputs,
    log_scale = temp_param$log_scale,
    permutation_analysis = temp_param$permutation_analysis,
    threshold_min_cells = temp_param$threshold_min_cells,
    threshold_pct = temp_param$threshold_pct,
    threshold_quantile_score = temp_param$threshold_quantile_score,
    threshold_p_value_specificity = temp_param$threshold_p_value_specificity,
    threshold_p_value_de = temp_param$threshold_p_value_de,
    threshold_logfc = temp_param$threshold_logfc,
    max_nL = temp_param$max_nL,
    max_nR = temp_param$max_nR,
    class_signature = class_signature
  )
  if (temp_param$permutation_analysis) {
    mes <- paste0(
      "Returning ",
      nrow(cci_dt),
      " detected CCIs."
    )
    if (verbose) message(mes)
  } else {
    mes <- paste0(
      "Returning ",
      nrow(cci_dt),
      " detected CCIs. (Note: incomplete filtering as statistical",
      " test has not been performed!)"
    )
    if (verbose) message(mes)
  }
  cci_dt <- clean_colnames(
    cci_dt = cci_dt,
    max_nL = temp_param$max_nL,
    max_nR = temp_param$max_nR,
    condition_inputs = condition_inputs,
    permutation_analysis = temp_param$permutation_analysis
  )
  if (nrow(cci_dt) == 0) {
    object@cci_table_detected <- list()
    object@ora_table <- list()
  } else {
    object@cci_table_detected <- cci_dt
    methods::validObject(object)
    if (skip_ora) {
      object@ora_table <- list()
    } else {
      object <- run_ora(
        object = object,
        categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
        overwrite = TRUE,
        stringent_or_default = "default",
        stringent_logfc_threshold = NULL,
        verbose = verbose,
        class_signature = class_signature,
        global = FALSE
      )
      if (class_signature == "scDiffComCombined") {
        if (verbose) message("Running global ORA.")
        object <- run_ora(
          object = object,
          categories = c(
            "ER_CELLTYPES",
            "LRI",
            "GO_TERMS",
            "KEGG_PWS",
            "ID"),
          overwrite = TRUE,
          stringent_or_default = "default",
          stringent_logfc_threshold = NULL,
          verbose = verbose,
          class_signature = class_signature,
          global = TRUE
        )
      }
    }
  }
  return(object)
}

process_cci_raw <- function(
  cci_dt,
  condition_inputs,
  log_scale,
  permutation_analysis,
  threshold_min_cells,
  threshold_pct,
  threshold_quantile_score,
  threshold_p_value_specificity,
  threshold_p_value_de,
  threshold_logfc,
  max_nL,
  max_nR,
  class_signature
) {
  BH_P_VALUE_DE <- P_VALUE_DE <- BH_P_VALUE <- P_VALUE <-
    IS_CCI_EXPRESSED <- ER_CELLTYPES <- EMITTER_CELLTYPE <-
    RECEIVER_CELLTYPE <- LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <-
    RECEPTOR_3 <- LOGFC <- LOGFC_ABS <- REGULATION <- IS_DE_LOGFC <-
    IS_DE_SIGNIFICANT <- DE_DIRECTION <- CCI_SCORE <- CCI <- LRI <- NULL
  if (class_signature == "scDiffCom") {
    temp_by <- NULL
  }
  if (class_signature == "scDiffComCombined") {
    temp_by <- "ID"
  }
  if (condition_inputs$is_cond) {
    cci_dt <- cci_dt[
      get(paste0("IS_CCI_EXPRESSED_", condition_inputs$cond1)) == TRUE |
        get(paste0("IS_CCI_EXPRESSED_", condition_inputs$cond2)) == TRUE
      ]
    if (permutation_analysis) {
      cci_dt[, c(
        paste0("BH_P_VALUE_", condition_inputs$cond1),
        paste0("BH_P_VALUE_", condition_inputs$cond2),
        "BH_P_VALUE_DE"
      ) := list(
        stats::p.adjust(
          get(paste0("P_VALUE_", condition_inputs$cond1)
              ),
          method = "BH"),
        stats::p.adjust(
          get(paste0("P_VALUE_", condition_inputs$cond2)
              ),
          method = "BH"),
        stats::p.adjust(P_VALUE_DE, method = "BH")
      ),
      by = temp_by
      ]
      cci_dt[, c(
        paste0("IS_CCI_SCORE_", condition_inputs$cond1),
        paste0("IS_CCI_SCORE_", condition_inputs$cond2),
        paste0("IS_CCI_SPECIFIC_", condition_inputs$cond1),
        paste0("IS_CCI_SPECIFIC_", condition_inputs$cond2),
        "IS_DE_LOGFC",
        "IS_DE_SIGNIFICANT",
        "DE_DIRECTION",
        paste0("IS_CCI_DETECTED_", condition_inputs$cond1),
        paste0("IS_CCI_DETECTED_", condition_inputs$cond2),
        "IS_CCI_DE"
      ) := {
        threshold_score_temp <- stats::quantile(
          x = c(
            .SD[get(paste0(
              "IS_CCI_EXPRESSED_",
              condition_inputs$cond1)) == TRUE][[paste0(
                "CCI_SCORE_",
                condition_inputs$cond1
                )]],
            .SD[get(paste0(
              "IS_CCI_EXPRESSED_",
              condition_inputs$cond2)
              ) == TRUE][[paste0(
                "CCI_SCORE_",
                condition_inputs$cond2
                )]]
          ),
          probs = threshold_quantile_score
        )
        is_cci_score_1 <- (get(paste0(
          "CCI_SCORE_",
          condition_inputs$cond1
          )) >= threshold_score_temp)
        is_cci_score_2 <- (get(paste0(
          "CCI_SCORE_",
          condition_inputs$cond2
          )) >= threshold_score_temp)
        is_cci_specific_1 <- get(paste0(
          "BH_P_VALUE_",
          condition_inputs$cond1
          )) <= threshold_p_value_specificity
        is_cci_specific_2 <- get(paste0(
          "BH_P_VALUE_",
          condition_inputs$cond2
          )) <= threshold_p_value_specificity
        is_de_logfc <- LOGFC_ABS >= threshold_logfc
        is_de_significant <- BH_P_VALUE_DE <= threshold_p_value_de
        de_direction <- fifelse(LOGFC > 0, "UP", "DOWN")
        is_cci_detected_1 <- (get(paste0(
          "IS_CCI_EXPRESSED_",
          condition_inputs$cond1
          )) == TRUE) &
          is_cci_score_1 & is_cci_specific_1
        is_cci_detected_2 <- (get(paste0(
          "IS_CCI_EXPRESSED_",
          condition_inputs$cond2
          )) == TRUE) &
          is_cci_score_2 & is_cci_specific_2
        is_cci_de <- is_de_logfc & is_de_significant
        list(
          is_cci_score_1, is_cci_score_2,
          is_cci_specific_1, is_cci_specific_2,
          is_de_logfc, is_de_significant, de_direction,
          is_cci_detected_1, is_cci_detected_2,
          is_cci_de
        )
      },
      by = temp_by
      ]
      cci_dt[, REGULATION :=
               ifelse(
                 !get(paste0("IS_CCI_DETECTED_", condition_inputs$cond1)) &
                   !get(paste0("IS_CCI_DETECTED_", condition_inputs$cond2)),
                 "NOT_DETECTED",
                 ifelse(
                   IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "UP",
                   "UP",
                   ifelse(
                     IS_DE_LOGFC & IS_DE_SIGNIFICANT & DE_DIRECTION == "DOWN",
                     "DOWN",
                     ifelse(
                       !IS_DE_LOGFC,
                       "FLAT",
                       ifelse(
                         IS_DE_LOGFC & !IS_DE_SIGNIFICANT,
                         "NSC",
                         "There is a problem here!"
                       )
                     )
                   )
                 )
               )
             ]
      if ("There is a problem here!" %in% cci_dt[["REGULATION"]]) {
        stop("Error when assigning regulation to CCIs.")
      }
      cci_dt <- cci_dt[REGULATION != "NOT_DETECTED"]
    }
  } else {
    cci_dt <- cci_dt[IS_CCI_EXPRESSED == TRUE]
    if (permutation_analysis) {
      cci_dt[
        ,
        BH_P_VALUE := stats::p.adjust(P_VALUE, method = "BH"),
        by = temp_by]
      cci_dt[, c(
        "IS_CCI_SCORE",
        "IS_CCI_SPECIFIC",
        "IS_CCI_DETECTED"
      ) := {
        threshold_score_temp <- stats::quantile(
          x = .SD[["CCI_SCORE"]],
          probs = threshold_quantile_score
        )
        is_cci_score <- CCI_SCORE >= threshold_score_temp
        is_cci_specific <- BH_P_VALUE <= threshold_p_value_specificity
        is_cci_detected <- IS_CCI_EXPRESSED & is_cci_score & is_cci_specific
        list(
          is_cci_score,
          is_cci_specific,
          is_cci_detected
        )
      },
      by = temp_by
      ]
    }
  }
  cci_dt[, ER_CELLTYPES := paste(
    EMITTER_CELLTYPE,
    RECEIVER_CELLTYPE,
    sep = "_")]
  cci_dt[, CCI := paste(ER_CELLTYPES, LRI, sep = "_")]
  return(cci_dt)
}

clean_colnames <- function(
  cci_dt,
  condition_inputs,
  max_nL,
  max_nR,
  permutation_analysis
) {
  first_cols <- c(
    "CCI",
    "ER_CELLTYPES",
    "EMITTER_CELLTYPE",
    "RECEIVER_CELLTYPE",
    "LRI")
  LR_COLNAMES <- c(
    paste0("LIGAND_", 1:max_nL),
    paste0("RECEPTOR_", 1:max_nR)
  )
  first_cols <- c(first_cols, LR_COLNAMES)
  if (!condition_inputs$is_cond) {
    second_cols <- c(
      "NCELLS_EMITTER",
      "NCELLS_RECEIVER",
      paste0("L", 1:max_nL, "_DETECTION_RATE"),
      paste0("R", 1:max_nR, "_DETECTION_RATE"),
      paste0("L", 1:max_nL, "_EXPRESSION"),
      paste0("R", 1:max_nR, "_EXPRESSION")
    )
    if (!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        second_cols,
        "CCI_SCORE", "IS_CCI_EXPRESSED"
      )
    } else {
      ordered_cols <- c(
        first_cols,
        second_cols,
        "CCI_SCORE", "P_VALUE", "BH_P_VALUE",
        "IS_CCI_EXPRESSED", "IS_CCI_SCORE", "IS_CCI_SPECIFIC",
        "IS_CCI_DETECTED"
      )
    }
  } else {
    second_cols <- c(
      paste0("NCELLS_EMITTER_", condition_inputs$cond1),
      paste0("NCELLS_EMITTER_", condition_inputs$cond2),
      paste0("NCELLS_RECEIVER_", condition_inputs$cond1),
      paste0("NCELLS_RECEIVER_", condition_inputs$cond2),
      paste0("L", 1:max_nL, "_DETECTION_RATE_", condition_inputs$cond1),
      paste0("L", 1:max_nL, "_DETECTION_RATE_", condition_inputs$cond2),
      paste0("R", 1:max_nR, "_DETECTION_RATE_", condition_inputs$cond1),
      paste0("R", 1:max_nR, "_DETECTION_RATE_", condition_inputs$cond2),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_inputs$cond1),
      paste0("L", 1:max_nL, "_EXPRESSION_", condition_inputs$cond2),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_inputs$cond1),
      paste0("R", 1:max_nR, "_EXPRESSION_", condition_inputs$cond2)
    )
    if (!permutation_analysis) {
      ordered_cols <- c(
        first_cols,
        second_cols,
        paste0("CCI_SCORE_", condition_inputs$cond1),
        paste0("CCI_SCORE_", condition_inputs$cond2),
        "LOGFC",
        "LOGFC_ABS",
        paste0("IS_CCI_EXPRESSED_", condition_inputs$cond1),
        paste0("IS_CCI_EXPRESSED_", condition_inputs$cond2)
      )
    } else {
      ordered_cols <- c(
        first_cols,
        second_cols,
        paste0("CCI_SCORE_", condition_inputs$cond1),
        paste0("CCI_SCORE_", condition_inputs$cond2),
        paste0("P_VALUE_", condition_inputs$cond1),
        paste0("BH_P_VALUE_", condition_inputs$cond1),
        paste0("P_VALUE_", condition_inputs$cond2),
        paste0("BH_P_VALUE_", condition_inputs$cond2),
        "LOGFC",
        "LOGFC_ABS",
        "P_VALUE_DE",
        "BH_P_VALUE_DE",
        paste0("IS_CCI_EXPRESSED_", condition_inputs$cond1),
        paste0("IS_CCI_EXPRESSED_", condition_inputs$cond2),
        paste0("IS_CCI_SCORE_", condition_inputs$cond1),
        paste0("IS_CCI_SCORE_", condition_inputs$cond2),
        paste0("IS_CCI_SPECIFIC_", condition_inputs$cond1),
        paste0("IS_CCI_SPECIFIC_", condition_inputs$cond2),
        paste0("IS_CCI_DETECTED_", condition_inputs$cond1),
        paste0("IS_CCI_DETECTED_", condition_inputs$cond2),
        "IS_DE_LOGFC",
        "IS_DE_SIGNIFICANT",
        "IS_CCI_DE",
        "DE_DIRECTION",
        "REGULATION"
      )
    }
  }
  setcolorder(
    x = cci_dt,
    neworder = ordered_cols
  )
  return(cci_dt)
}

get_table_cci <- function(
  object,
  type,
  simplified,
  class_signature
) {
  if (type == "raw") {
    table <- copy(object@cci_table_raw)
  }
  if (type == "detected") {
    table <- copy(object@cci_table_detected)
    condition_inputs <- list(
      is_cond = object@parameters$conditional_analysis,
      cond1 = object@parameters$seurat_condition_id$cond1_name,
      cond2 = object@parameters$seurat_condition_id$cond2_name
    )
    if (simplified) {
      if (class_signature == "scDiffComCombined") {
        first_cols <- c(
          "ID",
          "CCI",
          "ER_CELLTYPES",
          "LRI"
        )
      } else {
        first_cols <- c(
          "CCI",
          "ER_CELLTYPES",
          "LRI"
        )
      }
      if (!condition_inputs$is_cond) {
        second_cols <- c(
          "NCELLS_EMITTER",
          "NCELLS_RECEIVER"
        )
        if (!object@parameters$permutation_analysis) {
          ordered_cols <- c(
            first_cols,
            second_cols,
            "CCI_SCORE"
          )
        } else {
          ordered_cols <- c(
            first_cols,
            second_cols,
            "CCI_SCORE", "P_VALUE", "BH_P_VALUE",
            "IS_CCI_DETECTED"
          )
        }
      } else {
        second_cols <- c(
          paste0("NCELLS_EMITTER_", condition_inputs$cond1),
          paste0("NCELLS_EMITTER_", condition_inputs$cond2),
          paste0("NCELLS_RECEIVER_", condition_inputs$cond1),
          paste0("NCELLS_RECEIVER_", condition_inputs$cond2)
        )
        if (!object@parameters$permutation_analysis) {
          ordered_cols <- c(
            first_cols,
            second_cols,
            paste0("CCI_SCORE_", condition_inputs$cond1),
            paste0("CCI_SCORE_", condition_inputs$cond2),
            "LOGFC"
          )
        } else {
          ordered_cols <- c(
            first_cols,
            second_cols,
            paste0("CCI_SCORE_", condition_inputs$cond1),
            paste0("CCI_SCORE_", condition_inputs$cond2),
            paste0("P_VALUE_", condition_inputs$cond1),
            paste0("BH_P_VALUE_", condition_inputs$cond1),
            paste0("P_VALUE_", condition_inputs$cond2),
            paste0("BH_P_VALUE_", condition_inputs$cond2),
            "LOGFC",
            "P_VALUE_DE",
            "BH_P_VALUE_DE",
            paste0("IS_CCI_DETECTED_", condition_inputs$cond1),
            paste0("IS_CCI_DETECTED_", condition_inputs$cond2),
            "REGULATION"
          )
        }
      }
      table <- table[, ordered_cols, with = FALSE]
    }
  }
  table
}

