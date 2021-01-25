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
    temp_param$threshold_quantile_score <- new_threshold_quantile_score
  }
  if (!is.null(new_threshold_p_value_specificity)) {
    temp_param$threshold_p_value_specificity <- new_threshold_p_value_specificity
  }
  if (!is.null(new_threshold_p_value_de)) {
    temp_param$threshold_p_value_de <- new_threshold_p_value_de
  }
  if (!is.null(new_threshold_logfc)) {
    temp_param$threshold_logfc <- new_threshold_logfc
  }
  condition_inputs <- list(
    is_cond = temp_param$conditional_analysis,
    cond1 = temp_param$cond1_name,
    cond2 = temp_param$cond2_name
  )
  object@parameters <- temp_param
  #validate object?
  cci_dt <- copy(x = object@cci_raw)
  if (verbose) message("Filtering and cleaning `raw` CCIs.")
  cci_dt <- preprocess_cci_raw(
    cci_dt = cci_dt,
    condition_inputs = condition_inputs,
    log_scale = temp_param$log_scale,
    permutation_analysis = temp_param$permutation_analysis,
    threshold_min_cells = temp_param$threshold_min_cells,
    threshold_pct = temp_param$threshold_pct,
    max_nL = temp_param$max_nL,
    max_nR = temp_param$max_nR,
    class_signature = class_signature
  )
  if (temp_param$permutation_analysis) {
    cci_dt <- find_detected_cci(
      cci_dt = cci_dt,
      condition_inputs = condition_inputs,
      threshold_quantile_score = temp_param$threshold_quantile_score,
      threshold_p_value_specificity = temp_param$threshold_p_value_specificity,
      threshold_p_value_de = temp_param$threshold_p_value_de,
      threshold_logfc = temp_param$threshold_logfc,
      class_signature = class_signature
    )

    if (condition_inputs$is_cond) {
      #if (verbose) message("Assigning precise differential regulation.")
      # cci_dt <- assign_regulation(
      #   cci_dt = cci_dt,
      #   condition_inputs = condition_inputs,
      #   threshold_quantile_score = temp_param$threshold_quantile_score,
      #   threshold_p_value_specificity = temp_param$threshold_p_value_specificity,
      #   class_signature = class_signature
      # )
      ##
      #cci_dt <- cci_dt[REGULATION_SIMPLE != "NON_DETECTED"]
    }
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
      " detected CCIs. (Note: incomplete filtering as statistical test has not been performed!)"
    )
    if (verbose) message(mes)
  }
  cci_dt <- clean_colnames(
    cci_dt = cci_dt,
    max_nL = temp_param$max_nL,
    max_nR = temp_param$max_nR,
    condition_inputs = condition_inputs,
    permutation_analysis = temp_param$permutation_analysis#,
    #class_signature = class_signature
  )
  skip_ora <- TRUE
  if (nrow(cci_dt) == 0) {
    object@cci_detected <- list()
    object@ora_default <- list()
  } else {
    object@cci_detected <- cci_dt
    #validate object?
    if (skip_ora) {
      object@ora_default <- list()
    } else {
      object <- run_ora(
        object = object,
        categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS"),
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
          categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID"),
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

preprocess_cci_raw <- function(
  cci_dt,
  condition_inputs,
  log_scale,
  permutation_analysis,
  threshold_min_cells,
  threshold_pct,
  max_nL,
  max_nR,
  class_signature
) {
  BH_P_VALUE_DE <- P_VALUE_DE <-
    BH_P_VALUE <- P_VALUE <-
    CCI_DETECTED <- ER_CELLTYPES <- EMITTER_CELLTYPE <- RECEIVER_CELLTYPE <-
    LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  if (class_signature == "scDiffCom") {
    temp_by <- NULL
  }
  if (class_signature == "scDiffComCombined") {
    temp_by <- "ID"
  }
  if (condition_inputs$is_cond) {
    cci_dt <- cci_dt[
      get(paste0("CCI_DETECTED_", condition_inputs$cond1)) == TRUE |
        get(paste0("CCI_DETECTED_", condition_inputs$cond2)) == TRUE
      ]

    ####
    #additional columns, to remove?


    cci_dt[, paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond1) :=
             get(paste0("EMITTER_NCELLS_", condition_inputs$cond1))
           *
             do.call(
               pmin,
               c(lapply(1:max_nL, function(i) {
                 get(paste0("L", i, "_DETECTED_", condition_inputs$cond1))
               }), na.rm = TRUE)
             )
          ]
    cci_dt[, paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond2) :=
             get(paste0("EMITTER_NCELLS_", condition_inputs$cond2))
           *
             do.call(
               pmin,
               c(lapply(1:max_nL, function(i) {
                 get(paste0("L", i, "_DETECTED_", condition_inputs$cond2))
               }), na.rm = TRUE)
             )
           ]
    cci_dt[, paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond1) :=
             get(paste0("RECEIVER_NCELLS_", condition_inputs$cond1))
           *
             do.call(
               pmin,
               c(lapply(1:max_nR, function(i) {
                 get(paste0("R", i, "_DETECTED_", condition_inputs$cond1))
               }), na.rm = TRUE)
             )
           ]
    cci_dt[, paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond2) :=
             get(paste0("RECEIVER_NCELLS_", condition_inputs$cond2))
           *
             do.call(
               pmin,
               c(lapply(1:max_nR, function(i) {
                 get(paste0("R", i, "_DETECTED_", condition_inputs$cond2))
               }), na.rm = TRUE)
             )
           ]
    cci_dt[,  paste0("IS_EXPRESSED_", condition_inputs$cond1) :=
             ifelse(
               (
                 get(paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond1)) < threshold_min_cells &
                   (do.call(
                     pmin,
                     c(lapply(1:max_nL, function(i) {
                       get(paste0("L", i, "_DETECTED_", condition_inputs$cond1))
                     }), na.rm = TRUE)
                   ) >= 0.1)
               ) |
                 (
                   get(paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond1)) < threshold_min_cells &
                     (do.call(
                       pmin,
                       c(lapply(1:max_nR, function(i) {
                         get(paste0("R", i, "_DETECTED_", condition_inputs$cond1))
                       }), na.rm = TRUE)
                     ) >= 0.1)
                 ),
               "UN",
               ifelse(
                 (
                   get(paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond1)) >= threshold_min_cells &
                     (do.call(
                       pmin,
                       c(lapply(1:max_nL, function(i) {
                         get(paste0("L", i, "_DETECTED_", condition_inputs$cond1))
                       }), na.rm = TRUE)
                     ) >= 0.1)
                 ) &
                   (
                     get(paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond1)) >= threshold_min_cells &
                       (do.call(
                         pmin,
                         c(lapply(1:max_nR, function(i) {
                           get(paste0("R", i, "_DETECTED_", condition_inputs$cond1))
                         }), na.rm = TRUE)
                       ) >= 0.1)
                   ),
                 "YES",
                 "NO"
               )
             )]
    cci_dt[,  paste0("IS_EXPRESSED_", condition_inputs$cond2) :=
             ifelse(
               (
                 get(paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond2)) < threshold_min_cells &
                   (do.call(
                     pmin,
                     c(lapply(1:max_nL, function(i) {
                       get(paste0("L", i, "_DETECTED_", condition_inputs$cond2))
                     }), na.rm = TRUE)
                   ) >= 0.1)
               ) |
                 (
                   get(paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond2)) < threshold_min_cells &
                     (do.call(
                       pmin,
                       c(lapply(1:max_nR, function(i) {
                         get(paste0("R", i, "_DETECTED_", condition_inputs$cond2))
                       }), na.rm = TRUE)
                     ) >= 0.1)
                 ),
               "UN",
               ifelse(
                 (
                   get(paste0("EMITTER_NCELLS_EXPR_", condition_inputs$cond2)) >= threshold_min_cells &
                     (do.call(
                       pmin,
                       c(lapply(1:max_nL, function(i) {
                         get(paste0("L", i, "_DETECTED_", condition_inputs$cond2))
                       }), na.rm = TRUE)
                     ) >= 0.1)
                 ) &
                   (
                     get(paste0("RECEIVER_NCELLS_EXPR_", condition_inputs$cond2)) >= threshold_min_cells &
                       (do.call(
                         pmin,
                         c(lapply(1:max_nR, function(i) {
                           get(paste0("R", i, "_DETECTED_", condition_inputs$cond2))
                         }), na.rm = TRUE)
                       ) >= 0.1)
                   ),
                 "YES",
                 "NO"
               )
             )]
    ####

    if (permutation_analysis) {
      cci_dt[, BH_P_VALUE_DE := stats::p.adjust(P_VALUE_DE, method = "BH"), by = temp_by]
      cci_dt[, paste0("BH_P_VALUE_", condition_inputs$cond1) :=
               stats::p.adjust(get(paste0("P_VALUE_", condition_inputs$cond1)), method = "BH"), by = temp_by]
      cci_dt[, paste0("BH_P_VALUE_", condition_inputs$cond2) :=
               stats::p.adjust(get(paste0("P_VALUE_", condition_inputs$cond2)), method = "BH"), by = temp_by]
    }
  } else {
    cci_dt <- cci_dt[CCI_DETECTED == TRUE]
    if (permutation_analysis) {
      cci_dt[, BH_P_VALUE := stats::p.adjust(P_VALUE, method = "BH"), by = temp_by]
    }
  }
  cci_dt[, ER_CELLTYPES := paste(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, sep = "_")]
  return(cci_dt)
}

find_detected_cci <- function(
  cci_dt,
  condition_inputs,
  threshold_quantile_score,
  threshold_p_value_specificity,
  threshold_p_value_de,
  threshold_logfc,
  class_signature
) {
  DIFFERENTIALLY_EXPRESSED <- BH_P_VALUE_DE <- LOGFC_ABS <- DIFFERENTIAL_DIRECTION <-
    LOGFC <- CCI_DETECTED <- CCI_SCORE <- BH_P_VALUE <- NULL
  if (class_signature == "scDiffCom") {
    temp_by <- NULL
  }
  if (class_signature == "scDiffComCombined") {
    temp_by <- "ID"
  }
  if (condition_inputs$is_cond) {
    cci_dt[, c(
      paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1),
      paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2),
      paste0("IS_SCORE_", condition_inputs$cond1),
      paste0("IS_SCORE_", condition_inputs$cond2),
      paste0("IS_SPECIFIC_", condition_inputs$cond1),
      paste0("IS_SPECIFIC_", condition_inputs$cond2)
    ) := {
      threshold_score_temp <- stats::quantile(
        x = c(
          .SD[get(paste0("CCI_DETECTED_", condition_inputs$cond1)) == TRUE][[paste0("CCI_SCORE_", condition_inputs$cond1)]],
          .SD[get(paste0("CCI_DETECTED_", condition_inputs$cond2)) == TRUE][[paste0("CCI_SCORE_", condition_inputs$cond2)]]
        ),
        probs = threshold_quantile_score
      )
      res1 <- (get(paste0("CCI_DETECTED_", condition_inputs$cond1)) == TRUE) &
        (get(paste0("CCI_SCORE_", condition_inputs$cond1)) >= threshold_score_temp) &
        (get(paste0("BH_P_VALUE_", condition_inputs$cond1)) <= threshold_p_value_specificity)
      res2 <- (get(paste0("CCI_DETECTED_", condition_inputs$cond2)) == TRUE) &
        (get(paste0("CCI_SCORE_", condition_inputs$cond2)) >= threshold_score_temp) &
        (get(paste0("BH_P_VALUE_", condition_inputs$cond2)) <= threshold_p_value_specificity)
      res3 <- (get(paste0("CCI_SCORE_", condition_inputs$cond1)) >= threshold_score_temp)
      res4 <- (get(paste0("CCI_SCORE_", condition_inputs$cond2)) >= threshold_score_temp)
      res5 <- (get(paste0("BH_P_VALUE_", condition_inputs$cond1)) <= threshold_p_value_specificity)
      res6 <- (get(paste0("BH_P_VALUE_", condition_inputs$cond2)) <= threshold_p_value_specificity)
      list(
        res1,
        res2,
        res3,
        res4,
        res5,
        res6
      )
    },
    by = temp_by
    ]
    cci_dt[, IS_LOGFC := (LOGFC_ABS >= threshold_logfc)]
    cci_dt[, IS_DE := (BH_P_VALUE_DE <= threshold_p_value_de)]
    cci_dt[, DIFFERENTIAL_DIRECTION := fifelse(LOGFC > 0, "UP", "DOWN")]
    cci_dt[, DIFFERENTIALLY_EXPRESSED :=
             (get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1))
              | get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)))
           & (BH_P_VALUE_DE <= threshold_p_value_de)
           & (LOGFC_ABS >= threshold_logfc)]
    } else {
    cci_dt[, "CCI_DETECTED_AND_SIGNIFICANT" := {
      threshold_score_temp <- stats::quantile(
        x = .SD[["CCI_SCORE"]],
        probs = threshold_quantile_score
      )
      res <- (CCI_DETECTED == TRUE) &
        (CCI_SCORE >= threshold_score_temp) &
        (BH_P_VALUE <= threshold_p_value_specificity)
    },
    by = temp_by
    ]
  }
  return(cci_dt)
}

assign_regulation <- function(
  cci_dt,
  condition_inputs,
  threshold_quantile_score,
  threshold_p_value_specificity,
  class_signature
) {
  REGULATION <- REGULATION_SIMPLE <- DIFFERENTIALLY_EXPRESSED <- DIFFERENTIAL_DIRECTION <- NULL
  if (class_signature == "scDiffCom") {
    temp_by <- NULL
  }
  if (class_signature == "scDiffComCombined") {
    temp_by <- "ID"
  }
  cci_dt[,
         REGULATION := {
           threshold_score_temp <- stats::quantile(
             x = c(
               .SD[[paste0("CCI_SCORE_", condition_inputs$cond1)]],
               .SD[[paste0("CCI_SCORE_", condition_inputs$cond2)]]
             ),
             probs = threshold_quantile_score
           )
           res <- ifelse(
             get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
               get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
               DIFFERENTIALLY_EXPRESSED,
             ifelse(
               DIFFERENTIAL_DIRECTION == "UP",
               "UP",
               "DOWN"
             ),
             ifelse(
               get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                 get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                 !DIFFERENTIALLY_EXPRESSED,
               ifelse(
                 DIFFERENTIAL_DIRECTION == "UP",
                 "FLAT",
                 "FLAT" # TTFD"
               ),
               ifelse(
                 get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                   !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                   DIFFERENTIALLY_EXPRESSED,
                 ifelse(
                   DIFFERENTIAL_DIRECTION == "DOWN",
                   "DOWN_DISAPPEARS",
                   ifelse(
                       (get(paste0("BH_P_VALUE_", condition_inputs$cond2)) > threshold_p_value_specificity +
                          (!get(paste0("CCI_DETECTED_", condition_inputs$cond2))) +
                          (get(paste0("CCI_SCORE_", condition_inputs$cond2)) < threshold_score_temp))
                        == 1,
                     "UP",
                     "NON_DETECTED"
                   )
                 ),
                 ifelse(
                   get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                     !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                     !DIFFERENTIALLY_EXPRESSED,
                   ifelse(
                     ((get(paste0("BH_P_VALUE_", condition_inputs$cond2)) > threshold_p_value_specificity) +
                        (!get(paste0("CCI_DETECTED_", condition_inputs$cond2))) +
                        (get(paste0("CCI_SCORE_", condition_inputs$cond2)) < threshold_score_temp))
                      == 1,
                     ifelse(
                       DIFFERENTIAL_DIRECTION == "UP",
                       "FLAT", # "TTFU",
                       "FLAT" # "TTFD"
                     ),
                     "NON_DETECTED"
                   ),
                   ifelse(
                     !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                       get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                       DIFFERENTIALLY_EXPRESSED,
                     ifelse(
                       DIFFERENTIAL_DIRECTION == "UP",
                       "UP_APPEARS",
                       ifelse(
                         ((get(paste0("BH_P_VALUE_", condition_inputs$cond1)) > threshold_p_value_specificity) +
                            (!get(paste0("CCI_DETECTED_", condition_inputs$cond1))) +
                            (get(paste0("CCI_SCORE_", condition_inputs$cond1)) < threshold_score_temp))
                         == 1,
                         "DOWN",
                         "NON_DETECTED"
                       )
                     ),
                     ifelse(
                       !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                         get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                         !DIFFERENTIALLY_EXPRESSED,
                       ifelse(
                         ((get(paste0("BH_P_VALUE_", condition_inputs$cond1)) > threshold_p_value_specificity) +
                            (!get(paste0("CCI_DETECTED_", condition_inputs$cond1))) +
                            (get(paste0("CCI_SCORE_", condition_inputs$cond1)) < threshold_score_temp))
                         == 1,
                         ifelse(
                           DIFFERENTIAL_DIRECTION == "UP",
                           "FLAT", # "TTFU",
                           "FLAT" # "TTFD"
                         ),
                         "NON_DETECTED"
                       ),
                       "NON_DETECTED"
                     )
                   )
                 )
               )
             )
           )
         },
         by = temp_by
         ]
  cci_dt[, REGULATION_SIMPLE := ifelse(
    REGULATION %in% c("UP", "UP_APPEARS"),
    "UP",
    ifelse(
      REGULATION %in% c("DOWN", "DOWN_DISAPPEARS"),
      "DOWN",
      ifelse(
        REGULATION == "FLAT",
        "FLAT",
        ifelse(
          REGULATION == "NON_DETECTED",
          "NON_DETECTED",
          "OTHER"
        )
      )
    )
  )]
  if ("OTHER" %in% cci_dt$REGULATION_SIMPLE) stop("Problem of classification.")
  return(cci_dt)
}
