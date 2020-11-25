preprocess_cci_raw <- function(
  cci_dt,
  condition_inputs,
  log_scale,
  permutation_analysis
) {
    BH_P_VALUE_DE <- P_VALUE_DE <-
    BH_P_VALUE <- P_VALUE <-
    LR_GENES <- CCI_DETECTED <- ER_CELLTYPES <- EMITTER_CELLTYPE <- RECEIVER_CELLTYPE <-
    LIGAND_1 <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  if (condition_inputs$is_cond) {
    cci_dt <- cci_dt[
      get(paste0("CCI_DETECTED_", condition_inputs$cond1)) == TRUE |
        get(paste0("CCI_DETECTED_", condition_inputs$cond2)) == TRUE
      ]
    if (permutation_analysis) {
      cci_dt[, BH_P_VALUE_DE := stats::p.adjust(P_VALUE_DE, method = "BH")]
      cci_dt[, paste0("BH_P_VALUE_", condition_inputs$cond1) :=
               stats::p.adjust(get(paste0("P_VALUE_", condition_inputs$cond1)), method = "BH")]
      cci_dt[, paste0("BH_P_VALUE_", condition_inputs$cond2) :=
               stats::p.adjust(get(paste0("P_VALUE_", condition_inputs$cond2)), method = "BH")]
    }
  } else {
    cci_dt <- cci_dt[CCI_DETECTED == TRUE]
    if (permutation_analysis) {
      cci_dt[, BH_P_VALUE := stats::p.adjust(P_VALUE, method = "BH")]
    }
  }
  cci_dt[, LR_GENES := list(sapply(1:nrow(.SD), function(i) {
    temp1 <- c(LIGAND_1[[i]], LIGAND_2[[i]])
    temp1 <- temp1[!is.na(temp1)]
    temp1 <- paste0(temp1, collapse = "_")
    temp2 <- c(RECEPTOR_1[[i]], RECEPTOR_2[[i]], RECEPTOR_3[[i]])
    temp2 <- temp2[!is.na(temp2)]
    temp2 <- paste0(temp2, collapse = "_")
    return(paste(temp1, temp2, sep = ":"))
  }))]
  cci_dt[, ER_CELLTYPES := paste(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, sep = "_")]
  return(cci_dt)
}

find_detected_cci <- function(
  cci_dt,
  condition_inputs,
  threshold_quantile_score,
  threshold_p_value_specificity,
  threshold_p_value_de,
  threshold_logfc
) {
  DIFFERENTIALLY_EXPRESSED <- BH_P_VALUE_DE <- LOGFC_ABS <- DIFFERENTIAL_DIRECTION <-
    LOGFC <- CCI_DETECTED <- CCI_SCORE <- BH_P_VALUE <- NULL
  if (condition_inputs$is_cond) {
    threshold_score <- stats::quantile(
      x = c(
        cci_dt[[paste0("CCI_SCORE_", condition_inputs$cond1)]],
        cci_dt[[paste0("CCI_SCORE_", condition_inputs$cond2)]]
      ),
      probs = threshold_quantile_score
    )
    cci_dt[, paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1) :=
             (get(paste0("CCI_DETECTED_", condition_inputs$cond1)) == TRUE)
           & (get(paste0("CCI_SCORE_", condition_inputs$cond1)) >= threshold_score)
           & (get(paste0("BH_P_VALUE_", condition_inputs$cond1)) <= threshold_p_value_specificity)]
    cci_dt[, paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2) :=
             (get(paste0("CCI_DETECTED_", condition_inputs$cond2)) == TRUE)
           & (get(paste0("CCI_SCORE_", condition_inputs$cond2)) >= threshold_score)
           & (get(paste0("BH_P_VALUE_", condition_inputs$cond2)) <= threshold_p_value_specificity)]
    cci_dt[, DIFFERENTIALLY_EXPRESSED :=
             (get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1))
              | get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)))
           & (BH_P_VALUE_DE <= threshold_p_value_de)
           & (LOGFC_ABS >= threshold_logfc)]
    cci_dt[, DIFFERENTIAL_DIRECTION := fifelse(LOGFC > 0, "UP", "DOWN")]
  } else {
    threshold_score <- stats::quantile(
      x = cci_dt[["CCI_SCORE"]],
      probs = threshold_quantile_score
    )
    cci_dt[, "CCI_DETECTED_AND_SIGNIFICANT" :=
             (CCI_DETECTED == TRUE)
           & (CCI_SCORE >= threshold_score)
           & (BH_P_VALUE <= threshold_p_value_specificity)]
  }
  return(cci_dt)
}

assign_regulation <- function(
  cci_dt,
  condition_inputs,
  threshold_quantile_score,
  threshold_p_value_specificity
) {
  REGULATION <- REGULATION_SIMPLE <- DIFFERENTIALLY_EXPRESSED <- DIFFERENTIAL_DIRECTION <- NULL
  threshold_score <- stats::quantile(
    x = c(
      cci_dt[[paste0("CCI_SCORE_", condition_inputs$cond1)]],
      cci_dt[[paste0("CCI_SCORE_", condition_inputs$cond2)]]
    ),
    probs = threshold_quantile_score
  )
  cci_dt[, REGULATION := ifelse(
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
            sum(c(
              get(paste0("BH_P_VALUE_", condition_inputs$cond2)) > threshold_p_value_specificity,
              !get(paste0("CCI_DETECTED_", condition_inputs$cond2)),
              get(paste0("CCI_SCORE_", condition_inputs$cond2)) < threshold_score
            )) == 1,
            "UP",
            "NON_DETECTED"
          )
        ),
        ifelse(
          get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
            !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
            !DIFFERENTIALLY_EXPRESSED,
          ifelse(
            sum(c(
              get(paste0("BH_P_VALUE_", condition_inputs$cond2)) > threshold_p_value_specificity,
              !get(paste0("CCI_DETECTED_", condition_inputs$cond2)),
              get(paste0("CCI_SCORE_", condition_inputs$cond2)) < threshold_score
            )) == 1,
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
                sum(c(
                  get(paste0("BH_P_VALUE_", condition_inputs$cond1)) > threshold_p_value_specificity,
                  !get(paste0("CCI_DETECTED_", condition_inputs$cond1)),
                  get(paste0("CCI_SCORE_", condition_inputs$cond1)) < threshold_score
                )) == 1,
                "DOWN",
                "NON_DETECTED"
              )
            ),
            ifelse(
              !get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond1)) &
                get(paste0("CCI_DETECTED_AND_SIGNIFICANT_IN_", condition_inputs$cond2)) &
                !DIFFERENTIALLY_EXPRESSED,
              ifelse(
                sum(c(
                  get(paste0("BH_P_VALUE_", condition_inputs$cond1)) > threshold_p_value_specificity,
                  !get(paste0("CCI_DETECTED_", condition_inputs$cond1)),
                  get(paste0("CCI_SCORE_", condition_inputs$cond1)) < threshold_score
                )) == 1,
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
  )]
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
          # stop("Problem of classification.")
        )
      )
    )
  )]
  if ("OTHER" %in% cci_dt$REGULATION_SIMPLE) stop("Problem of classification.")
  return(cci_dt)
}
