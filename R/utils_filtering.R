add_convenience_cols <- function(
  cci_dt,
  condition_info,
  log_scale,
  permutation_analysis,
  pre_filtering
) {
  LOGFC <- LOGFC_ABS <-
    BH_PVAL_DIFF <- PVAL_DIFF <-
    BH_PVAL <- PVAL <-
    LR_NAME <- LR_DETECTED <- LR_CELLTYPE <- L_CELLTYPE <- R_CELLTYPE <-
    LIGAND_1  <- LIGAND_2 <- RECEPTOR_1 <- RECEPTOR_2 <- RECEPTOR_3 <- NULL
  if(condition_info$is_cond) {
    if(pre_filtering) {
      cci_dt <- cci_dt[get(paste0("LR_DETECTED_", condition_info$cond1)) == TRUE | get(paste0("LR_DETECTED_", condition_info$cond2)) == TRUE]
    }
    if(log_scale) {
      cci_dt[, LOGFC := get(paste0("LR_SCORE_", condition_info$cond2)) - get(paste0("LR_SCORE_", condition_info$cond1))]
    } else {
      cci_dt[, LOGFC := log(get(paste0("LR_SCORE_", condition_info$cond2)) / get(paste0("LR_SCORE_", condition_info$cond1)))]
    }
    cci_dt[, LOGFC_ABS := abs(LOGFC)]
    if(permutation_analysis) {
      cci_dt[, BH_PVAL_DIFF := stats::p.adjust(PVAL_DIFF, method = "BH")]
      cci_dt[, paste0("BH_PVAL_", condition_info$cond1) := stats::p.adjust(get(paste0("PVAL_", condition_info$cond1)), method = "BH")]
      cci_dt[, paste0("BH_PVAL_", condition_info$cond2) := stats::p.adjust(get(paste0("PVAL_", condition_info$cond2)), method = "BH")]
    }
  } else {
    if(pre_filtering) {
      cci_dt <- cci_dt[LR_DETECTED == TRUE]
    }
    if(permutation_analysis) {
      cci_dt[, BH_PVAL := stats::p.adjust(PVAL, method = "BH")]
    }
  }
  cci_dt[, LR_NAME := list(sapply(1:nrow(.SD), function(i) {
    temp1 <- c(LIGAND_1[[i]], LIGAND_2[[i]])
    temp1 <- temp1[!is.na(temp1)]
    temp1 <- paste0(temp1, collapse = "_")
    temp2 <- c(RECEPTOR_1[[i]], RECEPTOR_2[[i]], RECEPTOR_3[[i]])
    temp2 <- temp2[!is.na(temp2)]
    temp2 <- paste0(temp2, collapse = "_")
    return(paste(temp1, temp2, sep = ":"))
  }))]
  cci_dt[, LR_CELLTYPE := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]
  return(cci_dt)
}

find_detected_cci <- function(
  cci_dt,
  condition_info,
  cutoff_quantile_score,
  cutoff_pval_specificity,
  cutoff_pval_de,
  cutoff_logfc
) {
  DIFFERENTIALLY_EXPRESSED <- BH_PVAL_DIFF <- LOGFC_ABS <- DIFFERENTIAL_DIRECTION <-
    LOGFC <- LR_DETECTED <- LR_SCORE <- BH_PVAL <- NULL
  if(condition_info$is_cond) {
    cutoff_score <- stats::quantile(
      x = c(
        cci_dt[[paste0("LR_SCORE_", condition_info$cond1)]],
        cci_dt[[paste0("LR_SCORE_", condition_info$cond2)]]
        ),
      probs = cutoff_quantile_score
      )
    cci_dt[, paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1) :=
             (get(paste0("LR_DETECTED_", condition_info$cond1)) == TRUE)
           & (get(paste0("LR_SCORE_", condition_info$cond1)) >=  cutoff_score)
           & (get(paste0("BH_PVAL_", condition_info$cond1)) <= cutoff_pval_specificity)
           ]
    cci_dt[, paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2) :=
             (get(paste0("LR_DETECTED_", condition_info$cond2)) == TRUE)
           & (get(paste0("LR_SCORE_", condition_info$cond2)) >=  cutoff_score)
           & (get(paste0("BH_PVAL_", condition_info$cond2)) <= cutoff_pval_specificity)
           ]
    cci_dt[, DIFFERENTIALLY_EXPRESSED :=
           (get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1))
            | get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)))
         & (BH_PVAL_DIFF <= cutoff_pval_de)
         & (LOGFC_ABS >= cutoff_logfc)
         ]
    cci_dt[, DIFFERENTIAL_DIRECTION := fifelse(LOGFC > 0, "UP", "DOWN")]
  } else {
    cutoff_score <- stats::quantile(
      x = cci_dt[["LR_SCORE"]],
      probs = cutoff_quantile_score
    )
    cci_dt[, "LR_DETECTED_AND_SIGNIFICANT" :=
             (LR_DETECTED == TRUE)
           & (LR_SCORE >=  cutoff_score)
           & (BH_PVAL <= cutoff_pval_specificity)
           ]
  }
  return(cci_dt)
}

assign_regulation <- function(
  cci_dt,
  condition_info,
  cutoff_quantile_score,
  cutoff_pval_specificity
) {
  REGULATION <- REGULATION_SIMPLE <- DIFFERENTIALLY_EXPRESSED <- DIFFERENTIAL_DIRECTION <- NULL
  cutoff_score <- stats::quantile(
    x = c(
      cci_dt[[paste0("LR_SCORE_", condition_info$cond1)]],
      cci_dt[[paste0("LR_SCORE_", condition_info$cond2)]]
    ),
    probs = cutoff_quantile_score
  )
  cci_dt[, REGULATION := ifelse(
    get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
      get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
      DIFFERENTIALLY_EXPRESSED,
    ifelse(
      DIFFERENTIAL_DIRECTION == "UP",
      "UP",
      "DOWN"
    ),
    ifelse(
      get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
        get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
        !DIFFERENTIALLY_EXPRESSED,
      ifelse(
        DIFFERENTIAL_DIRECTION == "UP",
        "FLAT",
        "FLAT"#TTFD"
      ),
      ifelse(
        get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
          !get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
          DIFFERENTIALLY_EXPRESSED,
        ifelse(
          DIFFERENTIAL_DIRECTION == "DOWN",
          "DOWN_DISAPPEARS",
          ifelse(
            sum(c(get(paste0("BH_PVAL_", condition_info$cond2)) > cutoff_pval_specificity,
                  !get(paste0("LR_DETECTED_", condition_info$cond2)),
                  get(paste0("LR_SCORE_", condition_info$cond2)) < cutoff_score
            )
            ) == 1,
            "UP",
            "NON_DETECTED"
          )
        ),
        ifelse(
          get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
            !get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
            !DIFFERENTIALLY_EXPRESSED,
          ifelse(
            sum(c(get(paste0("BH_PVAL_", condition_info$cond2)) > cutoff_pval_specificity,
                  !get(paste0("LR_DETECTED_", condition_info$cond2)),
                  get(paste0("LR_SCORE_", condition_info$cond2)) < cutoff_score
            )
            ) == 1,
            ifelse(
              DIFFERENTIAL_DIRECTION == "UP",
              "FLAT", #"TTFU",
              "FLAT" #"TTFD"
            ),
            "NON_DETECTED"
          ),
          ifelse(
            !get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
              get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
              DIFFERENTIALLY_EXPRESSED,
            ifelse(
              DIFFERENTIAL_DIRECTION == "UP",
              "UP_APPEARS",
              ifelse(
                sum(c(get(paste0("BH_PVAL_", condition_info$cond1)) > cutoff_pval_specificity,
                      !get(paste0("LR_DETECTED_", condition_info$cond1)),
                      get(paste0("LR_SCORE_", condition_info$cond1)) < cutoff_score
                )
                ) == 1,
                "DOWN",
                "NON_DETECTED"
              )
            ),
            ifelse(
              !get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond1)) &
                get(paste0("LR_DETECTED_AND_SIGNIFICANT_IN_", condition_info$cond2)) &
                !DIFFERENTIALLY_EXPRESSED,
              ifelse(
                sum(c(get(paste0("BH_PVAL_", condition_info$cond1)) > cutoff_pval_specificity,
                      !get(paste0("LR_DETECTED_", condition_info$cond1)),
                      get(paste0("LR_SCORE_", condition_info$cond1)) < cutoff_score
                )
                ) == 1,
                ifelse(
                  DIFFERENTIAL_DIRECTION == "UP",
                  "FLAT", #"TTFU",
                  "FLAT" #"TTFD"
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
          #stop("Problem of classification.")
        )
      )
    )
  )]
  if("OTHER" %in% cci_dt$REGULATION_SIMPLE) stop("Problem of classification.")
  return(cci_dt)
}
