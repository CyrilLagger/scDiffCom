build_ora_dt <- function(
  cci_detected,
  logfc_threshold,
  regulation,
  category,
  species
) {
  COUNTS_VALUE_REGULATED <- COUNTS_NOTVALUE_REGULATED <-
    COUNTS_NOTVALUE_NOTREGULATED <- COUNTS_VALUE_NOTREGULATED <-
    COUNTS_VALUE_REGULATED_temp <- COUNTS_VALUE_NOTREGULATED_temp<- CATEGORY <- NULL
  cci_dt <- data.table::copy(
    cci_detected
  )
  ora_tables <- lapply(
    X = regulation,
    FUN = function(reg) {
      if ("GO_TERMS" == category) {
        counts_dt <- extract_category_counts(
          cci_detected = cci_dt,
          category = "LR_SORTED",
          reg = reg,
          logfc_threshold = logfc_threshold
        )
        if (species == "mouse") {
          go_intersection_dt <- scDiffCom::LRdb_mouse$LRdb_curated_GO
        }
        if (species == "human") {
          go_intersection_dt <- scDiffCom::LRdb_human$LRdb_curated_GO
        }
        counts_intersection_dt <- data.table::merge.data.table(
          go_intersection_dt,
          counts_dt,
          by.x = "LR_SORTED",
          by.y = "VALUE"
        )
        counts_intersection_dt[, c("COUNTS_VALUE_REGULATED_temp", "COUNTS_VALUE_NOTREGULATED_temp") :=
                                 list(
                                   sum(COUNTS_VALUE_REGULATED),
                                   sum(COUNTS_VALUE_NOTREGULATED)
                                 ),
                               by = "GO_ID"
                               ]
        counts_intersection_dt[, c("COUNTS_NOTVALUE_REGULATED_temp", "COUNTS_NOTVALUE_NOTREGULATED_temp") := list(
          COUNTS_NOTVALUE_REGULATED + COUNTS_VALUE_REGULATED - COUNTS_VALUE_REGULATED_temp,
          COUNTS_NOTVALUE_NOTREGULATED + COUNTS_VALUE_NOTREGULATED - COUNTS_VALUE_NOTREGULATED_temp
        )]
        counts_intersection_dt[, CATEGORY := "GO_intersection"]
        counts_intersection_dt[, c(
          "COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED",
          "COUNTS_NOTVALUE_REGULATED", "COUNTS_NOTVALUE_NOTREGULATED",
          "LR_SORTED"
        ) := NULL]
        counts_intersection_dt <- unique(counts_intersection_dt)
        data.table::setnames(
          counts_intersection_dt,
          old = c("GO_ID", "GO_NAME", "COUNTS_VALUE_REGULATED_temp", "COUNTS_VALUE_NOTREGULATED_temp",
                  "COUNTS_NOTVALUE_REGULATED_temp", "COUNTS_NOTVALUE_NOTREGULATED_temp"),
          new = c("VALUE_BIS", "VALUE", "COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED",
                  "COUNTS_NOTVALUE_REGULATED", "COUNTS_NOTVALUE_NOTREGULATED")
        )
        counts_dt <- counts_intersection_dt
      } else {
        counts_dt <- extract_category_counts(
          cci_detected = cci_dt,
          category = category,
          reg = reg,
          logfc_threshold = logfc_threshold
        )
      }
      counts_dt <- perform_ora_from_counts(
        counts_dt = counts_dt
      )
      cols_to_rename <- colnames(counts_dt)
      cols_to_rename <- cols_to_rename[!(cols_to_rename %in% c("VALUE", "VALUE_BIS", "CATEGORY"))]
      data.table::setnames(
        x = counts_dt,
        old = cols_to_rename,
        new = paste0(cols_to_rename, "_", reg)
      )
      return(counts_dt)
    }
  )
  ora_full <- Reduce(merge, ora_tables)
  return(ora_full)
}

extract_category_counts <- function(
  cci_detected,
  category,
  reg,
  logfc_threshold
) {
  REGULATION_SIMPLE <- COUNTS_VALUE_REGULATED <- COUNTS_VALUE_NOTREGULATED <- LOGFC <- NULL
  if (reg == "UP") {
    dt_regulated <- cci_detected[REGULATION_SIMPLE == "UP" & LOGFC >= logfc_threshold]
    dt_notregulated <- cci_detected[!(REGULATION_SIMPLE == "UP" & LOGFC >= logfc_threshold)]
  } else if (reg == "DOWN") {
    dt_regulated <- cci_detected[REGULATION_SIMPLE == "DOWN" & LOGFC <= -logfc_threshold]
    dt_notregulated <- cci_detected[!(REGULATION_SIMPLE == "DOWN" & LOGFC <= -logfc_threshold)]
  } else if (reg == "DIFF") {
    dt_regulated <- cci_detected[REGULATION_SIMPLE %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold]
    dt_notregulated <- cci_detected[!(REGULATION_SIMPLE %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold)]
  } else if (reg == "FLAT") {
    dt_regulated <- cci_detected[REGULATION_SIMPLE == "FLAT"]
    dt_notregulated <- cci_detected[REGULATION_SIMPLE != "FLAT"]
  } else {
    stop("Type of ORA not supported.")
  }
  dt_counts <- data.table::rbindlist(
    l = lapply(
      X = category,
      FUN = function(categ) {
        data.table::merge.data.table(
          x = dt_regulated[, list(CATEGORY = categ, COUNTS_VALUE_REGULATED = .N), by = list(VALUE = get(categ))],
          y = dt_notregulated[, list(CATEGORY = categ, COUNTS_VALUE_NOTREGULATED = .N), by = list(VALUE = get(categ))],
          by = c("VALUE", "CATEGORY"),
          all = TRUE
        )
      }
    ),
    use.names = TRUE
  )
  data.table::setnafill(
    x = dt_counts,
    type = "const",
    fill = 0,
    cols = c("COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED")
  )
  COUNTS_REGULATED <- dt_regulated[, .N]
  COUNTS_NOTREGULATED <- dt_notregulated[, .N]
  dt_counts <- dt_counts[, c("COUNTS_NOTVALUE_REGULATED", "COUNTS_NOTVALUE_NOTREGULATED") :=
                           list(
                             COUNTS_REGULATED - COUNTS_VALUE_REGULATED,
                             COUNTS_NOTREGULATED - COUNTS_VALUE_NOTREGULATED
                           )]

  return(dt_counts)
}

perform_ora_from_counts <- function(
  counts_dt
) {
  P_VALUE <- BH_P_VALUE <- OR <-  COUNTS_VALUE_REGULATED <- COUNTS_VALUE_NOTREGULATED <-
    COUNTS_NOTVALUE_REGULATED <- COUNTS_NOTVALUE_NOTREGULATED <- ORA_SCORE  <- NULL
  counts_dt[, c("OR", "P_VALUE") :=
              vfisher_2sided(
                COUNTS_VALUE_REGULATED,
                COUNTS_VALUE_NOTREGULATED,
                COUNTS_NOTVALUE_REGULATED,
                COUNTS_NOTVALUE_NOTREGULATED
              )]
  counts_dt[, c("KULC_DISTANCE", "IMBALANCE_RATIO") := list(
    kulc(
      COUNTS_VALUE_REGULATED,
      COUNTS_VALUE_NOTREGULATED,
      COUNTS_NOTVALUE_REGULATED,
      COUNTS_NOTVALUE_NOTREGULATED
    ),
    imbalance_ratio(
      COUNTS_VALUE_REGULATED,
      COUNTS_VALUE_NOTREGULATED,
      COUNTS_NOTVALUE_REGULATED,
      COUNTS_NOTVALUE_NOTREGULATED
    )
  )]
  counts_dt[, "BH_P_VALUE" := list(stats::p.adjust(P_VALUE, method = "BH"))]
  counts_dt <- clip_infinite_OR(
    ORA_dt = counts_dt
  )
  counts_dt[, ORA_SCORE := -log10(BH_P_VALUE) * log2(OR)]
  return(counts_dt)
}

clip_infinite_OR <- function(
  ORA_dt
) {
  OR <- COUNTS_VALUE_REGULATED <- COUNTS_NOTVALUE_NOTREGULATED <-
    COUNTS_VALUE_NOTREGULATED <- COUNTS_NOTVALUE_REGULATED <- NULL
  ORA_dt[, OR := ifelse(
    is.infinite(OR),
    (COUNTS_VALUE_REGULATED + 1) * (COUNTS_NOTVALUE_NOTREGULATED + 1) / ((COUNTS_VALUE_NOTREGULATED + 1) * (COUNTS_NOTVALUE_REGULATED + 1)),
    OR
  )]
  return(ORA_dt)
}

fisher_2sided <- function(
  COUNTS_VALUE_REGULATED,
  COUNTS_VALUE_NOTREGULATED,
  COUNTS_NOTVALUE_REGULATED,
  COUNTS_NOTVALUE_NOTREGULATED
) {
  m <- matrix(
    data = c(
      COUNTS_VALUE_REGULATED, COUNTS_VALUE_NOTREGULATED,
      COUNTS_NOTVALUE_REGULATED, COUNTS_NOTVALUE_NOTREGULATED
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  test <- stats::fisher.test(
    x = m,
    alternative = "two.sided"
  )
  res <- c(test$estimate, test$p.value)
  return(res)
}

vfisher_2sided <- function(
  COUNTS_VALUE_REGULATED,
  COUNTS_VALUE_NOTREGULATED,
  COUNTS_NOTVALUE_REGULATED,
  COUNTS_NOTVALUE_NOTREGULATED
) {
  v <- mapply(
    FUN = fisher_2sided,
    COUNTS_VALUE_REGULATED,
    COUNTS_VALUE_NOTREGULATED,
    COUNTS_NOTVALUE_REGULATED,
    COUNTS_NOTVALUE_NOTREGULATED
  )
  l <- list(OR = v[1, ], P_VALUE = v[2, ])
  return(l)
}

kulc <- function(
  COUNTS_VALUE_REGULATED,
  COUNTS_VALUE_NOTREGULATED,
  COUNTS_NOTVALUE_REGULATED,
  COUNTS_NOTVALUE_NOTREGULATED
) {
  P_AB <- COUNTS_VALUE_REGULATED / (COUNTS_VALUE_REGULATED + COUNTS_NOTVALUE_REGULATED)
  P_BA <- COUNTS_VALUE_REGULATED / (COUNTS_VALUE_REGULATED + COUNTS_VALUE_NOTREGULATED)
  avg <- (P_AB + P_BA) / 2
  return(avg)
}

imbalance_ratio <- function(
  COUNTS_VALUE_REGULATED,
  COUNTS_VALUE_NOTREGULATED,
  COUNTS_NOTVALUE_REGULATED,
  COUNTS_NOTVALUE_NOTREGULATED
) {
  numerator <- abs(COUNTS_VALUE_NOTREGULATED - COUNTS_NOTVALUE_REGULATED)
  denominator <- (COUNTS_VALUE_REGULATED
                  + COUNTS_VALUE_NOTREGULATED
                  + COUNTS_NOTVALUE_REGULATED)
  return(numerator / denominator)
}
