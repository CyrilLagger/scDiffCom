extract_category_counts <- function(
  scdiffcom_dt,
  categories,
  ora_type,
  logfc_threshold
  ) {
  REGULATION_SIMPLE <- Counts_value_significant <- Counts_value_notsignificant <- LOGFC <- NULL
  if(ora_type == "UP") {
    dt_significant <- scdiffcom_dt[REGULATION_SIMPLE == "UP" & LOGFC >= logfc_threshold]
    dt_notsignificant <- scdiffcom_dt[!(REGULATION_SIMPLE == "UP" & LOGFC >= logfc_threshold)]
  } else if(ora_type == "DOWN") {
    dt_significant <- scdiffcom_dt[REGULATION_SIMPLE == "DOWN" & LOGFC <= -logfc_threshold]
    dt_notsignificant <- scdiffcom_dt[!(REGULATION_SIMPLE == "DOWN" & LOGFC <= -logfc_threshold)]
  } else if(ora_type == "DIFF") {
    dt_significant <- scdiffcom_dt[REGULATION_SIMPLE %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold ]
    dt_notsignificant <- scdiffcom_dt[!(REGULATION_SIMPLE %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold)]
  } else if(ora_type == "FLAT") {
    dt_significant <- scdiffcom_dt[REGULATION_SIMPLE == "FLAT"]
    dt_notsignificant <- scdiffcom_dt[REGULATION_SIMPLE != "FLAT"]
  } else {
    stop("Type of ORA not supported.")
  }
  dt_counts <- data.table::rbindlist(
    l = lapply(
      X = categories,
      FUN = function(category) {
        data.table::merge.data.table(
          x = dt_significant[, list(Category = category, Counts_value_significant = .N), by = list(Value = get(category))],
          y = dt_notsignificant[, list(Category = category, Counts_value_notsignificant = .N), by = list(Value = get(category))],
          by = c("Value", "Category"),
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
    cols = c("Counts_value_significant", "Counts_value_notsignificant")
  )
  Counts_significant <- dt_significant[, .N]
  Counts_notsignificant <- dt_notsignificant[, .N]
  dt_counts <- dt_counts[, c("Counts_notvalue_significant", "Counts_notvalue_notsignificant") :=
                           list(
                             Counts_significant - Counts_value_significant,
                             Counts_notsignificant - Counts_value_notsignificant
                           )
                         ]

  return(dt_counts)
}

perform_ora_from_counts <- function(
  counts_dt
) {
  pval <- Counts_value_significant <- Counts_value_notsignificant <-
    Counts_notvalue_significant <- Counts_notvalue_notsignificant <- NULL
  counts_dt[, c("OR", "pval") :=
              vfisher_2sided(
                Counts_value_significant,
                Counts_value_notsignificant,
                Counts_notvalue_significant,
                Counts_notvalue_notsignificant
              )
            ]
  counts_dt[, c("Kulc_distance", "Imbalance_ratio") := list(
    kulc(
      Counts_value_significant,
      Counts_value_notsignificant,
      Counts_notvalue_significant,
      Counts_notvalue_notsignificant
    ),
    imbalance_ratio(
      Counts_value_significant,
      Counts_value_notsignificant,
      Counts_notvalue_significant,
      Counts_notvalue_notsignificant
    )
  )]
  counts_dt[, "pval_adjusted" := list(stats::p.adjust(pval, method="BH"))]
  return(counts_dt)
}

fisher_2sided <- function(
  Counts_value_significant,
  Counts_value_notsignificant,
  Counts_notvalue_significant,
  Counts_notvalue_notsignificant
) {
  m <- matrix(
    data = c(Counts_value_significant, Counts_value_notsignificant,
             Counts_notvalue_significant, Counts_notvalue_notsignificant),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  test <- stats::fisher.test(
    x = m,
    alternative = "two.sided"
    )
  res = c(test$estimate, test$p.value)
  return(res)
}

vfisher_2sided <- function(
  Counts_value_significant,
  Counts_value_notsignificant,
  Counts_notvalue_significant,
  Counts_notvalue_notsignificant
) {
  v <- mapply(
    FUN = fisher_2sided,
    Counts_value_significant,
    Counts_value_notsignificant,
    Counts_notvalue_significant,
    Counts_notvalue_notsignificant
  )
  l = list(OR = v[1, ], pval = v[2, ])
  return(l)
}

kulc <- function(
  Counts_value_significant,
  Counts_value_notsignificant,
  Counts_notvalue_significant,
  Counts_notvalue_notsignificant
) {
  P_AB <- Counts_value_significant / (Counts_value_significant + Counts_notvalue_significant)
  P_BA <- Counts_value_significant / (Counts_value_significant + Counts_value_notsignificant)
  avg <- (P_AB + P_BA) / 2
  return(avg)
}


imbalance_ratio <- function(
  Counts_value_significant,
  Counts_value_notsignificant,
  Counts_notvalue_significant,
  Counts_notvalue_notsignificant
) {
  numerator <- abs(Counts_value_notsignificant - Counts_notvalue_significant)
  denominator <- (Counts_value_significant
                 + Counts_value_notsignificant
                 + Counts_notvalue_significant)
  return(numerator / denominator)
}
