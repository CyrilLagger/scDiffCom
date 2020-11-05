build_ora_dt <- function(
  scdiffcom_dt,
  logfc_threshold,
  ora_types,
  categories
) {
  Counts_value_significant <- Counts_notvalue_significant <-
    Counts_notvalue_notsignificant <- Counts_value_notsignificant <-
    Val_sig <- Val_notsig <- Category <-  NULL
  cci_dt <- data.table::copy(
    scdiffcom_dt
  )
  ora_tables <- lapply(
    X = ora_types,
    FUN = function(ora_type) {
      if("GO" %in% categories) {
        counts_go_dt <- extract_category_counts(
          scdiffcom_dt = cci_dt,
          categories = "LR_SORTED",
          ora_type = ora_type,
          logfc_threshold = logfc_threshold
        )
        # go_union_dt <- LR6db$LR6db_GO$LR_GO_union
        # counts_go_union_dt <- data.table::merge.data.table(
        #   go_union_dt,
        #   counts_go_dt,
        #   by.x = "LR_SORTED",
        #   by.y = "Value"
        # )
        # counts_go_union_dt[, c("Val_sig", "Val_notsig") :=
        #                      list(
        #                        sum(Counts_value_significant),
        #                        sum(Counts_value_notsignificant)),
        #                    by = "GO_union"]
        # counts_go_union_dt[, c("NotVal_sig", "NotVal_notsig") := list(
        #   Counts_notvalue_significant + Counts_value_significant - Val_sig,
        #   Counts_notvalue_notsignificant + Counts_value_notsignificant - Val_notsig
        # )]
        # counts_go_union_dt[, Category := "GO_UNION" ]
        # counts_go_union_dt[, c("Counts_value_significant", "Counts_value_notsignificant",
        #                        "Counts_notvalue_significant", "Counts_notvalue_notsignificant",
        #                        "LR_SORTED") := NULL]
        # counts_go_union_dt <- unique(counts_go_union_dt)
        # data.table::setnames(
        #   counts_go_union_dt,
        #   old = c("GO_union", "Val_sig", "Val_notsig", "NotVal_sig", "NotVal_notsig"),
        #   new = c("Value", "Counts_value_significant", "Counts_value_notsignificant",
        #           "Counts_notvalue_significant", "Counts_notvalue_notsignificant")
        # )
        go_intersection_dt <- LR6db$LR6db_GO$LR_GO_intersection
        counts_go_intersection_dt <- data.table::merge.data.table(
          go_intersection_dt,
          counts_go_dt,
          by.x = "LR_SORTED",
          by.y = "Value"
        )
        counts_go_intersection_dt[, c("Val_sig", "Val_notsig") :=
                                    list(
                                      sum(Counts_value_significant),
                                      sum(Counts_value_notsignificant)),
                                  by = "GO_intersection"]
        counts_go_intersection_dt[, c("NotVal_sig", "NotVal_notsig") := list(
          Counts_notvalue_significant + Counts_value_significant - Val_sig,
          Counts_notvalue_notsignificant + Counts_value_notsignificant - Val_notsig
        )]
        counts_go_intersection_dt[, Category := "GO_intersection" ]
        counts_go_intersection_dt[, c("Counts_value_significant", "Counts_value_notsignificant",
                                      "Counts_notvalue_significant", "Counts_notvalue_notsignificant",
                                      "LR_SORTED") := NULL]
        counts_go_intersection_dt <- unique(counts_go_intersection_dt)
        data.table::setnames(
          counts_go_intersection_dt,
          old = c("GO_intersection", "Val_sig", "Val_notsig", "NotVal_sig", "NotVal_notsig"),
          new = c("Value", "Counts_value_significant", "Counts_value_notsignificant",
                  "Counts_notvalue_significant", "Counts_notvalue_notsignificant")
        )
        # counts_go_dt <- data.table::rbindlist(
        #   list(counts_go_union_dt, counts_go_intersection_dt)
        # )
        counts_go_dt <- counts_go_intersection_dt
      } else {
        counts_go_dt <- NULL
      }
      categories <- categories[categories != "GO"]
      if(length(categories) >=1 ) {
        counts_dt <- extract_category_counts(
          scdiffcom_dt = cci_dt,
          categories = categories,
          ora_type = ora_type,
          logfc_threshold = logfc_threshold
        )
      } else {
        counts_dt <- NULL
      }
      counts_dt <- rbindlist(
        list(counts_dt, counts_go_dt)
      )
      counts_dt <- perform_ora_from_counts(
        counts_dt = counts_dt
      )
      cols_to_rename <- colnames(counts_dt)
      cols_to_rename <- cols_to_rename[!(cols_to_rename %in% c("Value", "Category"))]
      data.table::setnames(
        x = counts_dt,
        old = cols_to_rename,
        new = paste0(cols_to_rename, "_", ora_type)
      )
      return(counts_dt)
    }
  )
  ora_full <- Reduce(merge, ora_tables)
  return(ora_full)
}

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
  counts_dt <- clipp_infinite_OR(
    ORA_dt = counts_dt
  )
  return(counts_dt)
}

clipp_infinite_OR <- function(
  ORA_dt
) {
  if(sum(is.finite(ORA_dt$OR)) == 0) {
    ORA_dt$OR <- 10
  } else {
    ORA_dt$OR[is.infinite(ORA_dt$OR)] <- ceiling(max(ORA_dt$OR[is.finite(ORA_dt$OR)]))
  }
  if(sum(ORA_dt$OR > 0) == 0) {
    ORA_dt$OR <- 0.1
  } else {
    ORA_dt$OR[ORA_dt$OR == 0] <- min(ORA_dt$OR[ORA_dt$OR > 0])
  }
  return(ORA_dt)
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
