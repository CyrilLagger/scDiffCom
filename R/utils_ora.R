run_ora <- function(
  object,
  categories,
  overwrite,
  stringent_or_default,
  stringent_logfc_threshold,
  verbose,
  class_signature,
  global
) {
  if (global == TRUE) {
    stop(
      paste0(
        "Global ORA analyis is not supported anymore.",
        " Use 'global == FALSE'"
      )
    )
  }
  if (stringent_or_default == "stringent") {
    stop(
      paste0(
        "Stringent ORA analysis is not supported anymore.",
        " Use 'stringent_or_default == 'default''"
      )
    )
  }
  regulation <- c("UP", "DOWN", "FLAT", "DIFF")
  temp_param <- object@parameters
  condition_inputs <- list(
    is_cond = temp_param$conditional_analysis,
    cond1 = temp_param$seurat_condition_id$cond1_name,
    cond2 = temp_param$seurat_condition_id$cond2_name
  )
  if (class_signature == "scDiffCom") {
    temp_by <- NULL
  }
  if (class_signature == "scDiffComCombined") {
    if (global) {
      temp_by <- NULL
    } else {
      temp_by <- "ID"
    }
  }
  if (temp_param$permutation_analysis &
      condition_inputs$is_cond) {
    logfc_threshold <- temp_param$threshold_logfc
    if (stringent_or_default == "default") {
      if (!global) {
        temp_ora <- object@ora_table
      } else {
        temp_ora <- object@ora_combined_default
      }
      categories_old <- names(temp_ora)
      if (is.null(categories_old)) {
        mes <- paste0(
          "Performing over-representation analysis on the categories: ",
          paste0(categories, collapse = ", "),
          "."
        )
        if (verbose) message(mes)
        categories_to_run <- categories
      } else {
        if (overwrite) {
          mes <- paste0(
            "Performing over-representation analysis on the categories: ",
            paste0(categories, collapse = ", "),
            ".\n",
            "Erasing all previous ORA results: ",
            paste0(categories_old, collapse = ", "),
            "."
          )
          if (verbose) message(mes)
          categories_to_run <- categories
        } else {
          categories_to_run <- setdiff(categories, categories_old)
          mes <- paste0(
            "Performing over-representation analysis on the categories: ",
            paste0(categories_to_run, collapse = ", "),
            ".\n",
            "Keeping previous ORA results: ",
            paste0(setdiff(categories_old, categories_to_run), collapse = ", "),
            "."
          )
          if (verbose) message(mes)
        }
      }
      ora_default <- sapply(
        categories_to_run,
        function(category) {
          temp_dt <- object@cci_table_detected
          res <- temp_dt[
            ,
            build_ora_dt(
              cci_table_detected = .SD,
              logfc_threshold = logfc_threshold,
              regulation = regulation,
              category = category,
              species = temp_param$LRI_species,
              global = global
            ),
            by = temp_by
          ]
        },
        USE.NAMES = TRUE,
        simplify = FALSE
      )
      if (is.null(categories_old)) {
        res_ora <- ora_default
      } else {
        if (overwrite) {
          res_ora <- ora_default
        } else {
          res_ora <- c(temp_ora, ora_default)
        }
      }
      if (global) {
        if (class_signature == "scDiffComCombined") {
          object@ora_combined_default <- res_ora
        } else {
          stop(
            paste0(
              "No ORA global analysis allowed for object",
              " of class 'scDiffComCombined'"
            )
          )
        }
      } else {
        object@ora_table <- res_ora
      }
    } else if (stringent_or_default == "stringent") {
      if (is.null(stringent_logfc_threshold)) {
        if (verbose) {
          message(
            paste0(
              "Choose a non-NULL 'stringent_logfc_threshold' to",
              " perform stringent over-representation analysis."
            )
          )
        }
      } else  {
        if(stringent_logfc_threshold > logfc_threshold) {
          if (!global) {
            temp_ora_stringent <- object@ora_stringent
          } else {
            temp_ora_stringent <- object@ora_combined_stringent
          }
          categories_old_stringent <- names(temp_ora_stringent)
          if (is.null(categories_old_stringent)) {
            mes <- paste0(
              "Performing stringent over-representation analysis ",
              "on the specified categories: ",
              paste0(categories, collapse = ", "),
              "."
            )
            if (verbose) message(mes)
            categories_stringent_to_run <- categories
          } else {
            if (overwrite) {
              mes <- paste0(
                "Performing stringent over-representation analysis ",
                "on the specified categories: ",
                paste0(categories, collapse = ", "),
                ".\n",
                "Erasing all previous ORA results: ",
                paste0(categories_old_stringent, collapse = ", "),
                "."
              )
              if (verbose) message(mes)
              categories_stringent_to_run <- categories
            } else {
              categories_stringent_to_run <- setdiff(
                categories,
                categories_old_stringent
              )
              mes <- paste0(
                "Performing stringent over-representation analysis ",
                "on the specified categories: ",
                paste0(categories_stringent_to_run, collapse = ", "),
                ".\n",
                "Keeping previous ORA results: ",
                paste0(
                  setdiff(
                    categories_old_stringent,
                    categories_stringent_to_run
                  ),
                  collapse = ", "
                ),
                "."
              )
              if (verbose) message(mes)
            }
          }
          ora_stringent <- sapply(
            categories_stringent_to_run,
            function(category) {
              temp_dt <- object@cci_table_detected
              res <- temp_dt[
                ,
                build_ora_dt(
                  cci_table_detected = .SD,
                  logfc_threshold = stringent_logfc_threshold,
                  regulation = regulation,
                  category = category,
                  species = temp_param$LRI_species,
                  global = global
                ),
                by = temp_by
              ]
            },
            USE.NAMES = TRUE,
            simplify = FALSE
          )
          if (is.null(categories_old_stringent)) {
            res_ora_stringent <- ora_stringent
          } else {
            if (overwrite) {
              res_ora_stringent <- ora_stringent
            } else {
              res_ora_stringent <- c(temp_ora_stringent, ora_stringent)
            }
          }
          if (global == TRUE) {
            if (class_signature == "scDiffComCombined") {
              object@ora_combined_stringent <- res_ora_stringent
            } else {
              stop(
                paste0(
                  "No ORA global analysis allowed for object of ",
                  "class 'scDiffComCombined'"
                )
              )
            }
          } else {
            object@ora_stringent <- res_ora_stringent
          }
        } else {
          if (verbose) {
            message(
              paste0(
                "The supposedly 'stringent' logfc is actually less ",
                "'stringent' than the default parameter. Choose a higher value."
              )
            )
          }
        }
      }
    } else {
      stop("Can't recognize parameter 'stringent_or_default")
    }
  } else {
    if (verbose) {
      message(
        "No over-representation analysis available for the selected parameters."
      )
    }
  }
  return(object)
}

build_ora_dt <- function(
  cci_table_detected,
  logfc_threshold,
  regulation,
  category,
  species,
  global
) {
  COUNTS_VALUE_REGULATED <- COUNTS_NOTVALUE_REGULATED <-
    COUNTS_NOTVALUE_NOTREGULATED <- COUNTS_VALUE_NOTREGULATED <-
    COUNTS_VALUE_REGULATED_temp <- COUNTS_VALUE_NOTREGULATED_temp<- CATEGORY <-
    ER_CELLTYPES <- ID <-  NULL
  cci_dt <- copy(
    cci_table_detected
  )
  if (global & (category == "ER_CELLTYPES")) {
    cci_dt[, ER_CELLTYPES := paste(ID, ER_CELLTYPES, sep = "_")]
  }
  ora_tables <- lapply(
    X = regulation,
    FUN = function(reg) {
      if (category %in% c("GO_TERMS", "KEGG_PWS")) {
        counts_dt <- extract_category_counts(
          cci_table_detected = cci_dt,
          category = "LRI",
          reg = reg,
          logfc_threshold = logfc_threshold
        )
        if (category == "GO_TERMS") {
          if (species == "mouse") {
            new_intersection_dt <- scDiffCom::LRI_mouse$LRI_curated_GO[
              ,
              c("LRI", "GO_ID", "GO_NAME")
            ]
          }
          if (species == "human") {
            new_intersection_dt <- scDiffCom::LRI_human$LRI_curated_GO[
              ,
              c("LRI", "GO_ID", "GO_NAME")
            ]
          }
          new_id <- "GO_ID"
          new_name <- "GO_NAME"
          new_category <- "GO_intersection"
        }
        if (category == "KEGG_PWS") {
          if (species == "mouse") {
            new_intersection_dt <- scDiffCom::LRI_mouse$LRI_curated_KEGG
          }
          if (species == "human") {
            new_intersection_dt <- scDiffCom::LRI_human$LRI_curated_KEGG
          }
          new_id <- "KEGG_ID"
          new_name <- "KEGG_NAME"
          new_category <- "KEGG_intersection"
        }
        counts_intersection_dt <- merge.data.table(
          new_intersection_dt,
          counts_dt,
          by.x = "LRI",
          by.y = "VALUE"
        )
        counts_intersection_dt[
          ,
          c(
            "COUNTS_VALUE_REGULATED_temp",
            "COUNTS_VALUE_NOTREGULATED_temp"
          ) :=
            list(
              sum(COUNTS_VALUE_REGULATED),
              sum(COUNTS_VALUE_NOTREGULATED)
            ),
          by = new_id
        ]
        counts_intersection_dt[
          ,
          c(
            "COUNTS_NOTVALUE_REGULATED_temp",
            "COUNTS_NOTVALUE_NOTREGULATED_temp"
          ) := list(
            COUNTS_NOTVALUE_REGULATED +
              COUNTS_VALUE_REGULATED -
              COUNTS_VALUE_REGULATED_temp,
            COUNTS_NOTVALUE_NOTREGULATED +
              COUNTS_VALUE_NOTREGULATED -
              COUNTS_VALUE_NOTREGULATED_temp
          )]
        counts_intersection_dt[, CATEGORY := new_category]
        counts_intersection_dt[, c(
          "COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED",
          "COUNTS_NOTVALUE_REGULATED", "COUNTS_NOTVALUE_NOTREGULATED",
          "LRI"
        ) := NULL]
        counts_intersection_dt <- unique(counts_intersection_dt)
        setnames(
          counts_intersection_dt,
          old = c(
            new_id, new_name,
            "COUNTS_VALUE_REGULATED_temp",
            "COUNTS_VALUE_NOTREGULATED_temp",
            "COUNTS_NOTVALUE_REGULATED_temp",
            "COUNTS_NOTVALUE_NOTREGULATED_temp"
          ),
          new = c(
            "VALUE_BIS", "VALUE",
            "COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED",
            "COUNTS_NOTVALUE_REGULATED", "COUNTS_NOTVALUE_NOTREGULATED"
          )
        )
        counts_dt <- counts_intersection_dt
      } else {
        counts_dt <- extract_category_counts(
          cci_table_detected = cci_dt,
          category = category,
          reg = reg,
          logfc_threshold = logfc_threshold
        )
      }
      counts_dt <- perform_ora_from_counts(
        counts_dt = counts_dt
      )
      cols_to_rename <- colnames(counts_dt)
      cols_to_rename <- cols_to_rename[
        !(cols_to_rename %in%
            c("VALUE", "VALUE_BIS", "CATEGORY"))
      ]
      setnames(
        x = counts_dt,
        old = cols_to_rename,
        new = paste0(cols_to_rename, "_", reg)
      )
      return(counts_dt)
    }
  )
  ora_full <- Reduce(merge, ora_tables)
  if (category == "GO_TERMS") {
    if (species == "mouse") {
      go_info <- scDiffCom::LRI_mouse$LRI_curated_GO
    }
    if (species == "human") {
      go_info <- scDiffCom::LRI_human$LRI_curated_GO
    }
    ora_full[
      unique(go_info[, c("GO_ID", "LEVEL", "ASPECT")]),
      on = "VALUE_BIS==GO_ID",
      c("LEVEL", "ASPECT") := mget(c("i.LEVEL", "i.ASPECT"))
    ]
  }
  return(ora_full)
}

extract_category_counts <- function(
  cci_table_detected,
  category,
  reg,
  logfc_threshold
) {
  REGULATION <- COUNTS_VALUE_REGULATED <-
    COUNTS_VALUE_NOTREGULATED <- LOGFC <- NULL
  if (reg == "UP") {
    dt_regulated <- cci_table_detected[
      REGULATION == "UP" & LOGFC >= logfc_threshold
    ]
    dt_notregulated <- cci_table_detected[
      !(REGULATION == "UP" & LOGFC >= logfc_threshold)
    ]
  } else if (reg == "DOWN") {
    dt_regulated <- cci_table_detected[
      REGULATION == "DOWN" & LOGFC <= -logfc_threshold
    ]
    dt_notregulated <- cci_table_detected[
      !(REGULATION == "DOWN" & LOGFC <= -logfc_threshold)
    ]
  } else if (reg == "DIFF") {
    dt_regulated <- cci_table_detected[
      REGULATION %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold
    ]
    dt_notregulated <- cci_table_detected[
      !(REGULATION %in% c("UP", "DOWN") & abs(LOGFC) >= logfc_threshold)
    ]
  } else if (reg == "FLAT") {
    dt_regulated <- cci_table_detected[REGULATION == "FLAT"]
    dt_notregulated <- cci_table_detected[REGULATION != "FLAT"]
  } else {
    stop("Type of ORA not supported.")
  }
  dt_counts <- rbindlist(
    l = lapply(
      X = category,
      FUN = function(categ) {
        merge.data.table(
          x = dt_regulated[
            ,
            list(
              CATEGORY = categ,
              COUNTS_VALUE_REGULATED = .N
            ),
            by = list(VALUE = get(categ))],
          y = dt_notregulated[
            ,
            list(
              CATEGORY = categ,
              COUNTS_VALUE_NOTREGULATED = .N
            ),
            by = list(VALUE = get(categ))],
          by = c("VALUE", "CATEGORY"),
          all = TRUE
        )
      }
    ),
    use.names = TRUE
  )
  setnafill(
    x = dt_counts,
    type = "const",
    fill = 0,
    cols = c("COUNTS_VALUE_REGULATED", "COUNTS_VALUE_NOTREGULATED")
  )
  COUNTS_REGULATED <- dt_regulated[, .N]
  COUNTS_NOTREGULATED <- dt_notregulated[, .N]
  dt_counts <- dt_counts[
    ,
    c(
      "COUNTS_NOTVALUE_REGULATED",
      "COUNTS_NOTVALUE_NOTREGULATED"
    ) :=
      list(
        COUNTS_REGULATED - COUNTS_VALUE_REGULATED,
        COUNTS_NOTREGULATED - COUNTS_VALUE_NOTREGULATED
      )]

  return(dt_counts)
}

perform_ora_from_counts <- function(
  counts_dt
) {
  P_VALUE <- BH_P_VALUE <- OR <-  COUNTS_VALUE_REGULATED <-
    COUNTS_VALUE_NOTREGULATED <-COUNTS_NOTVALUE_REGULATED <-
    COUNTS_NOTVALUE_NOTREGULATED <- ORA_SCORE <- i.BH_P_VALUE <- NULL
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
  #counts_dt[, "BH_P_VALUE" := list(stats::p.adjust(P_VALUE, method = "BH"))]
  temp_bh_dt <- counts_dt[COUNTS_VALUE_REGULATED > 0]
  temp_bh_dt[, "BH_P_VALUE" := list(stats::p.adjust(P_VALUE, method = "BH"))]
  counts_dt[temp_bh_dt, on = "VALUE", "BH_P_VALUE" := i.BH_P_VALUE]
  counts_dt <- clip_infinite_OR(
    ORA_dt = counts_dt
  )
  counts_dt[, ORA_SCORE := ifelse(
    OR <= 1,
    0,
    -log10(BH_P_VALUE) * log2(OR)
  )]
  return(counts_dt)
}

clip_infinite_OR <- function(
  ORA_dt
) {
  OR <- COUNTS_VALUE_REGULATED <- COUNTS_NOTVALUE_NOTREGULATED <-
    COUNTS_VALUE_NOTREGULATED <- COUNTS_NOTVALUE_REGULATED <- NULL
  ORA_dt[, OR := ifelse(
    is.infinite(OR),
    (COUNTS_VALUE_REGULATED + 1) *
      (COUNTS_NOTVALUE_NOTREGULATED + 1) /
      ((COUNTS_VALUE_NOTREGULATED + 1) *
         (COUNTS_NOTVALUE_REGULATED + 1)),
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
  P_AB <- COUNTS_VALUE_REGULATED /
    (COUNTS_VALUE_REGULATED + COUNTS_NOTVALUE_REGULATED)
  P_BA <- COUNTS_VALUE_REGULATED /
    (COUNTS_VALUE_REGULATED + COUNTS_VALUE_NOTREGULATED)
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

get_tables_ora <- function(
  object,
  categories,
  simplified,
  class_signature
) {
  if (!is.character(categories)) {
    stop("'categories' must be a character vector")
  }
  tables <- copy(object@ora_table)
  if (identical(tables, list())) {
    warning("The object does not contain ORA tables. Returning 'NULL'")
    return(NULL)
  }
  categories_in <- names(tables)
  if (identical(categories, "all")) {
    categories <- categories_in
  }
  if (!all(categories %in% categories_in)) {
    stop("All 'categories' must be present in the object.")
  }
  tables <- tables[categories]
  if (simplified) {
    if (class_signature == "scDiffComCombined") {
      first_col <- "ID"
    } else {
      first_col <- NULL
    }
    cols_to_keep <- c(
      first_col,
      "VALUE",
      "VALUE_BIS",
      "LEVEL",
      "ASPECT",
      as.vector(outer(
        c("OR_", "P_VALUE_", "BH_P_VALUE_", "ORA_SCORE_"),
        c("UP", "DOWN", "FLAT", "DIFF"),
        FUN = "paste0"
      ))
    )
    tables <- lapply(
      tables,
      function(table) {
        cols_temp <- intersect(cols_to_keep, colnames(table))
        res <- table[, cols_temp, with = FALSE]
      }
    )
  }
  tables
}

plot_ora <- function(
  ora_dt,
  category,
  regulation,
  max_terms_show,
  GO_aspect,
  OR_threshold,
  bh_p_value_threshold
) {
  VALUE <- ASPECT <- LEVEL <-  NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ggplot2\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if(regulation == "UP") {
    OR_ID <- "OR_UP"
    p_value_ID <- "BH_P_VALUE_UP"
    ORA_SCORE_ID <- "ORA_SCORE_UP"
  } else if(regulation == "DOWN") {
    OR_ID <- "OR_DOWN"
    p_value_ID <- "BH_P_VALUE_DOWN"
    ORA_SCORE_ID <- "ORA_SCORE_DOWN"
  } else if(regulation == "FLAT") {
    OR_ID <- "OR_FLAT"
    p_value_ID <- "BH_P_VALUE_FLAT"
    ORA_SCORE_ID <- "ORA_SCORE_FLAT"
  } else {
    stop("Can't find 'regulation' type")
  }
  ora_dt <- ora_dt[
    get(OR_ID) > OR_threshold &
      get(p_value_ID) <= bh_p_value_threshold
  ][order(-get(ORA_SCORE_ID))]
  if (category == "GO_TERMS") {
    ora_dt <- ora_dt[
      ASPECT == GO_aspect
    ]
    ora_dt[, VALUE := paste0(
      VALUE,
      " (L",
      LEVEL,
      ")"
    )]
  }
  if (nrow(ora_dt) == 0) {
    return("No over-represented results.")
  }
  n_row_tokeep <- min(max_terms_show, nrow(ora_dt))
  ora_dt <- ora_dt[1:n_row_tokeep]
  ora_dt$VALUE <- sapply(
    ora_dt$VALUE,
    function(i) {
      words <- strsplit(i, " ")[[1]]
      n_words <- length(words)
      if (n_words >= 5) {
        if (n_words%%2 == 0) {
          mid <- n_words/2
        } else {
          mid <- (n_words+1)/2
        }
        res <- paste0(
          paste0(words[1:mid], collapse = " "),
          "\n",
          paste0(words[(mid+1):length(words)], collapse = " ")
        )
      } else {
        res <- i
      }
      res
    }
  )
  category_label <- ifelse(
    category == "LRI",
    "Ligand-Receptor Interaction",
    ifelse(
      category == "ER_CELLTYPES",
      "Cell Type Emitter/Receiver",
      ifelse(
        category == "GO_TERMS",
        ifelse(
          GO_aspect == "biological_process",
          "GO Biological Process",
          ifelse(
            GO_aspect == "molecular_function",
            "GO Molecular Function",
            "GO Cellular Component"
          )
        ),
        ifelse(
          category == "KEGG_PWS",
          "KEGG Pathway",
          category
        )
      )
    )
  )
  ggplot2::ggplot(
    ora_dt,
    ggplot2::aes(
      get(ORA_SCORE_ID),
      stats::reorder(get("VALUE"), get(ORA_SCORE_ID))
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = -log10(get(p_value_ID)),
        color = log2(get(OR_ID))
      )
    ) +
    ggplot2::scale_color_gradient(low = "orange", high = "red") +
    ggplot2::xlab(paste0("ORA score ", regulation)) +
    ggplot2::ylab(category_label) +
    ggplot2::labs(
      size = "-log10(Adj. P-Value)",
      color = "log2(Odds Ratio)"
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::ggtitle(
      paste0(
        "Top ",
        n_row_tokeep,
        " keywords"
      )
    )
}

