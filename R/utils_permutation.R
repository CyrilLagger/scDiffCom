run_stat_analysis <- function(
  analysis_inputs,
  cci_dt_simple,
  iterations,
  return_distributions,
  score_type,
  verbose
) {
  P_VALUE <- P_VALUE_DE <- NULL
  internal_iter <- 1000
  if (!analysis_inputs$condition$is_cond) {
    sub_cci_template <- cci_dt_simple[get("IS_CCI_EXPRESSED") == TRUE]
    cci_score_actual <- sub_cci_template[["CCI_SCORE"]]
  } else {
    sub_cci_template <- cci_dt_simple[
      get(
        paste0(
          "IS_CCI_EXPRESSED_",
          analysis_inputs$condition$cond1
        )
      ) == TRUE |
        get(
          paste0(
            "IS_CCI_EXPRESSED_",
            analysis_inputs$condition$cond2
          )
        ) == TRUE
    ]
    cci_score_cond1_actual <- sub_cci_template[[
      paste0(
        "CCI_SCORE_",
        analysis_inputs$condition$cond1
      )
    ]]
    cci_score_cond2_actual <- sub_cci_template[[
      paste0(
        "CCI_SCORE_",
        analysis_inputs$condition$cond2
      )
    ]]
    cci_score_diff_actual <-  cci_score_cond2_actual - cci_score_cond1_actual
    ligand_score_diff_actual <- lapply(
      1:analysis_inputs$max_nL,
      function(i) {
        sub_cci_template[[
          paste0("L", i, "_EXPRESSION_", analysis_inputs$condition$cond2)
        ]] -
          sub_cci_template[[
            paste0("L", i, "_EXPRESSION_", analysis_inputs$condition$cond1)
          ]]
      }
    )
    receptor_score_diff_actual <- lapply(
      1:analysis_inputs$max_nR,
      function(i) {
        sub_cci_template[[
          paste0("R", i, "_EXPRESSION_", analysis_inputs$condition$cond2)
        ]] -
          sub_cci_template[[
            paste0("R", i, "_EXPRESSION_", analysis_inputs$condition$cond1)
          ]]
      }
    )
  }
  mes <- paste0(
    "Performing permutation analysis (",
    iterations,
    " iterations by batches of ",
    internal_iter,
    ") on ",
    nrow(sub_cci_template),
    " potential CCIs."
  )
  if (verbose) message(mes)
  LR_COLNAMES <- c(
    paste0("LIGAND_", 1:analysis_inputs$max_nL),
    paste0("RECEPTOR_", 1:analysis_inputs$max_nR)
  )
  sub_data_tr <- analysis_inputs$data_tr[
    ,
    colnames(analysis_inputs$data_tr) %in%
      unique(unlist(sub_cci_template[, LR_COLNAMES, with = FALSE]))
  ]
  analysis_inputs$data_tr <- sub_data_tr
  cols_keep <- c("LRI", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", LR_COLNAMES)
  lr_pvalue_names <- c(
    sapply(
      1:analysis_inputs$max_nL,
      function(i) {
        paste0("L", i, "_P_VALUE_DE")
      }
    ),
    sapply(
      1:analysis_inputs$max_nR,
      function(i) {
        paste0("R", i, "_P_VALUE_DE")
      }
    )
  )
  if (iterations <= internal_iter) {
    n_broad_iter <- 1
  } else {
    if(return_distributions) {
      stop(
        paste0(
          "Only a maximum of 1000 iterations is allowed when",
          "'return_distributions' is TRUE"
        )
      )
    } else {
      n_broad_iter <- floor(iterations/internal_iter)
    }
  }
  sub_cci_template_iter <- sub_cci_template[, cols_keep, with = FALSE]
  if(!return_distributions) {
    array_counts <- sapply(
      X = 1:n_broad_iter,
      FUN = function(i) {
        mes <- paste0("Performing batch ", i, " of ", n_broad_iter, ".")
        if (verbose) message(mes)
        if (i < n_broad_iter) {
          n_temp_iter <- internal_iter
        } else {
          n_temp_iter <- iterations %% internal_iter
          if(n_temp_iter == 0) {
            n_temp_iter <- internal_iter
          }
        }
        cci_perm <- future.apply::future_replicate(
          n = n_temp_iter,
          expr = run_stat_iteration(
            analysis_inputs = analysis_inputs,
            cci_template = sub_cci_template_iter,
            score_type = score_type
          ),
          simplify = "array",
          future.seed = TRUE
        )
        if (!analysis_inputs$condition$is_cond) {
          temp_distr <- cbind(cci_perm, cci_score_actual)
          temp_counts <- rowSums(
            temp_distr[, 1:n_temp_iter] >= temp_distr[, (n_temp_iter + 1)]
          )
          return(temp_counts)
        } else {
          temp_distr_diff <- cbind(cci_perm[, 1, ], cci_score_diff_actual)
          temp_distr_cond1 <- cbind(cci_perm[, 2, ], cci_score_cond1_actual)
          temp_distr_cond2 <- cbind(cci_perm[, 3, ], cci_score_cond2_actual)
          temp_distr_ligand <- lapply(
            1:analysis_inputs$max_nL,
            function(i) {
              cbind(cci_perm[, 3 + i, ], ligand_score_diff_actual[[i]])
            }
          )
          temp_distr_receptor <- lapply(
            1:analysis_inputs$max_nR,
            function(i) {
              cbind(
                cci_perm[, 3 + analysis_inputs$max_nL + i, ],
                receptor_score_diff_actual[[i]]
              )
            }
          )
          temp_counts_cci_diff <- rowSums(
            abs(
              temp_distr_diff[, 1:n_temp_iter]
            ) >= abs(
              temp_distr_diff[, (n_temp_iter + 1)]
            )
          )
          temp_counts_cond1 <- rowSums(
            temp_distr_cond1[, 1:n_temp_iter] >=
              temp_distr_cond1[, (n_temp_iter + 1)])
          temp_counts_cond2 <- rowSums(
            temp_distr_cond2[, 1:n_temp_iter] >=
              temp_distr_cond2[, (n_temp_iter + 1)])
          temp_counts_ligand_diff <- lapply(
            1:analysis_inputs$max_nL,
            function(i) {
              rowSums(
                abs(
                  temp_distr_ligand[[i]][, 1:n_temp_iter]
                ) >= abs(
                  temp_distr_ligand[[i]][, (n_temp_iter + 1)]
                )
              )
            }
          )
          temp_counts_receptor_diff <- lapply(
            1:analysis_inputs$max_nR,
            function(i) {
              rowSums(
                abs(
                  temp_distr_receptor[[i]][, 1:n_temp_iter]
                ) >= abs(
                  temp_distr_receptor[[i]][, (n_temp_iter + 1)]
                )
              )
            }
          )
          temp_counts_list <- c(
            list(
              temp_counts_cci_diff,
              temp_counts_cond1,
              temp_counts_cond2
            ),
            temp_counts_ligand_diff,
            temp_counts_receptor_diff
          )
          return(do.call(cbind, temp_counts_list))
        }
      },
      simplify = "array"
    )
    if (!analysis_inputs$condition$is_cond) {
      if (n_broad_iter == 1) {
        pvals <- array_counts / iterations
      } else {
        pvals <- rowSums(array_counts) / iterations
      }
      sub_cci_template[, P_VALUE := pvals]
      cols_new <- c("P_VALUE")
      cols_keep2 <- c(cols_keep, cols_new)
      sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
    } else {
      if (n_broad_iter == 1) {
        pvals_cci_diff <- array_counts[, 1, ] / iterations
        pvals_cond1 <-  array_counts[, 2, ] / iterations
        pvals_cond2 <-  array_counts[, 3, ] / iterations
        pvals_ligand_diff <- lapply(
          1:analysis_inputs$max_nL,
          function(i) {
            array_counts[, 3 + i, ] / iterations
          }
        )
        pvals_receptor_diff <- lapply(
          1:analysis_inputs$max_nR,
          function(i) {
            array_counts[, 3 + analysis_inputs$max_nL + i, ] / iterations
          }
        )
      } else {
        pvals_cci_diff <- rowSums(array_counts[, 1, ]) / iterations
        pvals_cond1 <-  rowSums(array_counts[, 2, ]) / iterations
        pvals_cond2 <-  rowSums(array_counts[, 3, ]) / iterations
        pvals_ligand_diff <- lapply(
          1:analysis_inputs$max_nL,
          function(i) {
            rowSums(array_counts[, 3 + i, ]) / iterations
          }
        )
        pvals_receptor_diff <- lapply(
          1:analysis_inputs$max_nR,
          function(i) {
            rowSums(array_counts[, 3 + analysis_inputs$max_nL + i, ]) /
              iterations
          }
        )
      }
      sub_cci_template[
        ,
        paste0(
          "P_VALUE_",
          c(
            analysis_inputs$condition$cond1,
            analysis_inputs$condition$cond2
          )
        ) := list(pvals_cond1, pvals_cond2)
      ]
      sub_cci_template[, P_VALUE_DE := pvals_cci_diff]
      sub_cci_template[
        ,
        c(lr_pvalue_names) :=
          c(
            lapply(
              1:analysis_inputs$max_nL,
              function(i) {
                pvals_ligand_diff[[i]]
              }
            ),
            lapply(
              1:analysis_inputs$max_nR,
              function(i) {
                pvals_receptor_diff[[i]]
              }
            )
          )
      ]
      cols_new <- c(
        paste0(
          "P_VALUE_",
          analysis_inputs$condition$cond1
        ),
        paste0(
          "P_VALUE_",
          analysis_inputs$condition$cond2
        ),
        "P_VALUE_DE",
        lr_pvalue_names
      )
      cols_keep2 <- c(cols_keep, cols_new)
      sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
    }
  } else {
    cci_perm <- future.apply::future_replicate(
      n = iterations,
      expr = run_stat_iteration(
        analysis_inputs = analysis_inputs,
        cci_template = sub_cci_template_iter,
        score_type = score_type
      ),
      simplify = "array",
      future.seed = TRUE
    )
    if (!analysis_inputs$condition$is_cond) {
      distr <- cbind(cci_perm, cci_score_actual)
      pvals <- rowSums(
        distr[, 1:iterations] >= distr[, (iterations + 1)]
      ) / iterations
      sub_cci_template[, P_VALUE := pvals]
      cols_new <- c("P_VALUE")
      cols_keep2 <- c(cols_keep, cols_new)
      sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
      rownames(distr) <- paste(
        sub_cci_template[["LRI"]],
        sub_cci_template[["EMITTER_CELLTYPE"]],
        sub_cci_template[["RECEIVER_CELLTYPE"]],
        sep = "_"
      )
    } else {
      distr_diff <- cbind(cci_perm[, 1, ], cci_score_diff_actual)
      distr_cond1 <- cbind(cci_perm[, 2, ], cci_score_cond1_actual)
      distr_cond2 <- cbind(cci_perm[, 3, ], cci_score_cond2_actual)
      distr_ligand <- lapply(
        1:analysis_inputs$max_nL,
        function(i) {
          cbind(cci_perm[, 3 + i, ], ligand_score_diff_actual[[i]])
        }
      )
      distr_receptor <- lapply(
        1:analysis_inputs$max_nR,
        function(i) {
          cbind(
            cci_perm[, 3 + analysis_inputs$max_nL + i, ],
            receptor_score_diff_actual[[i]]
          )
        }
      )
      pvals_cci_diff <- rowSums(
        abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations + 1)])
      ) / iterations
      pvals_cond1 <- rowSums(
        distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]
      ) / iterations
      pvals_cond2 <- rowSums(
        distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]
      ) / iterations
      pvals_ligand_diff <- lapply(
        1:analysis_inputs$max_nL,
        function(i) {
          rowSums(
            abs(distr_ligand[[i]][, 1:iterations]) >=
              abs(distr_ligand[[i]][, (iterations + 1)])
          ) / iterations
        }
      )
      pvals_receptor_diff <- lapply(
        1:analysis_inputs$max_nR,
        function(i) {
          rowSums(
            abs(distr_receptor[[i]][, 1:iterations]) >=
              abs(distr_receptor[[i]][, (iterations + 1)])
          ) / iterations
        }
      )
      sub_cci_template[
        , paste0(
          "P_VALUE_",
          c(
            analysis_inputs$condition$cond1,
            analysis_inputs$condition$cond2
          )
        ) := list(pvals_cond1, pvals_cond2)
      ]
      sub_cci_template[, P_VALUE_DE := pvals_cci_diff]
      sub_cci_template[
        ,
        c(lr_pvalue_names) :=
          c(
            lapply(
              1:analysis_inputs$max_nL,
              function(i) {
                pvals_ligand_diff[[i]]
              }
            ),
            lapply(
              1:analysis_inputs$max_nR,
              function(i) {
                pvals_receptor_diff[[i]]
              }
            )
          )
      ]
      cols_new <- c(
        paste0(
          "P_VALUE_",
          analysis_inputs$condition$cond1
        ),
        paste0(
          "P_VALUE_",
          analysis_inputs$condition$cond2
        ),
        "P_VALUE_DE",
        lr_pvalue_names
      )
      cols_keep2 <- c(cols_keep, cols_new)
      sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
      rownames(distr_diff) <- rownames(distr_cond1) <-
        rownames(distr_cond2) <- paste(
          sub_cci_template[["LRI"]],
          sub_cci_template[["EMITTER_CELLTYPE"]],
          sub_cci_template[["RECEIVER_CELLTYPE"]],
          sep = "_"
        )
    }
  }
  cci_dt <- data.table::merge.data.table(
    x = cci_dt_simple,
    y = sub_cci_template,
    by = cols_keep,
    all.x = TRUE,
    sort = FALSE
  )
  for (j in cols_new) {
    data.table::set(
      cci_dt,
      i = which(is.na(cci_dt[[j]])),
      j = j,
      value = 1
    )
  }
  if (analysis_inputs$condition$is_cond) {
    cci_dt[
      ,
      c(lr_pvalue_names) :=
        c(
          lapply(
            1:analysis_inputs$max_nL,
            function(i) {
              ifelse(
                is.na(get(paste0("LIGAND_", i))),
                NA,
                get(paste0("L", i, "_P_VALUE_DE"))
              )
            }
          ),
          lapply(
            1:analysis_inputs$max_nR,
            function(i) {
              ifelse(
                is.na(get(paste0("RECEPTOR_", i))),
                NA,
                get(paste0("R", i, "_P_VALUE_DE"))
              )
            }
          )
        )
    ]
  }
  if (!return_distributions) {
    return(list(
      cci_raw = cci_dt,
      distributions = list()
    ))
  } else {
    if (verbose) message("Returning distributions from the permutation test.")
    if (!analysis_inputs$condition$is_cond) {
      return(list(
        cci_raw = cci_dt,
        distributions = list("DISTRIBUTIONS" = distr)
      ))
    } else {
      res <- list(
        cci_raw = cci_dt,
        distributions = list(
          distr_diff,
          distr_cond1,
          distr_cond2
        )
      )
      names(res$distributions) <- c(
        "DISTRIBUTIONS_DE",
        paste0("DISTRIBUTIONS_", analysis_inputs$condition$cond1),
        paste0("DISTRIBUTIONS_", analysis_inputs$condition$cond2)
      )
      return(res)
    }
  }
}

run_stat_iteration <- function(
  analysis_inputs,
  cci_template,
  score_type
) {
  sample_temp <- sample_id <- condition <- NULL
  if (!analysis_inputs$condition$is_cond) {
    meta_ct <- copy(analysis_inputs$metadata)
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    analysis_inputs$metadata <- meta_ct
    permct <- run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      score_type = score_type,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    )
    return(permct)
  } else {
    if(analysis_inputs$condition$is_samp) {
      meta_cond <- copy(analysis_inputs$metadata)
      meta_cond <- meta_cond[, {
        temp <- unique(.SD[, c("sample_id", "condition")])
        temp[, sample_temp := sample(sample_id, replace = TRUE)]
        rbindlist(apply(temp, MARGIN = 1, function(i) {
          temp2 <- .SD[sample_id == i[["sample_temp"]]]
          temp2[, condition := i[["condition"]]]
          return(temp2)
        }))
      }, by = c("cell_type")]
      analysis_inputs$data_tr <- analysis_inputs$data_tr[meta_cond$cell_id, ]
      analysis_inputs$metadata <- meta_cond
    } else {
      meta_cond <- copy(analysis_inputs$metadata)
      for (x in unique(meta_cond$cell_type)) {
        meta_cond$condition[meta_cond$cell_type == x] <-
          sample(meta_cond$condition[meta_cond$cell_type == x])
      }
      analysis_inputs$metadata <- meta_cond
    }
    permcond_dt <- run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      score_type = score_type,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    )
    meta_ct <- copy(analysis_inputs$metadata)
    meta_ct$cell_type[
      meta_ct$condition == analysis_inputs$condition$cond1
    ] <-
      sample(
        meta_ct$cell_type[
          meta_ct$condition == analysis_inputs$condition$cond1
        ]
      )
    meta_ct$cell_type[
      meta_ct$condition == analysis_inputs$condition$cond2
    ] <-
      sample(
        meta_ct$cell_type[
          meta_ct$condition == analysis_inputs$condition$cond2
        ]
      )
    analysis_inputs$metadata <- meta_ct
    permct_dt <- run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      score_type = score_type,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    )
    permcond_dt_diff <- permcond_dt$cond2 - permcond_dt$cond1
    res_list <- c(
      list(
        permcond_dt_diff,
        permct_dt$cond1,
        permct_dt$cond2
      ),
      lapply(
        1:analysis_inputs$max_nL,
        function(i) {
          permcond_dt[[paste0("L", i, "_DIFF_EXPR")]]
        }
      ),
      lapply(
        1:analysis_inputs$max_nR,
        function(i) {
          permcond_dt[[paste0("R", i, "_DIFF_EXPR")]]
        }
      )
    )
    return(
      do.call(cbind, res_list)
    )
  }
}
