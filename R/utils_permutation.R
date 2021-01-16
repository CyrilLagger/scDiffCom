run_stat_analysis <- function(
  analysis_inputs,
  cci_dt_simple,
  iterations,
  return_distributions,
  verbose
) {
  P_VALUE <- P_VALUE_DE <- cols_keep <- cols_keep2 <- NULL
  if (!analysis_inputs$condition$is_cond) {
    sub_cci_template <- cci_dt_simple[get("CCI_DETECTED") == TRUE]
  } else {
    sub_cci_template <- cci_dt_simple[
      get(paste0("CCI_DETECTED_", analysis_inputs$condition$cond1)) == TRUE |
        get(paste0("CCI_DETECTED_", analysis_inputs$condition$cond2)) == TRUE
      ]
  }
  mes <- paste0(
    "Performing statistical analysis (",
    iterations,
    " iterations) on ",
    nrow(sub_cci_template),
    " CCIs."
  )
  if (verbose) message(mes)
  LR_COLNAMES <- c(
    paste0("LIGAND_", 1:analysis_inputs$max_nL),
    paste0("RECEPTOR_", 1:analysis_inputs$max_nR)
  )
  sub_data_tr <- analysis_inputs$data_tr[, colnames(analysis_inputs$data_tr) %in%
                           unique(unlist(sub_cci_template[, LR_COLNAMES, with = FALSE]))]
  analysis_inputs$data_tr <- sub_data_tr
  cols_keep <- c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", LR_COLNAMES)
  progressr::with_progress({
    prog <- progressr::progressor(steps = iterations)
    cci_perm <- future.apply::future_sapply(
      X = integer(iterations),
      FUN = function(iter) {
        if (iter %% 10 == 0) prog(sprintf("iter=%g", iter))
        run_stat_iteration(
          analysis_inputs = analysis_inputs,
          cci_template = sub_cci_template[, cols_keep, with = FALSE]
        )
      },
      simplify = "array",
      future.seed = TRUE,
      future.label = "future_replicate-%d"
    )
  })
  if (!analysis_inputs$condition$is_cond) {
    distr <- cbind(cci_perm, sub_cci_template[["CCI_SCORE"]])
    pvals <- rowSums(distr[, 1:iterations] >= distr[, (iterations + 1)]) / iterations
    sub_cci_template[, P_VALUE := pvals]
    cols_new <- c("P_VALUE")
    cols_keep2 <- c(cols_keep, cols_new)
    sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
    rownames(distr) <- paste(
      sub_cci_template[["LR_GENES"]],
      sub_cci_template[["EMITTER_CELLTYPE"]],
      sub_cci_template[["RECEIVER_CELLTYPE"]],
      sep = "_"
    )
  } else {
    distr_diff <- cbind(cci_perm[, 1, ], sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond2)]]
                        - sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond1)]])
    distr_cond1 <- cbind(cci_perm[, 2, ], sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond1)]])
    distr_cond2 <- cbind(cci_perm[, 3, ], sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond2)]])
    pvals_diff <- rowSums(abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations + 1)])) / iterations
    pvals_cond1 <- rowSums(distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]) / iterations
    pvals_cond2 <- rowSums(distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]) / iterations
    sub_cci_template[, paste0("P_VALUE_", c(analysis_inputs$condition$cond1, analysis_inputs$condition$cond2)) := list(pvals_cond1, pvals_cond2)]
    sub_cci_template[, P_VALUE_DE := pvals_diff]
    cols_new <- c(paste0("P_VALUE_", analysis_inputs$condition$cond1), paste0("P_VALUE_", analysis_inputs$condition$cond2), "P_VALUE_DE")
    cols_keep2 <- c(cols_keep, cols_new)
    sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
    rownames(distr_diff) <- rownames(distr_cond1) <- rownames(distr_cond2) <- paste(
      sub_cci_template[["LR_GENES"]],
      sub_cci_template[["EMITTER_CELLTYPE"]],
      sub_cci_template[["RECEIVER_CELLTYPE"]],
      sep = "_"
    )
  }
  cci_dt <- data.table::merge.data.table(
    x = cci_dt_simple,
    y = sub_cci_template,
    by = cols_keep,
    all.x = TRUE,
    sort = FALSE
  )
  for (j in cols_new) {
    data.table::set(cci_dt, i = which(is.na(cci_dt[[j]])), j = j, value = 1)
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
  cci_template
) {
  temp_md <- analysis_inputs$metadata
  if (!analysis_inputs$condition$is_cond) {
    meta_ct <- temp_md
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    analysis_inputs$metadata <- meta_ct
    return(run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    ))
  } else {
    meta_cond <- temp_md
    for (x in unique(meta_cond$cell_type)) {
      meta_cond$condition[meta_cond$cell_type == x] <- sample(meta_cond$condition[meta_cond$cell_type == x])
    }
    analysis_inputs$metadata <- meta_cond
    permcond_dt <- run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    )
    meta_ct <- temp_md
    meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond1] <- sample(meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond1])
    meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond2] <- sample(meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond2])
    analysis_inputs$metadata <- meta_ct
    permct_dt <- run_simple_cci_analysis(
      analysis_inputs = analysis_inputs,
      cci_template = cci_template,
      log_scale = FALSE,
      threshold_min_cells = NULL,
      threshold_pct = NULL,
      compute_fast = TRUE
    )
    return(cbind(permcond_dt$cond2 - permcond_dt$cond1, permct_dt$cond1, permct_dt$cond2))
  }
}
