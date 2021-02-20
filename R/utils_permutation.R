run_stat_analysis <- function(
  analysis_inputs,
  cci_dt_simple,
  iterations,
  return_distributions,
  score_type,
  verbose
) {
  P_VALUE <- P_VALUE_DE <- NULL
  if (!analysis_inputs$condition$is_cond) {
    sub_cci_template <- cci_dt_simple[get("IS_CCI_EXPRESSED") == TRUE]
    cci_score_actual <- sub_cci_template[["CCI_SCORE"]]
  } else {
    sub_cci_template <- cci_dt_simple[
      get(paste0("IS_CCI_EXPRESSED_", analysis_inputs$condition$cond1)) == TRUE |
        get(paste0("IS_CCI_EXPRESSED_", analysis_inputs$condition$cond2)) == TRUE
      ]
    cci_score_cond1_actual <- sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond1)]]
    cci_score_cond2_actual <- sub_cci_template[[paste0("CCI_SCORE_", analysis_inputs$condition$cond2)]]
    cci_score_diff_actual <-  cci_score_cond2_actual - cci_score_cond1_actual
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
  if (!return_distributions) {
    internal_iter <- 1000
    if (iterations <= internal_iter) {
      n_broad_iter <- 1
    } else {
      n_broad_iter <- floor(iterations/internal_iter)
    }
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
        sub_cci_template_iter <- sub_cci_template[, cols_keep, with = FALSE]

        # progressr::with_progress({
        #   prog <- progressr::progressor(steps = n_temp_iter)
        #   cci_perm <- future.apply::future_sapply(
        #     X = integer(n_temp_iter),
        #     FUN = function(iter) {
        #       if (iter %% 20 == 0) prog(sprintf("iter=%g", iter))
        #       run_stat_iteration(
        #         analysis_inputs = analysis_inputs,
        #         cci_template = sub_cci_template_iter,
        #         score_type = score_type
        #       )
        #     },
        #     simplify = "array",
        #     future.seed = TRUE,
        #     future.label = "future_replicate-%d"
        #   )
        # },
        # enable = verbose
        # )
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
          temp_counts <- rowSums(temp_distr[, 1:n_temp_iter] >= temp_distr[, (n_temp_iter + 1)])
          return(temp_counts)
        } else {
          temp_distr_diff <- cbind(cci_perm[, 1, ], cci_score_diff_actual)
          temp_distr_cond1 <- cbind(cci_perm[, 2, ], cci_score_cond1_actual)
          temp_distr_cond2 <- cbind(cci_perm[, 3, ], cci_score_cond2_actual)
          temp_counts_diff <- rowSums(abs(temp_distr_diff[, 1:n_temp_iter]) >= abs(temp_distr_diff[, (n_temp_iter + 1)]))
          temp_counts_cond1 <- rowSums(temp_distr_cond1[, 1:n_temp_iter] >= temp_distr_cond1[, (n_temp_iter + 1)])
          temp_counts_cond2 <- rowSums(temp_distr_cond2[, 1:n_temp_iter] >= temp_distr_cond2[, (n_temp_iter + 1)])
          return(cbind(temp_counts_diff, temp_counts_cond1, temp_counts_cond2))
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
        pvals_diff <- array_counts[, 1, ] / iterations
        pvals_cond1 <-  array_counts[, 2, ] / iterations
        pvals_cond2 <-  array_counts[, 3, ] / iterations
      } else {
        pvals_diff <- rowSums(array_counts[, 1, ]) / iterations
        pvals_cond1 <-  rowSums(array_counts[, 2, ]) / iterations
        pvals_cond2 <-  rowSums(array_counts[, 3, ]) / iterations
      }
      sub_cci_template[, paste0("P_VALUE_", c(analysis_inputs$condition$cond1, analysis_inputs$condition$cond2)) := list(pvals_cond1, pvals_cond2)]
      sub_cci_template[, P_VALUE_DE := pvals_diff]
      cols_new <- c(paste0("P_VALUE_", analysis_inputs$condition$cond1), paste0("P_VALUE_", analysis_inputs$condition$cond2), "P_VALUE_DE")
      cols_keep2 <- c(cols_keep, cols_new)
      sub_cci_template <- sub_cci_template[, cols_keep2, with = FALSE]
    }
  } else {
    sub_cci_template_iter <- sub_cci_template[, cols_keep, with = FALSE]

    # progressr::with_progress({
    #   prog <- progressr::progressor(steps = iterations)
    #   cci_perm <- future.apply::future_sapply(
    #     X = integer(iterations),
    #     FUN = function(iter) {
    #       if (iter %% 20 == 0) prog(sprintf("iter=%g", iter))
    #       run_stat_iteration(
    #         analysis_inputs = analysis_inputs,
    #         cci_template = sub_cci_template_iter,
    #         score_type = score_type
    #       )
    #     },
    #     simplify = "array",
    #     future.seed = TRUE,
    #     future.label = "future_replicate-%d"
    #   )
    # },
    # enable = verbose
    # )

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
      distr <- cbind(cci_perm, cci_score_actual)
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
      distr_diff <- cbind(cci_perm[, 1, ], cci_score_diff_actual)
      distr_cond1 <- cbind(cci_perm[, 2, ], cci_score_cond1_actual)
      distr_cond2 <- cbind(cci_perm[, 3, ], cci_score_cond2_actual)
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
        meta_cond$condition[meta_cond$cell_type == x] <- sample(meta_cond$condition[meta_cond$cell_type == x])
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
    meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond1] <-
      sample(meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond1])
    meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond2] <-
      sample(meta_ct$cell_type[meta_ct$condition == analysis_inputs$condition$cond2])
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
    return(cbind(permcond_dt_diff, permct_dt$cond1, permct_dt$cond2))
  }
}
