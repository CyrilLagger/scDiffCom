run_stat_analysis <- function(
                              cci_dt_simple,
                              expr_tr,
                              metadata,
                              condition_info,
                              pp_LR,
                              iterations,
                              return_distr,
                              verbose) {
  PVAL <- PVAL_DIFF <- NULL
  LR_names <- cols_keep <- cols_keep2 <- NULL
  if (!condition_info$is_cond) {
    if (verbose) message("Performing permutation analysis without conditions.")
    sub_template_cci_dt <- cci_dt_simple[get("LR_DETECTED") == TRUE]
  } else {
    if (verbose) message("Performing permutation analysis on the two conditions.")
    sub_template_cci_dt <- cci_dt_simple[get(paste0("LR_DETECTED_", condition_info$cond1)) == TRUE | get(paste0("LR_DETECTED_", condition_info$cond2)) == TRUE]
  }
  LR_names <- c(
    paste0("LIGAND_", 1:pp_LR$max_nL),
    paste0("RECEPTOR_", 1:pp_LR$max_nR)
  )
  sub_expr_tr <- expr_tr[, colnames(expr_tr) %in%
    unique(unlist(sub_template_cci_dt[, LR_names, with = FALSE]))]
  if (verbose) {
    message(paste0(
      "Initial number of genes: ",
      ncol(expr_tr),
      ". Number of detected genes for permutation test: ",
      ncol(sub_expr_tr)
    ))
  }
  cols_keep <- c("LR_SORTED", "L_CELLTYPE", "R_CELLTYPE", LR_names)
  progressr::with_progress({
    prog <- progressr::progressor(steps = iterations)
    # cci_perm <- pbapply::pbreplicate(
    cci_perm <- future.apply::future_sapply(
      X = integer(iterations),
      FUN = function(iter) {
        if (iter %% 10 == 0) prog(sprintf("iter=%g", iter))
        run_stat_iteration(
          expr_tr = sub_expr_tr,
          metadata = metadata,
          template_cci_dt = sub_template_cci_dt[, cols_keep, with = FALSE],
          pp_LR = pp_LR,
          condition_info = condition_info
        )
      },
      simplify = "array",
      future.seed = TRUE,
      future.label = "future_replicate-%d"
    )
    # cci_perm <- future.apply::future_replicate(
    #   n = iterations,
    #   expr = run_stat_iteration(
    #     expr_tr = sub_expr_tr,
    #     metadata = metadata,
    #     template_cci_dt = sub_template_cci_dt[, cols_keep, with = FALSE],
    #     pp_LR = pp_LR,
    #     condition_info = condition_info
    #   ),
    #   simplify = "array"
    # )
  })
  if (!condition_info$is_cond) {
    distr <- cbind(cci_perm, sub_template_cci_dt[["LR_SCORE"]])
    if (verbose) message("Computing p-values.")
    pvals <- rowSums(distr[, 1:iterations] >= distr[, (iterations + 1)]) / iterations
    sub_template_cci_dt[, PVAL := pvals]
    cols_new <- c("PVAL")
    cols_keep2 <- c(cols_keep, cols_new)
    sub_template_cci_dt <- sub_template_cci_dt[, cols_keep2, with = FALSE]
    rownames(distr) <- paste(
      sub_template_cci_dt[["LR_SORTED"]],
      sub_template_cci_dt[["L_CELLTYPE"]],
      sub_template_cci_dt[["R_CELLTYPE"]],
      sep = "_"
    )
  } else {
    distr_diff <- cbind(cci_perm[, 1, ], sub_template_cci_dt[[paste0("LR_SCORE_", condition_info$cond2)]]
    - sub_template_cci_dt[[paste0("LR_SCORE_", condition_info$cond1)]])
    distr_cond1 <- cbind(cci_perm[, 2, ], sub_template_cci_dt[[paste0("LR_SCORE_", condition_info$cond1)]])
    distr_cond2 <- cbind(cci_perm[, 3, ], sub_template_cci_dt[[paste0("LR_SCORE_", condition_info$cond2)]])
    if (verbose) message("Computing p-values.")
    pvals_diff <- rowSums(abs(distr_diff[, 1:iterations]) >= abs(distr_diff[, (iterations + 1)])) / iterations
    pvals_cond1 <- rowSums(distr_cond1[, 1:iterations] >= distr_cond1[, (iterations + 1)]) / iterations
    pvals_cond2 <- rowSums(distr_cond2[, 1:iterations] >= distr_cond2[, (iterations + 1)]) / iterations
    sub_template_cci_dt[, paste0("PVAL_", c(condition_info$cond1, condition_info$cond2)) := list(pvals_cond1, pvals_cond2)]
    sub_template_cci_dt[, PVAL_DIFF := pvals_diff]
    cols_new <- c(paste0("PVAL_", condition_info$cond1), paste0("PVAL_", condition_info$cond2), "PVAL_DIFF")
    cols_keep2 <- c(cols_keep, cols_new)
    sub_template_cci_dt <- sub_template_cci_dt[, cols_keep2, with = FALSE]
    rownames(distr_diff) <- rownames(distr_cond1) <- rownames(distr_cond2) <- paste(
      sub_template_cci_dt[["LR_SORTED"]],
      sub_template_cci_dt[["L_CELLTYPE"]],
      sub_template_cci_dt[["R_CELLTYPE"]],
      sep = "_"
    )
  }
  cci_dt <- data.table::merge.data.table(
    x = cci_dt_simple,
    y = sub_template_cci_dt,
    by = cols_keep,
    all.x = TRUE,
    sort = FALSE
  )
  for (j in cols_new) {
    data.table::set(cci_dt, i = which(is.na(cci_dt[[j]])), j = j, value = 1)
  }
  if (!return_distr) {
    return(list(
      cci_table_raw = cci_dt,
      distributions = list()
    ))
  } else {
    if (verbose) message("Storing the matrix of distributions from the permutation test.")
    if (!condition_info$is_cond) {
      return(list(
        cci_table_raw = cci_dt,
        distributions = list(distr)
      ))
    } else {
      res <- list(
        cci_table_raw = cci_dt,
        distributions = list(
          distr_diff,
          distr_cond1,
          distr_cond2
        )
      )
      names(res$distributions) <- c(
        "distr_diff",
        paste0("distr_", condition_info$cond1),
        paste0("distr_", condition_info$cond2)
      )
      return(res)
    }
  }
}

run_stat_iteration <- function(
                               expr_tr,
                               metadata,
                               template_cci_dt,
                               pp_LR,
                               condition_info) {
  if (!condition_info$is_cond) {
    meta_ct <- metadata
    meta_ct$cell_type <- sample(meta_ct$cell_type)
    return(run_simple_cci_analysis(
      expr_tr = expr_tr,
      metadata = meta_ct,
      template_cci_dt = template_cci_dt,
      pp_LR = pp_LR,
      condition_info = condition_info,
      pct_threshold = NULL,
      min_cells = NULL,
      compute_fast = TRUE
    ))
  } else {
    meta_cond <- metadata
    for (x in unique(meta_cond$cell_type)) {
      meta_cond$condition[meta_cond$cell_type == x] <- sample(meta_cond$condition[meta_cond$cell_type == x])
    }
    permcond_dt <- run_simple_cci_analysis(
      expr_tr = expr_tr,
      metadata = meta_cond,
      template_cci_dt = template_cci_dt,
      pp_LR = pp_LR,
      condition_info = condition_info,
      pct_threshold = NULL,
      min_cells = NULL,
      compute_fast = TRUE
    )
    meta_ct <- metadata
    meta_ct$cell_type[meta_ct$condition == condition_info$cond1] <- sample(meta_ct$cell_type[meta_ct$condition == condition_info$cond1])
    meta_ct$cell_type[meta_ct$condition == condition_info$cond2] <- sample(meta_ct$cell_type[meta_ct$condition == condition_info$cond2])
    permct_dt <- run_simple_cci_analysis(
      expr_tr = expr_tr,
      metadata = meta_ct,
      template_cci_dt = template_cci_dt,
      pp_LR = pp_LR,
      condition_info = condition_info,
      pct_threshold = NULL,
      min_cells = NULL,
      compute_fast = TRUE
    )
    return(cbind(permcond_dt$cond2 - permcond_dt$cond1, permct_dt$cond1, permct_dt$cond2))
  }
}
