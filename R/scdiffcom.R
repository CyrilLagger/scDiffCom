#' Run scDiffCom analysis
#'
#' @param seurat_object A Seurat object
#' @param LR_object A data.table with ligand-receptor interactions
#' @param celltype_col_id The column of the Seurat meta.data that contains the cell-type of each cell
#' @param condition_col_id The column of the Seurat meta.data that contains the group of each cell (only 2 groups allowed).
#'  Use NULL to run the analysis without conditions.
#' @param cond1_name Name of the first condition. Use NULL if no condition,
#'  otherwise it should corresponds to one of the groups contained in Celltype_col_id.
#'  During differential expression, LOGFC is positive when score(cond1) > score(cond1).
#' @param cond2_name Name of the second condition. Similar as cond1_name.
#' @param assay Seurat assay from which the data are extracted.
#' @param slot Seurat slot from which the data are extracted.
#' @param log_scale Logical indicating if working in log-space or not.
#' @param min_cells Minimal number of cells for a cell-type to be considered in the analysis.
#' @param pct_threshold Minimal percentage of cells required to express a gene, for the gene to be considered in the analysis.
#' @param permutation_analysis Logical indicating if performing the statistical analysis (specificity and differential expression).
#' @param iterations Number of permutations performed when permutation_analysis == TRUE
#' @param cutoff_quantile_score Quantile used to define the threshold of detection regarding the score of an interaction.
#' @param cutoff_pval_specificity P-value threshold that indicates when an interaction is specific
#' @param cutoff_pval_de P-value threshold that indicates when an interaction is differentially expressed
#' @param cutoff_logfc LOGFC threshold to filter interactions
#' @param return_distr Logical indicating if returning the distributions of the permutation test.
#' @param seed A seed for replicability
#' @param verbose Print messages
#' @param sparse Using sparse or dense data from Seurat
#'
#' @return A data.table with the CCI interactions.
#' @export
run_scdiffcom <- function(
  seurat_object,
  LR_object,
  celltype_col_id,
  condition_col_id,
  cond1_name,
  cond2_name,
  assay = "RNA",
  slot = "data",
  log_scale = FALSE,
  min_cells = 5,
  pct_threshold = 0.1,
  permutation_analysis = TRUE,
  iterations = 1000,
  cutoff_quantile_score = 0.25,
  cutoff_pval_specificity = 0.05,
  cutoff_pval_de = 0.05,
  cutoff_logfc = log(1.1),
  return_distr = FALSE,
  seed = 42,
  verbose = TRUE,
  sparse = TRUE
) {
  set.seed(seed)
  if(sparse) {
    return_type <- "sparse"
  } else {
    return_type <- "dense"
  }
  pp_seurat <- preprocess_seurat(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    return_type = return_type,
    celltype_col_id = celltype_col_id,
    condition_col_id = condition_col_id,
    min_cells = min_cells,
    verbose = verbose
  )
  pp_LR <- preprocess_LR(
    data = pp_seurat$data,
    LR_object = LR_object,
    verbose = verbose
  )
  condition_info <- preprocess_condition(
    condition_col_id = condition_col_id,
    cond1_name = cond1_name,
    cond2_name = cond2_name,
    metadata = pp_seurat$metadata,
    verbose = verbose
  )
  template_cci_dt <- create_template_cci(
    LR_db = pp_LR$LR_db,
    cell_types = pp_seurat$cell_types
  )
  template_cci_dt <- add_cell_number(
    cci_dt = template_cci_dt,
    condition_info = condition_info,
    metadata = pp_seurat$metadata
  )
  expr_tr <- DelayedArray::t(pp_LR$data)
  cci_dt_simple <- run_simple_cci_analysis(
    expr_tr = expr_tr,
    metadata = pp_seurat$metadata,
    template_cci_dt = template_cci_dt,
    pp_LR = pp_LR,
    condition_info = condition_info,
    pct_threshold = pct_threshold,
    min_cells = min_cells,
    compute_fast = FALSE
  )
  if (!permutation_analysis) {
    analysis_result <- list(
      scdiffcom_dt_raw = cci_dt_simple,
      scdiffcom_distributions = NA
    )
  } else {
    analysis_result <- run_stat_analysis(
      cci_dt_simple = cci_dt_simple,
      expr_tr = expr_tr,
      metadata = pp_seurat$metadata,
      pp_LR = pp_LR,
      condition_info = condition_info,
      iterations = iterations,
      return_distr = return_distr,
      verbose = verbose
    )
  }
  analysis_result[["parameters"]] <- list(
    celltype_col_id = celltype_col_id,
    LR_info = list(max_nL = pp_LR$max_nL, max_nR = pp_LR$max_nR),
    condition_info = condition_info,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    min_cells = min_cells,
    pct_threshold = pct_threshold,
    permutation_analysis = permutation_analysis,
    iterations = iterations,
    cutoff_quantile_score = cutoff_quantile_score,
    cutoff_pval_specificity = cutoff_pval_specificity,
    cutoff_pval_de = cutoff_pval_de,
    cutoff_logfc = cutoff_logfc,
    return_distr = return_distr
  )
  analysis_result <- run_filtering_and_ORA(
    scdiffcom_result = analysis_result,
    verbose = verbose
  )
  return(analysis_result)
}

#' Run filtering analysis on the result of scDiffCom
#'
#' @param scdiffcom_result The list of results returned by run_scdiffcom
#' @param verbose Print messages
#' @param new_cutoff_quantile_score New quantile used to define the threshold of detection regarding the score of an interaction.
#' @param new_cutoff_pval_specificity New P-value threshold that indicates when an interaction is specific
#' @param new_cutoff_pval_de New P-value threshold that indicates when an interaction is differentially expressed
#' @param new_cutoff_logfc New LOGFC threshold to filter interactions
#'
#' @return Return a list of result in the scDiffCom format with new filtering analysis.
#' @export
run_filtering_and_ORA <- function(
  scdiffcom_result,
  verbose = TRUE,
  new_cutoff_quantile_score = NULL,
  new_cutoff_pval_specificity = NULL,
  new_cutoff_pval_de = NULL,
  new_cutoff_logfc = NULL
) {
  REGULATION_SIMPLE <- NULL
  if(verbose) message("Filtering and cleaning results.")
  if(!is.null(new_cutoff_quantile_score)) {
    scdiffcom_result$parameters$cutoff_quantile_score <- new_cutoff_quantile_score
  }
  if(!is.null(new_cutoff_pval_specificity)) {
    scdiffcom_result$parameters$cutoff_pval_specificity <- new_cutoff_pval_specificity
  }
  if(!is.null(new_cutoff_pval_de)) {
    scdiffcom_result$parameters$cutoff_pval_de <- new_cutoff_pval_de
  }
  if(!is.null(new_cutoff_logfc)) {
    scdiffcom_result$parameters$cutoff_logfc <- new_cutoff_logfc
  }
  cci_dt <- data.table::copy(
    x = scdiffcom_result$scdiffcom_dt_raw
    )
  cci_dt <- add_convenience_cols(
    cci_dt = cci_dt,
    condition_info = scdiffcom_result$parameters$condition_info,
    log_scale = scdiffcom_result$parameters$log_scale,
    permutation_analysis = scdiffcom_result$parameters$permutation_analysis,
    pre_filtering = TRUE
  )
  if(scdiffcom_result$parameters$permutation_analysis) {
    cci_dt <- find_detected_cci(
      cci_dt = cci_dt,
      condition_info = scdiffcom_result$parameters$condition_info,
      cutoff_quantile_score = scdiffcom_result$parameters$cutoff_quantile_score,
      cutoff_pval_specificity = scdiffcom_result$parameters$cutoff_pval_specificity,
      cutoff_pval_de = scdiffcom_result$parameters$cutoff_pval_de,
      cutoff_logfc = scdiffcom_result$parameters$cutoff_logfc
    )
    if(scdiffcom_result$parameters$condition_info$is_cond) {
      cci_dt <- assign_regulation(
        cci_dt = cci_dt,
        condition_info = scdiffcom_result$parameters$condition_info,
        cutoff_quantile_score = scdiffcom_result$parameters$cutoff_quantile_score,
        cutoff_pval_specificity = scdiffcom_result$parameters$cutoff_pval_specificity
      )
      cci_dt <- cci_dt[REGULATION_SIMPLE != "NON_DETECTED"]
    }
  }
  cci_dt <- clean_colnames(
    cci_dt = cci_dt,
    max_nL = scdiffcom_result$parameters$LR_info$max_nL,
    max_nR = scdiffcom_result$parameters$LR_info$max_nR,
    condition_info = scdiffcom_result$parameters$condition_info,
    permutation_analysis = scdiffcom_result$parameters$permutation_analysis
  )
  if(nrow(cci_dt) == 0) {
    if (verbose) message("No detected interactions for this dataset.")
    scdiffcom_result[["scdiffcom_dt_filtered"]] <- NA
    scdiffcom_result[["ORA"]] <- NA
  } else {
    scdiffcom_result[["scdiffcom_dt_filtered"]] <- cci_dt
    scdiffcom_result <- run_ORA(
      scdiffcom_result = scdiffcom_result,
      verbose = TRUE
    )
  }
  return(scdiffcom_result)
}

#' Run over-representation analysis on the results from scDiffCom
#'
#' @param scdiffcom_result The list of results returned by run_scdiffcom
#' @param verbose Print messages
#' @param logfc_threshold Log fold-change threshold (in absolute value) to be used to define the categories of interest.
#' Default to cutoff_logfc used for the filtering analysis.
#' @param categories The categories over which the test is done
#' @param ora_types The type of regulation. Can be one or several of c("UP", "DOWN", "DIFF", "FLAT")
#'
#' @return A data.table
#' @export
run_ORA <- function(
  scdiffcom_result,
  verbose = TRUE,
  logfc_threshold = scdiffcom_result$parameters$cutoff_logfc,
  categories = c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE", "LR_NAME", "GO"),
  ora_types = c("UP", "DOWN", "DIFF", "FLAT")
) {
  if(scdiffcom_result$parameters$permutation_analysis &
     scdiffcom_result$parameters$condition_info$is_cond) {
    if(verbose) message("Performing ORA analysis.")
    ora_dt <- build_ora_dt(
      scdiffcom_dt = scdiffcom_result$scdiffcom_dt_filtered,
      logfc_threshold = logfc_threshold,
      ora_types = ora_types,
      categories = categories
    )
    scdiffcom_result[["ORA"]] <- ora_dt
  } else {
    if (verbose) message("No ORA analysis available for the selected parameters.")
    scdiffcom_result[["ORA"]] <- NA
  }
  return(scdiffcom_result)
}

#' Build a data.table of curated ligand-receptor interactions obtained from 6 databases.
#'
#' @return A data.table with ligands, receptors and some annotations (database of origin and source of curation).
#' @export
build_LR6db <- function(
) {
  LR6db_all <- combine_LR_db(
    one2one = FALSE,
    curated = FALSE
  )
  LR6db_curated <- combine_LR_db(
    one2one = FALSE,
    curated = TRUE
  )
  LR6db_GO <- get_GO_interactions(
    LR_db = LR6db_curated
  )
  return(list(
    LR6db_all = LR6db_all,
    LR6db_curated = LR6db_curated,
    LR6db_GO = LR6db_GO
  ))
}
