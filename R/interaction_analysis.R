#' Run (differential) intercellular communication analysis
#'
#' Perform (differential) cell-cell interaction analysis based on
#'  scRNA-seq data (as a Seurat object) and an internal database of
#'  ligand-receptor interactions (LRIs).
#'
#' TODO
#'
#' @param seurat_object \emph{Seurat} object that contains normalized data
#' and relevant \code{meta.data} annotations. Gene names must be MGI (mouse)
#' or HGNC (human).
#' @param LRI_species Either "mouse" or "human". It specifies which internal
#' LRI database to use.
#' @param seurat_celltype_id Column name of the \emph{Seurat} \code{meta.data}
#'  that indicates the cell-type of each cell (e.g.: \code{"CELL_TYPE"}).
#' @param seurat_condition_id For differential analysis, a list with three
#'  named items:
#'  \enumerate{
#'    \item column_name: Column name of the \emph{Seurat} \code{meta.data}
#'     that indicates to which of two groups each cell belongs (e.g. "AGE")
#'    \item cond1_name: name of the first condition (e.g. "YOUNG")
#'    \item cond2_name: name of the second condition (e.g. "OLD")
#'  }
#'  By convention, LOGFC are positive when \code{score(cond2) > score(cond1)}.
#'  Setting this parameter to \code{NULL} allows the user to not perform the
#'  differential analysis and only return detected cell-cell interactions.
#' @param iterations Number of permutations for the statistical analysis.
#'  Default is \code{1000}. Can be set to \code{0} to quickly return partial
#'  results without performing the statistical analysis.
#' @param scdiffcom_object_name Name of the \emph{scDiffCom} object that will
#'  be returned. Default is "scDiffCom_object". Having different names for
#'  different objects is useful to combine them later on in an
#'  \emph{scDiffComCombined} object.
#' @param seurat_assay \emph{Seurat} assay from which to extract the data.
#'  Default is "RNA".
#' @param seurat_slot \emph{Seurat} slot from which to extract the data.
#'  Default is "data".
#' @param log_scale Whether to use log-normalized data or non-log-normalized
#'  data. Default is \code{FALSE}.
#' @param score_type Either "geometric_mean" (default) or "arithmetic_mean".
#'  THe metric to compute the CCI score.
#' @param threshold_min_cells Used to filter out CCIs for which either the
#'  ligand(s) or receptor(s) are expressed in less cells than this threshold.
#'  Cell-types that have less cells that the threshold are removed initially.
#'  Default is \code{5}.
#' @param threshold_pct Used to filter out the CCIs for which either the
#'  ligand(s) or receptor(s) are expressed in less than this fraction of the
#'  cells. Default is \code{0.1}.
#' @param threshold_quantile_score Quantile value used to filter interactions
#' with low scores. Default is \code{0.2}. Can be modified at will after
#' statistical analysis.
#' @param threshold_p_value_specificity Maximal BH p-value to consider
#'  an interaction to be specific. Default is \code{0.05}. Can be modified at
#'  will after statistical analysis.
#' @param threshold_p_value_de Maximal BH p-value to consider an interaction
#'  to be differentially expressed. Default is \code{0.05}. Can be modified at
#'  will after statistical analysis.
#' @param threshold_logfc Minimal LOGFC (natural log scale) to consider an
#'  interaction to be differentially expressed. Default is \code{log(1.5)}.
#'  Can be modified at will after statistical analysis.
#' @param return_distributions Whether to return the null distributions
#' computed by the permutation test. For benchmarking and testing purpose.
#' Only available when computing less than or 1000 \code{iterations} for
#' memory concerns. Default is \code{FALSE}.
#' @param seed Set a random seed. Default is \code{42}.
#' @param verbose Whether to print messages. Default is \code{TRUE}.
#'
#' @return An S4 object of class \code{scDiffCom}.
#' @export
run_interaction_analysis <- function(
  seurat_object,
  LRI_species,
  seurat_celltype_id,
  seurat_condition_id,
  iterations = 1000,
  scdiffcom_object_name = "scDiffCom_object",
  seurat_assay = "RNA",
  seurat_slot = "data",
  log_scale = FALSE,
  score_type = "geometric_mean",
  threshold_min_cells = 5,
  threshold_pct = 0.1,
  threshold_quantile_score = 0.2,
  threshold_p_value_specificity = 0.05,
  threshold_p_value_de = 0.05,
  threshold_logfc = log(1.5),
  return_distributions = FALSE,
  seed = 42,
  verbose = TRUE
) {
  if (!methods::is(seurat_object, "Seurat")) {
    stop("`seurat_object` must be an object of class Seurat")
  }
  analysis_parameters <- list(
    LRI_species = LRI_species,
    seurat_celltype_id = seurat_celltype_id,
    seurat_condition_id = seurat_condition_id,
    object_name = scdiffcom_object_name,
    seurat_assay = seurat_assay,
    seurat_slot = seurat_slot,
    log_scale = log_scale,
    score_type = score_type,
    threshold_min_cells = threshold_min_cells,
    threshold_pct = threshold_pct,
    iterations = iterations,
    threshold_quantile_score = threshold_quantile_score,
    threshold_p_value_specificity = threshold_p_value_specificity,
    threshold_p_value_de = threshold_p_value_de,
    threshold_logfc = threshold_logfc,
    return_distributions = return_distributions,
    seed = seed,
    verbose = verbose
  )
  check_parameters <- validate_parameters(
    params = analysis_parameters,
    from_inputs = TRUE
  )
  if (!is.null(check_parameters$check)) {
    stop(paste0(
      "Invalid parameters: ",
      paste0(check_parameters$check, collapse = " and ")
      ))
  } else {
    analysis_parameters <- check_parameters$params
  }
  if (LRI_species == "human") {
    LRI_table <- scDiffCom::LRI_human$LRI_curated
  }
  if (LRI_species == "mouse") {
    LRI_table <- scDiffCom::LRI_mouse$LRI_curated
  }
  set.seed(seed)
  object <- run_internal_raw_analysis(
    seurat_object = seurat_object,
    LRI_table = LRI_table,
    params = analysis_parameters
  )
  object <- FilterCCI(
    object = object,
    verbose = verbose
  )
  if (verbose) message("Successfully returning final scDiffCom object.")
  return(object)
}

run_internal_raw_analysis <- function(
  seurat_object,
  LRI_table,
  params
) {
  verbose <- params$verbose
  analysis_inputs <- extract_analysis_inputs(
    seurat_object = seurat_object,
    celltype_column_id = params$seurat_celltype_id,
    sample_column_id = NULL,
    condition_column_id = params$seurat_condition_id$column_name,
    cond1_name = params$seurat_condition_id$cond1_name,
    cond2_name = params$seurat_condition_id$cond2_name,
    assay = params$seurat_assay,
    slot = params$seurat_slot,
    log_scale = params$log_scale,
    threshold_min_cells = as.integer(params$threshold_min_cells),
    LRI_table = LRI_table,
    verbose = params$verbose
  )
  cci_template <- create_cci_template(
    analysis_inputs = analysis_inputs
  )
  if (verbose) message("Building all cell-cell interactions.")
  cci_dt_simple <- run_simple_cci_analysis(
    analysis_inputs = analysis_inputs,
    cci_template = cci_template,
    log_scale = params$log_scale,
    score_type = params$score_type,
    threshold_min_cells = as.integer(params$threshold_min_cells),
    threshold_pct = params$threshold_pct,
    compute_fast = FALSE
  )
  mes <- paste0(
    "Total number of CCIs: ",
    nrow(cci_dt_simple),
    " (",
    cci_dt_simple[, uniqueN(get("EMITTER_CELLTYPE"))],
    " * ",
    cci_dt_simple[, uniqueN(get("EMITTER_CELLTYPE"))],
    " * ",
    cci_dt_simple[, uniqueN(get("LRI"))],
    ")."
  )
  if (verbose) message(mes)
  permutation_analysis <- ifelse(
    params$iterations == 0,
    FALSE,
    TRUE
  )
  if (!permutation_analysis) {
    cci_table_raw <- cci_dt_simple
    distributions <- list()
  } else {
    res_stat_analysis <- run_stat_analysis(
      analysis_inputs = analysis_inputs,
      cci_dt_simple = cci_dt_simple,
      iterations = as.integer(params$iterations),
      return_distributions = params$return_distributions,
      score_type = params$score_type,
      verbose = params$verbose
    )
    cci_table_raw <- res_stat_analysis$cci_raw
    distributions <- res_stat_analysis$distributions
  }
  params[["permutation_analysis"]] <- permutation_analysis
  params[["conditional_analysis"]] <- analysis_inputs$condition$is_cond
  params[["max_nL"]] <- analysis_inputs$max_nL
  params[["max_nR"]] <- analysis_inputs$max_nR
  if (verbose) message("Creating an `scDiffCom` object with all `raw` CCIs.")
  object <- methods::new(
    "scDiffCom",
    parameters = params,
    cci_table_raw = cci_table_raw,
    cci_table_detected = list(),
    ora_table = list(),
    distributions = distributions
  )
}
