#' Run intercellular communication analysis
#'
#' Main function of the scDiffCom package. From a Seurat object and a database of ligand-receptor interactions,
#'  it returns an S4 object of class scDiffCom that contains all possible cell-cell interactions (CCIs) with
#'  their expression scores. In default mode, it also performs a series of permutation tests to determine statistically
#'  significant CCIs, as well as differentially expressed CCIs in case the dataset contains two conditions of interest.
#'
#' @param seurat_object A Seurat object that contains pre-normalized data as well as cell-type annotations.
#' @param LRdb_species Either "mouse" or "human". It specifies which internal ligand-receptor database to use.
#' @param seurat_celltype_id The \code{meta.data} name of \code{seurat_object} that indicates the cell-type of each cell (e.g.: "CELL-TYPE")
#' @param seurat_condition_id The \code{meta.data} name of \code{seurat_object} that indicates the two groups of cells
#'  on which to perform the differential analysis (e.g: "AGE"). Set to \code{NULL} (default) to perform a detection analysis without
#'  testing for differential expression.
#' @param cond1_name The name of one of the two conditions contained in \code{seurat_condition_id} (e.g: "YOUNG"). Set to \code{NULL}
#'  (default) for no differential analysis. By convention, the returned LOGFC are positive when \code{score(cond2) > score(cond1)}.
#' @param cond2_name The name of the other condition contained in \code{seurat_condition_id} (e.g: "OLD"). Similar as \code{cond1_name}.
#' @param seurat_assay The assay of \code{seurat_object} from which the data are extracted. Set to "RNA" by default.
#' @param seurat_slot The slot of \code{seurat_object} from which the data are extracted. Set to "data" by default.
#' @param log_scale Should the data be log-transformed before the analysis? Set to "FALSE" by default (recommended).
#' @param threshold_min_cells The minimal number of cells that a cell-type need to contain to be considered in the analysis.
#'  Set to \code{5} by default.
#' @param threshold_pct The minimal percentage of cells that need to express a gene, for the gene to be considered in the analysis.
#'  Set to \code{0.1} by default.
#' @param object_name The name of the scDiffCom object that will be returned. Set to "scDiffCom_object" by default.
#' @param permutation_analysis Should the permutation analysis be performed? Set to "TRUE" by default. When "FALSE", only
#'  raw results will be returned (such as scores).
#' @param iterations The number of iterations to perform when \code{permutation_analysis == TRUE}. Set to \code{1000} by default.
#' @param threshold_quantile_score The quantile value used to define a detection threshold to filter lowly expressed interactions.
#'  Set to \code{0.25} by default (namely the 25% lower interactions are filtered out). Can be modified a posteriori
#'  (namely without the need to perform the permutation analysis again).
#' @param threshold_p_value_specificity The maximal BH p-value to consider an interaction as detected. Set to \code{0.05} by default.
#'  Can be modified a posteriori.
#' @param threshold_p_value_de The maximal BH p-value to consider an interaction as differentially expressed. Set to \code{0.05} by default.
#'  Can be modified a posteriori.
#' @param threshold_logfc The minimal LOGFC (natural log scale) to consider an interaction as differentially expressed.
#'  Set to \code{log(1.2)} by default. Can be modified a posteriori.
#' @param return_distributions Should the distributions from the permutation analysis be returned? Set to \code{FALSE} by default.
#' @param seed A seed for replicability. Set to \code{42} by default.
#' @param verbose Should messages be printed? Set to \code{TRUE} by default.
#'
#' @return An S4 object of class \code{scDiffCom}.
#' @export
run_interaction_analysis <- function(
  seurat_object,
  LRdb_species,
  seurat_celltype_id,
  seurat_condition_id = NULL,
  cond1_name = NULL,
  cond2_name = NULL,
  seurat_assay = "RNA",
  seurat_slot = "data",
  log_scale = FALSE,
  threshold_min_cells = 5,
  threshold_pct = 0.1,
  object_name = "scDiffCom_object",
  permutation_analysis = TRUE,
  iterations = 1000,
  threshold_quantile_score = 0.25,
  threshold_p_value_specificity = 0.05,
  threshold_p_value_de = 0.05,
  threshold_logfc = log(1.2),
  return_distributions = FALSE,
  seed = 42,
  verbose = TRUE
) {
  if (!methods::is(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object")
  }
  analysis_parameters <- as.list(match.call())[-1]
  formal_temp <- formals(run_interaction_analysis)
  for (arg_temp in names(formal_temp)) {
    if (!(arg_temp %in% names(analysis_parameters)))
      analysis_parameters <- append(analysis_parameters, formal_temp[arg_temp])
  }
  analysis_parameters$seurat_object <- NULL
  analysis_parameters <- lapply(
    analysis_parameters,
    eval
  )
  check_parameters <- validate_parameters(
    params = analysis_parameters,
    from_inputs = TRUE
  )
  if (!is.null(check_parameters$check)) {
    stop(paste0("Invalid parameters: ", paste0(check_parameters$check, collapse = " and ")))
  } else {
    analysis_parameters <- check_parameters$params
  }
  if (LRdb_species == "human") {
    LRdb_table <- scDiffCom::LRdb_human$LRdb_curated
  }
  if (LRdb_species == "mouse") {
    LRdb_table <- scDiffCom::LRdb_mouse$LRdb_curated
  }
  threshold_min_cells <- as.integer(threshold_min_cells)
  iterations <- as.integer(iterations)
  set.seed(seed)
  analysis_inputs <- extract_analysis_inputs(
    seurat_object = seurat_object,
    celltype_column_id = seurat_celltype_id,
    condition_column_id = seurat_condition_id,
    cond1_name = cond1_name,
    cond2_name = cond2_name,
    assay = seurat_assay,
    slot = seurat_slot,
    log_scale = log_scale,
    threshold_min_cells = threshold_min_cells,
    LRdb_table = LRdb_table,
    verbose = verbose
  )
  cci_template <- create_cci_template(
    analysis_inputs = analysis_inputs
  )
  if (verbose) message("Building all cell-cell interactions.")
  cci_dt_simple <- run_simple_cci_analysis(
    analysis_inputs = analysis_inputs,
    cci_template = cci_template,
    log_scale = log_scale,
    threshold_min_cells = threshold_min_cells,
    threshold_pct = threshold_pct,
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
    cci_dt_simple[, uniqueN(get("LR_SORTED"))],
    ")."
  )
  if (verbose) message(mes)
  if (!permutation_analysis) {
    cci_raw <- cci_dt_simple
    distributions <- list()
  } else {
    res_stat_analysis <- run_stat_analysis(
      analysis_inputs = analysis_inputs,
      cci_dt_simple = cci_dt_simple,
      iterations = iterations,
      return_distributions = return_distributions,
      verbose = verbose
    )
    cci_raw <- res_stat_analysis$cci_raw
    distributions <- res_stat_analysis$distributions
  }
  analysis_parameters[["conditional_analysis"]] <- analysis_inputs$condition$is_cond
  analysis_parameters[["max_nL"]] <- analysis_inputs$max_nL
  analysis_parameters[["max_nR"]] <- analysis_inputs$max_nR
  if (verbose) message("Creating an `scDiffCom` object with all `raw` CCIs.")
  object <- methods::new(
    "scDiffCom",
    parameters = analysis_parameters,
    cci_raw = cci_raw,
    cci_detected = list(),
    ora_default = list(),
    ora_stringent = list(),
    distributions = distributions
  )
  object <- FilterCCI(
    object = object,
    verbose = verbose
  )
  if (verbose) message("Successfully returning final scDiffCom object.")
  return(object)
}
