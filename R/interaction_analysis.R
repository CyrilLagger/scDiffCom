#' Run (differential) intercellular communication analysis
#'
#' Perform (differential) cell type to cell type communication analysis from a
#'  \href{https://satijalab.org/seurat/index.html}{Seurat object},
#'  using an internal database of ligand-receptor interactions (LRIs). It
#'  infers biologically relevant cell-cell interactions (CCIs) and how they
#'  change between two conditions of interest. Over-representation analysis
#'  is automatically performed to determine dominant differential signals at the
#'  level of the genes, cell types, GO Terms and KEGG Pathways.
#'
#'  The primary use of this function (and of the package) is to perform
#'  differential intercellular communication analysis. However, it is also
#'  possible to only perform a detection analysis (by setting
#'  \code{seurat_condition_id} to \code{NULL}), e.g. if one wants to
#'  infer cell-cell interactions from a dataset without having conditions on the cells.
#'
#'  By convention, when performing differential analysis, LOGFC are computed as
#'  \code{log(score(cond2_name)/score(cond1_name))}. In other words,
#'  "UP"-regulated CCIs have a larger score in \code{cond2_name}.
#'
#'  Parallel computing. If possible, it is recommended to
#'  run this function in parallel in order to speed up the analysis for large
#'  dataset and/or to obtain better accuracy on the p-values by setting a higher
#'  number of \code{iterations}. This is as simple as loading the
#'  \href{https://CRAN.R-project.org/package=future}{future}
#'  package and setting an appropriate \code{plan} (see also our
#'  \href{https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html}{vignette}).
#'
#'  Data extraction. The UMI or read counts matrix is extracted from
#'  the assay \code{seurat_assay} and the slot \code{seurat_slot}. By default,
#'  it is assumed that \code{seurat_object} contains log1p-transformed
#'  normalized data in the slot "data" of its assay "RNA". If \code{log_scale}
#'  is \code{FALSE} (as recommended), the data are \code{expm1()} transformed
#'  in order to recover normalized values not in log scale.
#'
#'  Modifying filtering parameters (differential analysis only). As long as
#'  the slot \code{cci_table_raw} of
#'  the returned scDiffCom object is not erased, filtering parameters can be
#'  modified to recompute the slots \code{cci_table_detected} and
#'  \code{ora_table}, without re-performing the time consuming permutation
#'  analysis. This may be useful if one wants a fast way to analyze how the
#'  results behave in function of, say, different LOGFC thresholds. In practice,
#'  this can be done by calling the functions \code{\link{FilterCCI}} or
#'  \code{\link{RunORA}} (see also our
#'  \href{https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html}{vignette}).
#'
#'
#' @param seurat_object Seurat object that must contain normalized
#'  data and relevant \code{meta.data} columns (see below). Gene names must be
#'   MGI (mouse) or HGNC (human) approved symbols.
#' @param LRI_species Either \code{"mouse"}, \code{"human"}, \code{"rat"} or \code{"custom"}.
#' Indicates which LRI database to use and corresponds to the species of the \code{seurat_object}.
#' Use \code{"custom"} at your own risk to use your own LRI table (see \code{custom_LRI_tables}).
#' @param seurat_celltype_id Name of the \code{meta.data} column in
#'  \code{seurat_object} that contains cell-type annotations
#'  (e.g.: \code{"CELL_TYPE"}).
#' @param seurat_condition_id List that contains information regarding the
#'  two conditions on which to perform differential analysis. Must contain
#'  the following three named items:
#'  \enumerate{
#'    \item \code{column_name}: name of the \code{meta.data} column in
#'    \code{seurat_object} that indicates the condition on each cell (e.g. \code{"AGE"})
#'    \item \code{cond1_name}: name of the first condition (e.g. \code{"YOUNG"})
#'    \item \code{cond2_name}: name of the second condition (e.g. \code{"OLD"})
#'  }
#'  Can also be set to \code{NULL} to only perform a detection analysis
#'  (see Details).
#' @param iterations Number of permutations to perform the statistical
#'  analysis. The default (\code{1000}) is a good compromise for an exploratory
#'  analysis and to obtain reasonably accurate p-values in a short time.
#'  Otherwise, we recommend using \code{10000} iterations and to run the
#'  analysis in parallel (see Details). Can also be set to \code{0} for
#'  debugging and quickly returning partial results without
#'  statistical significance.
#' @param scdiffcom_object_name Name of the \code{scDiffCom} S4 object that will
#'  be returned (\code{"scDiffCom_object"} by default).
#' @param seurat_assay Assay of \code{seurat_object} from which to extract data.
#'  See Details for an explanation on how data are extracted based on the three
#'  parameters \code{seurat_assay}, \code{seurat_slot} and \code{log_scale}.
#' @param seurat_slot Slot of \code{seurat_object} from which to extract data.
#'  See Details for an explanation on how data are extracted based on the three
#'  parameters \code{seurat_assay}, \code{seurat_slot} and \code{log_scale}.
#' @param log_scale When \code{FALSE} (the default, recommended), data are
#'  treated as normalized but not log1p-transformed. See Details for an
#'  explanation on how data are extracted based on the three
#'  parameters \code{seurat_assay}, \code{seurat_slot} and \code{log_scale}.
#' @param score_type Metric used to compute cell-cell interaction (CCI) scores.
#'  Can either be \code{"geometric_mean"} (default) or \code{"arithmetic_mean"}.
#'  It is strongly recommended to use the geometric mean, especially when
#'  performing differential analysis. The arithmetic mean might be used when
#'  uniquely doing a detection analysis or if the results want to be compared
#'  with those of another package.
#' @param threshold_min_cells Minimal number of cells - of a given cell type
#'  and condition - required to express a gene for this gene to be considered
#'  expressed in the corresponding cell type. Incidentally, cell types with
#'  less cells than this threshold are removed from the analysis.
#'  Set to \code{5} by default.
#' @param threshold_pct Minimal fraction of cells - of a given cell type
#'  and condition - required to express a gene for this gene to be considered
#'  expressed in the corresponding cell type. Set to  \code{0.1} by default.
#' @param threshold_quantile_score Threshold value used in conjunction with
#'  \code{threshold_p_value_specificity} to establish if a CCI is considered
#'  "detected". The default (\code{0.2}) indicates that CCIs with a score
#'  in the 20\% lowest-scores are not considered detected. Can be modified
#'  without the need to re-perform the permutation analysis (see Details).
#' @param threshold_p_value_specificity Threshold value used in conjunction
#'  with \code{threshold_quantile_score} to establish if a CCI is considered
#'  "detected". CCIs with a (BH-adjusted) specificity p-value above the
#'  threshold (\code{0.05} by default) are not considered detected. Can be
#'  modified without the need to re-perform the permutation analysis
#'  (see Details).
#' @param threshold_p_value_de Threshold value used in conjunction
#'  with \code{threshold_logfc} to establish how CCIs are differentially
#'  expressed between \code{cond1_name} and \code{cond2_name}. CCIs with a
#'  (BH-adjusted) differential p-value above the threshold (\code{0.05} by
#'  default) are not considered to change significantly. Can be modified
#'  without the need to re-perform the permutation analysis (see Details).
#' @param threshold_logfc Threshold value used in conjunction with
#'  \code{threshold_p_value_de} to establish how CCIs are differentially
#'  expressed between \code{cond1_name} and \code{cond2_name}. CCIs with an
#'  absolute logFC below the threshold (\code{log(1.5)} by default) are
#'  considered "FLAT". Can be modified without the need to
#'  re-perform the permutation analysis (see Details).
#' @param return_distributions \code{FALSE} by default. If \code{TRUE}, the
#'  distributions obtained from the permutation test are returned alongside
#'  the other results. May be used for testing or benchmarking purposes. Can
#'  only be enabled when \code{iterations} is less than \code{1000} in order
#'  to avoid out of memory issues.
#' @param seed Set a random seed (\code{42} by default) to obtain reproducible
#' results.
#' @param verbose If \code{TRUE} (default), print progress messages.
#' @param custom_LRI_tables A list containing a LRI table and, if known, 
#' tables with annotations supplied by the user. Overwrite
#' \code{LRI_species} and the corresponding internal LRI table. Use to
#' your own risk! Must contain at least the following named item:
#'  \enumerate{
#'    \item \code{LRI}: a data.table of LRIs
#'  }
#' The data.table of LRIs must be in the same format as the internal
#' LRI_tables, namely with the columns "LRI", "LIGAND_1", "LIGAND_2",
#' "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3". Other named data.tables
#' can be supplied for over-representation analysis (ORA) purposes.
#'
#' @return An S4 object of class \code{\link{scDiffCom-class}}.
#' @export
#' @examples
#' \dontrun{
#' run_interaction_analysis(
#'   seurat_object = seurat_sample_tms_liver,
#'   LRI_species = "mouse",
#'   seurat_celltype_id = "cell_type",
#'   seurat_condition_id = list(
#'     column_name = "age_group",
#'     cond1_name = "YOUNG",
#'     cond2_name = "OLD"
#'   )
#' )
#' }
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
  verbose = TRUE,
  custom_LRI_tables = NULL
) {
  if (!methods::is(seurat_object, "Seurat")) {
    stop(
      "'seurat_object' must be an object of class Seurat"
      )
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
    verbose = verbose,
    custom_LRI_tables = custom_LRI_tables
  )
  check_parameters <- validate_parameters(
    params = analysis_parameters,
    from_inputs = TRUE
  )
  if (!is.null(check_parameters$check)) {
    stop(
      paste0(
        "Invalid parameters: ",
        paste0(check_parameters$check, collapse = " and ")
      )
    )
  } else {
    analysis_parameters <- check_parameters$params
  }
  if (analysis_parameters$LRI_species != "custom") {
    if (!is.null(analysis_parameters$custom_LRI_tables)) {
      stop(
        "LRI_species is not 'custom' but a custom LRI table was supplied."
        )
    }
    if (LRI_species == "human") {
      LRI_table <- copy(scDiffCom::LRI_human$LRI_curated)
    }
    if (LRI_species == "mouse") {
      LRI_table <- copy(scDiffCom::LRI_mouse$LRI_curated)
    }
    if (LRI_species == "rat") {
      LRI_table <- copy(scDiffCom::LRI_rat$LRI_curated)
    }
  } else {
    if (is.null(analysis_parameters$custom_LRI_tables)) {
      stop(
        "LRI_species is 'custom' but no custom LRI table was supplied."
        )
    }
    LRI_table <- copy(analysis_parameters$custom_LRI_tables$LRI)
    if (verbose) message("Using custom LRI table. Use at your own risk!")
  }
  set.seed(seed)
  object <- run_internal_raw_analysis(
    seurat_object = seurat_object,
    LRI_table = LRI_table,
    LRI_species = LRI_species,
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
  LRI_species,
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
    LRI_species = LRI_species,
    verbose = params$verbose
  )
  cci_template <- create_cci_template(
    analysis_inputs = analysis_inputs
  )
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
    "Total number of potential cell-cell interactions (CCIs): ",
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
  object <- methods::new(
    "scDiffCom",
    parameters = params,
    cci_table_raw = cci_table_raw,
    cci_table_detected = list(),
    ora_table = list(),
    distributions = distributions
  )
}
