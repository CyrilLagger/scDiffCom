#' Run intercellular communication analysis
#'
#' Main function of the scDiffCom package. From a Seurat object and a database of ligand-receptor interactions,
#'  it returns an S4 object of class scDiffCom that contains all possible cell-cell interactions (CCIs) with
#'  their expression scores. In default mode, it also performs a series of permutation tests to determine statistically
#'  significant CCIs, as well as differentially expressed CCIs in case the dataset contains two conditions of interest.
#'
#' @param seurat_object A Seurat object that contains pre-normalized data as well as cell-type annotations.
#' @param LRdb_table A data.table with ligand-receptor interactions.
#' @param celltype_column_id The \code{meta.data} name of \code{seurat_object} that indicates the cell-type of each cell (e.g.: "CELL-TYPE")
#' @param condition_column_id The \code{meta.data} name of \code{seurat_object} that indicates the two groups of cells
#'  on which to perform the differential analysis (e.g: "AGE"). Set to \code{NULL} (default) to perform a detection analysis without
#'  testing for differential expression.
#' @param cond1_name The name of one of the two conditions contained in \code{condition_column_id} (e.g: "YOUNG"). Set to \code{NULL}
#'  (default) for no differential analysis. By convention, the returned LOGFC are positive when \code{score(cond2) > score(cond1)}.
#' @param cond2_name The name of the other condition contained in \code{condition_column_id} (e.g: "OLD"). Similar as \code{cond1_name}.
#' @param assay The assay of \code{seurat_object} from which the data are extracted. Set to "RNA" by default.
#' @param slot The slot of \code{seurat_object} from which the data are extracted. Set to "data" by default.
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
  LRdb_table,
  celltype_column_id,
  condition_column_id = NULL,
  cond1_name = NULL,
  cond2_name = NULL,
  assay = "RNA",
  slot = "data",
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
  if (!methods::is(LRdb_table, "data.table")) {
    stop("`LRdb_table must be a data.table")
  }
  analysis_parameters <- as.list(match.call())[-1]
  analysis_parameters$seurat_object <- NULL
  analysis_parameters$LRdb_table <- NULL
  analysis_parameters <- lapply(
    analysis_parameters,
    eval
  )
  check_parameters <- validate_parameters(
    params = analysis_parameters,
    from_inputs = TRUE
  )
  if (!is.null(check_parameters)) {
    stop(paste0("Invalid parameters: ", paste0(check_parameters, collapse = " and ")))
  }
  threshold_min_cells <- as.integer(threshold_min_cells)
  iterations <- as.integer(iterations)
  set.seed(seed)
  analysis_inputs <- extract_analysis_inputs(
    seurat_object = seurat_object,
    celltype_column_id = celltype_column_id,
    condition_column_id = condition_column_id,
    cond1_name = cond1_name,
    cond2_name = cond2_name,
    assay = assay,
    slot = slot,
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
    distributions = distributions,
    ora_default = list(),
    ora_stringent = list()
  )
  object <- run_filtering_and_ora(
    object = object,
    verbose = verbose
  )
  if (verbose) message("Successfully returning final scDiffCom object.")
  return(object)
}

#' Perform filtering and over-representation analysis on intercellular communication patterns
#'
#' This function is called internally by \code{run_interaction_analysis} after the permutation tests
#'  have been performed. Based on the threshold parameters, it returns detected and differentially expressed
#'  CCIs in the slot \code{cci_detected} and results of over-representation analysis in \code{cci_ora_default}.
#'  The function can be run with new threshold parameters on any scDiffCom object that already contain the slot
#'  \code{cci_raw}. This allows the user to test various filtering parameters without the need to rerun the
#'  potentially time-consuming permutation analysis. When new thresholds are defined the slot \code{parameters} is
#'  modified accordingly.
#'  Filtering and over-representation are not independent as the second depends on the first. Therefore, when
#'  running filtering with new parameters, we also need to update the ORA results.
#'
#' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' @param new_threshold_quantile_score A new threshold for the quantile score. Set to \code{NULL} by default.
#' @param new_threshold_p_value_specificity A new threshold for the specificity p-value specificity. Set to \code{NULL} by default.
#' @param new_threshold_p_value_de A new threshold for the differential expression p-value. Set to \code{NULL} by default.
#' @param new_threshold_logfc A new threshold for the differential expression logfc. Set to \code{NULL} by default.
#' @param skip_ora Should the over-representation analysis be skipped? Set to \code{FALSE} by default.
#'  Setting it to \code{FALSE} might be useful if we are interested in performing a rapid parameter scan only on the filtering
#'  parameters. Note that in such case, the slot \code{ora_default} is returned empty.
#' @param verbose Should messages be printed?
#'
#' @return An S4 object of class \code{scDiffCom}.
#' @export
run_filtering_and_ora <- function(
  object,
  new_threshold_quantile_score = NULL,
  new_threshold_p_value_specificity = NULL,
  new_threshold_p_value_de = NULL,
  new_threshold_logfc = NULL,
  skip_ora = FALSE,
  verbose = TRUE
) {
  REGULATION_SIMPLE <- NULL
  temp_param <- parameters(object)
  if (!is.null(new_threshold_quantile_score)) {
    temp_param$threshold_quantile_score <- new_threshold_quantile_score
  }
  if (!is.null(new_threshold_p_value_specificity)) {
    temp_param$threshold_p_value_specificity <- new_threshold_p_value_specificity
  }
  if (!is.null(new_threshold_p_value_de)) {
    temp_param$threshold_p_value_de <- new_threshold_p_value_de
  }
  if (!is.null(new_threshold_logfc)) {
    temp_param$threshold_logfc <- new_threshold_logfc
  }
  condition_inputs <- list(
    is_cond = temp_param$conditional_analysis,
    cond1 = temp_param$cond1_name,
    cond2 = temp_param$cond2_name
  )
  parameters(object) <- temp_param
  cci_dt <- copy(x = get_cci_raw(object))
  if (verbose) message("Filtering and cleaning `raw` CCIs.")
  cci_dt <- preprocess_cci_raw(
    cci_dt = cci_dt,
    condition_inputs = condition_inputs,
    log_scale = temp_param$log_scale,
    permutation_analysis = temp_param$permutation_analysis
  )
  if (temp_param$permutation_analysis) {
    cci_dt <- find_detected_cci(
      cci_dt = cci_dt,
      condition_inputs = condition_inputs,
      threshold_quantile_score = temp_param$threshold_quantile_score,
      threshold_p_value_specificity = temp_param$threshold_p_value_specificity,
      threshold_p_value_de = temp_param$threshold_p_value_de,
      threshold_logfc = temp_param$threshold_logfc
    )
    if (condition_inputs$is_cond) {
      if (verbose) message("Assigning precise differential regulation.")
      cci_dt <- assign_regulation(
        cci_dt = cci_dt,
        condition_inputs = condition_inputs,
        threshold_quantile_score = temp_param$threshold_quantile_score,
        threshold_p_value_specificity = temp_param$threshold_p_value_specificity
      )
      cci_dt <- cci_dt[REGULATION_SIMPLE != "NON_DETECTED"]
    }
    mes <- paste0(
      "Returning ",
      nrow(cci_dt),
      " detected CCIs."
    )
    if (verbose) message(mes)
  } else {
    mes <- paste0(
      "Returning ",
      nrow(cci_dt),
      " detected CCIs. (Note: incomplete filtering as statistical test has not been performed!)"
    )
    if (verbose) message(mes)
  }
  cci_dt <- clean_colnames(
    cci_dt = cci_dt,
    max_nL = temp_param$max_nL,
    max_nR = temp_param$max_nR,
    condition_inputs = condition_inputs,
    permutation_analysis = temp_param$permutation_analysis
  )
  if (nrow(cci_dt) == 0) {
    object <- set_cci_detected(object, list())
    object <- set_ora_default(object, list())
  } else {
    object <- set_cci_detected(object, cci_dt)
    if (skip_ora) {
      object <- set_ora_default(object, list())
    } else {
      object <- run_ora(
        object = object,
        categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS"),
        overwrite = TRUE,
        stringent_or_default = "default",
        stringent_logfc_threshold = NULL,
        verbose = verbose
      )
    }
  }
  return(object)
}

#' Perform over-representation analysis on detected intercellular communication patterns
#'
#' Over-representation analysis performed when an scDiffCom object contains results
#'  obtained from differential expression analysis.
#'
#' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#'  for each category. Set to \code{c("ER_CELLTYPES", "LR_GENES", "GO_TERMS")} by default.
#' @param overwrite Should existing categories be overwriten in case they match with new categories?
#' @param stringent_or_default Should the default or more stringent ORA be performed?
#' @param stringent_logfc_threshold A more stringent logfc threshold compared to the one stored in \code{parameters}.
#'  Set to \code{NULL} by default, namely no stringent ORA is done.
#' @param verbose Should messages be printed?
#'
#' @return An S4 object of class \code{scDiffCom} with updated slots \code{ora_defaut} (and potentially \code{ora_stringent}).
#' @export
run_ora <- function(
  object,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS"),
  overwrite = TRUE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  verbose = TRUE
) {
  regulation <- c("UP", "DOWN", "FLAT", "DIFF")
  temp_param <- parameters(object)
  condition_inputs <- list(
    is_cond = temp_param$conditional_analysis,
    cond1 = temp_param$cond1_name,
    cond2 = temp_param$cond2_name
  )
  if (temp_param$permutation_analysis &
      condition_inputs$is_cond) {
    logfc_threshold <- temp_param$threshold_logfc
    if (stringent_or_default == "default") {
      temp_ora <- get_ora_default(object)
      categories_old <- names(temp_ora)
      if (is.null(categories_old)) {
        mes <- paste0(
          "Performing over-representation analysis on the specified categories: ",
          paste0(categories, collapse = ", "),
          "."
        )
        if (verbose) message(mes)
        categories_to_run <- categories
      } else {
        if (overwrite) {
          mes <- paste0(
            "Performing over-representation analysis on the specified categories: ",
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
            "Performing over-representation analysis on the specified categories: ",
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
          build_ora_dt(
            cci_detected = get_cci_detected(object),
            logfc_threshold = logfc_threshold,
            regulation = regulation,
            category = category
          )
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
      object <- set_ora_default(object, res_ora)



    } else if (stringent_or_default == "stringent") {
      if (is.null(stringent_logfc_threshold)) {
        if (verbose) message("Choose a non-NULL `stringent_logfc_threshold` to perform stringent over-representation analysis.")
      } else  {
        if(stringent_logfc_threshold > logfc_threshold) {
          mes <- paste0(
            "Performing stringent over-representation analysis on all specified categories: ",
            paste0(categories, collapse = ", "),
            ".\n",
            "Erasing all previous stringent ORA results."
          )
          if (verbose) message(mes)
          ora_stringent <- sapply(
            categories,
            function(category) {
              build_ora_dt(
                cci_detected = get_cci_detected(object),
                logfc_threshold = stringent_logfc_threshold,
                regulation = regulation,
                category = category
              )
            },
            USE.NAMES = TRUE,
            simplify = FALSE
          )
          object <- set_ora_stringent(object, ora_stringent)
        } else {
          if (verbose) message("The supposedly `stringent` logfc is actually less `stringent` than the default parameter.
                               Choose a higher value.")
        }
      }
    } else {
      stop("Can't recognize parameter `stringent_or_default")
    }
  } else {
    if (verbose) message("No over-representation analysis analysis available for the selected parameters.")
  }
  return(object)
}


#' Plot over-represented terms of a given category
#'
#' @param object An scDiffCom object.
#' @param category The category over wich the ORA has been performed (e.g. "GO_TERMS").
#' @param regulation One from "UP", "DOWN", "DIFF" or "FLAT".
#' @param max_terms_show The maximal number of terms to show on the plot.
#' @param OR_threshold Only display terms with an odds ratio bigger than this threshold. Set to \code{1} by default.
#' @param p_value_threshold Only display terms with a BH-p-value smaller thatn this threshold. Set to \code{0.05} by default.
#' @param stringent Should we display the default ORA results or the stringent ORA results. Set to \code{FALSE} by default.
#'
#' @return A ggplot object.
#' @export
plot_ora <- function(
  object,
  category,
  regulation,
  max_terms_show,
  OR_threshold = 1,
  p_value_threshold = 0.05,
  stringent = FALSE
) {
  if (stringent) {
    ora_dt <- get_ora_stringent(object)
  } else {
    ora_dt <- get_ora_default(object)
  }
  if (identical(ora_dt, list())) {
    stop("No ORA data.table to exctract from scDiffCom object")
  }
  if (!(category %in% names(ora_dt))) {
    stop("Can't find the specified ORA category")
  }
  ora_dt <- ora_dt[[category]]
  VALUE_ID <- "VALUE"
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
    stop("Can't find `regulation` type")
  }
  ora_dt <- ora_dt[get(OR_ID) > OR_threshold & get(p_value_ID) <= p_value_threshold][order(-get(ORA_SCORE_ID))]
  if (nrow(ora_dt) == 0) {
    return("No significant ORA results for the selected parameters.")
  }
  n_row_tokeep <- min(max_terms_show, nrow(ora_dt))
  ora_dt <- ora_dt[1:n_row_tokeep]
  ggplot2::ggplot(ora_dt, aes(get(ORA_SCORE_ID), stats::reorder(get(VALUE_ID), get(ORA_SCORE_ID)))) +
    geom_point(aes(size = -log10(get(p_value_ID)), color = log2(get(OR_ID)))) +
    scale_color_gradient(low = "orange", high = "red") +
    xlab("ORA score") +
    ylab(category) +
    labs(size = "-log10(Adj. P-Value)", color = "log2(Odds Ratio)") +
    theme(text = element_text(size = 16))
}



#' Generate interactive networks with visNetwork
#'
#' @param object scDiffCom
#' @param network_type chr in "bipartite", "cells"
#'
#' @return visNetwork interactive network
#' @export
generate_interactive_network = function(
  object,
  network_type
) {
  return(build_interactive_network(object, network_type))
}

#'
#' #' Title
#' #'
#' #' @param object x
#' #' @param disperse x
#' #' @param dir x
#' #' @param from_shiny x
#' #'
#' #' @return x
#' #' @export
#' build_celltype_bipartite_graph <- function(
#'   object,
#'   disperse = FALSE,
#'   dir = NULL,
#'   from_shiny = FALSE
#' ) {
#'   ora_tables <- get_ora_tables(object)
#'   if ("ER_CELLTYPES" %in% names(ora_tables)) {
#'     ora_ct <- ora_tables[["ER_CELLTYPES"]]
#'   } else {
#'     temp_object <- run_ora(
#'       object = object,
#'       categories = "ER_CELLTYPES",
#'       overwrite = TRUE
#'     )
#'     ora_ct <- get_ora_tables(temp_object)[["ER_CELLTYPES"]]
#'   }
#'   graph_name <- parameters(object)$object_name
#'   G <- construct_graph(
#'     ora_ct = ora_ct,
#'     cci_table_filtered = get_cci_table_filtered(object),
#'     graph_name = graph_name
#'   )
#'   config <- define_graph_config()
#'   G <- setup_graph(
#'     G,
#'     config = config,
#'     use_adj_p_value = TRUE,
#'     disperse = disperse
#'   )
#'   plot_graph(
#'     G,
#'     config = config,
#'     path = NULL,
#'     show_legend = TRUE
#'   )
#'   # dir = NULL
#'   # analysis_name = NULL
#'   # ora_ct$Tissue <- tissue
#'   # if ( !is.null(dir) ) {
#'   #   if ( is.null(analysis_name) ) {
#'   #     stop('analyze_Graph: Analysis name not provided.')
#'   #   }
#'   #   subdirs = c('edge_tables', 'plots')
#'   #   create_analysis_dirs(dir, analysis_name, subdirs)
#'   #   write_as_edge_table(
#'   #     G,
#'   #     path = file.path(dir, analysis_name, 'edge_tables', paste0(tissue, '.tsv'))
#'   #   )
#'   #   plot_graph(
#'   #     G,
#'   #     path = file.path(dir, analysis_name, 'plots', paste0(tissue, '.pdf'))
#'   #   )
#'   # } else {
#'   #   if ( !is.null(analysis_name) ) {
#'   #     stop('analyze_Graph: Analysis name not null.')
#'   #   }
#'   #   plot_graph(
#'   #     G,
#'   #     path=NULL)
#'   # }
#'   # ora_has_tissue = 'Tissue' %in% names(dt_ora)
#'   # filt_has_tissue = 'Tissue' %in% names(dt_filtered)
#'   # if( !(ora_has_tissue & filt_has_tissue) ) {
#'   #   stop(paste0('analyze_Graph: Filtered and ora must have a tissue specified.',
#'   #               ' Insert a dummy tissue for compatibility with current code.'))
#'   # }
#'   # message('Solve the low statistical power issue by controlling num interactions
#'   #         in the BH adjustment.')
#'   # if( is.null(dt_filtered) ) {stop('analyze_Graph: dt_filtered is NULL.')}
#'
#'   # Process ora results and construct graph
#' }
