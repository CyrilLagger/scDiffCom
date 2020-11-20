#' Run scDiffCom analysis
#'
#' @param seurat_object A Seurat object
#' @param LR_object A data.table with ligand-receptor interactions
#' @param celltype_column_id The column of the Seurat meta.data that contains the cell-type of each cell
#' @param condition_column_id The column of the Seurat meta.data that contains the group of each cell (only 2 groups allowed).
#'  Use NULL to run the analysis without conditions.
#' @param cond1_name Name of the first condition. Use NULL if no condition,
#'  otherwise it should corresponds to one of the groups contained in celltype_column_id.
#'  During differential expression, LOGFC is positive when score(cond1) > score(cond1).
#' @param cond2_name Name of the second condition. Similar as cond1_name.
#' @param object_name Name of the scDiffCom object that will be returned.
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
run_interaction_analysis <- function(
  seurat_object,
  LR_object,
  celltype_column_id,
  condition_column_id,
  cond1_name,
  cond2_name,
  object_name = "scDiffCom_object",
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
  cutoff_logfc = log(1.2),
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
    celltype_column_id = celltype_column_id,
    condition_column_id = condition_column_id,
    min_cells = min_cells,
    verbose = verbose
  )
  pp_LR <- preprocess_LR(
    data = pp_seurat$data,
    LR_object = LR_object,
    verbose = verbose
  )
  condition_info <- preprocess_condition(
    condition_column_id = condition_column_id,
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
    cci_table_raw <- cci_dt_simple
    distributions <- list()
  } else {
    res_stat_analysis <- run_stat_analysis(
      cci_dt_simple = cci_dt_simple,
      expr_tr = expr_tr,
      metadata = pp_seurat$metadata,
      pp_LR = pp_LR,
      condition_info = condition_info,
      iterations = iterations,
      return_distr = return_distr,
      verbose = verbose
    )
    cci_table_raw <- res_stat_analysis$cci_table_raw
    distributions <- res_stat_analysis$distributions
  }
  parameters <- list(
    object_name = object_name,
    celltype_column_id = celltype_column_id,
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
  object <- methods::new(
    "scDiffCom",
    parameters = parameters,
    cci_table_raw = cci_table_raw,
    cci_table_filtered = list(),
    distributions = distributions,
    ora_tables = list()
  )
  object <- run_filtering_and_ora(
    object = object,
    verbose = verbose
  )
  return(object)
}

#' Run filtering analysis on the result of scDiffCom
#'
#' @param object The list of results returned by run_scdiffcom
#' @param verbose Print messages
#' @param new_cutoff_quantile_score New quantile used to define the threshold of detection regarding the score of an interaction.
#' @param new_cutoff_pval_specificity New P-value threshold that indicates when an interaction is specific
#' @param new_cutoff_pval_de New P-value threshold that indicates when an interaction is differentially expressed
#' @param new_cutoff_logfc New LOGFC threshold to filter interactions
#' @param skip_ora x
#'
#' @return Return a list of result in the scDiffCom format with new filtering analysis.
#' @export
run_filtering_and_ora <- function(
  object,
  verbose = TRUE,
  new_cutoff_quantile_score = NULL,
  new_cutoff_pval_specificity = NULL,
  new_cutoff_pval_de = NULL,
  new_cutoff_logfc = NULL,
  skip_ora = FALSE
) {
  REGULATION_SIMPLE <- NULL
  if(verbose) message("Filtering and cleaning results.")
  temp_param <- parameters(object)
  if(!is.null(new_cutoff_quantile_score)) {
    temp_param$cutoff_quantile_score <- new_cutoff_quantile_score
  }
  if(!is.null(new_cutoff_pval_specificity)) {
    temp_param$cutoff_pval_specificity <- new_cutoff_pval_specificity
  }
  if(!is.null(new_cutoff_pval_de)) {
    temp_param$cutoff_pval_de <- new_cutoff_pval_de
  }
  if(!is.null(new_cutoff_logfc)) {
    temp_param$cutoff_logfc <- new_cutoff_logfc
  }
  parameters(object) <- temp_param
  cci_dt <- data.table::copy(
    x = get_cci_table_raw(object)
    )
  cci_dt <- add_convenience_cols(
    cci_dt = cci_dt,
    condition_info = temp_param$condition_info,
    log_scale = temp_param$log_scale,
    permutation_analysis = temp_param$permutation_analysis,
    pre_filtering = TRUE
  )
  if(temp_param$permutation_analysis) {
    cci_dt <- find_detected_cci(
      cci_dt = cci_dt,
      condition_info = temp_param$condition_info,
      cutoff_quantile_score = temp_param$cutoff_quantile_score,
      cutoff_pval_specificity = temp_param$cutoff_pval_specificity,
      cutoff_pval_de = temp_param$cutoff_pval_de,
      cutoff_logfc = temp_param$cutoff_logfc
    )
    if(temp_param$condition_info$is_cond) {
      cci_dt <- assign_regulation(
        cci_dt = cci_dt,
        condition_info = temp_param$condition_info,
        cutoff_quantile_score = temp_param$cutoff_quantile_score,
        cutoff_pval_specificity = temp_param$cutoff_pval_specificity
      )
      cci_dt <- cci_dt[REGULATION_SIMPLE != "NON_DETECTED"]
    }
  }
  cci_dt <- clean_colnames(
    cci_dt = cci_dt,
    max_nL = temp_param$LR_info$max_nL,
    max_nR = temp_param$LR_info$max_nR,
    condition_info = temp_param$condition_info,
    permutation_analysis = temp_param$permutation_analysis
  )
  if(nrow(cci_dt) == 0) {
    if (verbose) message("No detected interactions for this dataset.")
    object <- set_cci_table_filtered(object, list())
    object <- set_ora_tables(object, list())
  } else {
    object <- set_cci_table_filtered(object, cci_dt)
    if(skip_ora) {
      object <- set_ora_tables(object, list())
    } else {
      object <- run_ora(
      object = object,
      verbose = TRUE,
      logfc_threshold = NULL,
      categories = c("LR_CELLTYPE", "LR_NAME", "GO"),
      overwrite = TRUE
    )
    }
  }
  return(object)
}

#' Run over-representation analysis on the results from scDiffCom
#'
#' @param object The list of results returned by run_scdiffcom
#' @param verbose Print messages
#' @param logfc_threshold Log fold-change threshold (in absolute value) to be used to define the categories of interest.
#' Default to cutoff_logfc used for the filtering analysis.
#' @param categories The categories over which the test is done
#' @param overwrite Should existing categories be overwriten in case they match with new categories?
#'
#' @return A data.table
#' @export
run_ora <- function(
  object,
  verbose = TRUE,
  logfc_threshold = NULL,
  categories = c("LR_CELLTYPE", "LR_NAME", "GO"),
  overwrite = FALSE
) {
  ora_types <- c("UP", "DOWN", "FLAT", "DIFF")
  temp_param <- parameters(object)
  if(is.null(logfc_threshold)) {
    logfc_threshold <- temp_param$cutoff_logfc
  }
  if(temp_param$permutation_analysis &
     temp_param$condition_info$is_cond) {
    temp_ora <- get_ora_tables(object)
    categories_old <- names(temp_ora)
    if(is.null(categories_old)) {
      if(verbose) message("Performing ORA analysis on the specified categories.")
      categories_to_keep <- categories
    } else {
      if(overwrite) {
        if(verbose) message("Performing ORA analysis on the specified categories (and erasing any previous ORA results).")
        categories_to_keep <- categories
      } else {
        if(verbose) message("Performing ORA analysis only on the new specified categories (and keeping previous results).")
        categories_to_keep <- setdiff(categories, categories_old)
      }
    }
    ora_tables <- sapply(
      categories_to_keep,
      function(category) {
        build_ora_dt(
          cci_table_filtered = get_cci_table_filtered(object),
          logfc_threshold = logfc_threshold,
          ora_types = ora_types,
          category = category
        )
      },
      USE.NAMES = TRUE,
      simplify = FALSE
    )
    if(is.null(categories_old)) {
      res_ora <- ora_tables
    } else {
      if(overwrite) {
        res_ora <- ora_tables
      } else {
        res_ora <- c(temp_ora, ora_tables)
      }
    }
    object <- set_ora_tables(object, res_ora)
  } else {
    if (verbose) message("No ORA analysis available for the selected parameters.")
    object <- set_ora_tables(object, list())
  }
  return(object)
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

#' Title
#'
#' @param object x
#' @param category x
#' @param OR_val x
#' @param pval_val x
#' @param ORA_score_val x
#' @param max_value x
#' @param OR_cutoff x
#' @param pval_cutoff x
#'
#' @return x
#' @export
plot_ora <- function(
  object,
  category,
  #ora_type,
  OR_val,
  pval_val,
  ORA_score_val,
  max_value,
  OR_cutoff = 1,
  pval_cutoff = 0.05
) {
  ora_tables <- get_ora_tables(object)
  ora_table <- ora_tables[[category]]
  if(category == "GO") {
    Value_val <- "Value_NAME"
  } else{
    Value_val <- "Value"
  }
  # if(ora_type == "UP") {
  #   OR_val <- "OR_UP"
  #   pval_val <- "pval_adjusted_UP"
  #   ORA_score_val <- "ORA_score_UP"
  # } else if(ora_type == "DOWN") {
  #   OR_val <- "OR_DOWN"
  #   pval_val <- "pval_adjusted_DOWN"
  #   ORA_score_val <- "ORA_score_DOWN"
  #
  # } else if(ora_type == "FLAT") {
  #   OR_val <- "OR_FLAT"
  #   pval_val <- "pval_adjusted_FLAT"
  #   ORA_score_val <- "ORA_score_FLAT"
  # } else {
  #   stop("ORA_type not supported.")
  # }
  dt <- ora_table[get(OR_val) > OR_cutoff & get(pval_val) <= pval_cutoff][order(-get(ORA_score_val))]
  n_row_tokeep <- min(max_value, nrow(dt))
  dt <- dt[1:n_row_tokeep]
  ggplot2::ggplot(dt, aes(get(ORA_score_val), reorder(get(Value_val), get(ORA_score_val)))) +
    geom_point(aes(size = -log10(get(pval_val)), color = log2(get(OR_val)))) +
    scale_color_gradient(low = "orange", high = "red") +
    xlab("ORA score") +
    ylab(category) +
    labs(size = "-log10(Adj. P-Value)", color = "log2(Odds Ratio)") +
    theme(text = element_text(size = 16))
}

#' Title
#'
#' @param object x
#' @param disperse x
#' @param dir x
#' @param from_shiny x
#'
#' @return x
#' @export
build_celltype_bipartite_graph <- function(
  object,
  disperse = FALSE,
  dir = NULL,
  from_shiny = FALSE
) {
  ora_tables <- get_ora_tables(object)
  if("LR_CELLTYPE" %in% names(ora_tables)) {
    ora_ct <- ora_tables[["LR_CELLTYPE"]]
  } else {
    temp_object <- run_ora(
      object = object,
      categories = "LR_CELLTYPE",
      overwrite = TRUE
    )
    ora_ct <- get_ora_tables(temp_object)[["LR_CELLTYPE"]]
  }
  graph_name <- parameters(object)$object_name
  G <- construct_graph(
    ora_ct = ora_ct,
    cci_table_filtered = get_cci_table_filtered(object),
    graph_name = graph_name
  )
  config <- define_graph_config()
  G <- setup_graph(
    G,
    config = config,
    use_adjpval = TRUE,
    disperse = disperse
  )
  plot_graph(
    G,
    config = config,
    path=NULL,
    show_legend = TRUE
  )
  #dir = NULL
  #analysis_name = NULL
  #ora_ct$Tissue <- tissue
  # if ( !is.null(dir) ) {
  #   if ( is.null(analysis_name) ) {
  #     stop('analyze_Graph: Analysis name not provided.')
  #   }
  #   subdirs = c('edge_tables', 'plots')
  #   create_analysis_dirs(dir, analysis_name, subdirs)
  #   write_as_edge_table(
  #     G,
  #     path = file.path(dir, analysis_name, 'edge_tables', paste0(tissue, '.tsv'))
  #   )
  #   plot_graph(
  #     G,
  #     path = file.path(dir, analysis_name, 'plots', paste0(tissue, '.pdf'))
  #   )
  # } else {
  #   if ( !is.null(analysis_name) ) {
  #     stop('analyze_Graph: Analysis name not null.')
  #   }
  #   plot_graph(
  #     G,
  #     path=NULL)
  # }
  # ora_has_tissue = 'Tissue' %in% names(dt_ora)
  # filt_has_tissue = 'Tissue' %in% names(dt_filtered)
  # if( !(ora_has_tissue & filt_has_tissue) ) {
  #   stop(paste0('analyze_Graph: Filtered and ora must have a tissue specified.',
  #               ' Insert a dummy tissue for compatibility with current code.'))
  # }
  # message('Solve the low statistical power issue by controlling num interactions
  #         in the BH adjustment.')
  # if( is.null(dt_filtered) ) {stop('analyze_Graph: dt_filtered is NULL.')}

  # Process ora results and construct graph
}
