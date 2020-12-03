#' @include generics.R
#' @import data.table
#' @import ggplot2
#' @importFrom methods new setClass setClassUnion setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####
# Class definitions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(
  name = "list_or_data.table",
  members = c("list", "data.table")
)

#' The scDiffComBase Class
#'
#' The scDiffComBase class is a virtual class that provides a template
#' for the scDiffCom and scDiffComCombined classes.
#'
#' @slot parameters A list of the parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table of all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \cdot m} rows,
#'  where n is the number of cell-types and m the number of ligand-receptor pairs.
#' @slot cci_detected A data.table of the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_default A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_detected} table.
#' @slot ora_stringent Same as \code{ora_default} but with a more stringent logfc threshold to specifically focus on strong signals.
#'
#' @name scDiffComBase-class
#' @rdname scDiffComBase-class
#'
setClass(
  Class = "scDiffComBase",
  contains = "VIRTUAL",
  slots = c(
    parameters = "list",
    cci_raw = "list_or_data.table",
    cci_detected = "list_or_data.table",
    ora_default = "list_or_data.table",
    ora_stringent = "list_or_data.table"
  ),
  prototype = list(
    parameters = list(),
    cci_raw = list(),
    cci_detected = list(),
    ora_default = list(),
    ora_stringent = list()
  )
)

#' The scDiffCom Class
#'
#'  A scDiffCom object stores the results of the intercellular communication analysis
#'  performed by \code{run_interaction_analysis} on a Seurat object.
#'
#' @slot parameters A list of the parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table of all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \cdot m} rows,
#'  where n is the number of cell-types and m the number of ligand-receptor pairs.
#' @slot cci_detected A data.table of the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_default A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_detected} table.
#' @slot ora_stringent Same as \code{ora_default} but with a more stringent logfc threshold to specifically focus on strong signals.
#' @slot distributions A list of matrices that contain the distributions over each CCI obtained from the permutation test(s).
#'  Storing such distributions will make the object very heavy for more than 1000 permutations. It is only intended to be used for
#'  benchmarking or debugging.
#'
#' @name scDiffCom-class
#' @rdname scDiffCom-class
#'
setClass(
  Class = "scDiffCom",
  slots = c(
    distributions = "list"
  ),
  prototype = list(
    distributions = list()
  ),
  contains = "scDiffComBase"
)

#' The scDiffComCombined Class
#'
#'  A scDiffComCombined object ...
#'
#' @slot parameters A list of the parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table of all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \cdot m} rows,
#'  where n is the number of cell-types and m the number of ligand-receptor pairs.
#' @slot cci_detected A data.table of the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_default A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_detected} table.
#' @slot ora_stringent Same as \code{ora_default} but with a more stringent logfc threshold to specifically focus on strong signals.
#' @slot combined_object xxx
#'
#' @name scDiffComCombined-class
#' @rdname scDiffComCombined-class
#'
setClass(
  Class = "scDiffComCombined",
  contains = "scDiffComBase",
  slots = c(
    combined_object = "logical",
    ora_combined_default = "list_or_data.table",
    ora_combined_stringent = "list_or_data.table"
  ),
  prototype = list(
    combined_object = TRUE,
    ora_combined_default = list(),
    ora_combined_stringent = list()
  )
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####
# Validity functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity(
  Class = "scDiffComBase",
  method = function(object) {
    validity_check <- c(
      validate_slot_parameters(
        parameters = object@parameters
      ),
      validate_slot_cci_raw(
        parameters = object@parameters,
        cci_raw = object@cci_raw
      ),
      validate_slot_cci_detected(
        parameters = object@parameters,
        cci_detected = object@cci_detected
      ),
      validate_slot_ora_default(
        parameters = object@parameters,
        ora_default = object@cci_ora_default
      ),
      validate_slot_ora_stringent(
        parameters = object@parameters,
        ora_default = object@cci_ora_default
      )
    )
    if (is.null(validity_check)) {
      TRUE
    } else {
      paste0(validity_check, collapse = " AND ")
    }
  }
)

setValidity(
  Class = "scDiffCom",
  method = function(object) {
    validity_check <- c(
      validate_slot_distributions(
        parameters = object@parameters,
        distributions = object@distributions
      )
    )
    if (is.null(validity_check)) {
      TRUE
    } else {
      paste0(validity_check, collapse = " AND ")
    }
  }
)

setValidity(
  Class = "scDiffComCombined",
  method = function(object) {
    validity_check <- NULL
    #validity_check <- c(
    #  validate_slot_combined_object()
    #)
    if (is.null(validity_check)) {
      TRUE
    } else {
      paste0(validity_check, collapse = " AND ")
    }
  }
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####
# Accessors
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Return scDiffCom \code{parameters}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "parameters",
  def = function(object) standardGeneric("parameters")
)

#' @describeIn scDiffComBase Return scDiffCom \code{parameters}.
setMethod(
  f = "parameters",
  signature = "scDiffComBase",
  definition = function(object) object@parameters
)

#' Return scDiffCom \code{cci_raw}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "get_cci_raw",
  def = function(object) standardGeneric("get_cci_raw")
)

#' @describeIn scDiffComBase Return scDiffCom \code{cci_raw}.
setMethod(
  f = "get_cci_raw",
  signature = "scDiffComBase",
  definition = function(object) object@cci_raw
)

#' Return scDiffCom \code{cci_detected}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "get_cci_detected",
  def = function(object) standardGeneric("get_cci_detected")
)

#' @describeIn scDiffComBase Return scDiffCom \code{cci_detected}.
setMethod(
  f = "get_cci_detected",
  signature = "scDiffComBase",
  definition = function(object) object@cci_detected
)

#' Return scDiffCom \code{ora_default}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "get_ora_default",
  def = function(object) standardGeneric("get_ora_default")
)

#' @describeIn scDiffComBase Return scDiffCom \code{ora_default}.
setMethod(
  f = "get_ora_default",
  signature = "scDiffComBase",
  definition = function(object) object@ora_default
)

#' Return scDiffCom \code{ora_stringent}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "get_ora_stringent",
  def = function(object) standardGeneric("get_ora_stringent")
)

#' @describeIn scDiffComBase Return scDiffCom \code{ora_stringent}.
setMethod(
  f = "get_ora_stringent",
  signature = "scDiffComBase",
  definition = function(object) object@ora_stringent
)

#' Return scDiffCom \code{distributions}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "distributions",
  def = function(object) standardGeneric("distributions")
)

#' @describeIn scDiffCom Return scDiffCom \code{distributions}.
setMethod(
  f = "distributions",
  signature = "scDiffCom",
  definition = function(object) object@distributions
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Combine_scDiffCom <- function(
  l,
  object_name,
  perform_combined_ora,
  verbose,
  ora_categories
) {
  list_of_names <- sapply(
    l,
    function(i) {
      if(!methods::is(object = i, "scDiffCom")) {
        stop("All objects to combine must be of class `scDiffCom`")
      }
      i@parameters$object_name
    }
  )
  list_of_params <-
    lapply(
    l,
    function(i) {
      temp <- i@parameters
      temp$object_name <- NULL
      temp$max_nL <- NULL
      temp$max_nR <- NULL
      temp$verbose <- NULL
      temp$seed <- NULL
      return(temp)
    }
  )
  if(length(unique(list_of_params)) != 1) {
    stop("Defining parameters from all scDiffCom objects must be identical.")
  }
  new_param <- l[[1]]@parameters
  new_param$object_name <- object_name
  new_param$max_nL <- max(sapply(l, function(i) i@parameters$max_nL))
  new_param$max_nR <- max(sapply(l, function(i) i@parameters$max_nR))
  new_param$verbose <- NA
  new_param$seed <- sapply(l, function(i) i@parameters$seed)
  list_of_cci_raw <- lapply(l, function(i) i@cci_raw)
  names(list_of_cci_raw) <- list_of_names
  new_cci_raw <- rbindlist(
    l = list_of_cci_raw,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  list_of_cci_detected <- lapply(l, function(i) i@cci_detected)
  names(list_of_cci_detected) <- list_of_names
  new_cci_detected <- rbindlist(
    l = list_of_cci_detected,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  new_object <- methods::new(
    "scDiffComCombined",
    parameters = new_param,
    cci_raw = new_cci_raw,
    cci_detected = new_cci_detected,
    ora_default = list(),
    ora_stringent = list(),
    ora_combined_default = list(),
    ora_combined_stringent = list(),
    combined_object = TRUE
  )
  if (perform_combined_ora) {
    new_object <- run_ora(
      object = new_object,
      categories = ora_categories,
      overwrite = TRUE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = verbose,
      class_signature = "scDiffComCombined"
    )
  }
  # list_of_old_ora_cat <- lapply(
  #   l,
  #   function(i)  names(i@ora_default)
  # )
  # if(length(unique(list_of_old_ora_cat)) != 1 | is.null(list_of_old_ora_cat[[1]])) {
  #   message("Not all ORA categories are the same. Returning empty ORA slots.")
  #   new_ora_default <- list()
  #   new_ora_stringent <- list()
  #   new_ora_combined_default <- list()
  #   new_ora_combined_stringent <- list()
  # } else {
  #   list_ora_default <-  lapply(
  #     l,
  #     function(i) {
  #       temp_ora_default <- i@ora_default
  #     }
  #   )
  #   names(list_ora_default) <- list_of_names
  #   new_ora_default <- rbindlist(
  #     l = list_ora_default,
  #     use.names = TRUE,
  #     fill = TRUE,
  #     idcol =
  #   )
  #   new_ora_stringent <- list()
  #   new_ora_combined_default <- list()
  #   new_ora_combined_stringent <- list()
  # }
  return(new_object)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####
# Methods for scDiffCom-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param new_threshold_quantile_score xxx
#' @param new_threshold_p_value_specificity xxx
#' @param new_threshold_p_value_de xxx
#' @param new_threshold_logfc xxx
#' @param skip_ora xxx
#' @param verbose Should messages be printed?
#'
#' @rdname FilterCCI
#' @export
#' @method FilterCCI scDiffCom
#'
FilterCCI.scDiffCom <- function(
  object,
  new_threshold_quantile_score = NULL,
  new_threshold_p_value_specificity = NULL,
  new_threshold_p_value_de = NULL,
  new_threshold_logfc = NULL,
  skip_ora = FALSE,
  verbose = TRUE,
  ...
) {
  run_filtering_and_ora(
    object = object,
    new_threshold_quantile_score = new_threshold_quantile_score,
    new_threshold_p_value_specificity = new_threshold_p_value_specificity,
    new_threshold_p_value_de = new_threshold_p_value_de,
    new_threshold_logfc = new_threshold_logfc,
    skip_ora = skip_ora,
    verbose = verbose,
    class_signature = "scDiffCom"
  )
}

#' @param new_threshold_quantile_score xxx
#' @param new_threshold_p_value_specificity xxx
#' @param new_threshold_p_value_de xxx
#' @param new_threshold_logfc xxx
#' @param skip_ora xxx
#' @param verbose Should messages be printed?
#'
#' @rdname FilterCCI
#' @export
#' @method FilterCCI scDiffComCombined
#'
FilterCCI.scDiffComCombined <- function(
  object,
  new_threshold_quantile_score = NULL,
  new_threshold_p_value_specificity = NULL,
  new_threshold_p_value_de = NULL,
  new_threshold_logfc = NULL,
  skip_ora = FALSE,
  verbose = TRUE,
  ...
) {
  run_filtering_and_ora(
    object = object,
    new_threshold_quantile_score = new_threshold_quantile_score,
    new_threshold_p_value_specificity = new_threshold_p_value_specificity,
    new_threshold_p_value_de = new_threshold_p_value_de,
    new_threshold_logfc = new_threshold_logfc,
    skip_ora = skip_ora,
    verbose = verbose,
    class_signature = "scDiffComCombined"
  )
}


#' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#'  for each category. Set to \code{c("ER_CELLTYPES", "LR_GENES", "GO_TERMS")} by default.
#' @param overwrite Should existing categories be overwriten in case they match with new categories?
#' @param stringent_or_default Should the default or more stringent ORA be performed?
#' @param stringent_logfc_threshold A more stringent logfc threshold compared to the one stored in \code{parameters}.
#'  Set to \code{NULL} by default, namely no stringent ORA is done.
#' @param verbose Should messages be printed?
#'
#' @rdname RunORA
#' @export
#' @method RunORA scDiffCom
RunORA.scDiffCom <- function(
  object,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS"),
  overwrite = TRUE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  verbose = TRUE,
  ...
) {
  run_ora(
    object = object,
    categories = categories,
    overwrite = overwrite,
    stringent_or_default = stringent_or_default,
    stringent_logfc_threshold = stringent_logfc_threshold,
    verbose = verbose,
    class_signature = "scDiffCom",
    global = FALSE
  )
}

#' @export
#' @method RunORA scDiffComCombined
RunORA.scDiffComCombined <- function(
  object,
  categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS"),
  overwrite = TRUE,
  stringent_or_default = "default",
  stringent_logfc_threshold = NULL,
  global = TRUE,
  verbose = TRUE,
  ...
) {
  run_ora(
    object = object,
    categories = categories,
    overwrite = overwrite,
    stringent_or_default = stringent_or_default,
    stringent_logfc_threshold = stringent_logfc_threshold,
    verbose = verbose,
    class_signature = "scDiffComCombined",
    global = global
  )

}


#' #' Perform filtering and over-representation analysis on intercellular communication patterns
#' #'
#' #' This function is called internally by \code{run_interaction_analysis} after the permutation tests
#' #'  have been performed. Based on the threshold parameters, it returns detected and differentially expressed
#' #'  CCIs in the slot \code{cci_detected} and results of over-representation analysis in \code{cci_ora_default}.
#' #'  The function can be run with new threshold parameters on any scDiffCom object that already contain the slot
#' #'  \code{cci_raw}. This allows the user to test various filtering parameters without the need to rerun the
#' #'  potentially time-consuming permutation analysis. When new thresholds are defined the slot \code{parameters} is
#' #'  modified accordingly.
#' #'  Filtering and over-representation are not independent as the second depends on the first. Therefore, when
#' #'  running filtering with new parameters, we also need to update the ORA results.
#' #'
#' #' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' #' @param new_threshold_quantile_score A new threshold for the quantile score. Set to \code{NULL} by default.
#' #' @param new_threshold_p_value_specificity A new threshold for the specificity p-value specificity. Set to \code{NULL} by default.
#' #' @param new_threshold_p_value_de A new threshold for the differential expression p-value. Set to \code{NULL} by default.
#' #' @param new_threshold_logfc A new threshold for the differential expression logfc. Set to \code{NULL} by default.
#' #' @param skip_ora Should the over-representation analysis be skipped? Set to \code{FALSE} by default.
#' #'  Setting it to \code{FALSE} might be useful if we are interested in performing a rapid parameter scan only on the filtering
#' #'  parameters. Note that in such case, the slot \code{ora_default} is returned empty.
#' #' @param verbose Should messages be printed?
#' #'
#' #' @return An S4 object of class \code{scDiffCom}.
#' #' @export



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
PlotORA <- function(
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
#'     temp_object <- RunORA(
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


