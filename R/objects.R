#' @import data.table
#' @importFrom methods new setClass setClassUnion
#'  setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#' @importFrom magrittr "%>%"
#'
NULL

utils::globalVariables(".")

################ Class definitions ################

setClassUnion(
  name = "list_or_data.table",
  members = c("list", "data.table")
)

setClass(
  Class = "scDiffComBase",
  contains = "VIRTUAL",
  slots = c(
    parameters = "list",
    cci_table_raw = "list_or_data.table",
    cci_table_detected = "list_or_data.table",
    ora_table = "list_or_data.table"#,
   #ora_stringent = "list_or_data.table"
  ),
  prototype = list(
    parameters = list(),
    cci_table_raw = list(),
    cci_table_detected = list(),
    ora_table = list()#,
    #ora_stringent = list()
  )
)

#' The scDiffCom Class
#'
#' An object of this class stores the intercellular communication results
#' obtained when calling \code{\link{run_interaction_analysis}}.
#'
#' @slot parameters List of parameters passed to
#'  \code{\link{run_interaction_analysis}} and used to build the object.
#' @slot cci_table_raw Data.table with all hypothetic CCIs induced from
#' the original Seurat object and the internal LRI database. Can be erased
#' with \code{\link{EraseRawCCI}} to obtain a lighter object, but might be
#' worth keeping if one intends to modify the filtering parameters
#' (see also our
#'  \href{https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html}{vignette}).
#' @slot cci_table_detected Data.table with only the detected CCIs. If
#'  \code{cci_table_raw} is not \code{NULL}, can be updated with new filtering
#'  parameters without running the full permutation analysis (see
#'  \code{\link{FilterCCI}})
#' @slot ora_table List of data.tables with the results of the
#'  over-representation analysis for each category. Results for additional
#'  categories can be added with \code{\link{RunORA}}.
#' @slot distributions List of matrices with the null distributions of each
#' CCI. \code{NULL} by default.
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

# #' The scDiffComCombined Class
# #'
# #' An object of this class stores the intercellular communication results
# #' from several scDiffCom objects.
# #'
# #' @slot parameters List of parameters used to build the object.
# #' @slot cci_table_raw Data.table with all potential CCIs induced from
# #' the original Seurat object and the internal LRI database.
# #' @slot cci_table_detected Data.table with only the detected CCIs.
# #' @slot ora_table List of data.tables with the results of the
# #'  over-representation analysis.
# #' @name scDiffComCombined-class
# #' @rdname scDiffComCombined-class
# #'
# setClass(
#   Class = "scDiffComCombined",
#   contains = "scDiffComBase"#,
#   # slots = c(
#   #   ora_combined_default = "list_or_data.table",
#   #   ora_combined_stringent = "list_or_data.table"
#   # ),
#   # prototype = list(
#   #   ora_combined_default = list(),
#   #   ora_combined_stringent = list()
#   # )
# )

################ Validity functions ################

setValidity(
  Class = "scDiffComBase",
  method = function(object) {
    validity_check <- c(
      validate_slot_parameters(
        parameters = object@parameters
      ),
      validate_slot_cci_table_raw(
        parameters = object@parameters,
        cci_table_raw = object@cci_table_raw
      ),
      validate_slot_cci_table_detected(
        parameters = object@parameters,
        cci_table_detected = object@cci_table_detected
      ),
      validate_slot_ora_table(
        parameters = object@parameters,
        ora_table = object@cci_ora_table
      )#,
      # validate_slot_ora_stringent(
      #   parameters = object@parameters,
      #   ora_table = object@cci_ora_table
      # )
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

# setValidity(
#   Class = "scDiffComCombined",
#   method = function(object) {
#     validity_check <- NULL
#     #validity_check <- c(
#     #  validate_slot_combined_object()
#     #)
#     if (is.null(validity_check)) {
#       TRUE
#     } else {
#       paste0(validity_check, collapse = " AND ")
#     }
#   }
# )

################ Accessors ################

#' Return the slot \code{parameters} from a scDiffCom object
#'
#' Return the parameters that have been passed to
#'  \code{\link{run_interaction_analysis}} as well as a few other
#'  parameters computed alongside the analysis.
#'
#' @param object \code{scDiffCom} object
#'
#' @return A list of parameters.
#'
#' @export
setGeneric(
  name = "GetParameters",
  def = function(object) standardGeneric("GetParameters")
)

#' @rdname GetParameters
setMethod(
  f = "GetParameters",
  signature = "scDiffComBase",
  definition = function(object) object@parameters
)

#' Return (a subset) of the slot \code{cci_table_raw} or
#'  \code{cci_table_detected} from a scDiffCom object
#'
#' @param object \code{scDiffCom} object
#' @param type Table to extract information from. Can be either
#' \code{"detected"} (default) or \code{"raw"}.
#' @param simplified If \code{TRUE} (default) only the most
#'  informative columns of the data.table are returned.
#'
#' @return A data.table.
#'
#' @export
setGeneric(
  name = "GetTableCCI",
  def = function(object, type, simplified) standardGeneric("GetTableCCI"),
  signature = "object"
)

#' @rdname GetTableCCI
setMethod(
  f = "GetTableCCI",
  signature = "scDiffCom",
  definition = function(
    object,
    type = c("detected", "raw"),
    simplified = TRUE
  ) {
    type <- match.arg(type)
    get_table_cci(
      object = object,
      type = type,
      simplified = simplified,
      class_signature = "scDiffCom"
    )
  }
)

# #' @rdname GetTableCCI
# setMethod(
#   f = "GetTableCCI",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     type = c("detected", "raw"),
#     simplified = TRUE
#   ) {
#     type <- match.arg(type)
#     get_table_cci(
#       object = object,
#       type = type,
#       simplified = simplified,
#       class_signature = "scDiffComCombined"
#     )
#   }
# )

#' Return some or all ORA tables from the slot \code{ora_table}
#' from a scDiffCom object
#'
#' @param object \code{scDiffCom} object
#' @param categories Names of the ORA categories to return. If \code{"all"}
#' (default), returns all of them.
#' @param simplified If \code{TRUE} (default) only the most
#'  informative columns of the data.table are returned.
#'
#' @return A list of data.tables.
#'
#' @export
setGeneric(
  name = "GetTableORA",
  def = function(object, categories, simplified) standardGeneric("GetTableORA"),
  signature = "object"
)

#' @rdname GetTableORA
setMethod(
  f = "GetTableORA",
  signature = "scDiffCom",
  definition = function(
    object,
    categories = "all",
    simplified = TRUE
    ) {
    get_tables_ora(
      object = object,
      categories = categories,
      simplified = simplified,
      class_signature = "scDiffCom"
    )
  }
)

# #' @rdname GetTableORA
# setMethod(
#   f = "GetTableORA",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     categories = "all",
#     simplified = TRUE
#   ) {
#     get_tables_ora(
#       object = object,
#       categories = categories,
#       simplified = simplified,
#       class_signature = "scDiffComCombined"
#     )
#   }
# )

#' Return the slot \code{distributions} from a scDiffCom object
#'
#' @return List of matrices with the null distributions of each
#' CCI.
#'
#' @param object \code{scDiffCom} object
#'
#' @export
setGeneric(
  name = "GetDistributions",
  def = function(object) standardGeneric("GetDistributions")
)

#' @rdname GetDistributions
setMethod(
  f = "GetDistributions",
  signature = "scDiffCom",
  definition = function(object) object@distributions
)

################ Generics and Methods ####

#' Display a scDiffCom object
#'
#' @param object \code{scDiffCom} object
setMethod(
  "show",
  signature = "scDiffCom",
  definition = function(
    object
  ) {
    line_1 <- paste0(
      "An object of class scDiffCom with name ",
      object@parameters$object_name
    )
    if (is.null(object@parameters$seurat_condition_id)) {
      if (!object@parameters$permutation_analysis) {
        line_2 <- paste0(
          "Analysis performed: ",
          "detection only (without statistical test!) "
        )
      } else {
        line_2 <- paste0(
          "Analysis performed: ",
          "detection only (with statistical test) "
        )
      }
    } else {
      if (!object@parameters$permutation_analysis) {
        line_2 <- paste0(
          "Analysis performed: ",
          "differential analysis (without statistical test!) between ",
          object@parameters$seurat_condition_id$cond1_name,
          " and ",
          object@parameters$seurat_condition_id$cond2_name,
          " cells"
        )
      } else {
        line_2 <- paste0(
          "Analysis performed: ",
          "differential analysis between ",
          object@parameters$seurat_condition_id$cond1_name,
          " and ",
          object@parameters$seurat_condition_id$cond2_name,
          " cells"

        )
      }
    }
    line_3 <- paste0(
      nrow(object@cci_table_detected),
      " detected CCIs across ",
      length(unique(object@cci_table_detected$EMITTER_CELLTYPE)),
      " cell types"
    )
    if (identical(object@ora_table, list())) {
      line_4 <- "No over-representation results"
    } else {
      is_ora <- TRUE
      ora_categories <- names(object@ora_table)
      line_4 <- paste0(
        "Over-representation results for ",
        paste0(ora_categories, collapse = ", ")
      )
    }
    cat(
      paste(
        line_1,
        line_2,
        line_3,
        line_4,
        sep = "\n"
      )
    )
  }
)

#' Create a copy of a scDiffCom object without \code{cci_table_raw}
#'
#' This function will replace \code{cci_table_raw} by an empty list. Useful to
#' save space for large datasets. However, after this operation,
#' no filtering can be re-run on the new object, meaning that obtaining
#' results for different filtering parameters will require the perform the full
#' analysis from scratch.
#'
#' @param object \code{scDiffCom} object
#'
#' @return A scDiffCom object with an empty list for \code{cci_table_raw}.
#'
#' @export
setGeneric(
  name = "EraseRawCCI",
  def = function(object) standardGeneric("EraseRawCCI"),
  signature = "object"
)

#' @rdname EraseRawCCI
setMethod(
  f = "EraseRawCCI",
  signature = "scDiffCom",
  definition = function(
    object
  ) {
    new_object <- object
    new_object@cci_table_raw <- list()
    new_object
  }
)

#' Filter a scDiffCom object with new filtering parameters
#'
#' Filtering (and ORA) is performed with new parameter on an existing
#' \code{scDiffCom} object. The slots \code{cci_table_detected} and
#' \code{ora_table} are updated accordingly.
#'
#' When \code{FilterCCI} is called with new parameters, both
#'  \code{cci_table_detected} and \code{ora_table} are updated. For
#'  ORA, a call to \code{RunORA} is automatically performed on all standard
#'  categories. Additional user-defined ORA categories can be added via the
#'  parameter \code{extra_annotations}. The data.frames or data.tables in this
#'  list must have exactly two columns that indicates a relationship between
#'  values from a standard category (first column) to values of the new
#'  category (second column). As a typical example, this
#'  \href{https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html}{vignette}
#'  shows how to perform ORA on cell type families attached to each cell type.
#'
#' @param object \code{scDiffCom} object
#' @param new_threshold_quantile_score New threshold value to update
#'  \code{threshold_quantile_score}. If \code{NULL} (default),
#'  the value is not updated.
#' @param new_threshold_p_value_specificity New threshold value to update
#'  \code{threshold_p_value_specificity}. If \code{NULL} (default),
#'  the value is not updated.
#' @param new_threshold_p_value_de New threshold value to update
#'  \code{threshold_p_value_de}. If \code{NULL} (default),
#'  the value is not updated.
#' @param new_threshold_logfc New threshold value to update
#'  \code{threshold_logfc}. If \code{NULL} (default),
#'  the value is not updated.
#' @param skip_ora Default is \code{FALSE}. If \code{TRUE}, ORA is not
#' performed with the new parameters and \code{ora_table} is set to an
#' empty list. May be useful if one wants to quickly test (loop-over) several
#' values of parameters and by-passing the ORA computing time.
#' @param extra_annotations Convenience parameter to perform ORA on user-defined
#' non-standard categories. If \code{NULL} (default), ORA is
#' performed on standard categories. Otherwise it must be a list of data.tables
#' or data.frames (see Details).
#' @param verbose If \code{TRUE} (default) progress messages are printed.
#'
#' @return A scDiffCom object with updated results in \code{cci_table_detected}
#' and \code{ora_table}.
#'
#' @export
setGeneric(
  name = "FilterCCI",
  def = function(
    object,
    new_threshold_quantile_score = NULL,
    new_threshold_p_value_specificity = NULL,
    new_threshold_p_value_de = NULL,
    new_threshold_logfc = NULL,
    skip_ora = FALSE,
    extra_annotations = NULL,
    verbose = TRUE
    ) standardGeneric("FilterCCI"),
  signature = "object"
)

#' @rdname FilterCCI
setMethod(
  f = "FilterCCI",
  signature = "scDiffCom",
  definition = function(
    object,
    new_threshold_quantile_score = NULL,
    new_threshold_p_value_specificity = NULL,
    new_threshold_p_value_de = NULL,
    new_threshold_logfc = NULL,
    skip_ora = FALSE,
    extra_annotations = NULL,
    verbose = TRUE
  ) {
    run_filtering_and_ora(
      object = object,
      new_threshold_quantile_score = new_threshold_quantile_score,
      new_threshold_p_value_specificity = new_threshold_p_value_specificity,
      new_threshold_p_value_de = new_threshold_p_value_de,
      new_threshold_logfc = new_threshold_logfc,
      skip_ora = skip_ora,
      extra_annotations = extra_annotations,
      verbose = verbose,
      class_signature = "scDiffCom"
    )
  }
)

#' Run over-representation analysis
#'
#' Perform over-representation analysis (ORA) on a scDiffCom object, with
#' the possibility to define new categories in addition to the standard
#' ones supported by default.
#'
#' Additional user-defined ORA categories can be added via the
#'  parameter \code{extra_annotations}. The data.frames or data.tables in this
#'  list must have exactly two columns that indicates a relationship between
#'  values from a standard category (first column) to values of the new
#'  category (second column). As a typical example, this
#'  \href{https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html}{vignette}
#'  shows how to perform ORA on cell type families attached to each cell type.
#'
#' @param object \code{scDiffCom} object
#' @param categories Names of the standard categories on which to perform ORA.
#'  Default is all standard categories, namely
#'  \code{c("LRI", "LIGAND_COMPLEX", "RECEPTOR_COMPLEX", "ER_CELLTYPES",
#'  "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "GO_TERMS", "KEGG_PWS")}
#' @param extra_annotations Convenience parameter to perform ORA on user-defined
#' non-standard categories. If \code{NULL} (default), ORA is
#' performed only on standard categories from \code{categories}. Otherwise it
#' must be a list of data.tables
#' or data.frames (see Details).
#' @param overwrite If \code{TRUE} (default), previous results are overwritten
#' in case they correspond to a category passed in \code{categories}.
#' @param verbose If \code{TRUE} (default), progress messages are printed.
#'
#' @return A scDiffCom object with updated slot \code{ora_table}.
#'
#' @export
setGeneric(
  name = "RunORA",
  def = function(
    object,
    categories = c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS"
    ),
    extra_annotations = NULL,
    overwrite = TRUE,
    verbose = TRUE
  ) standardGeneric("RunORA"),
  signature = "object"
)

#' @rdname RunORA
setMethod(
  f = "RunORA",
  signature = "scDiffCom",
  definition = function(
    object,
    categories =c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS"
    ),
    extra_annotations = NULL,
    overwrite = TRUE,
    #stringent_or_default = "default",
    #stringent_logfc_threshold = NULL,
    verbose = TRUE
  ) {
    run_ora(
      object = object,
      categories = categories,
      extra_annotations = extra_annotations,
      overwrite = overwrite,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = verbose,
      class_signature = "scDiffCom",
      global = FALSE
    )
  }
)


#' Display top over-represented keywords from a category of interest
#'
#' Plot a graph that shows the top over-represented terms of a given
#' category for a given regulation. Terms are ordered by their ORA scores,
#' computed from their odds ratios and adjusted p-values.
#'
#' The ORA score is computed as the product between \code{log2(odds ratio)} and
#' \code{-log10(adj. p-value)}.
#'
#' @param object \code{scDiffCom} object
#' @param category ORA category to display. Must be the name of one of the
#' category present in \code{ora_table}.
#' @param regulation ORA regulation to display. Can be either \code{UP}
#' (default), \code{DOWN} or \code{FLAT}.
#' @param max_terms_show Maximum number of terms to display. Default is
#'  \code{20}.
#' @param GO_aspect Name of the GO aspect to display when
#'  \code{category == "GO_TERMS"}. Can be either \code{biological_process} (
#'  default), \code{molecular_function} or \code{cellular_component}.
#' @param OR_threshold Only the terms with an odds ratio above this threshold
#' will be displayed. Default is \code{1}, meaning no filtering is performed.
#' @param bh_p_value_threshold Only the terms with an adjusted p-value below
#' this threshold (and always below 0.05) will be displayed. Default is
#'  \code{0.05}.
#'
#' @return A ggplot object.
#'
#' @export
setGeneric(
  name = "PlotORA",
  def = function(
    object,
    category,
    regulation = c("UP", "DOWN", "FLAT"),
    max_terms_show = 20,
    GO_aspect = c(
      "biological_process",
      "molecular_function",
      "cellular_component"
      ),
    OR_threshold = 1,
    bh_p_value_threshold = 0.05
  ) standardGeneric("PlotORA"),
  signature = "object"
)

#' @rdname PlotORA
setMethod(
  f = "PlotORA",
  signature = "scDiffCom",
  definition = function(
    object,
    category,
    regulation = c("UP", "DOWN", "FLAT"),
    max_terms_show = 20,
    GO_aspect = c(
      "biological_process",
      "molecular_function",
      "cellular_component"
    ),
    OR_threshold = 1,
    bh_p_value_threshold = 0.05
  ) {
    regulation <- match.arg(regulation)
    GO_aspect <- match.arg(GO_aspect)
    ora_dt <- object@ora_table
    if (identical(ora_dt, list())) {
      stop("No ORA data.table to extract from scDiffCom object")
    }
    if (!(category %in% names(ora_dt))) {
      stop("Can't find the specified ORA category")
    }
    ora_dt <- copy(ora_dt[[category]])
    plot_ora(
      ora_dt = ora_dt,
      category = category,
      regulation = regulation,
      max_terms_show = max_terms_show,
      GO_aspect = GO_aspect,
      OR_threshold = OR_threshold ,
      bh_p_value_threshold = bh_p_value_threshold
    )
  }
)

#' Display cell-type to cell-type interactive networks
#'
#' Create and plot an interactive network that summarize how
#' cell-types and their interactions are over-represented.
#'
#' @param object \code{scDiffCom} object
#' @param network_type Type of network to display. Currently, only
#' \code{ORA_network} (default) is supported.
#' @param layout_type Layout of the network to display. Can either be
#' \code{"bipartite"} (default) or \code{"conventional"}.
#' @param abbreviation_table Table with abbreviations
#' for the cell types present in the \code{object}. If \code{NULL} (default),
#' full names of the cell-types are displayed. Otherwise, it must be a
#' data.frame or data.table with exactly two columns with names
#' \code{ORIGINAL_CELLTYPE} and \code{ABBR_CELLTYPE}.
#'
#' @return A visNetwork object.
#'
#' @export
setGeneric(
  name = "BuildNetwork",
  def = function(
    object,
    network_type = c(
      #"condition1_network",
      #"condition2_network",
      #"difference_network",
      #"up_regulated_network",
      #"down_regulated_network",
      "ORA_network"
    ),
    layout_type = c(
      "bipartite",
      "conventional"
    ),
    abbreviation_table = NULL
  ) standardGeneric("BuildNetwork"),
  signature = "object"
)

#' @rdname BuildNetwork
setMethod(
  f = "BuildNetwork",
  signature = "scDiffCom",
  definition = function(
    object,
    network_type = c(
      #"condition1_network",
      #"condition2_network",
      #"difference_network",
      #"up_regulated_network",
      #"down_regulated_network",
      "ORA_network"
    ),
    layout_type = c(
      "bipartite",
      "conventional"
    ),
    abbreviation_table = NULL
  ) {
    network_type <- match.arg(network_type)
    layout_type <- match.arg(layout_type)
    build_interactive_network(
      object = object,
      network_type = network_type,
      layout_type = layout_type,
      class_signature = "scDiffCom",
      subobject_name = NULL,
      abbreviation_table = abbreviation_table#,
      #LRIs = LRIs
    )
  }
)

################ Not implemented generics/methods ################

# #' xxx
# #'
# #' @param object xxx
# setMethod(
#   "show",
#   signature = "scDiffComCombined",
#   definition = function(
#     object
#   ) {
#     line_1 <- paste0(
#       "An object of class scDiffComCombined"
#     )
#     cat(
#       line_1
#      )
#   }
# )

# #' @rdname FilterCCI
# setMethod(
#   f = "FilterCCI",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     new_threshold_quantile_score = NULL,
#     new_threshold_p_value_specificity = NULL,
#     new_threshold_p_value_de = NULL,
#     new_threshold_logfc = NULL,
#     skip_ora = FALSE,
#     verbose = TRUE
#   ) {
#     run_filtering_and_ora(
#       object = object,
#       new_threshold_quantile_score = new_threshold_quantile_score,
#       new_threshold_p_value_specificity = new_threshold_p_value_specificity,
#       new_threshold_p_value_de = new_threshold_p_value_de,
#       new_threshold_logfc = new_threshold_logfc,
#       skip_ora = skip_ora,
#       verbose = verbose,
#       class_signature = "scDiffComCombined"
#     )
#   }
# )

# #' @rdname RunORA
# setMethod(
#   f = "RunORA",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
#     overwrite = TRUE,
#     #stringent_or_default = "default",
#     #stringent_logfc_threshold = NULL,
#     #global = TRUE,
#     verbose = TRUE
#   ) {
#     run_ora(
#       object = object,
#       categories = categories,
#       overwrite = overwrite,
#       stringent_or_default = "default",
#       stringent_logfc_threshold = NULL,
#       verbose = verbose,
#       class_signature = "scDiffComCombined",
#       global = FALSE
#     )
#   }
# )

# #' @param subID xxx
# #' @rdname PlotORA
# setMethod(
#   f = "PlotORA",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     category,
#     regulation = c("UP", "DOWN", "FLAT"),
#     max_terms_show,
#     GO_aspect = c(
#       "biological_process",
#       "molecular_function",
#       "cellular_component"
#     ),
#     OR_threshold = 1,
#     bh_p_value_threshold = 0.05,
#     subID
#   ) {
#     ID <- NULL
#     regulation <- match.arg(regulation)
#     GO_aspect <- match.arg(GO_aspect)
#     ora_dt <- object@ora_table
#     if (identical(ora_dt, list())) {
#       stop("No ORA data.table to exctract from scDiffCom object")
#     }
#     if (!(category %in% names(ora_dt))) {
#       stop("Can't find the specified ORA category")
#     }
#     ora_dt <- copy(ora_dt[[category]])
#     if (!(subID %in% unique(ora_dt$ID))) {
#       stop("`subID` must be present in the scDiffComCombined object")
#     }
#     ora_dt <- ora_dt[ID == subID]
#     plot_ora(
#       ora_dt = ora_dt,
#       category = category,
#       regulation = regulation,
#       max_terms_show = max_terms_show,
#       GO_aspect = GO_aspect,
#       OR_threshold = OR_threshold ,
#       bh_p_value_threshold = bh_p_value_threshold
#     )
#   }
# )

# #' @param ID xxx
# #' @rdname BuildNetwork
# setMethod(
#   f = "BuildNetwork",
#   signature = "scDiffComCombined",
#   definition = function(
#     object,
#     network_type = c(
#       "condition1_network",
#       "condition2_network",
#       "difference_network",
#       "up_regulated_network",
#       "down_regulated_network",
#       "ORA_network"
#     ),
#     layout_type = c(
#       "conventional",
#       "bipartite"
#     ),
#     abbreviation_table,
#     ID
#   ) {
#     network_type <- match.arg(network_type)
#     layout_type <- match.arg(layout_type)
#     if (!(ID %in% unique(object@cci_table_detected$ID))) {
#       stop("`ID` must be present in the scDiffComCombined object")
#     }
#     build_interactive_network(
#       object = object,
#       network_type = network_type,
#       layout_type = layout_type,
#       class_signature = "scDiffComCombined",
#       subobject_name = ID,
#       abbreviation_table = abbreviation_table#,
#       #LRIs = LRIs
#     )
#   }
# )

################ Functions (not implemented) ################

# #' Combine several scDiffCom objects together.
# #'
# #' @param l A list of scDiffCom objects.
# #' @param object_name The name of the combined object that will be created.
# #' @param verbose Should messages be printed?
# #'
# #' @return An object of class scDiffComCombined
# #' @export
# Combine_scDiffCom <- function(
#   l,
#   object_name,
#   verbose
# ) {
#   list_of_names <- sapply(
#     l,
#     function(i) {
#       if(!methods::is(object = i, "scDiffCom")) {
#         stop("All objects to combine must be of class `scDiffCom`")
#       }
#       i@parameters$object_name
#     }
#   )
#   if(length(unique(list_of_names)) != length(list_of_names)) {
#     stop("All names of the objects to combine must be different")
#   }
#   list_of_params <-
#     lapply(
#     l,
#     function(i) {
#       temp <- i@parameters
#       temp$object_name <- NULL
#       temp$max_nL <- NULL
#       temp$max_nR <- NULL
#       temp$verbose <- NULL
#       temp$seed <- NULL
#       return(temp)
#     }
#   )
#   if(length(unique(list_of_params)) != 1) {
#     stop("Defining parameters from all scDiffCom objects must be identical.")
#   }
#   new_param <- l[[1]]@parameters
#   new_param$object_name <- object_name
#   new_param$max_nL <- max(sapply(l, function(i) i@parameters$max_nL))
#   new_param$max_nR <- max(sapply(l, function(i) i@parameters$max_nR))
#   new_param$verbose <- NA
#   new_param$seed <- sapply(l, function(i) i@parameters$seed)
#   if (verbose) message("Binding all raw CCIs.")
#   list_of_cci_table_raw <- lapply(l, function(i) i@cci_table_raw)
#   names(list_of_cci_table_raw) <- list_of_names
#   new_cci_table_raw <- rbindlist(
#     l = list_of_cci_table_raw,
#     use.names = TRUE,
#     fill = TRUE,
#     idcol = "ID"
#   )
#   if (verbose) message("Binding all detected CCIs.")
#   list_of_cci_table_detected <- lapply(l, function(i) i@cci_table_detected)
#   names(list_of_cci_table_detected) <- list_of_names
#   new_cci_table_detected <- rbindlist(
#     l = list_of_cci_table_detected,
#     use.names = TRUE,
#     fill = TRUE,
#     idcol = "ID"
#   )
#   list_of_ora_categories <- lapply(
#     l, function(i) names(i@ora_table)
#   )
#   if (is.null(list_of_ora_categories[[1]])) {
#     if (verbose) {
#       message(
#         paste0(
#           "Empty slot `ora_table` in input object.",
#           " Returning empty slot `ora_table` in combined object."
#           )
#       )
#     }
#     new_ora_table <- list()
#   } else  if (length(unique(list_of_ora_categories)) != 1) {
#     if (verbose) {
#       message(
#         paste0(
#           "Objects to bind don't have the same ORA (default) categories.",
#           " Returning empty slot `ora_table` in combined object."
#           )
#       )
#     }
#     new_ora_table <- list()
#   } else {
#     if (verbose) {
#       message("Binding `ora_table` results.")
#     }
#     new_ora_table <- sapply(
#       list_of_ora_categories[[1]],
#       function(categ) {
#         list_of_ora <- lapply(
#           l,
#           function(i) {
#             i@ora_table[[categ]]
#           }
#         )
#         names(list_of_ora) <- list_of_names
#         rbindlist(
#           l = list_of_ora ,
#           use.names = TRUE,
#           fill = TRUE,
#           idcol = "ID"
#         )
#       },
#       USE.NAMES = TRUE,
#       simplify = FALSE
#     )
#   }
#   # list_of_ora_categories <- lapply(
#   #   l, function(i) names(i@ora_stringent)
#   # )
#   # if (is.null(list_of_ora_categories[[1]])) {
#   #   if (verbose) {
#   #     message("Empty slot `ora_stringent` in input object. Returning empty slot `ora_stringent` in combined object.")
#   #   }
#   #   new_ora_stringent <- list()
#   # } else  if (length(unique(list_of_ora_categories)) != 1) {
#   #   if (verbose) {
#   #     message("Objects to bind don't have the same ORA (stringent) categories. Returning empty slot `ora_stringent` in combined object.")
#   #   }
#   #   new_ora_stringent <- list()
#   # } else {
#   #   if (verbose) {
#   #     message("Binding `ora_stringent` results.")
#   #   }
#   #   new_ora_stringent <- sapply(
#   #     list_of_ora_categories[[1]],
#   #     function(categ) {
#   #       list_of_ora <- lapply(
#   #         l,
#   #         function(i) {
#   #           i@ora_stringent[[categ]]
#   #         }
#   #       )
#   #       names(list_of_ora) <- list_of_names
#   #       rbindlist(
#   #         l = list_of_ora ,
#   #         use.names = TRUE,
#   #         fill = TRUE,
#   #         idcol = "ID"
#   #       )
#   #     },
#   #     USE.NAMES = TRUE,
#   #     simplify = FALSE
#   #   )
#   # }
#   new_object <- methods::new(
#     "scDiffComCombined",
#     parameters = new_param,
#     cci_table_raw = new_cci_table_raw,
#     cci_table_detected = new_cci_table_detected,
#     ora_table = new_ora_table#,
#     #ora_stringent = new_ora_stringent,
#     #ora_combined_default = list(),
#     #ora_combined_stringent = list()
#   )
#   # if (verbose) {
#   #   message("Performing combined ORA (default categories).")
#   # }
#   # new_object <- run_ora(
#   #   object = new_object,
#   #   categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS", "ID"),
#   #   overwrite = TRUE,
#   #   stringent_or_default = "default",
#   #   stringent_logfc_threshold = NULL,
#   #   verbose = verbose,
#   #   class_signature = "scDiffComCombined",
#   #   global = TRUE
#   # )
#   if (verbose) message("Returning combined object.")
#   return(new_object)
# }


