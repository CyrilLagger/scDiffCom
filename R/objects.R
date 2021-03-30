#' @import data.table
#' @import ggplot2
#' @importFrom methods new setClass setClassUnion setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#' @importFrom magrittr "%>%"
#'
NULL

utils::globalVariables(".")

################  Class definitions ################

setClassUnion(
  name = "list_or_data.table",
  members = c("list", "data.table")
)

#' The scDiffComBase Class
#'
#' The scDiffComBase class is a virtual class that provides a template
#' for the \code{scDiffCom} and \code{scDiffComCombined} classes.
#'
#' @slot parameters The list of parameters passed to \code{run_interaction_analysis}.
#' @slot cci_table_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_table_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_table A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_table_detected} table.
#' @name scDiffComBase-class
#' @rdname scDiffComBase-class
#'
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
#'  A scDiffCom object stores the results of the intercellular communication analysis
#'  performed by \code{run_interaction_analysis} on a Seurat object.
#'
#' @slot parameters The list of parameters passed to \code{run_interaction_analysis}.
#' @slot cci_table_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_table_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_table A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_table_detected} table.
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
#'  A scDiffComCombined object stores the results of multiple scDiffCom objects at a single place.
#'  It also allows to run a global over-representation analysis.
#'
#' @slot parameters The list of parameters passed to \code{run_interaction_analysis}.
#' @slot cci_table_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_table_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_table A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_table_detected} table.
#' @name scDiffComCombined-class
#' @rdname scDiffComCombined-class
#'
setClass(
  Class = "scDiffComCombined",
  contains = "scDiffComBase"#,
  # slots = c(
  #   ora_combined_default = "list_or_data.table",
  #   ora_combined_stringent = "list_or_data.table"
  # ),
  # prototype = list(
  #   ora_combined_default = list(),
  #   ora_combined_stringent = list()
  # )
)

################  Validity functions ################

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

################  Accessors ################

#' Return the slot \code{parameters} from an scDiffCom object
#'
#' Return the parameters that have been passed to
#'  \code{\link{run_interaction_analysis}} to compute
#'  the scDiffCom object.
#'
#' @param object An object of class `scDiffCom` or `scDiffComCombined`
#'
#' @return A list of parameters
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

#' Return one of the scDiffCom CCI data.tables
#'
#' @param object An object of class `scDiffCom` or `scDiffComCombined`
#' @param type Specifies the table to return: either "detected" or "raw".
#'  Default is \code{"detected"}.
#' @param simplified Logical indicating if all columns or only the most
#'  informative ones are returned. Default is \code{TRUE}.
#'
#' @return A data.table
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

#' @rdname GetTableCCI
setMethod(
  f = "GetTableCCI",
  signature = "scDiffComCombined",
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
      class_signature = "scDiffComCombined"
    )
  }
)

#' Return some of the scDiffCom ORA data.tables
#'
#' @param object An object of class `scDiffCom` or `scDiffComCombined`
#' @param categories A character vector indicating which ORA tables
#'  to return. Default (\code{"all"}) returns all of them.
#' @param simplified Logical indicating if all columns or only the most
#'  informative ones are returned for each data.table. Default is \code{TRUE}.
#'
#' @return A list of data.tables
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

#' @rdname GetTableORA
setMethod(
  f = "GetTableORA",
  signature = "scDiffComCombined",
  definition = function(
    object,
    categories = "all",
    simplified = TRUE
  ) {
    get_tables_ora(
      object = object,
      categories = categories,
      simplified = simplified,
      class_signature = "scDiffComCombined"
    )
  }
)

#' Return the slot \code{distributions} from an scDiffCom object
#'
#' @param object An scDiffCom object.
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

################  Functions ################

#' Combine several scDiffCom objects together.
#'
#' @param l A list of scDiffCom objects.
#' @param object_name The name of the combined object that will be created.
#' @param verbose Should messages be printed?
#'
#' @return An object of class scDiffComCombined
#' @export
Combine_scDiffCom <- function(
  l,
  object_name,
  verbose
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
  if(length(unique(list_of_names)) != length(list_of_names)) {
    stop("All names of the objects to combine must be different")
  }
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
  if (verbose) message("Binding all raw CCIs.")
  list_of_cci_table_raw <- lapply(l, function(i) i@cci_table_raw)
  names(list_of_cci_table_raw) <- list_of_names
  new_cci_table_raw <- rbindlist(
    l = list_of_cci_table_raw,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  if (verbose) message("Binding all detected CCIs.")
  list_of_cci_table_detected <- lapply(l, function(i) i@cci_table_detected)
  names(list_of_cci_table_detected) <- list_of_names
  new_cci_table_detected <- rbindlist(
    l = list_of_cci_table_detected,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  list_of_ora_categories <- lapply(
    l, function(i) names(i@ora_table)
  )
  if (is.null(list_of_ora_categories[[1]])) {
    if (verbose) {
      message("Empty slot `ora_table` in input object. Returning empty slot `ora_table` in combined object.")
    }
    new_ora_table <- list()
  } else  if (length(unique(list_of_ora_categories)) != 1) {
    if (verbose) {
      message("Objects to bind don't have the same ORA (default) categories. Returning empty slot `ora_table` in combined object.")
    }
    new_ora_table <- list()
  } else {
    if (verbose) {
      message("Binding `ora_table` results.")
    }
    new_ora_table <- sapply(
      list_of_ora_categories[[1]],
      function(categ) {
        list_of_ora <- lapply(
          l,
          function(i) {
            i@ora_table[[categ]]
          }
        )
        names(list_of_ora) <- list_of_names
        rbindlist(
          l = list_of_ora ,
          use.names = TRUE,
          fill = TRUE,
          idcol = "ID"
        )
      },
      USE.NAMES = TRUE,
      simplify = FALSE
    )
  }
  # list_of_ora_categories <- lapply(
  #   l, function(i) names(i@ora_stringent)
  # )
  # if (is.null(list_of_ora_categories[[1]])) {
  #   if (verbose) {
  #     message("Empty slot `ora_stringent` in input object. Returning empty slot `ora_stringent` in combined object.")
  #   }
  #   new_ora_stringent <- list()
  # } else  if (length(unique(list_of_ora_categories)) != 1) {
  #   if (verbose) {
  #     message("Objects to bind don't have the same ORA (stringent) categories. Returning empty slot `ora_stringent` in combined object.")
  #   }
  #   new_ora_stringent <- list()
  # } else {
  #   if (verbose) {
  #     message("Binding `ora_stringent` results.")
  #   }
  #   new_ora_stringent <- sapply(
  #     list_of_ora_categories[[1]],
  #     function(categ) {
  #       list_of_ora <- lapply(
  #         l,
  #         function(i) {
  #           i@ora_stringent[[categ]]
  #         }
  #       )
  #       names(list_of_ora) <- list_of_names
  #       rbindlist(
  #         l = list_of_ora ,
  #         use.names = TRUE,
  #         fill = TRUE,
  #         idcol = "ID"
  #       )
  #     },
  #     USE.NAMES = TRUE,
  #     simplify = FALSE
  #   )
  # }
  new_object <- methods::new(
    "scDiffComCombined",
    parameters = new_param,
    cci_table_raw = new_cci_table_raw,
    cci_table_detected = new_cci_table_detected,
    ora_table = new_ora_table#,
    #ora_stringent = new_ora_stringent,
    #ora_combined_default = list(),
    #ora_combined_stringent = list()
  )
  # if (verbose) {
  #   message("Performing combined ORA (default categories).")
  # }
  # new_object <- run_ora(
  #   object = new_object,
  #   categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS", "ID"),
  #   overwrite = TRUE,
  #   stringent_or_default = "default",
  #   stringent_logfc_threshold = NULL,
  #   verbose = verbose,
  #   class_signature = "scDiffComCombined",
  #   global = TRUE
  # )
  if (verbose) message("Returning combined object.")
  return(new_object)
}

################  Generics and Methods ####

#' xxx
#'
#' @param object xxx
#' @param new_threshold_quantile_score xxx
#' @param new_threshold_p_value_specificity xxx
#' @param new_threshold_p_value_de xxx
#' @param new_threshold_logfc xxx
#' @param skip_ora xxx
#' @param verbose Should messages be printed?
#'
#' @return xxx
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
    verbose = TRUE
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
)

#' @rdname FilterCCI
setMethod(
  f = "FilterCCI",
  signature = "scDiffComCombined",
  definition = function(
    object,
    new_threshold_quantile_score = NULL,
    new_threshold_p_value_specificity = NULL,
    new_threshold_p_value_de = NULL,
    new_threshold_logfc = NULL,
    skip_ora = FALSE,
    verbose = TRUE
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
)

#' xxx
#'
#' @param object xxx
#' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#'  for each category. Set to \code{c("ER_CELLTYPES", "LRI", "GO_TERMS")} by default.
#' @param overwrite Should existing categories be overwritten in case they match with new categories?
#' @param verbose Should messages be printed?
#'
#' @return xxx
#'
#' @export
setGeneric(
  name = "RunORA",
  def = function(
    object,
    categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
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
    categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
    overwrite = TRUE,
    #stringent_or_default = "default",
    #stringent_logfc_threshold = NULL,
    verbose = TRUE
  ) {
    run_ora(
      object = object,
      categories = categories,
      overwrite = overwrite,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = verbose,
      class_signature = "scDiffCom",
      global = FALSE
    )
  }
)

#' @rdname RunORA
setMethod(
  f = "RunORA",
  signature = "scDiffComCombined",
  definition = function(
    object,
    categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
    overwrite = TRUE,
    #stringent_or_default = "default",
    #stringent_logfc_threshold = NULL,
    #global = TRUE,
    verbose = TRUE
  ) {
    run_ora(
      object = object,
      categories = categories,
      overwrite = overwrite,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = verbose,
      class_signature = "scDiffComCombined",
      global = FALSE
    )
  }
)

#' xxx
#'
#' @param object xxx
#' @param category xxx
#' @param regulation xxx
#' @param max_terms_show xxx
#' @param OR_threshold xxx
#' @param p_value_threshold xxx
#' @param ... Extra named arguments passed to PlotORA
#'
#' @return xxx
#'
#' @export
setGeneric(
  name = "PlotORA",
  def = function(
    object,
    category,
    regulation,
    max_terms_show,
    OR_threshold = 1,
    p_value_threshold = 0.05,
    ...
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
    regulation,
    max_terms_show,
    OR_threshold = 1,
    p_value_threshold = 0.05
    #stringent = FALSE,
  ) {
    stringent <- FALSE
    if (stringent) {
      #ora_dt <- get_ora_stringent(object)
    } else {
      ora_dt <- object@ora_table
    }
    if (identical(ora_dt, list())) {
      stop("No ORA data.table to extract from scDiffCom object")
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
)

#' @param subID xxx
#' @rdname PlotORA
setMethod(
  f = "PlotORA",
  signature = "scDiffComCombined",
  definition = function(
    object,
    category,
    regulation,
    max_terms_show,
    OR_threshold = 1,
    p_value_threshold = 0.05,
    subID
    #global = FALSE,
  ) {
    global <- FALSE
    ID <- NULL
    stringent <- FALSE
    if (stringent) {
      if (global) {
        ora_dt <- object@ora_combined_stringent
      } else {
        ora_dt <- object@ora_stringent
      }
    } else {
      if (global) {
        ora_dt <- object@ora_combined_default
      } else {
        ora_dt <- object@ora_table
      }
    }
    if (identical(ora_dt, list())) {
      stop("No ORA data.table to exctract from scDiffCom object")
    }
    if (!(category %in% names(ora_dt))) {
      stop("Can't find the specified ORA category")
    }
    ora_dt <- ora_dt[[category]]
    if (!global) {
      if (!(subID %in% unique(ora_dt$ID))) {
        stop("`subID` must be present in the scDiffComCombined object")
      }
      ora_dt <- ora_dt[ID == subID]
    }

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
)

#' xxx
#'
#' @param object xxx
#' @param network_type xxx
#' @param layout_type xxx
#' @param abbreviation_table xxx
#' @param ... Extra named arguments passed to BuildNetwork
#'
#' @return xxx
#'
#' @export
setGeneric(
  name = "BuildNetwork",
  def = function(
    object,
    network_type = c(
      "condition1_network",
      "condition2_network",
      "difference_network",
      "up_regulated_network",
      "down_regulated_network",
      "ORA_network"
    ),
    layout_type = c(
      "conventional",
      "bipartite"
    ),
    abbreviation_table = NULL,
    ...
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
      "condition1_network",
      "condition2_network",
      "difference_network",
      "up_regulated_network",
      "down_regulated_network",
      "ORA_network"
    ),
    layout_type = c(
      "conventional",
      "bipartite"
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

#' @param ID xxx
#' @rdname BuildNetwork
setMethod(
  f = "BuildNetwork",
  signature = "scDiffComCombined",
  definition = function(
    object,
    network_type = c(
      "condition1_network",
      "condition2_network",
      "difference_network",
      "up_regulated_network",
      "down_regulated_network",
      "ORA_network"
    ),
    layout_type = c(
      "conventional",
      "bipartite"
    ),
    abbreviation_table,
    ID
  ) {
    network_type <- match.arg(network_type)
    layout_type <- match.arg(layout_type)
    if (!(ID %in% unique(object@cci_table_detected$ID))) {
      stop("`ID` must be present in the scDiffComCombined object")
    }
    build_interactive_network(
      object = object,
      network_type = network_type,
      layout_type = layout_type,
      class_signature = "scDiffComCombined",
      subobject_name = ID,
      abbreviation_table = abbreviation_table#,
      #LRIs = LRIs
    )
  }
)




