#' @include generics.R
#' @import data.table
#' @import ggplot2
#' @importFrom methods new setClass setClassUnion setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#'
NULL

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
#' @slot cci_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
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
#' @slot parameters The list of parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
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
#'  A scDiffComCombined object stores the results of multiple scDiffCom objects at a single place.
#'  It also allows to run a global over-representation analysis.
#'
#' @slot parameters The list of parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table with all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \times m} rows,
#'  where \eqn{n} is the number of cell-types and \eqn{m} the number of ligand-receptor pairs.
#' @slot cci_detected A data.table with the detected CCIs obtained from the raw data.table after filtering.
#' @slot ora_default A data.table storing the results of the over-representation analysis
#'  performed on the \code{cci_detected} table.
#' @slot ora_stringent Same as \code{ora_default} but with a more stringent logfc threshold to specifically focus on strong signals.
#' @slot ora_combined_default Same as \code{ora_default} but for the global analysis.
#' @slot ora_combined_stringent Same as \code{ora_stringent} but for the global analysis.
#' @name scDiffComCombined-class
#' @rdname scDiffComCombined-class
#'
setClass(
  Class = "scDiffComCombined",
  contains = "scDiffComBase",
  slots = c(
    ora_combined_default = "list_or_data.table",
    ora_combined_stringent = "list_or_data.table"
  ),
  prototype = list(
    ora_combined_default = list(),
    ora_combined_stringent = list()
  )
)


################  Validity functions ################

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

################  Accessors ####

#' Return scDiffCom \code{parameters}
#'
#' @param object An scDiffCom object.
#'
#' @export
setGeneric(
  name = "parameters",
  def = function(object) standardGeneric("parameters")
)

#' @param object xxx
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

#' @param object xxx
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

#' @param object xxx
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

#' @param object xxx
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

#' @param object xxx
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

#' @param object xxx
#' @describeIn scDiffCom Return scDiffCom \code{distributions}.
setMethod(
  f = "distributions",
  signature = "scDiffCom",
  definition = function(object) object@distributions
)

################  Functions ####

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
  list_of_cci_raw <- lapply(l, function(i) i@cci_raw)
  names(list_of_cci_raw) <- list_of_names
  new_cci_raw <- rbindlist(
    l = list_of_cci_raw,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  if (verbose) message("Binding all detected CCIs.")
  list_of_cci_detected <- lapply(l, function(i) i@cci_detected)
  names(list_of_cci_detected) <- list_of_names
  new_cci_detected <- rbindlist(
    l = list_of_cci_detected,
    use.names = TRUE,
    fill = TRUE,
    idcol = "ID"
  )
  list_of_ora_categories <- lapply(
    l, function(i) names(i@ora_default)
  )
  if (is.null(list_of_ora_categories[[1]])) {
    if (verbose) {
      message("Empty slot `ora_default`` in input object. Returning empty slot `ora_default` in combined object.")
    }
    new_ora_default <- list()
  } else  if (length(unique(list_of_ora_categories)) != 1) {
    if (verbose) {
      message("Objects to bind don't have the same ORA (default) categories. Returning empty slot `ora_default` in combined object.")
    }
    new_ora_default <- list()
  } else {
    if (verbose) {
      message("Binding `ora_default` results.")
    }
    new_ora_default <- sapply(
      list_of_ora_categories[[1]],
      function(categ) {
        list_of_ora <- lapply(
          l,
          function(i) {
            i@ora_default[[categ]]
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
  list_of_ora_categories <- lapply(
    l, function(i) names(i@ora_stringent)
  )
  if (is.null(list_of_ora_categories[[1]])) {
    if (verbose) {
      message("Empty slot `ora_stringent` in input object. Returning empty slot `ora_stringent` in combined object.")
    }
    new_ora_stringent <- list()
  } else  if (length(unique(list_of_ora_categories)) != 1) {
    if (verbose) {
      message("Objects to bind don't have the same ORA (stringent) categories. Returning empty slot `ora_stringent` in combined object.")
    }
    new_ora_stringent <- list()
  } else {
    if (verbose) {
      message("Binding `ora_stringent` results.")
    }
    new_ora_stringent <- sapply(
      list_of_ora_categories[[1]],
      function(categ) {
        list_of_ora <- lapply(
          l,
          function(i) {
            i@ora_stringent[[categ]]
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
  new_object <- methods::new(
    "scDiffComCombined",
    parameters = new_param,
    cci_raw = new_cci_raw,
    cci_detected = new_cci_detected,
    ora_default = new_ora_default,
    ora_stringent = new_ora_stringent,
    ora_combined_default = list(),
    ora_combined_stringent = list()
  )
  if (verbose) {
    message("Performing combined ORA (default categories).")
  }
  new_object <- run_ora(
    object = new_object,
    categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "ID"),
    overwrite = TRUE,
    stringent_or_default = "default",
    stringent_logfc_threshold = NULL,
    verbose = verbose,
    class_signature = "scDiffComCombined",
    global = TRUE
  )
  if (verbose) message("Returning combined object.")
  return(new_object)
}

################  Methods for scDiffCom-defined generics ####

#' @param object xxx
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

#' @param object xxx
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

#' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#'  for each category. Set to \code{c("ER_CELLTYPES", "LR_GENES", "GO_TERMS")} by default.
#' @param overwrite Should existing categories be overwriten in case they match with new categories?
#' @param stringent_or_default Should the default or more stringent ORA be performed?
#' @param stringent_logfc_threshold A more stringent logfc threshold compared to the one stored in \code{parameters}.
#'  Set to \code{NULL} by default, namely no stringent ORA is done.
#' @param global xxx
#' @param verbose Should messages be printed?
#'
#' @rdname RunORA
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

#' @param object xxx
#' @param category xxx
#' @param regulation xxx
#' @param max_terms_show xxx
#' @param OR_threshold xxx
#' @param p_value_threshold xxx
#' @param stringent xxx
#'
#' @rdname PlotORA
#' @export
#' @method PlotORA scDiffCom
PlotORA.scDiffCom <- function(
  object,
  category,
  regulation,
  max_terms_show,
  OR_threshold = 1,
  p_value_threshold = 0.05,
  stringent = FALSE,
  ...
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

#' @param object xxx
#' @param subID xxx
#' @param category xxx
#' @param regulation xxx
#' @param max_terms_show xxx
#' @param OR_threshold xxx
#' @param p_value_threshold xxx
#' @param stringent xxx
#' @param global xxx
#'
#' @rdname PlotORA
#' @export
#' @method PlotORA scDiffComCombined
PlotORA.scDiffComCombined <- function(
  object,
  subID,
  category,
  regulation,
  max_terms_show,
  global = FALSE,
  OR_threshold = 1,
  p_value_threshold = 0.05,
  stringent = FALSE,
  ...
) {
  ID <- NULL
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
      ora_dt <- object@ora_default
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

#' @param object xxx
#' @param network_type xxx
#' @param network_layout xxx
#'
#' @rdname BuildNetwork
#' @export
#' @method BuildNetwork scDiffCom
BuildNetwork.scDiffCom = function(
  object,
  network_type = c("ORA", "COUNTS_DIFF", "COUNTS_COND1", "COUNTS_COND2"),
  network_layout = c("bipartite", "celltypes"),
  ...
) {
  network_type <- match.arg(network_type)
  network_layout <- match.arg(network_layout)
  return(
    build_interactive_network(
      object = object,
      network_representation_type = network_type,
      network_layout_type = network_layout,
      class_signature = "scDiffCom",
      subobject_name = NULL
    )
  )
}

#' @param object xxx
#' @param network_type xxx
#' @param network_layout xxx
#' @param ID xxx
#'
#' @rdname BuildNetwork
#' @export
#' @method BuildNetwork scDiffComCombined
BuildNetwork.scDiffComCombined = function(
  object,
  network_type = c("ORA", "COUNTS_DIFF", "COUNTS_COND1", "COUNTS_COND2"),
  network_layout = c("bipartite", "celltypes"),
  ID,
  ...
) {
  network_type <- match.arg(network_type)
  network_layout <- match.arg(network_layout)
  if (!(ID %in% unique(object@cci_detected$ID))) {
    stop("`ID` must be present in the scDiffComCombined object")
  }
  return(
    build_interactive_network(
      object = object,
      network_representation_type = network_type,
      network_layout_type = network_layout,
      class_signature = "scDiffComCombined",
      subobject_name = ID
    )
  )
}
