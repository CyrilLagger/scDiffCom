#' @import data.table
#' @import ggplot2
#' @importFrom methods new setClass setClassUnion setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(
  name = "list_or_data.table",
  members = c("list", "data.table")
)

#' The scDiffCom Class
#'
#' The scDiffCom object stores the results of the intercellular communication analysis
#'  performed by \code{run_interaction_analysis} on a single-cell Seurat object.
#'
#' @slot parameters A list of the parameters passed to \code{run_interaction_analysis}.
#' @slot cci_raw A data.table of all possible cell-cell interactions (CCIs); namely with \eqn{n^2 \cdot m} rows,
#'  where n is the number of cell-types and m the number of ligand-receptor pairs.
#' @slot cci_detected A data.table of the detected CCIs obtained from the raw data.table after filtering.
#' @slot distributions A list of matrices that contain the distributions over each CCI obtained from the permutation test(s).
#'  Storing such distributions will make the object very heavy for more than 1000 permutations. It is only intended to be used for
#'  benchmarking or debugging.
#' @slot ora_default A list of data.tables (one per category) storing the results of the over-representation analysis
#'  performed on the \code{cci_detected} table.
#' @slot ora_stringent Same as \code{ora_default} but with a more stringent logfc threshold to specifically focus on strong signals.
#'
#' @name scDiffCom-class
#' @rdname scDiffCom-class
#' @exportClass scDiffCom
#'
scDiffCom <- setClass(
  Class = "scDiffCom",
  slots = c(
    parameters = "list",
    cci_raw = "list_or_data.table",
    cci_detected = "list_or_data.table",
    distributions = "list",
    ora_default = "list",
    ora_stringent = "list"
  ),
  prototype = list(
    parameters = list(),
    cci_raw = list(),
    cci_detected = list(),
    distributions = list(),
    ora_default = list(),
    ora_stringent = list()
  )
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Validity functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity(
  Class = "scDiffCom",
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
      validate_slot_distributions(
        parameters = object@parameters,
        distributions = object@distributions
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

#' Check the parameters
#'
#' @param params A list of parameters.
#'
#' @noRd
validate_parameters <- function(
  params,
  from_inputs
) {
  res <- NULL
  params <- lapply(
    params,
    eval
  )
  if (!is.character(params$celltype_column_id) | length(params$celltype_column_id) != 1) {
    res <- c(res, "`celltype_column_id` must be a character vector of length 1")
  }
  if (!is.null(params$condition_column_id)) {
    if (!is.character(params$condition_column_id) | length(params$condition_column_id) != 1) {
      res <- c(res, "`condition_column_id` must be NULL or a character vector of length 1")
    }
  }
  if (!is.null(params$cond1_name)) {
    if (!is.character(params$cond1_name) | length(params$cond1_name) != 1) {
      res <- c(res, "`cond1_name` must be NULL or a character vector of length 1")
    }
  }
  if (!is.null(params$cond2_name)) {
    if (!is.character(params$cond2_name) | length(params$cond2_name) != 1) {
      res <- c(res, "`cond2_name` must be NULL or a character vector of length 1")
    }
  }
  if (!is.character(params$assay) | length(params$assay) != 1) {
    res <- c(res, "`assay` must be NULL or a character vector of length 1")
  }
  if (!(params$slot %in% c("counts", "data"))) {
    res <- c(res, "`slot` must be either `data` or `counts`")
  }
  if(!is.logical(params$log_scale) | length(params$log_scale) != 1) {
    res <- c(res, "`log_scale` must be a logical vector of length 1")
  }
  if (!is.numeric(params$threshold_min_cells) | length(params$threshold_min_cells) > 1) {
    res <- c(res, "`threshold_min_cells` must be a numeric vector of length 1")
  } else if (params$threshold_min_cells < 0) {
    res <- c(res, "`threshold_min_cells` must be a non-negative numeric")
  }
  if (!is.numeric(params$threshold_pct) | length(params$threshold_pct) != 1) {
    res <- c(res, "`threshold_pct` must be a numeric vector of length 1")
  } else if(params$threshold_pct < 0 | params$threshold_pct >= 1) {
    res <- c(res, "`threshold_pct` must be a numeric in [0,1[")
  }
  if (!is.character(params$object_name) | length(params$object_name) != 1) {
    res <- c(res, "`object_name` must be a character vector of length 1")
  }
  if(!is.logical(params$permutation_analysis) | length(params$permutation_analysis) != 1) {
    res <- c(res, "`permutation_analysis` must be a logical vector of length 1")
  }
  if (!is.numeric(params$iterations) | length(params$iterations) > 1) {
    res <- c(res, "`iterations` must be a numeric vector of length 1")
  } else if (params$iterations <= 0) {
    res <- c(res, "`iterations` must be a positive integer")
  }
  if (params$iterations > 20000) {
    warning("performing more than 20'000 `iterations` will require a lot of memory")
  }
  if (!is.numeric(params$threshold_quantile_score) | length(params$threshold_quantile_score) != 1) {
    res <- c(res, "`threshold_quantile_score` must be a numeric vector of length 1")
  } else if (params$threshold_quantile_score < 0 |params$ threshold_quantile_score >= 1) {
    res <- c(res, "`threshold_quantile_score` must be a numeric in [0,1[")
  }
  if (!is.numeric(params$threshold_p_value_specificity) | length(params$threshold_p_value_specificity) != 1) {
    res <- c(res, "`threshold_p_value_specificity` must be a numeric vector of length 1")
  } else if (params$threshold_p_value_specificity <= 0 | params$threshold_p_value_specificity > 1) {
    res <- c(res, "`threshold_p_value_specificity` must be a numeric in ]0,1]")
  }
  if (!is.numeric(params$threshold_p_value_de) | length(params$threshold_p_value_de) != 1) {
    res <- c(res, "`threshold_p_value_de` must be a numeric vector of length 1")
  } else if(params$threshold_p_value_de <= 0 | params$threshold_p_value_de > 1) {
    res <- c(res, "`threshold_p_value_de` must be a numeric in ]0,1]")
  }
  if (!is.numeric(params$threshold_logfc) | length(params$threshold_logfc) != 1) {
    res <- c(res, "`threshold_logfc` must be a numeric vector of length 1")
  } else if(params$threshold_logfc <= 0) {
    res <- c(res, "`threshold_logfc` must be a positive numeric")
  }
  if(!is.logical(params$return_distributions) | length(params$return_distributions) != 1) {
    res <- c(res, "`return_distributions` must be a logical vector of length 1")
  }
  if(!is.logical(params$verbose) | length(params$verbose) != 1) {
    res <- c(res, "`verbose` must be a logical vector of length 1")
  }
  if(!from_inputs) {
    if(!is.logical(params$conditional_analysis) | length(params$conditional_analysis) != 1) {
      res <- c(res, "`conditional_analysis` must be a logical vector of length 1")
    }
    if (!is.numeric(params$max_nL) | length(params$max_nL) != 1) {
      res <- c(res, "`max_nL` must be a numeric vector of length 1")
    }
    if (!is.numeric(params$max_nR) | length(params$max_nR) != 1) {
      res <- c(res, "`max_nR` must be a numeric vector of length 1")
    }
  }
  res
}

#' Check the validity of \code{parameters}
#'
#' @param parameters A list of parameters to check the validity of.
#'
#' @noRd
validate_slot_parameters <- function(parameters) {
  res <- validate_parameters(
    params = parameters,
    from_inputs = FALSE
  )
  if(is.null(res)){
    NULL
  } else {
    "@parameters is not formatted the correct way"
  }
}

#' Check the validity of \code{cci_raw}
#'
#' @param parameters A list of parameters.
#' @param cci_raw A data.table to check the validity of.
#'
#' @noRd
validate_slot_cci_raw <- function(
  parameters,
  cci_raw
) {
  NULL
}

#' Check the validity of \code{cci_detected}
#'
#' @param parameters A list of parameters.
#' @param cci_detected A data.table to check the validity of.
#'
#' @noRd
validate_slot_cci_detected <- function(
  parameters,
  cci_detected
) {
  NULL
}

#' Check the validity of \code{distributions}
#'
#' @param parameters A list of parameters.
#' @param distributions A list of matrices to check the validity of.
#'
#' @noRd
validate_slot_distributions <- function(
  parameters,
  distributions
) {
  NULL
}

#' Check the validity of \code{ora_default}
#'
#' @param parameters A list of parameters.
#' @param ora_default A list of data.tables to check the validity of.
#'
#' @noRd
validate_slot_ora_default <- function(
  parameters,
  ora_default
) {
  NULL
}

#' Check the validity of \code{ora_stringent}
#'
#' @param parameters A list of parameters.
#' @param ora_default A list of data.tables to check the validity of.
#'
#' @noRd
validate_slot_ora_stringent <- function(
  parameters,
  ora_default
) {
  NULL
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Accessors
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Return scDiffCom \code{parameters}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("parameters", function(x) standardGeneric("parameters"))

#' @describeIn scDiffCom-class Return scDiffCom \code{parameters}.
setMethod("parameters", "scDiffCom", function(x) x@parameters)

#' Modify scDiffCom \code{parameters}
#'
#' @param x An scDiffCom object.
#' @param value A new list of parameters.
#'
#' @export
setGeneric("parameters<-", function(x, value) standardGeneric("parameters<-"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{parameters}.
#' @inheritParams parameters<-
setMethod("parameters<-", "scDiffCom", function(x, value) {
  x@parameters <- value
  validObject(x)
  return(x)
})

#' Return scDiffCom \code{cci_raw}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("get_cci_raw", function(x) standardGeneric("get_cci_raw"))

#' @describeIn scDiffCom-class Return scDiffCom \code{cci_raw}.
setMethod("get_cci_raw", "scDiffCom", function(x) x@cci_raw)

#' Modify scDiffCom \code{cci_raw}
#'
#' @param x An scDiffCom object.
#' @param new_cci_raw A \code{cci_raw} data.table.
#'
#' @export
setGeneric("set_cci_raw", function(x, new_cci_raw) standardGeneric("set_cci_raw"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{cci_raw}.
#' @inheritParams set_cci_raw
setMethod("set_cci_raw", "scDiffCom", function(x, new_cci_raw) {
  x@cci_raw <- new_cci_raw
  validObject(x)
  return(x)
})

#' Return scDiffCom \code{cci_detected}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("get_cci_detected", function(x) standardGeneric("get_cci_detected"))

#' @describeIn scDiffCom-class Return scDiffCom \code{cci_detected}.
setMethod("get_cci_detected", "scDiffCom", function(x) x@cci_detected)

#' Modify scDiffCom \code{cci_detected}
#'
#' @param x An scDiffCom object.
#' @param new_cci_detected A \code{cci_detected} data.table.
#'
#' @export
setGeneric("set_cci_detected", function(x, new_cci_detected) standardGeneric("set_cci_detected"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{cci_detected}.
#' @inheritParams set_cci_detected
setMethod("set_cci_detected", "scDiffCom", function(x, new_cci_detected) {
  x@cci_detected <- new_cci_detected
  validObject(x)
  return(x)
})

#' Return scDiffCom \code{distributions}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("distributions", function(x) standardGeneric("distributions"))

#' @describeIn scDiffCom-class Return scDiffCom \code{distributions}.
setMethod("distributions", "scDiffCom", function(x) x@distributions)

#' Modify scDiffCom \code{distributions}
#'
#' @param x An scDiffCom object.
#' @param value A new list of distributions.
#'
#' @export
setGeneric("distributions<-", function(x, value) standardGeneric("distributions<-"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{distributions}.
#' @inheritParams distributions<-
setMethod("distributions<-", "scDiffCom", function(x, value) {
  x@distributions <- value
  validObject(x)
  return(x)
})

#' Return scDiffCom \code{ora_default}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("get_ora_default", function(x) standardGeneric("get_ora_default"))

#' @describeIn scDiffCom-class Return scDiffCom \code{ora_default}.
setMethod("get_ora_default", "scDiffCom", function(x) x@ora_default)

#' Modify scDiffCom \code{ora_default}
#'
#' @param x An scDiffCom object.
#' @param new_ora_default A \code{ora_default} list of data.tables.
#'
#' @export
setGeneric("set_ora_default", function(x, new_ora_default) standardGeneric("set_ora_default"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{ora_default}.
#' @inheritParams set_ora_default
setMethod("set_ora_default", "scDiffCom", function(x, new_ora_default) {
  x@ora_default <- new_ora_default
  validObject(x)
  return(x)
})

#' Return scDiffCom \code{ora_stringent}
#'
#' @param x An scDiffCom object.
#'
#' @export
setGeneric("get_ora_stringent", function(x) standardGeneric("get_ora_stringent"))

#' @describeIn scDiffCom-class Return scDiffCom \code{ora_stringent}.
setMethod("get_ora_stringent", "scDiffCom", function(x) x@ora_stringent)

#' Modify scDiffCom \code{ora_stringent}
#'
#' @param x An scDiffCom object.
#' @param new_ora_stringent A \code{ora_stringent} list of data.tables.
#'
#' @export
setGeneric("set_ora_stringent", function(x, new_ora_stringent) standardGeneric("set_ora_stringent"))

#' @describeIn scDiffCom-class Modify scDiffCom \code{ora_stringent}.
#' @inheritParams set_ora_stringent
setMethod("set_ora_stringent", "scDiffCom", function(x, new_ora_stringent) {
  x@ora_stringent <- new_ora_stringent
  validObject(x)
  return(x)
})

