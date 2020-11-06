#' @import data.table
#' @importFrom methods new setClass setClassUnion setValidity setGeneric validObject
#' @importFrom DelayedArray rowsum
#'
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion("list_or_dt", c("data.table", "list"))

scDiffCom <- setClass(
  Class = "scDiffCom",
  slots = c(
    parameters = "list",
    cci_table_raw = "list_or_dt",
    cci_table_filtered = "list_or_dt",
    distributions = "list",
    ORA = "list"
  ),
  prototype = list(
    parameters = list(),
    cci_table_raw = list(),
    cci_table_filtered = list(),
    distributions = list(),
    ORA = list()
  )
)

setValidity(
  Class = "scDiffCom",
  method = function(
    object
  ) {
    ## add check about ORA duplicated
    parameters_expected <- c(
      object_name = "character",
      celltype_column_id = "character",
      LR_info = "list",
      condition_info = "list",
      assay = "character",
      slot = "character",
      log_scale = "logical",
      min_cells = "double",
      pct_threshold = "double",
      permutation_analysis = "logical",
      iterations = "double",
      cutoff_quantile_score = "double",
      cutoff_pval_specificity = "double",
      cutoff_pval_de = "double",
      cutoff_logfc = "double",
      return_distr = "logical"
    )
    parameters_actual <- sapply(object@parameters, typeof)
    if(!identical(parameters_expected, parameters_actual)) {
      "@parameters is not formatted the correct way"
    } else {
      TRUE
    }
  }
)

scDiffCom <- function(
  parameters,
  cci_table_raw,
  cci_table_filtered,
  distributions,
  ORA
) {
  new(
    "scDiffCom",
    parameters = parameters,
    cci_table_raw = cci_table_raw,
    cci_table_filtered = cci_table_filtered,
    distributions = distributions,
    ORA = ORA
  )
}

setGeneric("parameters", function(x) standardGeneric("parameters"))
setMethod("parameters", "scDiffCom", function(x) x@parameters)
setGeneric("parameters<-", function(x, value) standardGeneric("parameters<-"))
setMethod("parameters<-", "scDiffCom", function(x, value) {
  x@parameters <- value
  validObject(x)
  return(x)
})

setGeneric("get_cci_table_raw", function(x) standardGeneric("get_cci_table_raw"))
setMethod("get_cci_table_raw", "scDiffCom", function(x) x@cci_table_raw)
setGeneric("set_cci_table_raw", function(x, new_cci_table_raw) standardGeneric("set_cci_table_raw"))
setMethod("set_cci_table_raw", "scDiffCom", function(x, new_cci_table_raw) {
  x@cci_table_raw <- new_cci_table_raw
  validObject(x)
  return(x)
})

setGeneric("get_cci_table_filtered", function(x) standardGeneric("get_cci_table_filtered"))
setMethod("get_cci_table_filtered", "scDiffCom", function(x) x@cci_table_filtered)
setGeneric("set_cci_table_filtered", function(x, new_cci_table_filtered) standardGeneric("set_cci_table_filtered"))
setMethod("set_cci_table_filtered", "scDiffCom", function(x, new_cci_table_filtered) {
  x@cci_table_filtered <- new_cci_table_filtered
  validObject(x)
  return(x)
})

setGeneric("distributions", function(x) standardGeneric("distributions"))
setMethod("distributions", "scDiffCom", function(x) x@distributions)
setGeneric("distributions<-", function(x, value) standardGeneric("distributions<-"))
setMethod("distributions<-", "scDiffCom", function(x, value) {
  x@distributions <- value
  validObject(x)
  return(x)
})

setGeneric("get_ora_tables", function(x) standardGeneric("get_ora_tables"))
setMethod("get_ora_tables", "scDiffCom", function(x) x@ORA)
setGeneric("set_ora_tables", function(x, new_ora_tables) standardGeneric("set_ora_tables"))
setMethod("set_ora_tables", "scDiffCom", function(x, new_ora_tables) {
  x@ORA <- new_ora_tables
  validObject(x)
  return(x)
})

