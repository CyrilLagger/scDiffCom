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
    ## add various checks here
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


