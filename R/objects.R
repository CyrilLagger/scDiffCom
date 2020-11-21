#' @include generics.R
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
  name = "list_or_dt",
  members = c("list", "data.table")
)

#' The scDiffCom Class
#'
#' The scDiffCom Class is the class that stores the results of the
#' cell-cell interaction (CCI) analysis. Each object of this class is intended
#' to store the results of a single "tissue".
#'
#' @slot parameters List of the parameters used to perform the CCI analysis.
#' @slot cci_table_raw Data.table of all possible CCIs inferred from the single-cell data and Ligand-Receptor database.
#' @slot cci_table_filtered Data.table of all detected CCIs obtained from the raw data.table after filtering.
#' @slot distributions A list of matrices with the distributions returned by the permutation test(s).
#' @slot ora_tables A list of data.tables with the results of the over-representation analysis for each specified category.
#'
#' @name scDiffCom-class
#' @export scDiffCom
#'
scDiffCom <- setClass(
  Class = "scDiffCom",
  slots = c(
    parameters = "list",
    cci_table_raw = "list_or_dt",
    cci_table_filtered = "list_or_dt",
    distributions = "list",
    ora_tables = "list"
  ),
  prototype = list(
    parameters = list(),
    cci_table_raw = list(),
    cci_table_filtered = list(),
    distributions = list(),
    ora_tables = list()
  )
)

scDiffCom <- function(
                      parameters,
                      cci_table_raw,
                      cci_table_filtered,
                      distributions,
                      ora_tables) {
  new(
    "scDiffCom",
    parameters = parameters,
    cci_table_raw = cci_table_raw,
    cci_table_filtered = cci_table_filtered,
    distributions = distributions,
    ora_tables = ora_tables
  )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Validity functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity(
  Class = "scDiffCom",
  method = function(
                    object) {
    validity_check <- c(
      validity_parameters(
        parameters = object@parameters
      ),
      validity_cci_table_raw(
        parameters = object@parameters,
        cci_table_raw = object@cci_table_raw
      ),
      validity_cci_table_filtered(
        parameters = object@parameters,
        cci_table_filtered = object@cci_table_filtered
      ),
      validity_distributions(
        parameters = object@parameters,
        distributions = object@distributions
      ),
      validity_ora_tables(
        parameters = object@parameters,
        ora_tables = object@cci_ora_tables
      )
    )
    if (is.null(validity_check)) {
      TRUE
    } else {
      paste0(validity_check, collapse = " AND ")
    }
  }
)

validity_parameters <- function(
                                parameters) {
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
  parameters_actual <- sapply(parameters, typeof)
  if (identical(parameters_expected, parameters_actual)) {
    NULL
  } else {
    "@parameters is not formatted the correct way"
  }
}

validity_cci_table_raw <- function(
                                   parameters,
                                   cci_table_raw) {
  NULL
}

validity_cci_table_filtered <- function(
                                        parameters,
                                        cci_table_filtered) {
  NULL
}

validity_distributions <- function(
                                   parameters,
                                   distributions) {
  NULL
}

validity_ora_tables <- function(
                                parameters,
                                ora_tables) {
  NULL
}
