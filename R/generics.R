#' Perform filtering and over-representation analysis
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
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return An object with updated \code{cci_detected} and \code{ora_default}
#'
#' @rdname FilterCCI
#' @export FilterCCI
#'
FilterCCI <- function(object, ...) {
  UseMethod(generic = "FilterCCI", object = object)
}

#' Perform over-representation analysis on detected intercellular communication patterns
#'
#' Over-representation analysis performed when an scDiffCom object contains results
#'  obtained from differential expression analysis.
#'
#' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' @param ... Arguments passed to other methods
#'
#' @return An S4 object of class \code{scDiffCom} with updated slots \code{ora_defaut} (and potentially \code{ora_stringent}).
#'
#' @rdname RunORA
#' @export RunORA
#'
RunORA <- function(object, ...) {
  UseMethod(generic = "RunORA", object = object)
}

PlotORA <- function(object, ...) {
  UseMethod(generic = "PlotORA", object = object)
}

#' xxx
#'
#' xxx
#'
#' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' @param ... xxx
#'
#' @return xxx
#'
#' @rdname BuildNetwork
#' @export BuildNetwork
#'
BuildNetwork <- function(object, ...) {
  UseMethod(generic = "BuildNetwork", object = object)
}
