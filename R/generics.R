#' #' Perform filtering and over-representation analysis
#' #'
#' #' This function is called internally by \code{run_interaction_analysis} after the permutation tests
#' #'  have been performed. Based on the threshold parameters, it returns detected and differentially expressed
#' #'  CCIs in the slot \code{cci_table_detected} and results of over-representation analysis in \code{cci_ora_table}.
#' #'  The function can be run with new threshold parameters on any scDiffCom object that already contain the slot
#' #'  \code{cci_table_raw}. This allows the user to test various filtering parameters without the need to rerun the
#' #'  potentially time-consuming permutation analysis. When new thresholds are defined the slot \code{parameters} is
#' #'  modified accordingly.
#' #'  Filtering and over-representation are not independent as the second depends on the first. Therefore, when
#' #'  running filtering with new parameters, we also need to update the ORA results.
#' #'
#' #' @param object An object
#' #' @param ... Arguments passed to other methods
#' #'
#' #' @return An object with updated \code{cci_table_detected} and \code{ora_table}
#' #'
#' #' @rdname FilterCCI
#' #' @export FilterCCI
#' #'
#' FilterCCI <- function(object, ...) {
#'   UseMethod(generic = "FilterCCI", object = object)
#' }
#'
#' #' @param object xxx
#' #' @param new_threshold_quantile_score xxx
#' #' @param new_threshold_p_value_specificity xxx
#' #' @param new_threshold_p_value_de xxx
#' #' @param new_threshold_logfc xxx
#' #' @param skip_ora xxx
#' #' @param verbose Should messages be printed?
#' #'
#' #' @rdname FilterCCI
#' #' @export
#' #' @method FilterCCI scDiffCom
#' #'
#' FilterCCI.scDiffCom <- function(
#'   object,
#'   new_threshold_quantile_score = NULL,
#'   new_threshold_p_value_specificity = NULL,
#'   new_threshold_p_value_de = NULL,
#'   new_threshold_logfc = NULL,
#'   skip_ora = FALSE,
#'   verbose = TRUE,
#'   ...
#' ) {
#'   run_filtering_and_ora(
#'     object = object,
#'     new_threshold_quantile_score = new_threshold_quantile_score,
#'     new_threshold_p_value_specificity = new_threshold_p_value_specificity,
#'     new_threshold_p_value_de = new_threshold_p_value_de,
#'     new_threshold_logfc = new_threshold_logfc,
#'     skip_ora = skip_ora,
#'     verbose = verbose,
#'     class_signature = "scDiffCom"
#'   )
#' }
#'
#' #' @param object xxx
#' #' @param new_threshold_quantile_score xxx
#' #' @param new_threshold_p_value_specificity xxx
#' #' @param new_threshold_p_value_de xxx
#' #' @param new_threshold_logfc xxx
#' #' @param skip_ora xxx
#' #' @param verbose Should messages be printed?
#' #'
#' #' @rdname FilterCCI
#' #' @export
#' #' @method FilterCCI scDiffComCombined
#' #'
#' FilterCCI.scDiffComCombined <- function(
#'   object,
#'   new_threshold_quantile_score = NULL,
#'   new_threshold_p_value_specificity = NULL,
#'   new_threshold_p_value_de = NULL,
#'   new_threshold_logfc = NULL,
#'   skip_ora = FALSE,
#'   verbose = TRUE,
#'   ...
#' ) {
#'   run_filtering_and_ora(
#'     object = object,
#'     new_threshold_quantile_score = new_threshold_quantile_score,
#'     new_threshold_p_value_specificity = new_threshold_p_value_specificity,
#'     new_threshold_p_value_de = new_threshold_p_value_de,
#'     new_threshold_logfc = new_threshold_logfc,
#'     skip_ora = skip_ora,
#'     verbose = verbose,
#'     class_signature = "scDiffComCombined"
#'   )
#' }
#'
#'
#' #' Perform over-representation analysis on detected intercellular communication patterns
#' #'
#' #' Over-representation analysis performed when an scDiffCom object contains results
#' #'  obtained from differential expression analysis.
#' #'
#' #' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' #' @param ... Arguments passed to other methods
#' #'
#' #' @return An S4 object of class \code{scDiffCom} with updated slots \code{ora_defaut} (and potentially \code{ora_stringent}).
#' #'
#' #' @rdname RunORA
#' #' @export RunORA
#' #'
#' RunORA <- function(object, ...) {
#'   UseMethod(generic = "RunORA", object = object)
#' }
#'
#'
#' #' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#' #'  for each category. Set to \code{c("ER_CELLTYPES", "LRI", "GO_TERMS")} by default.
#' #' @param overwrite Should existing categories be overwriten in case they match with new categories?
#' #' @param verbose Should messages be printed?
#' #'
#' #' @rdname RunORA
#' #' @export
#' #' @method RunORA scDiffCom
#' RunORA.scDiffCom <- function(
#'   object,
#'   categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
#'   overwrite = TRUE,
#'   #stringent_or_default = "default",
#'   #stringent_logfc_threshold = NULL,
#'   verbose = TRUE,
#'   ...
#' ) {
#'   run_ora(
#'     object = object,
#'     categories = categories,
#'     overwrite = overwrite,
#'     stringent_or_default = "default",
#'     stringent_logfc_threshold = NULL,
#'     verbose = verbose,
#'     class_signature = "scDiffCom",
#'     global = FALSE
#'   )
#' }
#'
#' #' @param categories A character vector specifying the categories on which to perform the analysis. One data.table is returned
#' #'  for each category. Set to \code{c("ER_CELLTYPES", "LRI", "GO_TERMS")} by default.
#' #' @param overwrite Should existing categories be overwritten in case they match with new categories?
#' #' @param verbose Should messages be printed?
#' #'
#' #' @rdname RunORA
#' #' @export
#' #' @method RunORA scDiffComCombined
#' RunORA.scDiffComCombined <- function(
#'   object,
#'   categories = c("ER_CELLTYPES", "LRI", "GO_TERMS", "KEGG_PWS"),
#'   overwrite = TRUE,
#'   #stringent_or_default = "default",
#'   #stringent_logfc_threshold = NULL,
#'   #global = TRUE,
#'   verbose = TRUE,
#'   ...
#' ) {
#'   run_ora(
#'     object = object,
#'     categories = categories,
#'     overwrite = overwrite,
#'     stringent_or_default = "default",
#'     stringent_logfc_threshold = NULL,
#'     verbose = verbose,
#'     class_signature = "scDiffComCombined",
#'     global = FALSE
#'   )
#' }
#'
#' #' xxx
#' #'
#' #' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' #' @param ... Arguments passed to other methods
#' #'
#' #' @return xxx
#' #'
#' #' @rdname PlotORA
#' #' @export PlotORA
#' #'
#' PlotORA <- function(object, ...) {
#'   UseMethod(generic = "PlotORA", object = object)
#' }
#'
#' #' @param object xxx
#' #' @param category xxx
#' #' @param regulation xxx
#' #' @param max_terms_show xxx
#' #' @param OR_threshold xxx
#' #' @param p_value_threshold xxx
#' #'
#' #' @rdname PlotORA
#' #' @export
#' #' @method PlotORA scDiffCom
#' PlotORA.scDiffCom <- function(
#'   object,
#'   category,
#'   regulation,
#'   max_terms_show,
#'   OR_threshold = 1,
#'   p_value_threshold = 0.05,
#'   #stringent = FALSE,
#'   ...
#' ) {
#'   stringent <- FALSE
#'   if (stringent) {
#'     #ora_dt <- get_ora_stringent(object)
#'   } else {
#'     ora_dt <- object@ora_table
#'   }
#'   if (identical(ora_dt, list())) {
#'     stop("No ORA data.table to extract from scDiffCom object")
#'   }
#'   if (!(category %in% names(ora_dt))) {
#'     stop("Can't find the specified ORA category")
#'   }
#'   ora_dt <- ora_dt[[category]]
#'   VALUE_ID <- "VALUE"
#'   if(regulation == "UP") {
#'     OR_ID <- "OR_UP"
#'     p_value_ID <- "BH_P_VALUE_UP"
#'     ORA_SCORE_ID <- "ORA_SCORE_UP"
#'   } else if(regulation == "DOWN") {
#'     OR_ID <- "OR_DOWN"
#'     p_value_ID <- "BH_P_VALUE_DOWN"
#'     ORA_SCORE_ID <- "ORA_SCORE_DOWN"
#'
#'   } else if(regulation == "FLAT") {
#'     OR_ID <- "OR_FLAT"
#'     p_value_ID <- "BH_P_VALUE_FLAT"
#'     ORA_SCORE_ID <- "ORA_SCORE_FLAT"
#'   } else {
#'     stop("Can't find `regulation` type")
#'   }
#'   ora_dt <- ora_dt[get(OR_ID) > OR_threshold & get(p_value_ID) <= p_value_threshold][order(-get(ORA_SCORE_ID))]
#'   if (nrow(ora_dt) == 0) {
#'     return("No significant ORA results for the selected parameters.")
#'   }
#'   n_row_tokeep <- min(max_terms_show, nrow(ora_dt))
#'   ora_dt <- ora_dt[1:n_row_tokeep]
#'   ggplot2::ggplot(ora_dt, aes(get(ORA_SCORE_ID), stats::reorder(get(VALUE_ID), get(ORA_SCORE_ID)))) +
#'     geom_point(aes(size = -log10(get(p_value_ID)), color = log2(get(OR_ID)))) +
#'     scale_color_gradient(low = "orange", high = "red") +
#'     xlab("ORA score") +
#'     ylab(category) +
#'     labs(size = "-log10(Adj. P-Value)", color = "log2(Odds Ratio)") +
#'     theme(text = element_text(size = 16))
#' }
#'
#' #' @param object xxx
#' #' @param subID xxx
#' #' @param category xxx
#' #' @param regulation xxx
#' #' @param max_terms_show xxx
#' #' @param OR_threshold xxx
#' #' @param p_value_threshold xxx
#' #' @param global xxx
#' #'
#' #' @rdname PlotORA
#' #' @export
#' #' @method PlotORA scDiffComCombined
#' PlotORA.scDiffComCombined <- function(
#'   object,
#'   subID,
#'   category,
#'   regulation,
#'   max_terms_show,
#'   global = FALSE,
#'   OR_threshold = 1,
#'   p_value_threshold = 0.05,
#'   #stringent = FALSE,
#'   ...
#' ) {
#'   ID <- NULL
#'   stringent <- FALSE
#'   if (stringent) {
#'     if (global) {
#'       ora_dt <- object@ora_combined_stringent
#'     } else {
#'       ora_dt <- object@ora_stringent
#'     }
#'   } else {
#'     if (global) {
#'       ora_dt <- object@ora_combined_default
#'     } else {
#'       ora_dt <- object@ora_table
#'     }
#'   }
#'   if (identical(ora_dt, list())) {
#'     stop("No ORA data.table to exctract from scDiffCom object")
#'   }
#'   if (!(category %in% names(ora_dt))) {
#'     stop("Can't find the specified ORA category")
#'   }
#'   ora_dt <- ora_dt[[category]]
#'   if (!global) {
#'     if (!(subID %in% unique(ora_dt$ID))) {
#'       stop("`subID` must be present in the scDiffComCombined object")
#'     }
#'     ora_dt <- ora_dt[ID == subID]
#'   }
#'
#'   VALUE_ID <- "VALUE"
#'   if(regulation == "UP") {
#'     OR_ID <- "OR_UP"
#'     p_value_ID <- "BH_P_VALUE_UP"
#'     ORA_SCORE_ID <- "ORA_SCORE_UP"
#'   } else if(regulation == "DOWN") {
#'     OR_ID <- "OR_DOWN"
#'     p_value_ID <- "BH_P_VALUE_DOWN"
#'     ORA_SCORE_ID <- "ORA_SCORE_DOWN"
#'
#'   } else if(regulation == "FLAT") {
#'     OR_ID <- "OR_FLAT"
#'     p_value_ID <- "BH_P_VALUE_FLAT"
#'     ORA_SCORE_ID <- "ORA_SCORE_FLAT"
#'   } else {
#'     stop("Can't find `regulation` type")
#'   }
#'   ora_dt <- ora_dt[get(OR_ID) > OR_threshold & get(p_value_ID) <= p_value_threshold][order(-get(ORA_SCORE_ID))]
#'   if (nrow(ora_dt) == 0) {
#'     return("No significant ORA results for the selected parameters.")
#'   }
#'   n_row_tokeep <- min(max_terms_show, nrow(ora_dt))
#'   ora_dt <- ora_dt[1:n_row_tokeep]
#'   ggplot2::ggplot(ora_dt, aes(get(ORA_SCORE_ID), stats::reorder(get(VALUE_ID), get(ORA_SCORE_ID)))) +
#'     geom_point(aes(size = -log10(get(p_value_ID)), color = log2(get(OR_ID)))) +
#'     scale_color_gradient(low = "orange", high = "red") +
#'     xlab("ORA score") +
#'     ylab(category) +
#'     labs(size = "-log10(Adj. P-Value)", color = "log2(Odds Ratio)") +
#'     theme(text = element_text(size = 16))
#' }
#'
#' #' xxx
#' #'
#' #' xxx
#' #'
#' #' @param object An scDiffCom object previously returned by \code{run_interaction_analysis}.
#' #' @param ... xxx
#' #'
#' #' @return xxx
#' #'
#' #' @rdname BuildNetwork
#' #' @export BuildNetwork
#' #'
#' BuildNetwork <- function(object, ...) {
#'   UseMethod(generic = "BuildNetwork", object = object)
#' }
#'
#' #' @param object xxx
#' #' @param network_type xxx
#' #' @param layout_type xxx
#' #' @param abbreviation_table xxx
#' #'
#' #' @rdname BuildNetwork
#' #' @export
#' #' @method BuildNetwork scDiffCom
#' BuildNetwork.scDiffCom = function(
#'   object,
#'   network_type = c(
#'     "condition1_network",
#'     "condition2_network",
#'     "difference_network",
#'     "up_regulated_network",
#'     "down_regulated_network",
#'     "ORA_network"
#'   ),
#'   layout_type = c(
#'     "conventional",
#'     "bipartite"
#'   ),
#'   abbreviation_table = NULL,
#'   #LRIs = NULL,
#'   ...
#' ) {
#'   network_type <- match.arg(network_type)
#'   layout_type <- match.arg(layout_type)
#'   return(
#'     build_interactive_network(
#'       object = object,
#'       network_type = network_type,
#'       layout_type = layout_type,
#'       class_signature = "scDiffCom",
#'       subobject_name = NULL,
#'       abbreviation_table = abbreviation_table#,
#'       #LRIs = LRIs
#'     )
#'   )
#' }
#'
#' #' @param object xxx
#' #' @param network_type xxx
#' #' @param layout_type xxx
#' #' @param ID xxx
#' #'
#' #' @rdname BuildNetwork
#' #' @export
#' #' @method BuildNetwork scDiffComCombined
#' BuildNetwork.scDiffComCombined = function(
#'   object,
#'   network_type = c(
#'     "condition1_network",
#'     "condition2_network",
#'     "difference_network",
#'     "up_regulated_network",
#'     "down_regulated_network",
#'     "ORA_network"
#'   ),
#'   layout_type = c(
#'     "conventional",
#'     "bipartite"
#'   ),
#'   ID,
#'   abbreviation_table,
#'   #LRIs = NULL,
#'   ...
#' ) {
#'   network_type <- match.arg(network_type)
#'   layout_type <- match.arg(layout_type)
#'   if (!(ID %in% unique(object@cci_table_detected$ID))) {
#'     stop("`ID` must be present in the scDiffComCombined object")
#'   }
#'   return(
#'     build_interactive_network(
#'       object = object,
#'       network_type = network_type,
#'       layout_type = layout_type,
#'       class_signature = "scDiffComCombined",
#'       subobject_name = ID,
#'       abbreviation_table = abbreviation_table#,
#'       #LRIs = LRIs
#'     )
#'   )
#' }
#'
