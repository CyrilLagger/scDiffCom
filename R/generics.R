#' @include objects.R
#'
NULL

#' Return the slot 'parameters'
#'
#' @param x scDiffCom object
#'
#' @return
#' @export
setGeneric("parameters", function(x) standardGeneric("parameters"))

#' @rdname parameters
setMethod("parameters", "scDiffCom", function(x) x@parameters)

#' Title
#'
#' @param x x
#' @param value x
#' @name get_param
#'
#' @return x
#' @export
setGeneric("parameters<-", function(x, value) standardGeneric("parameters<-"))

#' @rdname get_param
setMethod("parameters<-", "scDiffCom", function(x, value) {
  x@parameters <- value
  validObject(x)
  return(x)
})

#' Title
#'
#' @param x x
#'
#' @return x
#' @export
setGeneric("get_cci_table_raw", function(x) standardGeneric("get_cci_table_raw"))

#' @rdname get_cci_table_raw
setMethod("get_cci_table_raw", "scDiffCom", function(x) x@cci_table_raw)

#' Title
#'
#' @param x x
#' @param new_cci_table_raw x
#'
#' @return x
#' @export

setGeneric("set_cci_table_raw", function(x, new_cci_table_raw) standardGeneric("set_cci_table_raw"))

#' @rdname set_cci_table_raw
setMethod("set_cci_table_raw", "scDiffCom", function(x, new_cci_table_raw) {
  x@cci_table_raw <- new_cci_table_raw
  validObject(x)
  return(x)
})

#' Title
#'
#' @param x x
#'
#' @return x
#' @export
setGeneric("get_cci_table_filtered", function(x) standardGeneric("get_cci_table_filtered"))

#' @rdname get_cci_table_filtered
setMethod("get_cci_table_filtered", "scDiffCom", function(x) x@cci_table_filtered)

#' Title
#'
#' @param x x
#' @param new_cci_table_filtered x
#'
#' @return x
#' @export
setGeneric("set_cci_table_filtered", function(x, new_cci_table_filtered) standardGeneric("set_cci_table_filtered"))

#' @rdname set_cci_table_filtered
setMethod("set_cci_table_filtered", "scDiffCom", function(x, new_cci_table_filtered) {
  x@cci_table_filtered <- new_cci_table_filtered
  validObject(x)
  return(x)
})

#' Title
#'
#' @param x x
#'
#' @return x
#' @export
setGeneric("distributions", function(x) standardGeneric("distributions"))

#' @rdname distributions
setMethod("distributions", "scDiffCom", function(x) x@distributions)

#' Title
#'
#' @param x x
#' @param value x
#' @name get_distr
#'
#' @return
#' @export
setGeneric("distributions<-", function(x, value) standardGeneric("distributions<-"))

#' @rdname get_distr
setMethod("distributions<-", "scDiffCom", function(x, value) {
  x@distributions <- value
  validObject(x)
  return(x)
})

#' Title
#'
#' @param x x
#'
#' @return x
#' @export
setGeneric("get_ora_tables", function(x) standardGeneric("get_ora_tables"))

#' @rdname get_ora_tables
setMethod("get_ora_tables", "scDiffCom", function(x) x@ora_tables)

#' Title
#'
#' @param x x
#' @param new_ora_tables x
#'
#' @return x
#' @export
setGeneric("set_ora_tables", function(x, new_ora_tables) standardGeneric("set_ora_tables"))

#' @rdname set_ora_tables
setMethod("set_ora_tables", "scDiffCom", function(x, new_ora_tables) {
  x@ora_tables <- new_ora_tables
  validObject(x)
  return(x)
})

