% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{EraseRawCCI}
\alias{EraseRawCCI}
\alias{EraseRawCCI,scDiffCom-method}
\title{Create a copy of a scDiffCom object without \code{cci_table_raw}}
\usage{
EraseRawCCI(object)

\S4method{EraseRawCCI}{scDiffCom}(object)
}
\arguments{
\item{object}{\code{scDiffCom} object}
}
\value{
A scDiffCom object with an empty list for \code{cci_table_raw}.
}
\description{
This function will replace \code{cci_table_raw} by an empty list. Useful to
save space for large datasets. However, after this operation,
no filtering can be re-run on the new object, meaning that obtaining
results for different filtering parameters will require the perform the full
analysis from scratch.
}
