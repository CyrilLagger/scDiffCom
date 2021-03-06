% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{BuildNetwork}
\alias{BuildNetwork}
\alias{BuildNetwork,scDiffCom-method}
\title{Display cell-type to cell-type interactive networks}
\usage{
BuildNetwork(
  object,
  network_type = c("ORA_network"),
  layout_type = c("bipartite", "conventional"),
  abbreviation_table = NULL
)

\S4method{BuildNetwork}{scDiffCom}(
  object,
  network_type = c("ORA_network"),
  layout_type = c("bipartite", "conventional"),
  abbreviation_table = NULL
)
}
\arguments{
\item{object}{\code{scDiffCom} object}

\item{network_type}{Type of network to display. Currently, only
\code{ORA_network} (default) is supported.}

\item{layout_type}{Layout of the network to display. Can either be
\code{"bipartite"} (default) or \code{"conventional"}.}

\item{abbreviation_table}{Table with abbreviations
for the cell types present in the \code{object}. If \code{NULL} (default),
full names of the cell-types are displayed. Otherwise, it must be a
data.frame or data.table with exactly two columns with names
\code{ORIGINAL_CELLTYPE} and \code{ABBR_CELLTYPE}.}
}
\value{
A visNetwork object.
}
\description{
Create and plot an interactive network that summarize how
cell-types and their interactions are over-represented.
}
