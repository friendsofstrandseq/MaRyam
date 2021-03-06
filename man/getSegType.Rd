% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDispParAndSegType.R
\name{getSegType}
\alias{getSegType}
\title{Compute the segment type given a cell (majority) type and the segment status and returns a string containing a number of w characters followed by a number of c characters.}
\usage{
getSegType(cellType, status)
}
\arguments{
\item{cellType}{The (majority) type of a cell that can have one these possible values: "ww","cc","wc","cw", or "?".}

\item{status}{A \code{vector} of length 4 containing {CN in hap1, inv CN in hap1, CN in hap2, inv CN in hap2} respectively.}
}
\description{
Compute the segment type given a cell (majority) type and the segment status and returns a string containing a number of w characters followed by a number of c characters.
}
\author{
Maryam Ghareghani
}
