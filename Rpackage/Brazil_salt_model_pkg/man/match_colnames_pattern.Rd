% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{match_colnames_pattern}
\alias{match_colnames_pattern}
\title{Obtain matching names corresponding to patterns}
\usage{
match_colnames_pattern(dt, ...)
}
\description{
`match_colnames_pattern` returns the matching names of the argument `dt`
(i.e. \code{names(dt)}) corresponding to the regular expression patterns
provided. The patterns must be supported by \code{\link{grep}}.
This is based on `data.table:::patterns`
}
