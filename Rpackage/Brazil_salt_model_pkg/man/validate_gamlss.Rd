% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{validate_gamlss}
\alias{validate_gamlss}
\title{Stochastic prediction from a gamlss object}
\usage{
validate_gamlss(dt, gamlss_obj, mc = 10L, orig_data = dt)
}
\description{
`validate_gamlss` returns a data.table with the observed and predicted
 variable. If \code{mc > 1} multiple predictions are drawn from the predicted
 distributions. Useful for plotting with ggplot
}
