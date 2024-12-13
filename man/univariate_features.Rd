% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/characterization.R
\name{univariate_features}
\alias{univariate_features}
\title{Univariate group difference tests for feature selection}
\usage{
univariate_features(dat, group, remove_zero_var = TRUE, parallel = 1)
}
\arguments{
\item{dat}{molecular data, samples on columns, rows must be named}

\item{group}{categorical variable matching \code{dat} columns}

\item{remove_zero_var}{whether to remove zero variance genes}

\item{parallel}{number of threads}
}
\value{
\code{list} of test results
}
\description{
Computes various statistics for molecular differences between groups.
}
