% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{pareto_fronts}
\alias{pareto_fronts}
\alias{get_metric_comparators}
\title{Pareto based multi-objective optimization}
\usage{
pareto_fronts(
  scores,
  metrics,
  metric_comparators = get_metric_comparators(metrics)
)

get_metric_comparators(metrics)
}
\arguments{
\item{scores}{\code{data.frame} of scores}

\item{metrics}{column names in \code{scores} to plot}

\item{metric_comparators}{list of comparison functions for each metric, 
should be \code{.Primitive(">")} for maximization and \code{.Primitive("<")} 
for minimization.}
}
\value{
integer vector of Pareto front number for each row of \code{scores}
}
\description{
Pareto based multi-objective optimization
}
\section{Functions}{
\itemize{
\item \code{get_metric_comparators()}: Default metric comparators.

}}
