% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{subsample_clustering_evaluation}
\alias{subsample_clustering_evaluation}
\title{Clustering analysis on cross-validated data sets}
\usage{
subsample_clustering_evaluation(
  dat_embedded,
  parallel = 1,
  by = c("datname", "drname", "run", "fold"),
  silhouette_dissimilarity = NULL,
  dat_list = NULL,
  ...
)
}
\arguments{
\item{dat_embedded}{list of \code{data.frame}s}

\item{parallel}{number of threads}

\item{by}{variables to split input data by}

\item{silhouette_dissimilarity}{dissimilarity matrix used for silhouette evaluation}

\item{dat_list}{list of input data matrices used for calculating clustering indices}

\item{...}{extra arguments are passed through to \code{\link{clustering_dissimilarity_from_data}}, \code{\link{clustering_analysis}} and \code{\link{clustering_metrics}}}
}
\value{
Returns a \code{list} of \code{data.frames} containing \code{\link{clustering_analysis}} and \code{\link{clustering_metrics}} outputs for every
        combination of CV run, CV fold, clustering method, number of clusters as well as all combinations of
        data sets and dimensionality reduction techniques found in the input \code{data.frame}.
}
\description{
Performs clustering analysis on each fold of an external cross validation.
}
\details{
Produces clusterings using multiple methods and settings while computing internal validation metrics 
such as Connectivity, Dunn and Silhouette scores. Also computes chi-squared tests with respect to 
a batch label if one is provided.
}
