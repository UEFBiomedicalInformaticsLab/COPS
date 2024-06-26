% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{PIK}
\alias{PIK}
\alias{PIK_KEGG}
\alias{PIK_from_networks}
\alias{PIK_GNGS}
\title{Pathway induced kernel (Manica et al. 2019)}
\usage{
PIK(x, L, rwr_smoothing = FALSE, rwr_restart_prob = 0.75, ...)

PIK_KEGG(x, gene_key = "SYMBOL", ...)

PIK_from_networks(x, networks, normalized_laplacian = TRUE, parallel = 1, ...)

PIK_GNGS(x, gene_network, gene_sets, ...)
}
\arguments{
\item{x}{gene feature matrix}

\item{L}{pathway graph Laplacian matrix}

\item{rwr_smoothing}{apply feature smoothing}

\item{rwr_restart_prob}{restart probability for smoothing}

\item{...}{passed on to \code{\link{PIK}}}

\item{gene_key}{column in \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} that KEGG IDs should be translated to}

\item{networks}{list of igraph objects corresponding to pathway graphs}

\item{normalized_laplacian}{normalize Laplacian for calculations}

\item{parallel}{number of parallel threads}

\item{gene_network}{\code{igraph} object}

\item{gene_sets}{\code{list} of gene sets}
}
\value{
kernel \code{matrix}

\code{list} of KEGG-based pathway-kernel matrices

list of kernel matrices

list of kernel matrices
}
\description{
Pathway induced kernel (Manica et al. 2019)
}
\section{Functions}{
\itemize{
\item \code{PIK_KEGG()}: Use KEGG pathways to form pathway induced kernels

\item \code{PIK_from_networks()}: Pathway induced kernel from networks

\item \code{PIK_GNGS()}: Extract subnetworks based on pathway genesets and compute PIKs

}}
