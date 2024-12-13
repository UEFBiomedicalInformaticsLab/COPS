% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{mkkm_mr}
\alias{mkkm_mr}
\title{Multiple kernel k-means with matrix induced regularization (Liu et al. 2016)}
\usage{
mkkm_mr(
  K_list,
  k,
  lambda,
  tolerance = 1e-06,
  parallel = 0,
  use_mosek = FALSE,
  mosek_verbosity = 0L,
  mkkm_mr_maxiter = 10,
  no_stop = FALSE
)
}
\arguments{
\item{K_list}{\code{list} of kernel matrices}

\item{k}{number of clusters}

\item{lambda}{\code{numeric} regularization parameter}

\item{tolerance}{\code{numeric} stopping criterion value}

\item{parallel}{number of parallel threads used by the quadratic solver}

\item{use_mosek}{If \code{TRUE}, the optimization will be run with 
\code{Rmosek} instead of \code{CVXR}.}

\item{mosek_verbosity}{MOSEK logging parameter.}

\item{mkkm_mr_maxiter}{maximum number of iterations, usually the algorithm 
converges soon after 2 iterations}

\item{no_stop}{If \code{TRUE}, always runs \code{mkkm_mr_maxiter} iterations}
}
\value{
a kernel \code{matrix}
}
\description{
Multiple kernel k-means with matrix induced regularization (Liu et al. 2016)
}
