\name{CovInit}
\alias{CovInit}
\title{Alternative sampling method...}
\usage{CovInit(pars, nseq, muX, qcov, bound.handling)}
\description{Alternative initial sampling method
used in examples, may not be useful in practice}
\value{
matrix dim nseq x ndim. parameter wise range [xmin,xmax] assured by handleBounds
When used as input to dream, muX, qcov and bound.handling must be passed as extra parameters}
\arguments{
  \item{pars}{list of parameter vectors}
  \item{nseq}{Number of chains to initialise. scalar}
  \item{muX}{vector of length ndim}
  \item{qcov}{Covariance matrix}
  \item{bound.handling}{one of: "none","reflect","bound","fold","rand"}
}
