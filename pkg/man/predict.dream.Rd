\name{predict.dream}
\alias{predict.dream}
\title{Predict values using dream object}
\usage{
\S3method{predict}{dream}(object, newdata=NULL,newFUN=NULL,method="maxLik", level=0.99, last.prop=0.5, use.thinned=TRUE,\dots)
}
\description{
Predict values using function calibrated by dream, optionally with new data,
using various methods of summarising the posterior parameter and output distributions}
\value{whatever newFUN returns (either numeric, ts or list). For CI,
  either a matrix with upper and lower bound or list of matrices.

  N.B. newFUN can return a list, though FUN for use in dream must return
a numeric}
\arguments{
  \item{object}{dream object}
  \item{newdata}{new FUN.pars list. If NULL, use object's.}
  \item{newFUN}{A new function run, for the case where FUN returned a posterior density, but the model result is now required. newFUN must have exactly the same parameters as FUN. If NULL, don't use a different function.}
  \item{method}{CI or a \code{method} of \code{\link{coef.dream}}}
  \item{level}{Requested two-sided level of confidence. For CI method.}
  \item{last.prop}{Proportion of MCMC chains to keep}
  \item{use.thinned}{Whether to use existing thinned chains}
  \item{...}{Unused - matches generic}
}
