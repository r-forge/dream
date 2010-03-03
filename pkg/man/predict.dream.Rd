\name{predict.dream}
\alias{predict.dream}
\title{Predict values using dream object}
\usage{
\S3method{predict}{dream}(object, newdata=NULL,method="maxLik", level=0.99, last.prop=0.5, use.thinned=TRUE,\dots)
}
\description{
Predict values using function calibrated by dream, optionally with new data,
using various methods of summarising the posterior parameter and output distributions}
\value{data frame with each column corresponding to a returned vector}
\arguments{
  \item{object}{dream object}
  \item{newdata}{new FUN.pars list. If NULL, use object's.}
  \item{method}{CI or a \code{method} of \code{\link{coef.dream}}}
  \item{level}{Requested two-sided level of confidence. For CI method.}
  \item{last.prop}{Proportion of MCMC chains to keep}
  \item{use.thinned}{Whether to use existing thinned chains}
  \item{...}{Unused - matches generic}
}
