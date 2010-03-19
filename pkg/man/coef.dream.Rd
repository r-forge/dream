\name{coef.dream}
\alias{coef.dream}
\title{Extract parameter values from dream object}
\usage{
\S3method{coef}{dream}(object, last.prop=0.5, use.thinned=FALSE, method=c("maxLik", "mean", "median"),\dots)
}
\description{Extract parameter values using a choice of methods (or an
  arbitrary function)}
\value{named vector of parameter values}
\arguments{
  \item{object}{dream object}
  \item{last.prop}{proportion of total sequence to use (0,1]
    if 1, use whole sequence}
  \item{use.thinned}{Use thinned MCMC chains}
  \item{method}{either one of maxLik,mean,median or a function parameter
    vector=f(dream object)}
  \item{...}{Unused. Matches generic}
}
\details{
  e.g. of using arbitrary function for method: 20\% quantile.
  \code{
    coef(object,method=function(sss) apply(as.matrix(sss),2,quantile,0.2)
  }
}