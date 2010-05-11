\name{coef.dream}
\alias{coef.dream}
\title{Extract parameter values from dream object}
\usage{
\method{coef}{dream}(object, last.prop = 0.5, use.thinned = FALSE,
     method = c("uni.mode", "mean", "median", "sample.ml"), \dots)
}
\description{Extract parameter values using a choice of methods (or an
  arbitrary function)}
\value{named vector of parameter values}
\arguments{
  \item{object}{dream object}
  \item{last.prop}{proportion of total sequence to use (0,1]
    if 1, use whole sequence}
  \item{use.thinned}{Use thinned MCMC chains}
  \item{method}{method for extracting a parameter set from the MCMC
    chains. One of:
    \describe{
      \item{\code{"uni.mode"}}{
        mode of the univariate density estimate
        for each parameter, using settings as in
        \code{\link{densityplot.mcmc}}.
      }
      \item{\code{"mean"}}{
        mean of each univariate parameter distribution.
      }
      \item{\code{"median"}}{
        median of each univariate parameter distribution.
      }
      \item{\code{"sample.ml"}}{
        parameter set with maximum likelihood (according to the chosen
        likelihood function) from the generated MCMC chains. 
      }
      \item{\code{function(object)}}{
        a function of the dream object which returns a parameter vector.
      }
    }
  }
  \item{...}{Unused. Matches generic}
}
\details{
  e.g. of using arbitrary function for method: 20\% quantile.
  \code{
    coef(object,method=function(sss) apply(as.matrix(sss),2,quantile,0.2)
  }
}