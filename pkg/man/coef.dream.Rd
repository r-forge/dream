\name{coef.dream}
\alias{coef.dream}
\alias{coef.dream_model}
\alias{maxLikCoda}
\title{Extract parameter values from dream or dream_model object}
\usage{
\method{coef}{dream}(object, method = c("sample.ml","uni.mode", "mean",
"median"), \dots)
maxLikCoda(x)
}
\description{Extract parameter values using a choice of methods (or an
  arbitrary function)}
\value{named vector of parameter values}
\arguments{
  \item{x}{for maxLikCoda, usually mcmc or mcmc.list object, or any representation of
    MCMC chains that can be converted with \code{as.matrix}}
  \item{object}{dream object}
  \item{method}{method for extracting a parameter set from the MCMC
    chains. One of:
    \describe{
      \item{\code{"uni.mode"}}{ using maxLikCoda,
        mode of the univariate density estimate
        for each parameter, using settings as in
        \code{\link{densityplot.mcmc}}.
	May not find the optimal parameter combination of the
	distribution is multi-modal. 
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
  \item{...}{Passed to \code{\link{window.dream}}}
}
\details{
  maxLikCoda re-uses code from \code{\link{densityplot.mcmc}}
  
  e.g. of using arbitrary function for method: 20\% quantile.
  \code{
    coef(object,method=function(sss) apply(as.matrix(sss),2,quantile,0.2)
  }
}