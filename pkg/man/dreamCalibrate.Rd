\name{dreamCalibrate}
\alias{dreamCalibrate}
\title{Utility to calibrate a function using dream}
\usage{

dreamCalibrate(FUN, pars, obs, lik.fun=calc.loglik,lik.control=NULL,FUN.pars = list(), ...)

}
\description{
  Calibrate a function using \code{\link{dream}}, a specified
  likelihood function \code{lik.fun} and observed values \code{obs}
}
\arguments{
  \item{FUN}{
    model function with first argument a vector of parameter values of
    length ndim.
  }
  \item{pars}{
    a list of variable ranges. Any names will be propagated to output.
  }
  \item{obs}{
    a numeric vector of observed values, corresponding to the output of \code{FUN}
  }
  \item{lik.fun}{
    A function that returns the log likelihood of model predictions
    matching observed values. \code{log.lik=f(predicted,observed,control)} 
  }
  \item{lik.control}{
    A list of any extra arguments to be passed to \code{lik.fun}
  }
  \item{FUN.pars}{
    A list of any extra arguments to be passed to \code{FUN}.
  }
  \item{...}{
    Extra arguments to be passed to dream, e.g. \code{control}
  }
}
\details{

  There are S3 methods for:
  \code{\link[=predict.dream_model]{predict}},
  \code{\link[=coef.dream]{coef}}.

}
\value{
  An object inheriting from \code{\link{dream}}, i.e. with the same elements and:
  \item{FUN}{The function calibrated}
  \item{FUN.pars}{The extra arguments originally passed to that function}
}

\seealso{
  See \code{\link{dream}} for details on the calibration method,
  visualisation of its results and diagnostics.

  Example in demo folder:
  \itemize{
    \item{FME non linear model: }{Calibrating the non-linear model shown
      in the FME vignette}
  }
}
