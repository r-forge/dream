\name{simulate.dream}
\alias{simulate.dream}
\title{Simulate values from a distribution using converged MCMC chains
  from a DREAM object}
\usage{\S3method{simulate}{dream}(object,nsim=1000,seed=NULL,\dots)}
\arguments{
  \item{object}{\code{\link{dream}} object}
  \item{nsim}{Approximate number of function evaluations (length of
    sequence will be lower due to thinning)}
  \item{seed}{passed to \code{\link{set.seed}} before continuing}
  \item{...}{ignored, for compatibility with S3 generic}
}
\description{
  Continues a converged MCMC chain to provide an additional
  sample.

  Outlier detection, reporting and convergence diagnostics are
  turned off.

  If no thinning was previously used, the sample is thinned with \code{t=10}. This value can be overridden by specifying a value for \code{thin.t} in the \code{control} parameter to \code{\link{dream}}.
}
\value{
  a dream object with approximately the requested number of function evaluations
}
