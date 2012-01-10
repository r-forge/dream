\name{plotMCMCQQ}
\alias{plotMCMCQQ}
\title{
  QQ plot of MCMC chains
}
\description{
  Compare two MCMC objects using their QQ plot
}
\usage{
plotMCMCQQ(mat.m, R.m)
}
\arguments{
  \item{mat.m}{
    matrix of MCMC chain parameter values
}
\item{R.m}{
  matrix of MCMC chain parameter values
}
}
\value{
  lattice \code{\link{xyplot}} showing quantile-quantile plots for each parameter
}
