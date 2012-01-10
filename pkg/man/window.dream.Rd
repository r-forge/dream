\name{window.dream}
\alias{window.dream}
\title{Extract MCMC chains from a DREAM object}
\usage{
\method{window}{dream}(x, start = 1 + (end(x$Sequences) - 1) *
(1 - fraction),fraction = 0.5,\dots)
}
\arguments{
  \item{x}{dream object}
  \item{start}{the first iteration of interest}
  \item{fraction}{The fraction of the MCMC chains to keep}
  \item{\dots}{extra arguments for \code{\link{window.mcmc}}}
  }
\value{\code{\link{mcmc.list}} object}
\description{Extract part or all of the MCMC chains from a DREAM object,
specifying a burn-in period and allowing thinning}