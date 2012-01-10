\name{plot.dream}
\alias{plot.dream}
\title{Visualising DREAM results}
\usage{
\S3method{plot}{dream}(x,interactive=TRUE,...)
}
\arguments{
  \item{x}{A \code{\link{dream}} object}
  \item{interactive}{whether to use devAskNewPage (the default) or print
    all at once}
  \item{...}{Placeholder for S3method}
}
\description{
  Plots various visualisation of MCMC chain results
}
\details{
  Extracts the second half of the MCMC chains as a coda object

  Calls:
  \itemize{
    \item coda's \code{\link[=plot.mcmc]{plot}}
    \item coda's \code{\link[=xyplot.mcmc]{xyplot}}
    \item coda's \code{\link[=densityplot.mcmc]{densityplot}}
    \item Plots the distribution of acceptance rates
    \item coda's \code{\link{gelman.plot}} (may fail depending on characteristics of
    MCMC chains)
    \item Plots the evolution of the Gelman R statistic, using previously
    internally compueted results
    \item Plots the multi-variate density of the first chain
  }
}
\examples{
\dontrun{
  ##To run the equivalent coda visualisations:

    ss <- window(dream.obj$Sequences, start = end(dream.obj$Sequences)/2 + 1)
    plot(ss)
    print(xyplot(ss))

    ## Lattice graphics parameter density plots
    ## Each chain is superposed
    print(densityplot(ss))

    ## All chains combined (treated as one)
    print(densityplot(as.mcmc(as.matrix(ss))))

    ## Alternatively, specify burn-in within the call
    print(densityplot(dream.obj$Sequences, start = end(dream.obj$Sequences)/2 + 1))

    ## Base graphics parameter density plots
    plot(ss,trace=F)

    ## Uses densplot behind the scenes, e.g.
    densplot(ss[,1])

  }
}

