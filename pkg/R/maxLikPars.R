##' Naive maximum density parameter selection
##' @param ss MCMC chains. mcmc or mcmclist object
##' @return named character vector of parameter values
##'
##' Uses which.max and density function with default parameters
maxLikPars <- function(ss){
  xx <- as.matrix(ss)
  pars.maxp <- rep(NA,ncol(xx))
  names(pars.maxp) <- colnames(xx)
  maxp.res <- rep(NA,ncol(xx))
  names(maxp.res) <- colnames(xx)
  for (n in colnames(xx)){
    den <- density(xx[,n])
    ii <- which.max(den$y)
    pars.maxp[n] <- den$x[ii]
    maxp.res[n] <- mean(diff(den$x))
  }
  return(pars.maxp)
} ##maxLikPars