##' Computes the density of each x value
##'
##' @param x matrix nseq x ndim
##' @param control. list containing gamma,Wb,Cb
##' @param FUN - the model to run
##'   R function with first argument a vector of length ndim.
##'   returning a single scalar corresponding to one of the options:
##' @param option unambiguous abbrev of:
##'   posterior.density, calc.loglik, calc.rmse, logposterior.density,calc.weighted.rmse
##' @param measurement list containing TODO: not sure
##'   data: vector of observations corresponding to model output
##'   sigma: scalar
##' @param ... additional arguments to FUN
##' @return ... list with components
##'   p vector of length nseq
##'   logp vector of length nseq
##
## TODO: p may be erroneously equal to logp?
CompDensity <- function(x,control,FUN,option,
                        measurement,...){

  ## dimensions:
  ##  i. iter 1:nseq
  ##  modpred. scalar or vector commensurate to measurement$data
  ##  err. vector of same length as modpred
  ##  SSR scalar

  p <- rep(NA,nseq)
  logp <- rep(NA,nseq)
  
  ## Sequential evaluation
  for (ii in 1:nrow(x)){
    ## Call model to generate simulated data
    ## TODO: correct use of optional pars?
    modpred <- FUN(x[ii,],...)

    switch(option,
           ## Model directly computes posterior density
           posterior.density={
             p[ii] <- modpred
             logp[ii] <- log(modpred)
           },
           ## Model computes output simulation           
           calc.loglik={
             err <- measurement$data-modpred
                 
             ## Compute the number of measurement data
             N <- length(measurement$data)
             
             ## Derive the log likelihood
             logp[ii] <- N*log(control$Wb/measurement$sigma)-
               control$Cb*(sum((abs(err/measurement$sigma))^(2/(1+control$gamma))))
             ## And retain in memory
             p[ii] <- logp[ii]
           },
           ## Model computes output simulation
           calc.rmse={
             err <- measurement$data-modpred
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p[ii] <- -SSR
             logp[ii] <- -0.5*SSR
             
           },
           ## Model directly computes log posterior density
           logposterior.density={
             p[ii] <- modpred
             logp[ii] <- modpred
           },
           ## Similar as 3, but now weights with the Measurement Sigma
           ## TODO: appears to be no difference to calc.rmse
           calc.weighted.rmse={
             ## Define the error
             err <- measurement$data-modpred
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p[ii] <- -SSR
             logp[ii] <- -0.5*SSR
           }) ##switch
  }           ## for rows
  return(list(p=p,logp=logp))
} ##CompDensity
