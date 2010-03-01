##' Computes the density of each pars value
##'
##' @param pars matrix nseq x ndim
##' @param control. list containing gamma,Wb,Cb,nseq
##' @param FUN - the model to run
##'   R function with first argument a vector of length ndim.
##'   returns a scalar or vector corresponding to one of the options below.
##' @param func.type Type of function output. One of:
##'   posterior.density, logposterior.density,
##'   calc.loglik. requires optional parameter measurement with elements data & sigma
##'   calc.rmse, calc.weighted.rmse.  requires measurement$data
##' @param measurement list containing TODO: not sure
##'   data: vector of observations corresponding to model output
##'   sigma: scalar
##' @param ... additional arguments to FUN
##' @return list with elements
##'   p vector of length nseq
##'   logp vector of length nseq
##
## TODO: p may be erroneously equal to logp?
## TODO: more appropriate naming of options?
## TODO: allow shortenings of option?
CompDensity <- function(pars,control,FUN,func.type,
                        measurement=NULL,...){

  ## Should be guaranteed by dream
  ## stopifnot(!is.null(measurement) || func.type%in% c("posterior.density","logposterior.density"))
  
  stopifnot(!any(is.na(pars)))
  
  ## dimensions:
  ##  i. iter 1:nseq
  ##  modpred. scalar or vector commensurate to measurement$data
  ##  err. vector of same length as modpred
  ##  SSR scalar

  p <- rep(NA,control$nseq)
  logp <- rep(NA,control$nseq)
  
  ## Sequential evaluation
  for (ii in 1:nrow(pars)){
    ## Call model to generate simulated data
    ## TODO: correct use of optional pars?
    modpred <- FUN(pars[ii,],...)

    switch(func.type,
           ## Model directly computes posterior density
           posterior.density={
             p[ii] <- modpred
             logp[ii] <- log(modpred)
           },
           ## Model computes output simulation           
           calc.loglik={
             err <- as.numeric(measurement$data-modpred)
                 
             ## Derive the log likelihood
             logp[ii] <- measurement$N*log(control$Wb/measurement$sigma)-
               control$Cb*(sum((abs(err/measurement$sigma))^(2/(1+control$gamma))))
             ## And retain in memory
             p[ii] <- logp[ii]
           },
           ## Model computes output simulation
           ## TODO: may need as.numeric
           calc.rmse={
             
             err <- as.numeric(measurement$data-modpred)
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
           ## TODO: identical to rmse because difference is in metrop
           calc.weighted.rmse={
             ## Define the error
             err <- as.numeric(measurement$data-modpred)
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p[ii] <- -SSR
             logp[ii] <- -0.5*SSR
           }) ##switch
  }           ## for rows

  stopifnot(!any(is.na(p)))
  ##stopifnot(!any(is.na(logp))) ##Not used anyway
  return(list(p=p,logp=logp))
} ##CompDensity
