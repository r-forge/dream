## control requires gamma, (N)
calc.rmse <- function(predicted,observed,control){
  ## TODO: more general case with multidim data?
  if (! "N" %in% names(control)) control$N <- length(observed)
  
  err <- as.numeric(observed-predicted)
  
  ## Derive the sum of squared error
  SSR <- sum(abs(err)^(2/(1+control$gamma)))

  ## TODO: Correctness: is this a p or logp?
  p <- -SSR
  logp <- -0.5*SSR

  ## Computational issues
  (-control$N*(1+control$gamma)/2)*log(p)
}## calc.rmse

## TODO: calc.weighted.rmse

## control: Wb,Cb,gamma,measurement.sigma
calc.loglik <- function(predicted,observed,control){

  if (!all(c("Cb","Wb") %in% names(control))) control <- modifyList(control,CalcCbWb(control$gamma))
  if (! "sigma" %in% names(control)) control$sigma <- sd(observed)
  if (! "N" %in% names(control)) control$N <- length(observed)

  err <- as.numeric(observed-predicted)
  
  ## Derive the log likelihood
  logp <- control$N*log(control$Wb/control$sigma)-
    control$Cb*(sum((abs(err/control$sigma))^(2/(1+control$gamma))))
  ## And retain in memory
  p <- logp
}## calc.loglik

## Design decisions:
##  lik.fun is log.likelihood=f(predicted,observed,control)
dreamCalibrate <- function(fun,
                           pars,
                           obs,
                           lik.fun=calc.rmse,
                           lik.control=list(),
                           ... ##Extra arguments to dream
                           ){
  
  FUN <- function(pars,...) lik.fun(fun(pars,...),obs,lik.control)
  dd <- dream(FUN=fun,
              pars=pars,
              func.type="logposterior.density",
              ... ##INIT,control,FUN.pars
              )
  
  class(dd) <- c("dream-model",class(dd))
  dd
} ##dreamCalibrate
