##' Continue an existing dream MCMC set of chains
##' @param object. dream object
##' @param nsim. approximate number of function evaluations. default 1000
##' @param seed passed to set.seed before continuing
##' @return a dream object with approximately the requested number of function evaluations
## TODO. extra parameters to set in control?
## TODO: does not seem to yield stationary distribution?
simulate.dream <- function(object,nsim=1000,seed=NULL,...){
  
  ## Generate more results from converged chains
  object$control$REPORT <- 0
  object$control$Rthres <- 0
  object$control$ndraw <- nsim
  object$control$burnin.length <- 0
  if (is.na(object$control$thin.t)) object$control$thin.t <- 10
  object$call$control <- object$control
  object$call$INIT <- function(pars,nseq) object$X[,1:object$control$ndim]

  print(sprintf("Will require %d function evaluations",nsim))

  if(!is.null(seed)) set.seed(seed)
  
  ee <- eval(object$call)
  return(ee)
  
}##simulate.dream

##' Predict values using dream object
##' Predict values using function calibrated by dream, optionally with new data,
##' using various methods of summarising the posterior parameter and output distributions
##' @param object dream object
##' @param newdata. new FUN.pars list. If NULL, use object's.
##' @param method CI or a \code{\link{method}} of coef
##' @param level. Requested two-sided level of confidence. For CI method.
##' @param last.prop Proportion of MCMC chains to keep
##' @param use.thinned Whether to use existing thinned chains
##' @return  data frame with each column corresponding to a returned vector
predict.dream <- function(object,newdata=NULL,
                          method="maxLik",level=0.99,
                          last.prop=0.5,use.thinned=TRUE,...
                          ){

  ## Check and initialise parameters
  stopifnot(is.null(newdata) || is.list(newdata))
  stopifnot(!"CI" %in% method || !is.null(level))
  stopifnot(last.prop>0)
  
###
  ## Fetch function and parameters from dream object
  
  if (is.null(newdata)) newdata <- eval(object$call$FUN.pars)
  par.name <- names(formals(eval(object$call$FUN)))[1]
  
  wrap <- function(p) {
    newdata[[par.name]] <- p
    as.numeric(do.call(eval(object$call$FUN),newdata))
  }
###
  
  ## Predict for desired method(s)

  if (method=="CI"){

    if (use.thinned) sss <- object$Reduced.Seq
    else sss <- object$Sequences
    if (last.prop<1) sss <- window(sss, start = end(sss)*(1-last.prop) + 1)
    
    ff <- apply(as.matrix(sss),1,wrap)
    
    return(t(apply(ff,1,quantile,c((1-level)/2,1-(1-level)/2))))
    
  } else {
    return(wrap(coef(object,method=method,
                     use.thinned=use.thinned,last.prop=last.prop)
                ))
  }

} ##predict.dream
