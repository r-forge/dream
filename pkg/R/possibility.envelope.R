
##' Obtain prediction confidence intervals for model, around new input values
##' @param dd dream object
##' @param FUN.pars new par values. if missing, use those from dream object
##' @param ndraw number of new iterations
##' @param conf % two-sided confidence interval
## TODO: use caching rather than a second application of FUN to parameter sets
possibility.envelope <- function(dd,FUN.pars=NULL,ndraw=1000,conf=99){
  stopifnot(is.null(FUN.pars) || is.list(FUN.pars))
  
  ## Generate more results from converged chains
  dd$control$REPORT <- 0
  dd$control$Rthres <- 0
  dd$control$ndraw <- ndraw
  if (is.na(dd$control$thin.t)) dd$control$thin.t <- 10
  dd$call$control <- dd$control
  dd$call$INIT <- function(pars,nseq) dd$X[,1:dd$control$ndim]

  print(sprintf("Will require %d function evaluations",ndraw+ndraw/dd$control$thin.t))
  
  ee <- eval(dd$call)

  if (is.null(FUN.pars)) FUN.pars <- eval(dd$call$FUN.pars)
  par.name <- names(formals(eval(dd$call$FUN)))[1]

  ff <- apply(as.matrix(ee$Reduced.Seq),1,
              function(p) {
                FUN.pars[[par.name]] <- p
                do.call(eval(dd$call$FUN),FUN.pars)
              })

  gg <- t(apply(ff,1,quantile,c((100-conf)/200,1-(100-conf)/200)))

  return(gg)
}##possibility.envelope
