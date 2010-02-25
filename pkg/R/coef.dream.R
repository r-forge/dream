##' Extract maximum likelihood parameter values
##' @param dream object
##' @last.prop proportion of total sequence to use (0,1]
##'  if 1, use whole sequence
##' @return named vector of parameter values
## TODO: potentially use reduced.seq instead
coef.dream <- function(dream.obj,last.prop=.5){
  stopifnot(last.prop>0)
  if (last.prop==1) return(maxLikPars(dream.obj$Sequences))
  else {
    ss <- window(dream.obj$Sequences, start = end(dd$Sequences)*(1-last.prop) + 1)
    return(maxLikPars(ss))
  }
}##coef.dream
