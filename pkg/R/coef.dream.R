##' Extract maximum likelihood parameter values
##' @param dream object
##' @last.prop proportion of total sequence to use (0,1]
##'  if 1, use whole sequence
##' @return named vector of parameter values
coef.dream <- function(dream.obj,last.prop=.5,use.thinned=FALSE){
  stopifnot(last.prop>0)
  if (use.thinned) sss <- dream.obj$Reduced.Seq
  else sss <- dream.obj$Sequences
  
  if (last.prop==1) return(maxLikPars(sss))
  else {
    ss <- window(dream.obj$Sequences, start = end(sss)*(1-last.prop) + 1)
    return(maxLikPars(ss))
  }
}##coef.dream
