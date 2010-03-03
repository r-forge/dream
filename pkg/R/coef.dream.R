##' Extract maximum likelihood parameter values
##' @param dream object
##' @param last.prop proportion of total sequence to use (0,1]
##'  if 1, use whole sequence
##' @param method. either a function or one of maxLik,mean,median
##' @return named vector of parameter values
coef.dream <- function(object,last.prop=.5,use.thinned=FALSE,
                       method=c("maxLik","mean","median"),...){
  stopifnot(last.prop>0)
  if (use.thinned) sss <- object$Reduced.Seq
  else sss <- object$Sequences

  if (class(method)!="function") {
    method <- switch(
                     match.arg(method),
                     "maxLik"=maxLikPars,
                     "mean"=function(sss) colMeans(as.matrix(sss)),
                     "median"=function(sss) apply(as.matrix(sss),2,median)
                     )
  }
  
  if (last.prop==1) return(method(sss))
  else {
    ss <- window(object$Sequences, start = end(sss)*(1-last.prop) + 1)
    return(method(ss))
  }
}##coef.dream
