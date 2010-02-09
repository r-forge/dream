
## Checks the bounds of the parameters
## in:
##  lower, upper: vectors of length ndim
##  x: matrix nseq x ndim
##  bound.handling: one of reflect, bound, fold, none
handleBounds <- function(x, lower, upper, bound.handling)
{
  if (is.vector(x)) x<-t(x)
  stopifnot(is.matrix(x))
  for (p in 1:ncol(x)){
    ##At each pass, different parameter
    ## Modify points that are below or above bounds
    too.low<-which(x[,p]<lower[p])
    too.high<-which(x[,p]>upper[p])
    switch(bound.handling,
           reflect = {
             x[too.low,p] <- 2*lower[p]-x[too.low,p]
             x[too.high,p] <- 2*upper[p]-x[too.high,p]
           },
           bound = {
             x[too.low,p] <- lower[p]
             x[too.high,p] <- upper[p]
           },
           fold = {
             ## ------- New approach that maintains detailed balance ----------
             x[too.low,p] <- upper[p]-(lower[p]-x[too.low,p])
             x[too.high,p] <- lower[p]+(x[too.high,p]-upper[p])
           },
           none = x,
           stop("unrecognised value of 'bound.handling'")
           )#switch
    if (bound.handling!="none") stopifnot(all(x[,p]>=lower[p] & x[,p]<=upper[p]))
  } ##for p
  return(x)
}#handleBounds
