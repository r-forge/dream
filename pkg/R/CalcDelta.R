## Calculate total normalized Euclidean distance for each crossover value
## in:
##   delta.tot	vector of length nCR ? TODO
##   delta.normX	vector of length ndim ? TODO
##   CR		vector of length nseq ? TODO
CalcDelta <- function(control,delta.tot,delta.normX,CR){

  ## Derive sum_p2 for each different CR value
  for (zz in 1:control$nCR){
    ## Find which chains are updated with zz/MCMCPar.nCR
    ## TODO: possible that floating point error prevents exact comparison?
    idx <- which(CR==zz/control$NCR)
    
    ## Add the normalized squared distance tot the current delta_tot;
    delta.tot[zz] <- delta.tot[zz]+sum(delta.normX[idx])
  } ## for CRs
  return(delta.tot)
} ##CalcDelta
