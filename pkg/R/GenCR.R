GenCR <- function(pCR,control){
  ## Generates CR values based on current probabilities
  
  ## How many candidate points for each crossover value?
                                        # TODO: rmultinom may not be equivalent to multrnd
  L <- rmultinom(1,size=control$nseq*control$steps,p=pCR)
  L2 <- c(0,cumsum(L))
  
  ## Then select which candidate points are selected with what CR
  r <- sample(control$nseq*control$steps)
  
  ## Then generate CR values for each chain
  for (zz in 1:control$nCR){
    ## Define start and end
    i.start <- L2[zz]+1
    i.end <- L2[zz+1]
    
    ## Take the appropriate elements of r
    idx <- r[i.start:i.end]
    
    ## Assign these indices MCMCPar.CR(zz)
    CR[idx] <- zz/control$nCR
    
    ## Now reshape CR
    CR <- array(CR,control$nseq,control$steps)\
  } ## for nCR
return(CR)
}   ## GenCR
