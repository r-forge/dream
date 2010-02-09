## Updates the probabilities of the various crossover values
## in:
## lCR.old vector of length nCR
AdaptpCR <- function(CR,delta.tot,lCR.old,control){
  ## Make CR to be a single vector
  CR <- c(CR)
  
  ## Determine lCR
  lCR <- rep(NA,control$nCR)
  for (zz in 1:control$nCR){
    
    ## Determine how many times a particular CR value is used
    ## TODO: shouldn't this be which(CR==CR[z]])?
    idx <- which(CR==zz/control$nCR)
    
    ## This is used to weight delta.tot
    lCR[zz] <- lCR.old[zz]+length(idx)
  }                                     #for CRs
  
  ## Adapt pCR using information from averaged normalized jumping distance
  pCR <- control$nseq*delta.tot/lCR / sum(delta.tot)
    
  ## Normalize pCR
  pCR <- pCR/sum(pCR)

  return(list(pCR=pCR,lCR=lCR))
} ##AdaptpCR
