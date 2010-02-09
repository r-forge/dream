DEStrategy<-function(control){
  ## Determine which sequences to evolve with what DE strategy
  ## in: control$ DEpairs, nseq
  ## out: DEversion vector of length control$nseq
  
  ## Determine probability of selecting a given number of pairs
  p.pair <- (1/control$DEpairs)*rep(1,control$DEpairs)
  p.pair <- c(0,cumsum(p.pair))
  
  ## Generate a random number between 0 and 1
  Z <- runif(control$nseq)
  
  ## Select number of pairs
  DEversion<-rep(0,control$nseq)
  for (qq in 1:control$nseq){
    z <- which(Z[qq]>p.pair)
    DEversion[qq]<-tail(z,1)
  }
  return(DEversion)
}#DEStrategy
