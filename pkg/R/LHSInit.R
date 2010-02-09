## Latin Hypercube sampling
LHSInit <- function(pars, nseq){
  ##TODO: R method?
  ## lapply(pars, function(r)
  ##        sample(seq(min(r), max(r), length = nseq))
  ##        )

  xmin <- sapply(pars, function(x) min(x[[1]]))
  xmax <- sapply(pars, function(x) max(x[[1]]))

  ## Define the size of xmin
  nvar <- length(xmin)
  ## Initialize array ran with random numbers
  ran <- rand(nseq,nvar)
  
  ## Initialize array s with zeros
  s <- matrix(0,nseq,nvar)
  
  ## Now fill s
  for (j in 1:nvar){
    ## Random permutation
    idx <- sample(nseq)
    P <- idx-ran[,j]/nseq
    s[,j] <- xmin[j]+P*(xmax[j]-xmin[j])
  } ##for pars

  return(x)
} ## LHSInit