## Generate offspring using METROPOLIS HASTINGS monte-carlo markov chain
rand<-function(m,n) matrix(runif(m*n),m,n)
randn<-function(m,n) matrix(rnorm(m*n),m,n)

offde<-function(x.old,control,CR,
                lower,upper,bound.handling, #For bound handling
		Table.JumpRate #TODO: may be better in another accessible scope rather than as a parameter
                ){

  #Generate ergodicity term
  eps <- 1e-6 * randn(control$nseq,control$ndim)

  #Not a delayed rejection step -> generate proposal with DE
  ## Determine which sequences to evolve with what DE strategy
  DEversion <- DEStrategy(control)
  
  ## Generate series of permutations of chains
  # TODO: sort collapses matrix !!!
  #tt<-sort(rand(control$nseq-1,control$nseq))
  m<-control$nseq-1
  n<-control$nseq
  tt<- matrix(sort(runif(m*n)),m,n)
  
  ## Generate uniform random numbers for each chain to determine which dimension to update
  D <-rand(control$nseq,control$ndim)

  ## Ergodicity for each individual chain
  noise.x<-control$eps*(2*rand(control$nseq,control$ndim)-1)
    
  ## Initialize the delta update to zero
  delta.x<-matrix(0,control$nseq,control$ndim)
  
  ## Each chain evolves using information from other chains to create offspring
  for (qq in 1:control$nseq){
    
    ## Define ii and remove current member as an option
    ii <-rep(1,control$nseq)
    ii [qq] <-0
    idx <- which(ii>0)
    
    ## randomly select two members of ii that have value == 1
    rr <- idx[tt[1:2*DEversion[qq,1],qq]]
    
    ## --- WHICH DIMENSIONS TO UPDATE? DO SOMETHING WITH CROSSOVER ----
    i <- which(D[qq,]>(1-CR[qq,1]))

    ## Update at least one dimension
    if (length(i)==0) i <- sample(control$ndim,1)
    
    ## ----------------------------------------------------------------
    ## Determine the number of dimensions that are going to be updated
    NrDim <- length(i)
    
    ## Determine the associated JumpRate and compute the jump
    if (runif(1)<4/5){
      ## Lookup Table
      JumpRate <- Table.JumpRate(NrDim,DEversion[qq,1])
      
      ## Produce the difference of the pairs used for population evolution
      delta <- colSums(x.old[rr[1:DEversion[qq,1]],]-x.old[rr[DEversion[qq,1]+1:2*DEversion[qq,1]]])
      
      ## Then fill update the dimension
      delta.x[qq,i] <- (1+noise.x[qq,i])*JumpRate*delta[1,i]
      
    } else {
      ## Set the JumpRate to 1 and overwrite CR and DEversion
      JumpRate <- 1
      CR[qq,1] <- -1
      
      ## Compute delta from one pair
      delta <- x.old[rr[1],]-x.old[rr[2],]
                       
      ## Now jumprate to facilitate jumping from one mode to the other in all dimensions
      delta.x[qq,] <- JumpRate*delta
      
      ## Check this line to avoid that jump = 0 and xnew is similar to xold
      if (rowSums(delta.x[qq,]^2)==0){
        ## Compute the Cholesky Decomposition of x.old
        R <- (2.38/sqrt(control$nseq))*chol(cov(x.old)+1e-5*diag(control$nseq))
        
        ## Generate jump using multinormal distribution
        delta.x[qq,] <- rnorm(control$ndim)*R
      }
    }#runif
  }#for qq
  
  ## TODO?:
  ## If delayed rejection step --> generate proposal with DR
  ## Loop over all chains -- all dimensions are updated
  ## Generate a new proposal distance using standard procedure

  
  ## Update x_old with delta_x and eps;
  x.new <- x.old+delta.x+eps
  x.new <- handleBounds(x,lower,upper,bound.handling)

  return(list(x.new=x.new,CR=CR))
}#function offde
