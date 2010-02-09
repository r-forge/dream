

dreamDefaults <- function()
    list(nseq = NA,              ## Number of Markov Chains / sequences (defaults to N)
         nCR = 3,                ## Crossover values used to generate proposals (geometric series)
         gamma = 0,              ## Kurtosis parameter Bayesian Inference Scheme
         DEpairs = 3,            ## Number of DEpairs
         steps = 10,             ## Number of steps in sem
         eps = 5e-2,             ## Random error for ergodicity
         outlierTest = 'IQR_test', ## What kind of test to detect outlier chains?
         pCR.Update = TRUE,      ## Adaptive tuning of crossover values
         boundHandling = 'reflect', ## Boundary handling: "reflect", "bound", "fold", "none"
         ndraw = Inf,            ## maximum number of iterations
         maxeval = Inf,          ## maximum number of function evaluations
         maxtime = 60,           ## maximum duration of optimization in seconds
         trace = 0,              ## level of user feedback
         REPORT = 10,            ## number of iterations between reports when trace >= 1
         thin=FALSE,              ## do reduced sample collection
         thin.t=NA               ## parameter for reduced sample collection
         )            

## MATLAB function:
# function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option)

dream <- function(FUN, pars = list(x = range(0, 1e6)),
                  ...,
                  INIT = LHSInit,
                  control = list()) #, do.SSE = FALSE)
{
  If (is.character(FUN))
  FUN <- get(FUN, mode = "function")
  stopifnot(is.function(FUN))
  stopifnot(is.list(par))
  stopifnot(length(par) > 0)

  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(dreamDefaults(), control)
  isValid <- names(control) %in% names(dreamDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))

  ## Initialize variables
  
  ## determine number of variables to be optimized
  control$ndim<-length(par)

  ## for each iteration...
  Iter <- nseq                          #? 1
  counter <- 2
  iloc <- 1
  teller <- 2
  new_teller <- 1
  
  ## Calculate the parameters in the exponential power density function of Box and Tiao (1973)
  cbwb <- CalcCbWb(control$gamma)
  control$Cb<-cbwb$cb
  control$Wb <- cbwb$wb
  
  ## Initialize the array that contains the history of the log_density of each chain
  hist.logp<-matrix(NA,floor(control$ndraw/control$nseq),control$nseq)
  
  ## Derive the number of elements in the output file
  n.elem<-floor(control$ndraw/control$seq)+1

## the output object
  obj <- list()
  class(obj) <- c("dream", class(obj))
  obj$call <- match.call()
  obj$control <- control

  EXITFLAG <- NA
  EXITMSG <- NULL
 
  ## Initialize output information -- AR
  obj$AR<-matrix(0,n.elem,2)
  obj$AR[1,1:2]<-control$nseq-1
  
  ## Initialize output information -- Outlier chains
  obj$outlier<-c()
  
  ## Initialize output information -- R statistic
  obj$R.stat<-matrix(floor(n.elem/control$steps),control$ndim+1)

  if (control$pCR.Update){
    ## Calculate multinomial probabilities of each of the nCR CR values
    pCR <- (1/control$nCR)*rep(1,control$nCR)
    
    ## Calculate the actual CR values based on p
    CR <- GenCR(control,pCR)
    lCR <- rep(0,control$nCR)
  } else {
    pCR <- 1/control$NCR
    ## Define
    CR <- pCR*matrix(1,control$nseq,control$steps)
    lCR <- control$nseq*control$steps
  } ##pCR.Update
  
  ## Initialize output information -- N_CR
  obj$CR <- matrix(0,floor(n.elem/control$steps),length(pCR)+1)

  if (isTRUE(control$saveInMemory)) {
    ## TODO: support this case?
    Sequences <- array(NA, c(floor(1.25*n.elem),control$ndim+1,control$nseq))
    ## TODO: embellishments. accurate?
    if (!is.null(names(pars))) colnames(Sequences) <- names(pars)
    Sequences[1,] <- sapply(pars, mean)
  } else {
    Sequences<-NULL
  }

  ## Generate the Table with JumpRates (dependent on number of dimensions and number of pairs)
  Table.JumpRate<-matrix(NA,control$ndim,control$DEpairs)
  for (zz in 1:control$DEpairs) Table.JumpRate[,zz] <- 2.38/sqrt(2*zz*1:control$ndim)

  ## Check whether will save a reduced sample
  # TODO
  
  ## Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
  control$steps<-control$steps-1
  
  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0

  
  ## Step 1: Sample s points in the parameter space
  x <- INIT(pars, nseq,...)

  ## make each element of pars a list and extract lower / upper
  pars <- lapply(pars, function(x) if (is.list(x)) x else list(x))
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  #Step 2: Calculate posterior density associated with each value in x
  X<-CompDensity(x, control = control, FUN = FUN, ...)
  #Save the initial population, density and log density in one list X
  X<-cbind(x=x,p=X$p[,1],logp=X$logp)
    
  #Initialise the sequences
  for (qq in 1:control$nseq){
    Sequences[1,1:(control$ndim+2),qq] <- X[qq,]
  }

  #Save N_CR in memory and initialize delta_tot
  obj$CR[1,1:(length(pCR)+1)] <- c(Iter,pCR)
  delta.tot <- rep(NA,control$nCR)
  
  #Save history log density of individual chains
  hist.logp[1,1:control$nseq] <- c(Iter, X$logp)
  
  #Compute R-statistic. Using coda package
  # TODO: gelman.diag takes mcmc.list as input
  obj$R.stat[1,1:(control$ndim+1)] <- c(Iter,gelman.diag(Sequences(1:iloc,1:control$ndim,1:control$nseq)))
  
  ##Start iteration
  while (Iter < control$ndraw) {

    for (gen.number in 1:control$steps) {

      ## Initialize DR properties
      new_teller <- new_teller + 1 ## counter for thinning

      ## Define the current locations and associated posterior densities
      x.old <- X[,1:control$ndim]
      p.old <- X[,1:(control$ndim+1)]
      logp.old <- X[,1:(control$ndim+2)]

      ## Now generate candidate in each sequence using current point and members of X
      tmp <- offde(x.old, control = control,
                   CR=CR,
                   lower = lower, upper = upper,bound.handling=bound.handling,
                   Table.JumpRate=Table.JumpRate)
      x.new <- tmp$x.new
      CR[,gen.number] <- tmp$CR

      ## Now compute the likelihood of the new points
      tmp <- CompDensity(xnew, control = control, FUN = FUN, ...)
      p.new <- tmp$p
      logp.new <- tmp$logp

      ## Now apply the acceptance/rejectance rule for the chain itself
      tmp <- metrop(x.new,p.new,logp.new,
                    x.old,p.old,logp.old,
                    N,sigma,control,option
                    )
      newgen <- tmp$newgen
      alpha12 <- tmp$alpha
      accept <- tmp$accept
      

      ## NOTE: original MATLAB code had option for DR Delayed Rejection here)
      ## accept2,ItExtra not required

      ## Update location in sequence and update the locations of the Sequences with the current locations
      if (isTRUE(control$saveInMemory)) iloc <- iloc + 1
      ##Sequences[iloc, 1:(NDIM+2), 1:nseq] <- matrix(newgen, NDIM+2, nseq)
      Sequences[1:(NDIM+2), 1:nseq] <- matrix(newgen, NDIM+2, nseq)

      ## Check reduced sample collection
      if (thin && new_teller == thin.t){
        ## Update iloc_2 and new_teller
        iloc.2 <- iloc.2+1
        new_teller <- 0
        ## Reduced sample collection
        Reduced.Seq[iloc.2,1:(control$ndim+2),1:control$nseq] <-
          array(newgen,c(1,control$ndim+2,control$nseq))
      }

      ## And update X using current members of Sequences
      X <- newgen; rm(newgen)

      if (control$pCR.Update) {
        ## Calculate the standard deviation of each dimension of X
        ## TODO: matlab syntax is unclear - seems to be elementwise?
        r <- sd(c(X[,1:control$ndim]))
        ## Compute the Euclidean distance between new X and old X
        delta.normX <- colSums(((x.old-X[,1:control$ndim])/r)^2,2)
        ## Use this information to update sum_p2 to update N_CR
        delta.tot <- CalcDelta(control,delta.tot,delta.normX,CR[,gen.number])
      }

      ## Update hist.logp
      hist.logp[counter,1:(control$nseq+1)] <- c(Iter,X[,control$ndim+2])
      
      ## Save some important output -- Acceptance Rate
      obj$AR[counter] <- 100 * sum(accept) / nseq

      ## Update Iteration and counter
      Iter <- Iter + control$nseq
      counter <- counter + 1
    } ##for gen.number steps

    ## ---------------------------------------------------------------------

    ## Store Important Diagnostic information -- Probability of individual crossover values
    obj$CR[teller, ] <- pCR

    ## Do this to get rounded iteration numbers
    if (teller == 2) steps <- steps + 1

    ## Check whether to update individual pCR values
    if (Iter <= 0.1 * control$ndraw) {
      if (control$pCR.Update) {
        ## Update pCR values
        tmp <- AdaptpCR(CR, delta_tot, lCR, control)
        pCR <- tmp$pCR
        lCR <- tmp$lCR
      }
    } else {
      ## See whether there are any outlier chains, and remove them to current best value of X
      tmp <- RemOutlierChains(X,
                              Sequences[iloc,,drop=FALSE],
                              hist.logp[1:(counter-1),2:(control$nseq+1)],Iter,outlier,control
                              )
      Sequences[iloc,,] <- tmp$Sequences
      X <- tmp$X
      hist.logp[1:(counter-1),2:(control$nseq+1)] <- tmp$hist.logp
      obj$outlier <- tmp$outlier
    }

    if (control$pCR.Update) {
      ## Generate CR values based on current pCR values
      CR <- GenCR(pCR, control = control)
    } else {
      CR <- matrix(pCR, nrow = nseq, ncol = steps)
    }

    ## Calculate Gelman and Rubin convergence diagnostic
    ## Compute the R-statistic using 50% burn-in from Sequences
    ## TODO: NQR
    obj$R.stat <- gelman.diag(Sequences[seq(iloc %/% 2, iloc), control$ndim, control$nseq))

    ## break if maximum time exceeded
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      EXITMSG <- 'Exceeded maximum time.'
      EXITFLAG <- 2
      break
    }

    ## Update the teller
    teller = teller + 1
  } ##while

  ## Postprocess output from DREAM before returning arguments
  ## Remove extra rows from Sequences
  i <- which(rowSums(Sequences[,,1])==0)
  if (length(i)>0) {
    i <- i[1]-1
    obj$Sequences <- Sequences[1:i,,]
  }
  ## Remove extra rows from Reduced.Seq
  i <- which(rowSums(Reduced.Seq)==0)
  if (length(i)>0) {
    i <- i[1]-1
    obj$Reduced.Seq <- Reduced.Seq[1:i,,]
  }
  ##Remove extra rows R.stat
  i <- which(rowSums(obj$R.stat[,,1])==0)
  if (length(i)>0) {
    i <- i[1]-1
    obj$R.stat <- obj$R.stat[1:i,,]
  }
  ##Remove extra rows AR
  i <- which(rowSums(obj$AR)==0)
  if (length(i)>0) {
    i <- i[1]-1
    obj$AR <- obj$AR[1:i,]
  }
  ## Remove extra rows from CR
  i <- which(rowSums(obj$CR)==0)
  if (length(i)>0) {
    i <- i[1]-1
    obj$CR <- obj$CR[1:i,]
  }
  
  ## store number of function evaluations
  obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc

  obj
}##dream


