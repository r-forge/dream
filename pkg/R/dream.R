

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
         thin=FALSE,             ## do reduced sample collection
         thin.t=NA,               ## parameter for reduced sample collection
         metrop.opt=NA            ## ?? which option to use for metropolis acceptance/rejection. Some also need measurement$sigma
         )            

library(coda)

## MATLAB function:
## function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option)

dream <- function(FUN, func.type,
                  pars = list(x = range(0, 1e6)),
                  ...,
                  INIT = LHSInit,
                  control = list(),
                  measurement=NULL
                  )
{

  
  ## dimensions
  ##  x points in parameter space. matrix nseq x ndim
  ##  hist.logp matrix. ndraw/nseq x nseq. length nearly ndraw.
  ##    TODO: removed Iter for simplicity. should have been kept?
  ##  CR nseq x steps
  ##  pCR length nCR or scalar
  ##  lCR length nCR or scalar
  ##  Table.JumpRate ndim x DEpairs. range (0,~1.683]
  ##  delta.tot vector of length nCR

  
  if (is.character(FUN))
    FUN <- get(FUN, mode = "function")
  stopifnot(is.function(FUN))
  stopifnot(is.list(par))
  stopifnot(length(par) > 0)
  stopifnot(!is.null(measurement) || func.type %in% c("posterior.density","logposterior.density"))


  pars <- lapply(pars, function(x) if (is.list(x)) x else list(x))

############################
  ## Initialize variables
  
  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(dreamDefaults(), control)
  isValid <- names(control) %in% names(dreamDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))

  ## determine number of variables to be optimized
  control$ndim<-length(par)
  if (is.na(control$nseq)) control$nseq <- control$ndim

  NDIM <- control$ndim
  NCR <- control$ncr
  NSEQ <- control$nseq
  
  ## for each iteration...
  Iter <- NSEQ                          #? 1
  counter <- 2
  iloc <- 1
  teller <- 2
  new_teller <- 1
  
  ## Calculate the parameters in the exponential power density function of Box and Tiao (1973)
  cbwb <- CalcCbWb(control$gamma)
  control$Cb <- cbwb$cb
  control$Wb <- cbwb$wb

  ## Generate the Table with JumpRates (dependent on number of dimensions and number of pairs)
  Table.JumpRate<-matrix(NA,NDIM,control$DEpairs)
  for (zz in 1:control$DEpairs) Table.JumpRate[,zz] <- 2.38/sqrt(2*zz*1:NDIM)
  
  ## Initialize the array that contains the history of the log_density of each chain
  hist.logp<-matrix(NA,floor(control$ndraw/NSEQ),NSEQ)
  
  if (control$pCR.Update){
    ## Calculate multinomial probabilities of each of the nCR CR values
    pCR <- rep(1/NCR,NCR)
    
    ## Calculate the actual CR values based on p
    CR <- GenCR(pCR,control)
    lCR <- rep(0,NCR)
  } else {
    pCR <- 1/NCR
    ## Define
    CR <- matrix(pCR,NSEQ,control$steps)
    lCR <- NSEQ*control$steps
  } ##pCR.Update


############################
  ## Initialise output object

  ## dimensions:
  ## AR matrix n.elem x 2
  ## outlier vector of variable length
  ## R.stat matrix n.elem/steps x 1+ndim
  ## CR matrix n.elem/steps x 1+length(pCR)
  ## Sequences array n.elem*1.125 x ndim+2 x nseq
  ## Reduced.Seq array n.elem*1.125 x ndim+2 x nseq
  
  obj <- list()
  class(obj) <- c("dream", class(obj))
  obj$call <- match.call()
  obj$control <- control

  EXITFLAG <- NA
  EXITMSG <- NULL

  ## Derive the number of elements in the output file
  n.elem<-floor(control$ndraw/NSEQ)+1

  ## Iter + AR at each step
  obj$AR<-matrix(NA,n.elem,2)
  obj$AR[1,1]<-NSEQ-1
  
  obj$outlier<-NULL

  ##Iter + R statistic for each variable at each step
  obj$R.stat<-matrix(NA,floor(n.elem/control$steps),NDIM+1)
  
  ##Iter + pCR for each CR
  obj$CR <- matrix(NA,floor(n.elem/control$steps),length(pCR)+1)
  
  Sequences <- array(NA, c(floor(1.25*n.elem),NDIM+2,NSEQ))
  if (!is.null(names(pars))) colnames(Sequences) <- c(names(pars),"p","logp")
  ## Sequences[1,] <- sapply(pars, mean) ## TODO: include?

  ## Check whether will save a reduced sample
  if (control$thin){
    Reduced.Seq <- array(NA,c(floor(n.elem/control$thin.t),NDIM+2,NSEQ))
  } else Reduced.Seq <- NULL

############################
  
  ## Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
  control$steps<-control$steps-1
  
  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0

################################
  
  ## Step 1: Sample s points in the parameter space
  x <- INIT(pars, NSEQ,...)

  ## make each element of pars a list and extract lower / upper
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  ##Step 2: Calculate posterior density associated with each value in x
  tmp<-CompDensity(x, control = control, FUN = FUN, func.type=func.type,measurement=measurement...)
  ##Save the initial population, density and log density in one list X
  X<-cbind(x=x,p=tmp$p,logp=tmp$logp)
  if (!is.null(names(pars))) colnames(X) <- c(names(pars),"p","logp")
    
  ##Initialise the sequences
  for (qq in 1:NSEQ){
    Sequences[1,,qq] <- X[qq,]
  }

  ##Save N_CR in memory and initialize delta.tot
  obj$CR[1,] <- c(Iter,pCR)
  delta.tot <- rep(NA,NCR)
  
  ##Save history log density of individual chains
  hist.logp[1,] <- X[,"logp"]
  
  ##Compute R-statistic. Using coda package
  ## TODO: more elegant way of using coda. And check for correctness
  ## TODO: alternatively, convert matlab implementation
  obj$R.stat[1,] <- c(Iter,gelman.diag(as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[1:iloc,1:NDIM,i])))))

################################
  ##Start iteration
  while (Iter < control$ndraw) {

    for (gen.number in 1:control$steps) {

      ## Initialize DR properties
      new_teller <- new_teller + 1 ## counter for thinning

      ## Define the current locations and associated posterior densities
      x.old <- X[,1:NDIM]
      p.old <- X[,1:(NDIM+1)]
      logp.old <- X[,1:(NDIM+2)]

      ## Now generate candidate in each sequence using current point and members of X
      tmp <- offde(x.old, control = control,
                   CR=CR[,gen.number],
                   lower = lower, upper = upper,
                   Table.JumpRate=Table.JumpRate)
      x.new <- tmp$x.new
      CR[,gen.number] <- tmp$CR

      ## Now compute the likelihood of the new points
      tmp <- CompDensity(x.new, control = control, FUN = FUN, ...)
      p.new <- tmp$p
      logp.new <- tmp$logp

      ## Now apply the acceptance/rejectance rule for the chain itself
      tmp <- metrop(x.new,p.new,logp.new,
                    x.old,p.old,logp.old,
                    measurement,control
                    )
      newgen <- tmp$newgen
      alpha12 <- tmp$alpha
      accept <- tmp$accept
      

      ## NOTE: original MATLAB code had option for DR Delayed Rejection here)
      ## accept2,ItExtra not required

      ## Update location in sequence and update the locations of the Sequences with the current locations
      iloc <- iloc + 1
      Sequences[iloc,,] <- t(newgen)

      ## Check reduced sample collection
      if (control$thin && new_teller == control$thin.t){
        ## Update iloc_2 and new_teller
        iloc.2 <- iloc.2+1
        new_teller <- 0
        ## Reduced sample collection
        Reduced.Seq[iloc.2,,] <- t(newgen)
      }

      ## And update X using current members of Sequences
      X <- newgen; rm(newgen)

      if (control$pCR.Update) {
        ## Calculate the standard deviation of each dimension of X
        ## TODO: matlab syntax is unclear - seems to be elementwise?
        r <- sd(c(X[,1:NDIM]))
        ## Compute the Euclidean distance between new X and old X
        delta.normX <- rowSums(((x.old-X[,1:NDIM])/r)^2)
        ## Use this information to update sum_p2 to update N_CR
        delta.tot <- CalcDelta(NCR,delta.tot,delta.normX,CR[,gen.number])
      }

      ## Update hist.logp
      hist.logp[counter,] <- X[,NDIM+2]
      
      ## Save some important output -- Acceptance Rate
      obj$AR[counter] <- 100 * sum(accept) / NSEQ

      ## Update Iteration and counter
      Iter <- Iter + NSEQ
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
        tmp <- AdaptpCR(CR, delta.tot, lCR, control)
        pCR <- tmp$pCR
        lCR <- tmp$lCR
      }
    } else {
      ## See whether there are any outlier chains, and remove them to current best value of X
      tmp <- RemOutlierChains(X,hist.logp[1:(counter-1),],control)
      ## Loop over each outlier chain (if length>0)
      for (out.id in tmp$chain.id){
        ## Draw random other chain -- cannot be the same as current chain
        r.idx <- which.max(tmp$mean.hist.logp)
        ## Added -- update hist_logp -- chain will not be considered as an outlier chain then
        hist.logp[1:(counter-1),out.id] <- hist.logp[,r.idx]
        ## Jump outlier chain to r_idx -- Sequences
        Sequences[iloc,1:(NSEQ+2),out.id] <- X[r.idx,]
        ## Jump outlier chain to r_idx -- X
        X[out.id,1:(NSEQ+2)] <- X[r.idx,]
        ## Add to chainoutlier
        outlier <- rbind(outlier,c(Iter,out.id))
      } ##for remove outliers
    }   ##else

    if (control$pCR.Update) {
      ## Generate CR values based on current pCR values
      CR <- GenCR(pCR, control)
    } else {
      CR <- matrix(pCR, nrow = NSEQ, ncol = control$steps)
    }

    ## Calculate Gelman and Rubin convergence diagnostic
    ## Compute the R-statistic using 50% burn-in from Sequences
    obj$R.stat <- gelman.diag(as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[,1:NDIM,i]))),autoburnin=TRUE)

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
  if (control$thin){
    i <- which(rowSums(Reduced.Seq)==0)
    if (length(i)>0) {
      i <- i[1]-1
      obj$Reduced.Seq <- Reduced.Seq[1:i,,]
    }
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
  ## TODO: funevals is not calculated atm
  ## obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc

  obj
} ##dream


