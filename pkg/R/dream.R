

dreamDefaults <- function()
    list(
         nCR = 3,                ## Crossover values used to generate proposals (geometric series)
         gamma = 0,              ## Kurtosis parameter Bayesian Inference Scheme
         steps = 10,             ## Number of steps in sem
         eps = 5e-2,             ## Random error for ergodicity
         outlierTest = 'IQR_test', ## What kind of test to detect outlier chains?
         pCR.Update = TRUE,      ## Adaptive tuning of crossover values
         boundHandling = 'reflect', ## Boundary handling: "reflect", "bound", "fold", "none"
### Termination criteria. TODO: are the 2nd two valid, given that ndraw is used in adaptive pcr
         ndraw = 1e5,            ## maximum number of function evaluations
         maxtime = Inf,           ## maximum duration of optimization in seconds
         Rthres=1.01,            ## R value at which to stop. Vrugt suggests 1.2
### Thinning
         thin.t=NA,            ## parameter for reduced sample collection
### Reporting
         REPORT = 1000,            ## approximate number of function evaluations between reports. >0. 0=none  TODO: when trace >= 1
### Parameters with auto-set values
         ndim=NA,			 ## number of parameters (automatically set from length of pars)
         DEpairs = NA,          ## Number of DEpairs. defaults to max val floor((nseq-1)/2)
         nseq = NA,              ## Number of Markov Chains / sequences (defaults to N)
         Cb=NA,Wb=NA
         ## Currently unused parameters
##         trace = 0,              ## level of user feedback
         )

library(coda)



##' @param FUN model function with first argument a vector of parameter values of length ndim
##' @param func.type type of value FUN returns.
##'  one of: posterior.density, logposterior.density,calc.loglik,calc.rmse,calc.weighted.rmse
##' @param pars a list of variable ranges
##' @param INIT f(pars,nseq,...) returns nseq x ndim matrix of initial parameter values
##' @param control
##' @param measurement list. must be included unless func.type=posterior.density or logposterior.density is selected
##'  for calc.rmse: must have element data

##' @return ...
##'   TODO
##'   X converged nseq points in parameter space. matrix nseq x ndim
##'   Sequences mcmc.list object. nseq mcmc elements of ndim variables
##'   Reduced.Seq mcmc.list object. nseq mcmc elements of ndim variables
##'   AR acceptance rate for each draw. matrix max.counter x 2
##'   outlier vector of variable length
##'   R.stat Gelman.Diag statistic for each variable at each step. matrix max.counter/steps x 1+ndim
##'   CR. Probability of crossover. matrix max.counter/steps x 1+length(pCR)
##'

##' Terminates either when control$ndraw or control$maxtime is reached
##'
##' MATLAB function:
##' function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option)

dream <- function(FUN, func.type,pars,
                  FUN.pars=list(),
                  INIT = LHSInit,
                  INIT.pars=list(),
                  control = list(),
                  measurement=NULL
                  )
{

  
  ## dimensions
  ##  hist.logp matrix. ndraw/nseq x nseq. length nearly ndraw.
  ##    TODO: removed counter.fun.evals for simplicity. should have been kept?
  ##  CR nseq x steps
  ##  pCR length nCR or scalar
  ##  lCR length nCR or scalar
  ##  Table.JumpRate ndim x DEpairs. range (0,~1.683]
  ##  delta.tot vector of length nCR

  ## Sequences. array max.counter*1.125 x ndim+2 x nseq
  ## Reduced.Seq array max.counter*1.125 x ndim+2 x nseq

  ##  counter.outloop [2,ndraw/nseq]. count number of outside loops
  ## counter.fun.evals [nseq,ndraw(+steps*nseq)],
  ## counter [2,ndraw/nseq] . number of generations - iterations of inner loop
  ## iloc

############################
  ## Process parameters

  ## Check validity of parameters
  if (is.character(FUN))  FUN <- get(FUN, mode = "function")
  stopifnot(is.function(FUN))
  stopifnot(is.list(pars))
  stopifnot(length(pars) > 0)
  pars <- lapply(pars, function(x) if (is.list(x)) x else list(x))
  stopifnot(is.list(control))
  stopifnot(func.type %in% c("calc.rmse","calc.loglik","calc.weighted.rmse","posterior.density","logposterior.density"))
  stopifnot(!is.null(measurement) || func.type %in% c("posterior.density","logposterior.density"))
  stopifnot(!func.type %in% c("calc.rmse","calc.loglik","calc.weighted.rmse") || "data" %in% names(measurement))
  
  ## Check INIT and FUN have required extra parameters in INIT.pars & FUN.pars
  req.args.init <- names(formals(INIT))
  req.args.FUN <- names(formals(FUN))

  if(!all(req.args.init %in% c("pars","nseq",names(INIT.pars)))) stop(paste(c("INIT Missing extra arguments:",req.args.init[!req.args.init %in% c("pars","nseq",names(INIT.pars))]),sep=" "))
  if(!all(req.args.FUN[2:length(req.args.FUN)] %in% c(names(FUN.pars)))) stop(paste(c("FUN Missing extra arguments:",req.args.FUN[!req.args.FUN %in% c("x",names(FUN.pars))]),sep=" "))
  
  ## Update default settings with supplied settings

  control <- modifyList(dreamDefaults(), control)
  isValid <- names(control) %in% names(dreamDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))

  if (!is.null("measurement")){
    if (! "sigma" %in% names(measurement)) measurement$sigma <- sd(measurement$data)
    if (! "N" %in% names(measurement)) measurement$N <- length(measurement$data)
  }

  ## Set automatically determined values
  control$ndim<-length(pars)
  if (is.na(control$nseq)) control$nseq <- control$ndim
  if (is.na(control$DEpairs)) control$DEpairs <- floor((control$nseq-1)/2)

  ## Correct to match nseq
  control$REPORT <- (control$REPORT%/%control$nseq) * control$nseq


  ## Check validity of settings
  if (control$DEpairs==0) stop("control$DEpairs set to 0. Increase nseq?")
  stopifnot(control$DEpairs<=(control$nseq-1)/2) ## Requirement of offde
  stopifnot(control$boundHandling %in% c("reflect", "bound", "fold", "none")) 
  if (control$boundHandling == 'none') warning("No bound handling in use, parameters may cause errors elsewhere")
  stopifnot(control$REPORT>=0)
  
############################
  ## Initialize variables

  NDIM <- control$ndim
  NCR <- control$nCR
  NSEQ <- control$nseq
  
  ## Counters
  counter.fun.evals <- NSEQ
  counter <- 2
  iloc <- 1
  counter.outloop <- 2
  counter.thin <- 1

  ## Max number of times through loops
  max.counter <- ceiling((control$ndraw+control$steps*NSEQ)/NSEQ)+1
  max.counter.outloop <- ceiling((control$ndraw+control$steps*NSEQ)/NSEQ/control$steps)+1
  
  ## Calculate the parameters in the exponential power density function of Box and Tiao (1973)
  cbwb <- CalcCbWb(control$gamma)
  control$Cb <- cbwb$Cb
  control$Wb <- cbwb$Wb

  ## Generate the Table with JumpRates (dependent on number of dimensions and number of pairs)
  Table.JumpRate<-matrix(NA,NDIM,control$DEpairs)
  for (zz in 1:control$DEpairs) Table.JumpRate[,zz] <- 2.38/sqrt(2*zz*1:NDIM)
  
  ## Initialize the array that contains the history of the log_density of each chain
  hist.logp<-matrix(NA,max.counter,NSEQ)
  
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

  obj <- list()
  class(obj) <- c("dream", class(obj))
  obj$call <- match.call()
  obj$control <- control

  EXITFLAG <- NA
  EXITMSG <- NULL

  ## counter.fun.evals + AR at each step
  obj$AR<-matrix(NA,max.counter,2)
  obj$AR[1,2]<-NSEQ-1 ##Number if only one rejected
  colnames(obj$AR) <- c("fun.evals","AR")
  
  ##counter.fun.evals + R statistic for each variable at each step
  ## TODO: now using counter.report
  obj$R.stat<-matrix(NA,max.counter.outloop,NDIM+1)
  ##  n<10 matlab: -2 * ones(1,MCMCPar.n);
  obj$R.stat[1,] <- c(counter.fun.evals,rep(-2,NDIM))
  if (!is.null(names(pars))) colnames(obj$R.stat) <- c("fun.evals",names(pars))
  else   colnames(obj$R.stat) <- c("fun.evals",paste("p",1:length(pars),sep=""))

  ##counter.fun.evals + pCR for each CR
  obj$CR <- matrix(NA,max.counter.outloop,length(pCR)+1)
  colnames(obj$CR) <- c("fun.evals",paste("CR",1:length(pCR),sep=""))

  obj$outlier<-NULL
  
  Sequences <- array(NA, c(floor(1.25*max.counter),NDIM+2,NSEQ))
  if (!is.null(names(pars))) colnames(Sequences) <- c(names(pars),"p","logp")
  ## Sequences[1,] <- sapply(pars, mean) ## TODO: include?

  ## Check whether will save a reduced sample
  if (!is.na(control$thin.t)){
    iloc.2 <- 0
    Reduced.Seq <- array(NA,c(floor(max.counter/control$thin.t),NDIM+2,NSEQ))
  } else Reduced.Seq <- NULL

############################
  
  ## Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
  control$steps<-control$steps-1
  
  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0
  counter.report <- 1

################################
  
  ## Step 1: Sample s points in the parameter space

  x <- do.call(INIT,modifyList(INIT.pars,list(pars=pars,nseq=NSEQ)))

  ## make each element of pars a list and extract lower / upper
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  ##Step 2: Calculate posterior density associated with each value in x
  tmp<-do.call(CompDensity,modifyList(FUN.pars,list(pars=x,control=control,FUN=FUN,func.type=func.type,measurement=measurement)))

  ##Save the initial population, density and log density in one list X
  X<-cbind(x=x,p=tmp$p,logp=tmp$logp)
  if (!is.null(names(pars))) colnames(X) <- c(names(pars),"p","logp")
    
  ##Initialise the sequences
  for (qq in 1:NSEQ){
    Sequences[1,,qq] <- X[qq,]
  }

  ##Save N_CR in memory and initialize delta.tot
  obj$CR[1,] <- c(counter.fun.evals,pCR)
  delta.tot <- rep(0,NCR)
  
  ##Save history log density of individual chains
  hist.logp[1,] <- X[,"logp"]
  

################################
  ##Start iteration
  ## max times (ndraw+steps*NSEQ)/(NSEQ*steps)
  while (counter.fun.evals < control$ndraw) {

    ## max times ceiling(ndraw+steps*NSEQ)/NSEQ
    for (gen.number in 1:control$steps) {

      ## Initialize DR properties
      counter.thin <- counter.thin + 1 

      ## Define the current locations and associated posterior densities
      x.old <- X[,1:NDIM]
      p.old <- X[,NDIM+1]
      logp.old <- X[,NDIM+2]

      ## Now generate candidate in each sequence using current point and members of X
      ## Table.JumpRate appears to match matlab version
      tmp <- offde(x.old, control = control,
                   CR=CR[,gen.number],
                   lower = lower, upper = upper,
                   Table.JumpRate=Table.JumpRate)
      x.new <- tmp$x.new
      stopifnot(!identical(x.new,x.old))
      CR[,gen.number] <- tmp$CR

      ## Now compute the likelihood of the new points
      tmp<-do.call(CompDensity,modifyList(FUN.pars,list(pars=x.new,control=control,FUN=FUN,func.type=func.type,measurement=measurement)))
      p.new <- tmp$p
      logp.new <- tmp$logp

      ## Now apply the acceptance/rejectance rule for the chain itself
      tmp <- metrop(x.new,p.new,logp.new,
                    x.old,p.old,logp.old,
                    func.type,control,
                    measurement
                    )
      newgen <- tmp$newgen
      alpha12 <- tmp$alpha
      accept <- tmp$accept
      ## stopifnot(any(accept)) #Unlikely, but possible

      ## NOTE: original MATLAB code had option for DR Delayed Rejection here)
      ## accept2,ItExtra not required

      ## Update location in sequence and update the locations of the Sequences with the current locations
      iloc <- iloc + 1
      Sequences[iloc,,] <- t(newgen)

      ## Check reduced sample collection
      if (!is.na(control$thin.t) && counter.thin == control$thin.t){
        ## Update iloc_2 and counter.thin
        iloc.2 <- iloc.2+1
        counter.thin <- 0
        ## Reduced sample collection
        Reduced.Seq[iloc.2,,] <- t(newgen)
      }

      ## And update X using current members of Sequences
      X <- newgen; rm(newgen)

      
      if (control$pCR.Update) {
        ## Calculate the standard deviation of each dimension of X
        ## TODO: matlab syntax is unclear - seems to be columnwise
        ## element-wise: sd(c(X[,1:NDIM]))
        r <- apply(X[,1:NDIM],2,sd)
        ## Compute the Euclidean distance between new X and old X
        delta.normX <- rowSums(((x.old-X[,1:NDIM])/r)^2)
        ## Use this information to update sum_p2 to update N_CR
        delta.tot <- CalcDelta(NCR,delta.tot,delta.normX,CR[,gen.number])

        ##0s in delta.tot, delta.normX -> pCR has NaN in pCR.Update
        stopifnot(any(delta.tot!=0) | any(delta.normX!=0)) 

      }

      ## Update hist.logp
      hist.logp[counter,] <- X[,NDIM+2]
      
      ## Save Acceptance Rate
      obj$AR[counter,] <- c(counter.fun.evals,100 * sum(accept) / NSEQ)

      ## CompDensity executes function NSEQ times per loop
      counter.fun.evals <- counter.fun.evals + NSEQ
      counter <- counter + 1
    } ##for gen.number steps

    ## ---------------------------------------------------------------------

    ## Store Important Diagnostic information -- Probability of individual crossover values
    obj$CR[counter.outloop, ] <- c(counter.fun.evals,pCR)

    ## Do this to get rounded iteration numbers
    if (counter.outloop == 2) control$steps <- control$steps + 1

    ## Check whether to update individual pCR values
    if (counter.fun.evals <= 0.1 * control$ndraw) {
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
        hist.logp[1:(counter-1),out.id] <- hist.logp[1:(counter-1),r.idx]
        ## Jump outlier chain to r_idx -- Sequences
        Sequences[iloc,1:(NDIM+2),out.id] <- X[r.idx,]
        ## Jump outlier chain to r_idx -- X
        X[out.id,1:(NDIM+2)] <- X[r.idx,]
        ## Add to chainoutlier
        obj$outlier <- rbind(obj$outlier,c(counter.fun.evals,out.id))
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
    ## TODO: alternatively, convert matlab implementation
    if (control$REPORT>0 && counter.fun.evals %% control$REPORT==0) {

      counter.report <- counter.report+1
            
      try(
          obj$R.stat[counter.report,] <- c(counter.fun.evals,gelman.diag(
                    as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[1:iloc,1:NDIM,i]))),
                       autoburnin=TRUE)$psrf[,1])
        )
      
      if (all(!is.na(obj$R.stat[counter.report,])) &&
          all(obj$R.stat[counter.report,-1]<control$Rthres)) {
        obj$EXITMSG <- 'Convergence criteria reached'
        break
        ## obj$EXITFLAG <- 3
      }

    }##counter.report
    
    ## Update the counter.outloop
    counter.outloop = counter.outloop + 1

      
    ## break if maximum time exceeded
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      obj$EXITMSG <- 'Exceeded maximum time.'
      obj$EXITFLAG <- 2
      break
    }

  } ##while

  toc <- as.numeric(Sys.time()) - tic

  if (counter.fun.evals>= control$ndraw){
    obj$EXITMSG <- "Maximum function evaluations reached"
   ## obj$EXITFLAG <- 4
  }

  obj$X <- X

  ## Trim outputs to collected data - remove extra rows
  ## Convert sequences to mcmc objects
  Sequences <- Sequences[1:iloc,,]
  obj$Sequences <- as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[,1:NDIM,i])))
  if (!is.na(control$thin.t)){
    Reduced.Seq <- Reduced.Seq[1:iloc.2,,]
    obj$Reduced.Seq <- as.mcmc.list(lapply(1:NSEQ,function(i) mcmc(
                                                               Reduced.Seq[,1:NDIM,i],
                                                                   start=1,
                                                                   end=iloc,
                                                                   thin=control$thin.t)
                                           ))
  }

  obj$R.stat <- obj$R.stat[1:counter.report,,drop=FALSE]
  obj$AR <- obj$AR[1:(counter-1),]
  obj$CR <- obj$CR[1:(counter.outloop-1),]
  
  ## store number of function evaluations
  ## store number of iterations
  obj$fun.evals <- counter.fun.evals
  ## store the amount of time taken
  obj$time <- toc

  obj
} ##dream


