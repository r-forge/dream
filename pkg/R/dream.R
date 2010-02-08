

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
         maxit = Inf,            ## maximum number of iterations
         maxeval = Inf,          ## maximum number of function evaluations
         maxtime = 60,           ## maximum duration of optimization in seconds
         trace = 0,              ## level of user feedback
         REPORT = 10)            ## number of iterations between reports when trace >= 1


LHSInit <- function(pars, nseq)
{
    lapply(pars, function(r)
           sample(seq(min(r), max(r), length = nseq))
           )
}
CovInit <- function(pars, nseq)
{

}

## MATLAB function:
# function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option)

dream <- function(FUN, pars = list(x = range(0, 1e6)),
                  ...,
                  INIT = LHSInit,
                  thin = 0, control = list()) #, do.SSE = FALSE)
{
    if (is.character(FUN))
        FUN <- get(FUN, mode = "function")
    stopifnot(is.function(FUN))
    stopifnot(is.list(par))
    stopifnot(length(par) > 0)
    stopifnot(is.numeric(thin))

    ## determine number of variables to be optimized
    NDIM <- length(par)

    ## update default options with supplied options
    stopifnot(is.list(control))
    control <- modifyList(dreamDefaults(), control)
    isValid <- names(control) %in% names(dreamDefaults())
    if (any(!isValid))
        stop("unrecognised options: ",
             toString(names(control)[!isValid]))

    nseq <- control$nseq
    nCR <- control$nCR

    ## Sample s points in the parameter space
    x <- INIT(pars, nseq)

    ## make each element of pars a list and extract lower / upper
    pars <- lapply(pars, function(x) if (is.list(x)) x else list(x))
    lower <- sapply(pars, function(x) min(x[[1]]))
    upper <- sapply(pars, function(x) max(x[[1]]))

    x <- handleBounds(x, lower, upper, control$bound.handling)


    ## initialize variables
    if (isTRUE(control$saveInMemory)) {
        ## support this case?
    } else {
        Sequences <- matrix(as.numeric(NA), nrow = nseq, ncol = NDIM+2)
    }
    if (!is.null(names(pars)))
        colnames(Sequences) <- names(pars)
    Sequences[1,] <- sapply(pars, mean)

    ## the output object
    obj <- list()
    class(obj) <- c("dream", class(obj))
    obj$call <- match.call()
    obj$control <- control

    EXITFLAG <- NA
    EXITMSG <- NULL

    ## initialize timer
    tic <- as.numeric(Sys.time())
    toc <- 0

    ## for each iteration...
    Iter <- nseq #? 1
    counter <- 2
    iloc <- 1
    teller <- 2
    new_teller <- 1
    
    while (Iter < MAXIT) {

        for (gen_number in 1:control$steps) {
            new_teller <- new_teller + 1 ## counter for thinning

            ## Define the current locations and associated posterior densities
            xold <- GetLocation(X, control)

            ## Now generate candidate in each sequence using current point and members of X
            xnew <- offde(xold, X, control = control, lower = lower, upper = upper)

            ## Now compute the likelihood of the new points
            xnew <- CompDensity(xnew, control = control, FUN = FUN, ...)

            ## Now apply the acceptance/rejectance rule for the chain itself
            tmp <- metrop(xnew, xold, ETC)

            ## (NOTE: original MATLAB code had option for DR Delayed Rejection here)

            ## Now update the locations of the Sequences with the current locations
            if (isTRUE(control$saveInMemory))
                iloc <- iloc + 1
            #Sequences[iloc, 1:(NDIM+2), 1:nseq] <- matrix(newgen, NDIM+2, nseq)
            Sequences[1:(NDIM+2), 1:nseq] <- matrix(newgen, NDIM+2, nseq)

            ## Check reduced sample collection
            if (thin > 0) {
                ## TODO
            }

            ## And update X using current members of Sequences
            X <- newgen; rm(newgen)

            if (control$pCR.Update) {
                ## Calculate the standard deviation of each dimension of X

                ## Compute the Euclidean distance between new X and old X

                ## Use this information to update sum_p2 to update N_CR

            }

            ## Update hist_logp

            ## Save some important output -- Acceptance Rate
            obj$AR[counter] <- 100 * sum(accept) / nseq

            ## Update Iteration and counter
            Iter <- Iter + nseq
            counter <- counter + 1
        }

        ## ---------------------------------------------------------------------

        ## Store Important Diagnostic information -- Probability of individual crossover values
        obj$CR[teller, ] <- pCR

        ## Do this to get rounded iteration numbers
        if (teller == 2)
            steps <- steps + 1

        ## Check whether to update individual pCR values
        if (Iter <= 0.1 * MAXIT) {
            if (control$pCR.Update) {
                ## Update pCR values
                tmp <- AdaptpCR(CR, delta_tot, lCR, control)
                pCR <- tmp$pCR
                lCR <- tmp$lCR
            }
        } else {
            ## See whether there are any outlier chains, and remove them to current best value of X
            outliers <- OutlierChains(X, Sequences[iloc,,drop=FALSE],
                                      ..., ETC)
            Sequences[iloc,,] <- Sequences[iloc,,] outliers
            ETC
        }

        if (control$pCR.Update) {
            ## Generate CR values based on current pCR values
            CR <- GenCR(pCR, control = control)
        } else {
            CR <- matrix(pCR, nrow = nseq, ncol = steps)
        }

        ## Calculate Gelman and Rubin convergence diagnostic
        ## Compute the R-statistic using 50% burn-in from Sequences
        obj$Rstat <- Gelman(Sequences[seq(iloc %/% 2, iloc), NDIM, nseq, control))

        ## break if maximum time exceeded
        toc <- as.numeric(Sys.time()) - tic
        if (toc > control$maxtime) {
            EXITMSG <- 'Exceeded maximum time.'
            EXITFLAG <- 2
            break
        }

        ## Update the teller
        teller = teller + 1
    }

    ## store number of function evaluations
    obj$counts <- funevals
    ## store number of iterations
    obj$iterations <- i
    ## store the amount of time taken
    obj$time <- toc

    obj
}

handleBounds <- function(x, lower, upper, bound.handling)
{
    switch(bound.handling,
           reflect = {
           },
           bound = {
               x <- lapply(x, function(r)
                           r <- pmax(pmin(r, upper), lower)
           },
           fold = {
           },
           none = x,
           stop("unrecognised value of 'bound.handling'")
}
