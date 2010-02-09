CompDensity <- function(x,control,FUN,option,
                        measurement,...){
  ## This function computes the density of each x value
  ## option is unambiguous abbrev of:
  ##   posterior.density, calc.loglik, calc.rmse, logposterior.density,calc.weighted.rmse
  
  ## Sequential evaluation
  for (ii in 1:nrow(x)){
    ## Call model to generate simulated data
    ## TODO: correct use of optional pars?
    modpred <- FUN(...)

    switch(option,
           ## Model directly computes posterior density
           posterior.density={
             p[ii,1:2] <- c(ModPred,ii)
             logp <- log(p[ii,1])
           },
           ## Model computes output simulation           
           calc.loglik={
             err <- measurement$data-modpred
                 
             ## Compute the number of measurement data
             N <- nrow(measurement)
             
             ## Derive the log likelihood
             logp[ii,1] <- N*log(control$Wb/measurement$sigma)-
               control$CB*(sum((abs(Err/measurement$sigma))^(2/(1+control$gamma))))
             ## And retain in memory
             p[ii,1:2] <- c(logp[ii,1],ii)
           },
           ## Model computes output simulation
           calc.rmse={
             err <- measurement$data-modpred
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p[ii,1:2] <- c(-SSR,ii)
             logp[ii,1] <- -0.5*SSR
             
           },
           ## Model directly computes log posterior density
           logposterior.density={
             p[ii,1:2] <- [modpred ii]
             logp[ii,1] <- p[ii,1]
           },
           ## Similar as 3, but now weights with the Measurement Sigma
           ## TODO: appears to be no difference to calc.rmse
           calc.weighted.rmse={
             ## Define the error
             err <- measurement$data-modpred
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p[ii,2] <- c(-SSR,ii)
             logp[ii,1] <- -0.5*SSR
           }
           )
  }## for rows
return(list(p=p,logp=logp))
}##CompDensity
