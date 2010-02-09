## Metropolis rule for acceptance or rejection
metrop<-function(x,p.x,logp.x,
                 x.old,p.old,logp.old,
                 N,sigma,control,option
                 ){

  ## Calculate the number of Chains
  nr.chains <- nrow(x)
  
  ## First set newgen to the old positions in X
  newgen <- cbind(x.old,p.old,logp.old)
  
  ## And initialize accept with false
  accept <- rep(FALSE,nr.chains)
  
  switch(option,
         1={
           alpha <- min(p.x/p.old,1)
         },
         2={ ## Lnp probability evaluation
           alpha <- min(exp(p.x-p.old),1)
         },
         3={ ## SSE probability evaluation
           alpha <- min((p.x/p.old)^(-N*(1+control$gamma)/2),1)
         },
         4={ ## Lnp probability evaluation
           alpha <- min(exp(p.x-p.old),1)
         },
         5={ ## Similar to 3 but now weighted with Measurement.Sigma
           ## signs are different because we write -SSR
           alpha <- min(exp(-0.5*(-p.x + p.old)./sigma^2),1);
         })

  ## Generate random numbers
  Z <- runif(nr.chains)
  ## Find which alpha's are greater than Z
  idx <- which(alpha>Z)
  ## And update these chains
  newgen[idx,] <- cbind(x[idx,],p.x[idx],logp.x[idx])
         
  ## And indicate that these chains have been accepted
  accept[idx] <- TRUE

return(list(newgen=newgen,alpha=alpha,accept=accept))
} ##metrop
