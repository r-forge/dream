## Finds outlier chains and removes them when needed
## in:
## hist.logp matrix dim ndraw/nseq x nseq = length ndraw
RemOutlierChains <- function(X,Sequences,hist.logp,Iter,outlier,control){
  ## Determine the number of elements of L_density
  idx.end <- nrow(hist.logp)
  idx.start <- floor(0.5*idx.end)
  
  ## Then determine the mean log density of the active chains
  ## vector length nseq
  mean.hist.logp <- colMeans(hist.logp[idx.start:idx.end,])
  
  ## Initialize chain_id and Nid
  chain.id <- NULL
  Nid <- 0
  
  ## Check whether any of these active chains are outlier chains
  switch(control$outlierTest,
         'IQR_test'={
           ## TODO: shouldn't upper.range be Q3+3*IQR?
           
           ## Derive the upper and lower quantile of the data
           q13<-quantile(mean.hist.logp,c(0.75,0.25))
           ## Derive the Inter quartile range
           iqr <- q13[1]-q12[2]
           ## Compute the upper range -- to detect outliers
           upper.range <- q13[2]-2*iqr
           ## See whether there are any outlier chains
           chain.id <- which(mean.hist.logp<upper.range)
           Nid <- length(chain.id)
         },
         'Grubbs_test'={
           ## Test whether minimum log_density is outlier
           G <- (mean(mean.hist.logp)-min(mean.hist.logp))/sd(mean.hist.logp)
           
           ## Determine t-value of one-sided interval
           t2 = tinv(1 - 0.01/control$nseq,control$nseq-2)^2; ## 95% interval
           
           ## Determine the critical value
           Gcrit <- ((control$nseq-1)/sqrt(control$nseq))*sqrt(t2/(control$nseq-2+t2))
           
           ## Then check this
           if (G > Gcrit) { ## Reject null-hypothesis
             chain.id <- which.min(mean.hist.logp)
             Nid <- 1
           }
         },
         'Mahal_test'={
           ## Use the Mahalanobis distance to find outlier chains
           alpha <- 0.01
           upper.range <- ACR(control$ndim,control$nseq-1,alpha)
           ## Find which chain has minimum log_density
           idx <- which.min(mean.hist.logp)
           ## Then check the Mahalanobis distance
           d1 <- mahalanobis(X[idx,control$ndim],X[-idx,control$ndim])
           ## Then see whether idx is an outlier in X
           if (d1>upper.range) {
             chain.id <- idx
             Nid <- 1
           }
         }
         )

  if (Nid>0){
    ## Loop over each outlier chain
    for (qq in 1:Nid){
      ## Draw random other chain -- cannot be the same as current chain
      r.idx <- which.max(mean.hist.logp)
      ## Added -- update hist_logp -- chain will not be considered as an outlier chain then
      hist.logp[,chain.id[qq]] <- hist.logp[,r.idx]
      ## Jump outlier chain to r_idx -- Sequences
      Sequences[1,1:(control$nseq+2),chain.id[qq]] <- X[r.idx,1:(control$nseq+2)]
      ## Jump outlier chain to r_idx -- X
      X[chain.id[qq],1:(control$nseq+2)]] <- X[r.idx,1:(control$nseq+2)]
    ## Add to chainoutlier
    outlier <- rbind(outlier,c(Iter,chain.id[qq]))
  }

return(list(
X=X,
Sequences=Sequences,
hist.logp=hist.logp,
outlier=outlier
))

} ##RemOutlierChains
