unloadNamespace("dream")
library(dream)

## n-dimensional banana shaped gaussian distribution
## Hyperbolic shaped posterior probability distribution
## @param x vector of length ndim
## @param bpar banana-ness. scalar.
## @param imat matrix ndim x ndim
## @return SSR. scalar
## Cursorily verified to match matlab version
Banshp <- function(x,bpar,imat){
  x[2] <- x[2]+ bpar * x[1]^2 - 100*bpar
  S <- -0.5 * (x %*% imat) %*% as.matrix(x)
  return(S)
}

control <- list(
                ndim=10,
                DEpairs=3,
                gamma=0,
                nCR=3,
                ndraw=1e5,
                steps=10,
                eps=5e-2,
                outlierTest='IQR_test',
                pCR.Update=TRUE,
                thin.t=10,
                boundHandling='none'
                )


pars <- lapply(1:control$ndim,function(x) c(-Inf,Inf))
names(pars) <- paste("b",1:control$ndim,sep="")
cmat <- diag(control$ndim)
cmat[1,1] <- 100
muX=rep(0,control$ndim)
qcov=diag(control$ndim)*5

FUN=Banshp
func.type="logposterior.density"
FUN.pars=list(
  imat=solve(cmat),
  bpar=0.1
  )


