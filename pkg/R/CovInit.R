## in:
##  muX. vector of length ndim
## out:
##  x. dim nseq x ndim
CovInit <- function(pars, nseq,muX,qcov,bound.handling)
{

  ##[x] = repmat(Extra.muX,MCMCPar.seq,1) + randn(MCMCPar.seq,MCMCPar.n) * chol(Extra.qcov);
  x <- t(matrix(rep(a,nseq),length(a)))+randn(control$nseq,length(pars)) %*% chol(qcov)
  
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  x <- handleBounds(x,lower,upper,bound.handling)
  return(x)
}