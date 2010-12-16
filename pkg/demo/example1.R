library(dream)

## n-dimensional banana shaped gaussian distribution
## Hyperbolic shaped posterior probability distribution
## @param x vector of length ndim
## @param bpar banana-ness. scalar.
## @param imat matrix ndim x ndim
## @return SSR. scalar
Banshp <- function(x,bpar,imat){
  x[2] <- x[2]+ bpar * x[1]^2 - 100*bpar
  S <- -0.5 * (x %*% imat) %*% as.matrix(x)
  return(S)
}
## Output manually verified against matlab version
##Banshp(1:10,bpar,imat) ## Banshp(1:10,Extra)
##Banshp(rep(1,10),bpar,imat) ## Banshp(ones(1,10),Extra)

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
                boundHandling='none',
                burnin.length=Inf, ##for compatibility with matlab code
                REPORT=1e4 ##reduce frequency of progress reports
                )


pars <- replicate(control$ndim,c(-Inf,Inf),simplify=F)
cmat <- diag(control$ndim)
cmat[1,1] <- 100
bpar <- 0.1
imat <- solve(cmat)
muX=rep(0,control$ndim)
qcov=diag(control$ndim)*5

set.seed(11)

dd <- dream(
            FUN=Banshp, func.type="logposterior.density",
            pars = pars,
            FUN.pars=list(
              imat=imat,
              bpar=bpar),
            INIT = CovInit,
            INIT.pars=list(
              muX=muX,
              qcov=qcov,
              bound.handling=control$boundHandling
              ),
            control = control
            )

summary(dd)


## Show bananity
library(lattice)
ddm <- as.matrix(window(dd))
plot(ddm[,1],ddm[,2])
splom(ddm)

## Compare to two matlab runs
fn.example1a <- system.file("extdata/example1a.mat",package="dream")
fn.example1b <- system.file("extdata/example1b.mat",package="dream")

for (fn.example in c(fn.example1a,fn.example1b)){
  compareToMatlab(fn.example,dd)
  mat <- readMat(fn.example)
  all.equal(muX,as.numeric(mat$Extra[,,1]$muX))
  all.equal(qcov,mat$Extra[,,1]$qcov)
  all.equal(imat,mat$Extra[,,1]$imat)
  all.equal(bpar,as.numeric(mat$Extra[,,1]$bpar))
}

## While banana doesn't appear to match, it appears this is because it is a difficult problem
## Separate matlab results do not match either

## Compare matlab runs
matb <- getMatlabSeq(readMat(fn.example1b))
matb <- as.matrix(window(matb,start=1+(end(matb)-1)*(1-0.5)))

mata <- getMatlabSeq(readMat(fn.example1a))
mata <- as.matrix(window(mata,start=1+(end(mata)-1)*(1-0.5)))

pvals <- sapply(1:ncol(mata), function(i) ks.test(mata[,i], matb[, i])$p.value)
round(pvals,2)

plotMCMCQQ(matb,mata)

## Results from matlab version
plot(mata[,1],mata[,2])

splom(mata)
