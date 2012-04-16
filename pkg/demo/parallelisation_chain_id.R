## Example of "snow.chains" parallelisation option
## Intended to allow running external model instances using separate batch files
## Uses FME data, see demo("FME.nonlinear.model")

obs.all <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour

## Define control settings
control <- list(
                ## Number of chains to run in parallel
                nseq=4
                )
## Define parameter ranges
pars <- list(p1=c(0,1),p2=c(0,100))

## Define function to return log posterior density for a parameter set
##  doing something different for each parallel chain
Model.fit <- function(id,pars){
  ##Select out this parameter set
  this.par.set=pars[id,]

  ##Do something to set the parameters for this instance
  ## Because they're running on other nodes, they don't have access to the data
  ##  or functions in this session of R
  x <- c(   28,  55,   83,  110,  138,  225,  375)
  Model.y <- function(p,x) {
    p[1]*x/(x+p[2])
  }

  ##Run the model
  ## If an external model, it could be e.g. a different batch script for each parameter chain
  ans <- Model.y(this.par.set,x)

  ##Do something to collect the answer as log posterior density
  Qs <- ans
  Qobs <- c(0.053,0.06,0.112,0.105,0.099,0.122,0.125)

  res<-Qobs-Qs
  sd.res<-sd(res)
  var.res<-var(res)
  logp<-sum(sapply(res,function(res_k)
                   log(1/sqrt(2*pi)*sd.res)*exp(-res_k^2/(2*var.res))))
  return(logp)
}

################################################################################
## Example here is for SNOW with a socket cluster, which doesn't require any
##   other software or setup
## Parallelisation incurs an overhead and is not worthwhile for fast functions
library(snow)
cl <- makeCluster(2, type = "SOCK")

## Using dream directly, with Model.fit, which returns a likelihood
set.seed(456)
control$parallel <- "snow.chains"
result.dream <- dream(FUN = Model.fit,
            pars = pars,
            func.type = "logposterior.density",
            control=control
            )
stopCluster(cl)

################################################################################
## Exploring results

## Summary of settings
print(result.dream)

## Summary of results of fit
summary(result.dream)

## Selection of key plots of fit
plot(result.dream)

## Return parameters having highest likelihood, see ?coef.dream
coef(result.dream)

################################################################################
## Extract parameters and likelihood

## Extract MCMC object, removing burn.in and thinning, see ?window.dream
mcmc <- window(result.dream)
## Convert likelihood to MCMC object
logp <- as.mcmc(window(result.dream$hist.logp, start = start(mcmc)))

## Calculate effective size of MCMC object, and thin
effsz <- effectiveSize(mcmc)
thin <- 2 * ceiling(max(nrow(mcmc[[1]])/effsz))
mcmc <- window(mcmc, thin = thin)
logp <- window(logp, thin = thin)
## Convert to more usable formats
psets <- as.matrix(mcmc)
objseq <- as.vector(logp)

summary(psets)
summary(objseq)

## Plot log likelihood as a function of parameter value
##   (also known as 'dotty' plot)
##  with vertical lines showing maximum likelihood parameter value
par(mfrow=c(1,2))
plot(psets[,1],objseq,ylim=c(-1,0))
abline(v=coef(result.dream)[1],col="red")
plot(psets[,2],objseq,ylim=c(-1,0))
abline(v=coef(result.dream)[2],col="red")

